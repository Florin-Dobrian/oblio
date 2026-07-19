// test_pipeline.cpp - the whole pipeline, ordering through solve, judged by the residual.
//
// This suite answers "does this combination work end to end", where a combination is an ordering,
// a factorization, a traversal and a storage. It is deliberately not the place to check any single
// phase: test_numfactor checks the factor against an oracle, test_solve checks the solve, and this
// file checks that the phases compose. When a residual here goes wrong, the focused suites are
// where to look next.
//
// **Ordering is Natural throughout, and that is a choice, not a default.** A fill-reducing ordering
// would make the test depend on AMD's tie-breaking, so instead the matrices are built already in a
// good order, banded or grid-structured, and the ordering step is asked to do nothing. That keeps
// the numerical behavior the only variable.
//
// Mostly real. There is one complex section at the end, covering the cells that reach complex
// input at all.
//
// The tiers, which are about how hard the matrix is to pivot rather than how large it is:
//
//   Tier 0   No pivoting required. Diagonally dominant, so dynamic LDL must delay nothing and
//            choose no 2x2 pivot. Every factorization should handle it, and dynamic LDL should
//            reduce to static LDL.
//   Tier 1   Mild pivoting. A few columns cannot be pivoted where they stand, so a handful of
//            delays and 2x2 pivots occur, but the structure barely moves.
//
// Tier 2, heavy pivoting, comes later. Singular matrices are excluded on purpose: they have no
// residual to hit, and asserting something weaker about them would only look like coverage.
//
// Every assertion here is listed in docs/TESTING_SPECIFICATION.md. The two are kept in sync: a
// change to one is a change to the other.

#include "oblio/DirectSolver.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/NumFactorStatic.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SolveEngine.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Oblio;

namespace {

int pass = 0, fail = 0;
void ck(bool ok, const std::string& what) {
    if (ok) { ++pass; std::cout << "  PASS  " << what << "\n"; }
    else    { ++fail; std::cout << "  FAIL  " << what << "\n"; }
}

// Dense to CSC, full storage. **The diagonal is stored even when it is numerically zero**: a
// direct solver needs it structurally present, since symbolic factorization builds a column's
// index set from A's column structure and a column missing its own diagonal never enters it. See
// the input-validation entry in docs/TODO.md; nothing enforces this yet.
SparseMatrix<double> toSparse(const std::vector<std::vector<double>>& A) {
    const std::size_t n = A.size();
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i)
            if (A[i][j] != 0.0 || i == j) {
                rowIdx.push_back(static_cast<std::int32_t>(i));
                val.push_back(A[i][j]);
            }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<double>(n, std::move(colPtr), std::move(rowIdx), std::move(val));
}

// ---------------------------------------------------------------------------------------------
// The matrices.
//
// **Random values are derived from the engine by hand, not through a distribution.** `std::mt19937`
// has its sequence fixed by the standard, so it produces the same 32-bit stream everywhere; the
// distribution templates do not, their algorithms being implementation-defined, so the same seed
// can yield different doubles under libstdc++ and libc++. That does not matter where a test asserts
// only that the residual is small, which holds for the whole family. It matters here, because the
// tier 1 assertions pin exact delay and pivot counts, and those are properties of the particular
// matrix. Deriving the doubles ourselves makes the matrix a pure function of the seed.
// ---------------------------------------------------------------------------------------------

// A double in [0, 1), and one in [-1, 1), from the engine's raw output.
double u01(std::mt19937& rng) { return static_cast<double>(rng()) / 4294967296.0; }
double sym(std::mt19937& rng) { return 2.0 * u01(rng) - 1.0; }

// Tier 0. The 5-point Laplacian on a g by g grid, numbered row-major, which is already a good
// ordering (bandwidth g). Symmetric positive definite and diagonally dominant, so no pivoting is
// needed and Cholesky applies.
std::vector<std::vector<double>> gridLaplacian(std::size_t g) {
    const std::size_t n = g * g;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    const auto id = [g](std::size_t r, std::size_t c) { return r * g + c; };
    for (std::size_t r = 0; r < g; ++r)
        for (std::size_t c = 0; c < g; ++c) {
            const std::size_t i = id(r, c);
            A[i][i] = 4.0;
            if (r + 1 < g) { A[i][id(r + 1, c)] = A[id(r + 1, c)][i] = -1.0; }
            if (c + 1 < g) { A[i][id(r, c + 1)] = A[id(r, c + 1)][i] = -1.0; }
        }
    return A;
}

// Tier 1. A banded matrix of half-bandwidth w with random off-diagonals, in which a fraction of
// the diagonal entries are zero. Banded, so Natural is a sensible ordering; indefinite, and the
// zero diagonals are what force a handful of delays and 2x2 pivots.
//
// The zero diagonals do not all delay: most fill in from the Schur complement before they are
// reached, which is why quasi-definite systems factor without pivoting at all. The ones that do
// delay are those reached while still small relative to their column, which is exactly the case
// the machinery exists for.
std::vector<std::vector<double>> bandIndefinite(std::size_t n, std::size_t w,
                                                double zeroFraction, unsigned seed) {
    std::mt19937 rng(seed);
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        A[i][i] = (u01(rng) < zeroFraction) ? 0.0 : 4.0;
        for (std::size_t k = 1; k <= w && i + k < n; ++k) {
            const double v = sym(rng);
            A[i][i + k] = A[i + k][i] = v;
        }
    }
    return A;
}

// Tier 2. A saddle point system, [[H, B^T], [B, 0]], with **both** blocks carrying a zero diagonal:
// H is tridiagonal with nothing on its diagonal, and the constraint block is exactly zero. The
// honest use case for an indefinite solver, and the family that delays hardest, since a constraint
// column has no diagonal to pivot on and no update can give it one.
//
// A nonzero H diagonal makes this tier 0 again: with `hdiag` positive nothing delays at all, which
// is worth knowing before reaching for this family and finding it quiet.
std::vector<std::vector<double>> saddlePoint(std::size_t m, std::size_t k, double hdiag,
                                             unsigned seed) {
    std::mt19937 rng(seed);
    const std::size_t n = m + k;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    for (std::size_t i = 0; i < m; ++i) {
        A[i][i] = hdiag * (1.0 + u01(rng));
        if (i + 1 < m) { const double v = 0.5 * sym(rng); A[i][i + 1] = A[i + 1][i] = v; }
    }
    for (std::size_t r = 0; r < k; ++r)
        for (std::size_t c = 0; c < m; ++c)
            if ((r + c) % 3 == 0) { const double v = sym(rng); A[m + r][c] = v; A[c][m + r] = v; }
    return A;
}

// Tier 2, the extreme. A tridiagonal matrix with nothing on its diagonal: no 1x1 pivot can ever be
// accepted, so every pivot is a 2x2 and the factorization is exact.
//
// **Even order only.** At odd order this matrix is exactly singular (condition number around 1e16),
// so it has no residual to hit and is not a test of anything. That trap cost an hour once; n = 10
// is a fine case and n = 25 is not a case at all.
std::vector<std::vector<double>> zeroDiagonalTridiagonal(std::size_t n) {
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i + 1 < n; ++i)
        A[i][i + 1] = A[i + 1][i] = 1.0 + 0.1 * static_cast<double>(i);
    return A;
}

// A complex HERMITIAN band: A = A^H, conjugate off-diagonals and a real diagonal. What LDL^H
// factors, and a genuinely different matrix from the symmetric one below.
std::vector<std::vector<std::complex<double>>> bandComplexHermitian(std::size_t n, std::size_t w,
                                                                    double zeroFraction,
                                                                    unsigned seed) {
    using C = std::complex<double>;
    std::mt19937 rng(seed);
    std::vector<std::vector<C>> A(n, std::vector<C>(n, C(0.0, 0.0)));
    for (std::size_t i = 0; i < n; ++i) {
        A[i][i] = (u01(rng) < zeroFraction) ? C(0.0, 0.0) : C(4.0, 0.0);   // real diagonal
        for (std::size_t k = 1; k <= w && i + k < n; ++k) {
            const C v(sym(rng), sym(rng));
            A[i][i + k] = v;
            A[i + k][i] = std::conj(v);
        }
    }
    return A;
}

// A complex SYMMETRIC band, A = A^T with complex entries on the diagonal. Not Hermitian: over the
// complex field LDL^T and LDL^H factor genuinely different matrices, and this is the one LDL^T
// means. Same zero-diagonal trick as the real tier 1 family, to make it delay.
std::vector<std::vector<std::complex<double>>> bandComplexSymmetric(std::size_t n, std::size_t w,
                                                                    double zeroFraction,
                                                                    unsigned seed) {
    using C = std::complex<double>;
    std::mt19937 rng(seed);
    std::vector<std::vector<C>> A(n, std::vector<C>(n, C(0.0, 0.0)));
    for (std::size_t i = 0; i < n; ++i) {
        A[i][i] = (u01(rng) < zeroFraction) ? C(0.0, 0.0) : C(4.0, 1.0);
        for (std::size_t k = 1; k <= w && i + k < n; ++k) {
            const C v(sym(rng), sym(rng));
            A[i][i + k] = A[i + k][i] = v;
        }
    }
    return A;
}

SparseMatrix<std::complex<double>> toSparseComplex(
        const std::vector<std::vector<std::complex<double>>>& A) {
    using C = std::complex<double>;
    const std::size_t n = A.size();
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<C>            val;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i)
            if (A[i][j] != C(0.0, 0.0) || i == j) {
                rowIdx.push_back(static_cast<std::int32_t>(i));
                val.push_back(A[i][j]);
            }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<C>(n, std::move(colPtr), std::move(rowIdx), std::move(val));
}

// ---------------------------------------------------------------------------------------------
// One run of the pipeline.
// ---------------------------------------------------------------------------------------------

struct Outcome {
    bool        ran            = false;   // the factorization was produced
    bool        solved         = false;
    double      residual       = -1.0;    // ||Ax - b|| / ||b||
    std::int32_t delayed       = 0;       // columns delayed, summed over supernodes
    std::int32_t snodesDelaying = 0;
    std::int32_t pivots1x1     = 0;
    std::int32_t pivots2x2     = 0;
    std::size_t snodeSize      = 0;
};

// Factor is NumFactorStatic<double> or NumFactorDynamic<double>. The pivot statistics exist only
// on the dynamic factor, so they are gathered under `if constexpr` and left zero otherwise, which
// is the truth: a static factor delays nothing.
template<class Val, class Factor>
Outcome run(const SparseMatrix<Val>& A, OrderMethod om, Factorization fz, Traversal tr) {
    Outcome o;
    const std::size_t n = A.size();

    OrderEngine ord(om);
    Permutation p;
    if (!ord.compute(A, p)) return o;

    ElmForest f;
    ElmForestEngine fe;
    if (!fe.compute(A, p, f)) return o;

    SymFactor s;
    SymFactorEngine se;
    if (!se.compute(A, p, f, s)) return o;

    Factor nf;
    NumFactorEngine ne(fz, tr);
    if (!ne.compute(A, p, s, nf)) return o;
    o.ran = true;
    o.snodeSize = nf.snodeSize();

    if constexpr (std::is_same_v<Factor, NumFactorDynamic<Val>>) {
        for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.snodeSize()); ++kk) {
            const std::int32_t d = nf.numberOfDelayedColumns(kk);
            o.delayed += d;
            if (d > 0) ++o.snodesDelaying;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (nf.pivotType()[i] == 1) ++o.pivots1x1;
            if (nf.pivotType()[i] == 2) ++o.pivots2x2;
        }
    }

    Vector<Val> b(n), x(n);
    for (std::size_t i = 0; i < n; ++i)
        b[i] = Val(1.0 + 0.3 * static_cast<double>(i % 5));

    SolveEngine sol;
    if (!sol.compute(p, nf, b, x)) return o;
    o.solved = true;

    MultiplyEngine mul;
    Vector<Val> r(n);
    if (!mul.residual(A, x, b, r)) return o;
    o.residual = r.norm() / b.norm();
    return o;
}

// The worst residual over both implemented traversals, for a factorization that supports both.
// Returns -1 if any of them failed to run or solve.
template<class Val, class Factor>
double worstOverTraversals(const SparseMatrix<Val>& A, OrderMethod om, Factorization fz) {
    double worst = 0.0;
    for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
        const Outcome o = run<Val, Factor>(A, om, fz, tr);
        if (!o.solved) return -1.0;
        worst = std::max(worst, o.residual);
    }
    return worst;
}

std::string with(const std::string& label, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.2e", v);
    return label + " (" + buf + ")";
}

std::string counts(const std::string& label, const Outcome& o) {
    char buf[128];
    std::snprintf(buf, sizeof buf, " (%d delayed in %d snodes, %d 1x1, %d 2x2, resid %.2e)",
                  o.delayed, o.snodesDelaying, o.pivots1x1, o.pivots2x2, o.residual);
    return label + buf;
}

} // namespace

int main() {
    const double tol = 1e-12;

    using FS = NumFactorStatic<double>;
    using FD = NumFactorDynamic<double>;

    // =============================================================================================
    // Tier 0: no pivoting required. A 6x6 grid Laplacian, 36 columns, natural (row-major) order.
    // =============================================================================================
    {
        const SparseMatrix<double> A = toSparse(gridLaplacian(6));

        // The three statically pivoted factorizations, each in both storages, each worst-cased over
        // both traversals. For real input LDLT and LDLH are the same computation, and both are run
        // rather than one being assumed to stand in for the other.
        const double chS = worstOverTraversals<double, FS>(A, OrderMethod::Natural, Factorization::Cholesky);
        const double chD = worstOverTraversals<double, FD>(A, OrderMethod::Natural, Factorization::Cholesky);
        ck(chS >= 0 && chS < tol, with("tier 0 Cholesky   : residual, flat storage, both traversals", chS));
        ck(chD >= 0 && chD < tol, with("tier 0 Cholesky   : residual, per-supernode storage, both traversals", chD));

        const double ltS = worstOverTraversals<double, FS>(A, OrderMethod::Natural, Factorization::StaticLDLT);
        const double ltD = worstOverTraversals<double, FD>(A, OrderMethod::Natural, Factorization::StaticLDLT);
        ck(ltS >= 0 && ltS < tol, with("tier 0 StaticLDLT : residual, flat storage, both traversals", ltS));
        ck(ltD >= 0 && ltD < tol, with("tier 0 StaticLDLT : residual, per-supernode storage, both traversals", ltD));

        const double lhS = worstOverTraversals<double, FS>(A, OrderMethod::Natural, Factorization::StaticLDLH);
        const double lhD = worstOverTraversals<double, FD>(A, OrderMethod::Natural, Factorization::StaticLDLH);
        ck(lhS >= 0 && lhS < tol, with("tier 0 StaticLDLH : residual, flat storage, both traversals", lhS));
        ck(lhD >= 0 && lhD < tol, with("tier 0 StaticLDLH : residual, per-supernode storage, both traversals", lhD));

        // Dynamic LDL on an input that needs no pivoting. Two separate claims: the answer is right,
        // and the machinery correctly decided to do nothing. The second is the one that would catch
        // a pivot search that delays out of confusion rather than necessity.
        const Outcome dynL = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                     Traversal::LeftLooking);
        const Outcome dynR = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                     Traversal::RightLooking);
        ck(dynL.solved && dynR.solved && dynL.residual < tol && dynR.residual < tol,
           with("tier 0 DynamicLDLT: residual, both traversals", std::max(dynL.residual, dynR.residual)));
        ck(dynL.delayed == 0 && dynL.pivots2x2 == 0 && dynR.delayed == 0 && dynR.pivots2x2 == 0,
           counts("tier 0 DynamicLDLT: nothing delayed, no 2x2 chosen, both traversals", dynL));

        // And the same for DynamicLDLH, which over the reals is the same computation: the dynamic
        // path never reads the Hermitian flag, and the solve's conjugate is the identity for
        // double. Run rather than assumed.
        const Outcome dynHL = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLH,
                                      Traversal::LeftLooking);
        const Outcome dynHR = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLH,
                                      Traversal::RightLooking);
        ck(dynHL.solved && dynHR.solved && dynHL.residual < tol && dynHR.residual < tol,
           with("tier 0 DynamicLDLH: residual, both traversals", std::max(dynHL.residual, dynHR.residual)));
        ck(dynHL.delayed == 0 && dynHL.pivots2x2 == 0 && dynHR.delayed == 0 && dynHR.pivots2x2 == 0,
           counts("tier 0 DynamicLDLH: nothing delayed, no 2x2 chosen, both traversals", dynHL));

        // Combinations that must be refused rather than answered. These are as much a part of the
        // specification as the ones that work: a cell that starts returning a plausible wrong
        // answer instead of false is exactly the failure a port invites.
        ck(!run<double, FS>(A, OrderMethod::Natural, Factorization::DynamicLDLT, Traversal::LeftLooking).ran,
           "tier 0 refusal    : dynamic pivoting into flat storage");
        ck(!run<double, FS>(A, OrderMethod::Natural, Factorization::Cholesky, Traversal::Multifrontal).ran,
           "tier 0 refusal    : multifrontal traversal (not ported yet)");
    }

    // =============================================================================================
    // Tier 1: mild pivoting. Banded, indefinite, with a fraction of zero diagonals. Two matrices,
    // so a single lucky seed cannot carry the tier.
    //
    // The counts are pinned exactly rather than bounded. The matrices are reproducible (see the
    // note on the generators above), so the counts are facts about a specific matrix, not about a
    // family, and an exact assertion says more: it fails if the pivot search changes what it
    // chooses, not merely if it stops choosing anything.
    //
    // They are therefore change detectors, deliberately. What they detect is a change in *pivoting
    // behavior*, which is precisely the thing no other assertion in the suite can see. Legitimate
    // causes exist, a different default pivot threshold or a different amalgamation, and when one
    // of those lands these numbers are expected to move and should be re-recorded here and in the
    // specification together.
    // =============================================================================================
    {
        const SparseMatrix<double> A = toSparse(bandIndefinite(40, 3, 0.50, 7));
        const Outcome o = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol, with("tier 1 band n=40  : residual", o.residual));
        ck(o.ran && o.delayed == 5 && o.snodesDelaying == 5 && o.pivots2x2 == 4,
           counts("tier 1 band n=40  : 5 delayed in 5 snodes, 4 2x2", o));

        // The claim that the two transposes coincide over the reals, tested where it could
        // plausibly fail: an input that actually pivots. Bit-identical, not merely close, because
        // the two select the same arithmetic rather than equivalent arithmetic.
        const Outcome oH = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLH,
                                   Traversal::LeftLooking);
        ck(oH.solved && oH.delayed == o.delayed && oH.snodesDelaying == o.snodesDelaying
                     && oH.pivots2x2 == o.pivots2x2 && oH.residual == o.residual,
           "tier 1 band n=40  : DynamicLDLH is bit-identical to DynamicLDLT over the reals");

        // The two traversals are two different drivers over the same two kernels, and they grow a
        // front by opposite means: left-looking discards an empty front and rebuilds it, while
        // right-looking must carry forward the values already accumulated in it. Agreement on a
        // matrix that actually delays is what says the second of those is right.
        const Outcome oR = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                   Traversal::RightLooking);
        ck(oR.solved && oR.delayed == o.delayed && oR.snodesDelaying == o.snodesDelaying
                     && oR.pivots2x2 == o.pivots2x2 && oR.residual == o.residual,
           "tier 1 band n=40  : right-looking is bit-identical to left-looking");
    }
    {
        const SparseMatrix<double> A = toSparse(bandIndefinite(24, 3, 0.50, 7));
        const Outcome o = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol, with("tier 1 band n=24  : residual", o.residual));
        ck(o.ran && o.delayed == 3 && o.snodesDelaying == 3 && o.pivots2x2 == 3,
           counts("tier 1 band n=24  : 3 delayed in 3 snodes, 3 2x2", o));
    }

    // =============================================================================================
    // Tier 2: heavy pivoting. Where tier 1 delays a handful of columns, these delay dozens, and the
    // second family delays or pairs every column there is.
    //
    // **Counts are bounded here, not pinned**, which is the opposite of tier 1 and for a reason.
    // A pivot is accepted on `|d| >= threshold * max`, both sides coming out of BLAS updates, so
    // the more pivot decisions a matrix forces the likelier one of them sits near enough the
    // threshold to be decided differently under different rounding. Tier 1 makes a handful of
    // decisions and was verified identical across two BLAS implementations; tier 2 makes dozens and
    // that verification would not be worth relying on.
    //
    // The exception is the all-2x2 claim below, which is pinned. It is structural rather than
    // numerical: with nothing on the diagonal no 1x1 is available at all, so the count follows from
    // the matrix rather than from a comparison that rounding could tip.
    // =============================================================================================
    {
        const double tol2 = 1e-11;   // looser: these are far worse conditioned than tiers 0 and 1
        const SparseMatrix<double> A = toSparse(saddlePoint(30, 12, 0.0, 13));

        const Outcome L = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);
        const Outcome R = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                  Traversal::RightLooking);

        ck(L.solved && R.solved && L.residual < tol2 && R.residual < tol2,
           with("tier 2 saddle 30+12: residual, both traversals", std::max(L.residual, R.residual)));
        ck(L.delayed >= 40 && L.pivots2x2 >= 5 && R.delayed == L.delayed && R.pivots2x2 == L.pivots2x2,
           counts("tier 2 saddle 30+12: heavy delaying, traversals agree", L));
    }
    {
        // Every pivot a 2x2, so half the columns are marked as a pair's first. Even order only.
        for (std::size_t n : {std::size_t{12}, std::size_t{24}}) {
            const SparseMatrix<double> A = toSparse(zeroDiagonalTridiagonal(n));
            const Outcome o = run<double, FD>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                      Traversal::LeftLooking);

            const std::int32_t half = static_cast<std::int32_t>(n / 2);
            ck(o.solved && o.residual < tol && o.pivots1x1 == 0 && o.pivots2x2 == half,
               counts("tier 2 zero-diag tri : every pivot a 2x2, none 1x1", o));
        }
    }

    // =============================================================================================
    // Complex. Nine of the ten (factorization, scalar type) cells are supported; the missing one is
    // complex DynamicLDLH, which is refused rather than answered and is asserted to be.
    //
    // The dynamic kernels needed nothing to accept complex input. 0.9's complex factorDynamicLDL_
    // differs from its real one in six lines, all declaring the pivot magnitudes real rather than
    // scalar, and this port declared them double from the start; updateDynamicLDL_ is byte-identical
    // between 0.9's two engines. So this section is checking a claim about the port's shape as much
    // as it is checking arithmetic.
    //
    // A complex *symmetric* matrix, deliberately: A = A^T with complex diagonal entries, which is
    // what LDL^T means over the complex field and is not the same thing as Hermitian.
    // =============================================================================================
    {
        using C = std::complex<double>;
        using FDC = NumFactorDynamic<C>;
        using FSC = NumFactorStatic<C>;

        // **Two matrices, because these factorizations cannot share one.** A static factorization
        // cannot pivot, so a zero diagonal is *perturbed* rather than delayed and the residual is
        // then honestly poor: it factored a slightly different matrix and said so. So StaticLDLT
        // gets the dominant matrix and DynamicLDLT gets the one with zero diagonals, which is the
        // whole point of dynamic pivoting.
        //
        // Cholesky is absent on purpose: it needs Hermitian positive definite input, which is a
        // third matrix and a different property, and it is covered for complex in test_numfactor
        // and test_solve already. A complex symmetric matrix is not a valid Cholesky input at all.
        const SparseMatrix<C> dominant = toSparseComplex(bandComplexSymmetric(32, 3, 0.0, 7));
        const SparseMatrix<C> A        = toSparseComplex(bandComplexSymmetric(32, 3, 0.5, 7));

        const double slt = worstOverTraversals<C, FSC>(dominant, OrderMethod::Natural,
                                                       Factorization::StaticLDLT);
        ck(slt >= 0 && slt < tol, with("complex StaticLDLT: residual, both traversals", slt));

        const Outcome dL = run<C, FDC>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                       Traversal::LeftLooking);
        const Outcome dR = run<C, FDC>(A, OrderMethod::Natural, Factorization::DynamicLDLT,
                                       Traversal::RightLooking);
        ck(dL.solved && dR.solved && dL.residual < tol && dR.residual < tol,
           with("complex DynamicLDLT: residual, both traversals", std::max(dL.residual, dR.residual)));
        ck(dL.ran && dL.delayed >= 1 && dL.pivots2x2 >= 1
                  && dR.delayed == dL.delayed && dR.pivots2x2 == dL.pivots2x2,
           counts("complex DynamicLDLT: delaying happened, traversals agree", dL));

        // And the Hermitian one, on a genuinely Hermitian matrix: conjugate off-diagonals, real
        // diagonal, half of it zeroed so the pivot search has work to do. This is the extension
        // rather than a port, 0.9's complex LDL being symmetric only, so nothing here was checked
        // against a reference; the oracle is the residual and, in test_numfactor, reconstruction.
        const SparseMatrix<C> H = toSparseComplex(bandComplexHermitian(32, 3, 0.5, 7));

        const Outcome hL = run<C, FDC>(H, OrderMethod::Natural, Factorization::DynamicLDLH,
                                       Traversal::LeftLooking);
        const Outcome hR = run<C, FDC>(H, OrderMethod::Natural, Factorization::DynamicLDLH,
                                       Traversal::RightLooking);
        ck(hL.solved && hR.solved && hL.residual < tol && hR.residual < tol,
           with("complex DynamicLDLH: residual, both traversals", std::max(hL.residual, hR.residual)));
        ck(hL.ran && hL.delayed >= 1 && hL.pivots2x2 >= 1
                  && hR.delayed == hL.delayed && hR.pivots2x2 == hL.pivots2x2,
           counts("complex DynamicLDLH: delaying happened, traversals agree", hL));
    }

    // =============================================================================================
    // The same tier 0 matrix through DirectSolver, the facade over everything above.
    //
    // This is not redundant with the by-hand sweep. The facade owns both factors and picks between
    // them with dynamicPivoting(), so it can reach a combination the by-hand caller cannot reach by
    // accident and, more to the point, can fail to reach one it should. That is exactly the defect
    // examples/pipeline.cpp had: it fixed the storage at NumFactorStatic, so every dynamic cell
    // reported "not implemented" long after it was implemented, and nothing noticed because
    // examples are built but never run.
    // =============================================================================================
    {
        const SparseMatrix<double> A = toSparse(gridLaplacian(6));
        const std::size_t n = A.size();

        Vector<double> b(n), x(n);
        for (std::size_t i = 0; i < n; ++i)
            b[i] = 1.0 + 0.3 * static_cast<double>(i % 5);

        double worst = 0.0;
        int    reached = 0;
        for (Factorization fz : {Factorization::Cholesky, Factorization::StaticLDLT,
                                 Factorization::StaticLDLH, Factorization::DynamicLDLT,
                                 Factorization::DynamicLDLH})
            for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
                DirectSolver<double> solver(OrderMethod::Natural, fz, tr);
                if (!solver.analyze(A) || !solver.factor(A) || !solver.solve(b, x))
                    continue;
                ++reached;
                worst = std::max(worst, solver.relativeResidual(A, b, x));
            }

        ck(reached == 10, "DirectSolver      : all five factorizations reached, both traversals");
        ck(reached == 10 && worst < tol,
           with("DirectSolver      : worst residual over all ten", worst));
    }

    std::cout << "\nPipeline tests: " << pass << "/" << (pass + fail) << " passed\n";
    return fail == 0 ? 0 : 1;
}
