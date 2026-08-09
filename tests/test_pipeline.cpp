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
#include <cmath>
#include <complex>
#include <cstdint>
#include <optional>
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

// The same grid with sigma taken off the diagonal, A - sigma I. The pattern is untouched, since the
// diagonal is stored whatever its value, so one analysis serves every shift. Used to drive the
// inertia assertions, where the shift is what makes eigenvalues negative.
std::vector<std::vector<double>> shiftedGridLaplacian(std::size_t g, double sigma) {
    std::vector<std::vector<double>> A = gridLaplacian(g);
    for (std::size_t i = 0; i < A.size(); ++i) A[i][i] -= sigma;
    return A;
}

// How many eigenvalues of that matrix are negative, from the closed form rather than from anything
// the library computes: the eigenvalues of the g by g five-point Laplacian are
// 4 - 2cos(pi i / (g+1)) - 2cos(pi j / (g+1)). This is the oracle the inertia is checked against,
// and it shares no code with the factorization.
std::size_t negativeEigenvalues(std::size_t g, double sigma) {
    std::size_t negative = 0;
    for (std::size_t i = 1; i <= g; ++i)
        for (std::size_t j = 1; j <= g; ++j) {
            const double lambda = 4.0
                - 2.0 * std::cos(M_PI * static_cast<double>(i) / static_cast<double>(g + 1))
                - 2.0 * std::cos(M_PI * static_cast<double>(j) / static_cast<double>(g + 1))
                - sigma;
            if (lambda < 0.0) ++negative;
        }
    return negative;
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
    std::size_t rank           = 0;       // dynamic only: full order less the zero pivots taken
};

// Factor is NumFactorStatic<double> or NumFactorDynamic<double>. The pivot statistics exist only
// on the dynamic factor, so they are gathered under `if constexpr` and left zero otherwise, which
// is the truth: a static factor delays nothing.
template<class Val, class Factor>
Outcome run(const SparseMatrix<Val>& A, Ordering om, Factorization fz, Traversal tr) {
    Outcome o;
    const std::size_t n = A.size();

    OrderEngine ord(om);
    Permutation P;
    if (!ord.compute(A, P)) return o;

    ElmForest f;
    ElmForestEngine fe;
    if (!fe.compute(A, P, f)) return o;

    SymFactor s;
    SymFactorEngine se;
    if (!se.compute(A, P, f, s)) return o;

    Factor nf;
    NumFactorEngine ne(fz, tr);
    if (!ne.compute(A, P, s, nf)) return o;
    o.ran = true;
    o.snodeSize = nf.snodeSize();

    if constexpr (std::is_same_v<Factor, NumFactorDynamic<Val>>) {
        o.rank = std::as_const(nf).rank();
        for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.snodeSize()); ++kk) {
            const std::int32_t d = static_cast<std::int32_t>(nf.delaySize(kk));
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
    if (!sol.compute(P, nf, b, x)) return o;
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
double worstOverTraversals(const SparseMatrix<Val>& A, Ordering om, Factorization fz) {
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
        const double chS = worstOverTraversals<double, FS>(A, Ordering::Natural, Factorization::Cholesky);
        const double chD = worstOverTraversals<double, FD>(A, Ordering::Natural, Factorization::Cholesky);
        ck(chS >= 0 && chS < tol, with("tier 0 Cholesky   : residual, flat storage, both traversals", chS));
        ck(chD >= 0 && chD < tol, with("tier 0 Cholesky   : residual, per-supernode storage, both traversals", chD));

        const double ltS = worstOverTraversals<double, FS>(A, Ordering::Natural, Factorization::StaticLDLT);
        const double ltD = worstOverTraversals<double, FD>(A, Ordering::Natural, Factorization::StaticLDLT);
        ck(ltS >= 0 && ltS < tol, with("tier 0 StaticLDLT : residual, flat storage, both traversals", ltS));
        ck(ltD >= 0 && ltD < tol, with("tier 0 StaticLDLT : residual, per-supernode storage, both traversals", ltD));

        const double lhS = worstOverTraversals<double, FS>(A, Ordering::Natural, Factorization::StaticLDLH);
        const double lhD = worstOverTraversals<double, FD>(A, Ordering::Natural, Factorization::StaticLDLH);
        ck(lhS >= 0 && lhS < tol, with("tier 0 StaticLDLH : residual, flat storage, both traversals", lhS));
        ck(lhD >= 0 && lhD < tol, with("tier 0 StaticLDLH : residual, per-supernode storage, both traversals", lhD));

        // Dynamic LDL on an input that needs no pivoting. Two separate claims: the answer is right,
        // and the machinery correctly decided to do nothing. The second is the one that would catch
        // a pivot search that delays out of confusion rather than necessity.
        const Outcome dynL = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                     Traversal::LeftLooking);
        const Outcome dynR = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                     Traversal::RightLooking);
        ck(dynL.solved && dynR.solved && dynL.residual < tol && dynR.residual < tol,
           with("tier 0 DynamicLDLT: residual, both traversals", std::max(dynL.residual, dynR.residual)));
        ck(dynL.delayed == 0 && dynL.pivots2x2 == 0 && dynR.delayed == 0 && dynR.pivots2x2 == 0,
           counts("tier 0 DynamicLDLT: nothing delayed, no 2x2 chosen, both traversals", dynL));

        // And the same for DynamicLDLH, which over the reals is the same computation: the dynamic
        // path never reads the Hermitian flag, and the solve's conjugate is the identity for
        // double. Run rather than assumed.
        const Outcome dynHL = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLH,
                                      Traversal::LeftLooking);
        const Outcome dynHR = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLH,
                                      Traversal::RightLooking);
        ck(dynHL.solved && dynHR.solved && dynHL.residual < tol && dynHR.residual < tol,
           with("tier 0 DynamicLDLH: residual, both traversals", std::max(dynHL.residual, dynHR.residual)));
        ck(dynHL.delayed == 0 && dynHL.pivots2x2 == 0 && dynHR.delayed == 0 && dynHR.pivots2x2 == 0,
           counts("tier 0 DynamicLDLH: nothing delayed, no 2x2 chosen, both traversals", dynHL));

        // Combinations that must be refused rather than answered. These are as much a part of the
        // specification as the ones that work: a cell that starts returning a plausible wrong
        // answer instead of false is exactly the failure a port invites.
        ck(!run<double, FS>(A, Ordering::Natural, Factorization::DynamicLDLT, Traversal::LeftLooking).ran,
           "tier 0 refusal    : dynamic pivoting into flat storage");

        // Dynamic multifrontal on an input that needs no pivoting: the basic path, before delayed
        // columns enter. Same answer as the other two dynamic traversals.
        const Outcome dynM = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                     Traversal::Multifrontal);
        ck(dynM.solved && dynM.residual < tol && dynM.delayed == 0,
           with("tier 0 DynamicLDLT: multifrontal residual, nothing delayed", dynM.residual));

        // Multifrontal, the third traversal, for every static factorization now: Cholesky and both
        // static LDLs, the same factor as left- and right-looking reached by a postorder pass that
        // carries each supernode's contribution block up a stack to its parent. Checked end to end by
        // residual (real here; complex multifrontal is checked directly against the dense oracle in
        // test_numfactor).
        const Outcome mfC = run<double, FS>(A, Ordering::Natural, Factorization::Cholesky,
                                    Traversal::Multifrontal);
        const Outcome mfT = run<double, FS>(A, Ordering::Natural, Factorization::StaticLDLT,
                                    Traversal::Multifrontal);
        const Outcome mfH = run<double, FS>(A, Ordering::Natural, Factorization::StaticLDLH,
                                    Traversal::Multifrontal);
        ck(mfC.solved && mfC.residual < tol, with("tier 0 Cholesky  : multifrontal residual, real", mfC.residual));
        ck(mfT.solved && mfT.residual < tol, with("tier 0 StaticLDLT: multifrontal residual, real", mfT.residual));
        ck(mfH.solved && mfH.residual < tol, with("tier 0 StaticLDLH: multifrontal residual, real", mfH.residual));

        // The same three into *dynamic* storage. A static factorization needs no pivoting and so
        // delays nothing, but the storage is agnostic about that and the engine's second compute
        // overload accepts it: the cell exists to keep the two storages interchangeable for the
        // factorizations both can hold. Left- and right-looking were already dispatched here;
        // multifrontal returned false as "not ported yet" until 2026-07-26.
        const Outcome mfCD = run<double, FD>(A, Ordering::Natural, Factorization::Cholesky,
                                     Traversal::Multifrontal);
        const Outcome mfTD = run<double, FD>(A, Ordering::Natural, Factorization::StaticLDLT,
                                     Traversal::Multifrontal);
        const Outcome mfHD = run<double, FD>(A, Ordering::Natural, Factorization::StaticLDLH,
                                     Traversal::Multifrontal);
        ck(mfCD.solved && mfCD.residual < tol && mfCD.delayed == 0,
           with("tier 0 Cholesky  : multifrontal into dynamic storage", mfCD.residual));
        ck(mfTD.solved && mfTD.residual < tol && mfTD.delayed == 0,
           with("tier 0 StaticLDLT: multifrontal into dynamic storage", mfTD.residual));
        ck(mfHD.solved && mfHD.residual < tol && mfHD.delayed == 0,
           with("tier 0 StaticLDLH: multifrontal into dynamic storage", mfHD.residual));
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
        const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol, with("tier 1 band n=40  : residual", o.residual));
        ck(o.ran && o.delayed == 5 && o.snodesDelaying == 5 && o.pivots2x2 == 4,
           counts("tier 1 band n=40  : 5 delayed in 5 snodes, 4 2x2", o));

        // The claim that the two transposes coincide over the reals, tested where it could
        // plausibly fail: an input that actually pivots. Bit-identical, not merely close, because
        // the two select the same arithmetic rather than equivalent arithmetic.
        const Outcome oH = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLH,
                                   Traversal::LeftLooking);
        ck(oH.solved && oH.delayed == o.delayed && oH.snodesDelaying == o.snodesDelaying
                     && oH.pivots2x2 == o.pivots2x2 && oH.residual == o.residual,
           "tier 1 band n=40  : DynamicLDLH is bit-identical to DynamicLDLT over the reals");

        // The two traversals are two different drivers over the same two kernels, and they grow a
        // front by opposite means: left-looking discards an empty front and rebuilds it, while
        // right-looking must carry forward the values already accumulated in it. Agreement on a
        // matrix that actually delays is what says the second of those is right.
        const Outcome oR = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                   Traversal::RightLooking);
        ck(oR.solved && oR.delayed == o.delayed && oR.snodesDelaying == o.snodesDelaying
                     && oR.pivots2x2 == o.pivots2x2 && oR.residual == o.residual,
           "tier 1 band n=40  : right-looking is bit-identical to left-looking");

        // Multifrontal is the third driver over the same two kernels, reaching the delayed columns by
        // yet another route: it carries each supernode's contribution block up the stack rather than
        // pulling or pushing per ancestor. This is the assertion that exercises delayed columns
        // meeting the stack. The pivoting decisions are the same (same delays, same 2x2 pivots); the
        // residual is checked to tolerance rather than bit-for-bit, since the assembly sums the
        // front in a different order than the pull.
        const Outcome oM = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                   Traversal::Multifrontal);
        ck(oM.solved && oM.residual < tol && oM.delayed == o.delayed
                     && oM.snodesDelaying == o.snodesDelaying && oM.pivots2x2 == o.pivots2x2,
           with("tier 1 band n=40  : multifrontal, same delays and 2x2, residual", oM.residual));
    }
    {
        const SparseMatrix<double> A = toSparse(bandIndefinite(24, 3, 0.50, 7));
        const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol, with("tier 1 band n=24  : residual", o.residual));
        ck(o.ran && o.delayed == 3 && o.snodesDelaying == 3 && o.pivots2x2 == 4,
           counts("tier 1 band n=24  : 3 delayed in 3 snodes, 4 2x2", o));
    }
    {
        // Roots pivot by bounded Bunch-Kaufman, and this is the family that separates it from what
        // 0.9 did there. The whole matrix is one root front. Column 0 has a zero diagonal, and its
        // largest off-diagonal is mutually maximal with column 1's, so the old acceptance test,
        // `max1 == max2` on the magnitudes alone, took a 2x2 whose partner diagonal is 1e6. The
        // determinant is then tiny beside the entries that divide by it, and L picks up an entry of
        // 1e6 against a matrix whose largest entry is 1e6 and whose smallest is 0.
        //
        // The chase cannot make that choice: it reaches a 2x2 only after both diagonals have failed
        // their own tests, and 1e6 passes, so it is taken as a 1x1 instead. max|L| is 1e-6 rather
        // than 1e6, and the residual moves from 1.45e-11 to 9.7e-17. The residual alone separates
        // them at this tolerance, which is why no factor-norm plumbing is needed to pin it.
        const SparseMatrix<double> A = toSparse({{0.0, 1.0, 1.0},
                                                 {1.0, 1.0e6, 0.0},
                                                 {1.0, 0.0, 0.0}});
        const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol,
           with("root pivoting     : zero diagonal beside a huge one, residual", o.residual));
        ck(o.ran && o.delayed == 0 && o.pivots1x1 == 1 && o.pivots2x2 == 1,
           counts("root pivoting     : no delay at a root, 1 1x1 and 1 2x2", o));
    }
    {
        // rank counts the zero pivots taken, and the last column of a root front is where it is
        // easiest to lose. That column's scan sees two empty ranges, so if the scan reported a
        // negative sentinel rather than zero, the "nothing to eliminate" test would not fire, the
        // column would be accepted as a 1x1 by a comparison against a negative bound, and a zero
        // diagonal would go uncounted. Two isolated columns, one of them zero, so both supernodes
        // are singleton roots and both take that path.
        std::vector<std::size_t>  colPtr{0, 1, 2};
        std::vector<std::int32_t> rowIdx{0, 1};
        std::vector<double>       val{0.0, 1.0};
        const SparseMatrix<double> A(2, std::move(colPtr), std::move(rowIdx), std::move(val));
        const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.ran && o.rank == 1 && o.delayed == 0 && o.pivots1x1 == 2,
           counts("root rank         : one zero diagonal drops rank to 1", o));
    }
    {
        // The non-root counterpart, and the reason the symmetric-maximum clause was deleted from
        // acceptPivot2x2 rather than repaired. Found by sweeping small banded matrices with zero
        // diagonals and a wide magnitude spread, which is what lets one diagonal be huge while its
        // neighbour is zero. The clause fires exactly once here, at a non-root front, on a block
        // the growth-bound test rejects: it took the 2x2, L picked up an entry of 1.98e5, and the
        // residual came out at 1.17e-11. With the clause gone the column is delayed instead,
        // max|L| is 2.76 and the residual is 4.8e-16.
        //
        // As with the root case the residual alone separates the two at this tolerance, so the
        // assertion needs no factor-norm plumbing. The delay count is pinned as well, since that is
        // where the refused block goes.
        const double h = 547866.5096, a = -0.3620555447, b = 0.9559790324, c = 0.07699174341;
        const double d = -0.8265131298, e = 0.05552958231, f = 0.7376029119;
        const double g = -0.8681273111, p1 = -0.001664438805, p2 = -0.00933406837;
        const double p3 = 0.002869110033;
        const SparseMatrix<double> A = toSparse({{0.0, 1.0, a,   0.0, 0.0, 0.0},
                                                 {1.0, h,   b,   c,   0.0, 0.0},
                                                 {a,   b,   0.0, d,   1.0, 0.0},
                                                 {0.0, c,   d,   p1,  e,   f  },
                                                 {0.0, 0.0, 1.0, e,   p2,  g  },
                                                 {0.0, 0.0, 0.0, f,   g,   p3 }});
        const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);

        ck(o.solved && o.residual < tol,
           with("non-root pivoting : zero diagonal beside a huge one, residual", o.residual));
        ck(o.ran && o.delayed == 2 && o.pivots1x1 == 2 && o.pivots2x2 == 2,
           counts("non-root pivoting : 2 delayed, 2 1x1 and 2 2x2", o));
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

        const Outcome L = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::LeftLooking);
        const Outcome R = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::RightLooking);

        ck(L.solved && R.solved && L.residual < tol2 && R.residual < tol2,
           with("tier 2 saddle 30+12: residual, both traversals", std::max(L.residual, R.residual)));
        ck(L.delayed >= 40 && L.pivots2x2 >= 5 && R.delayed == L.delayed && R.pivots2x2 == L.pivots2x2,
           counts("tier 2 saddle 30+12: heavy delaying, traversals agree", L));

        // Multifrontal under heavy delaying. Counts are bounded, not matched to left-looking: with
        // dozens of threshold decisions and a different summation order in the assembly, a decision
        // can tip, which is a different valid factorization, not an error. The residual is the oracle.
        const Outcome M = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                  Traversal::Multifrontal);
        ck(M.solved && M.residual < tol2 && M.delayed >= 40 && M.pivots2x2 >= 5,
           with("tier 2 saddle 30+12: multifrontal, heavy delaying, residual", M.residual));
    }
    {
        // Every pivot a 2x2, so half the columns are marked as a pair's first. Even order only.
        for (std::size_t n : {std::size_t{12}, std::size_t{24}}) {
            const SparseMatrix<double> A = toSparse(zeroDiagonalTridiagonal(n));
            const Outcome o = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                      Traversal::LeftLooking);

            const std::int32_t half = static_cast<std::int32_t>(n / 2);
            ck(o.solved && o.residual < tol && o.pivots1x1 == 0 && o.pivots2x2 == half,
               counts("tier 2 zero-diag tri : every pivot a 2x2, none 1x1", o));
        }
    }

    // =============================================================================================
    // Complex. Nine of the ten (factorization, scalar type) cells are supported; the missing one is
    // complex DynamicLDLH, the one cell with no 0.9 reference behind it.
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

        const double slt = worstOverTraversals<C, FSC>(dominant, Ordering::Natural,
                                                       Factorization::StaticLDLT);
        ck(slt >= 0 && slt < tol, with("complex StaticLDLT: residual, both traversals", slt));

        const Outcome dL = run<C, FDC>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                       Traversal::LeftLooking);
        const Outcome dR = run<C, FDC>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                       Traversal::RightLooking);
        ck(dL.solved && dR.solved && dL.residual < tol && dR.residual < tol,
           with("complex DynamicLDLT: residual, both traversals", std::max(dL.residual, dR.residual)));
        ck(dL.ran && dL.delayed >= 1 && dL.pivots2x2 >= 1
                  && dR.delayed == dL.delayed && dR.pivots2x2 == dL.pivots2x2,
           counts("complex DynamicLDLT: delaying happened, traversals agree", dL));

        const Outcome dM = run<C, FDC>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                       Traversal::Multifrontal);
        ck(dM.solved && dM.residual < tol && dM.delayed >= 1 && dM.pivots2x2 >= 1,
           with("complex DynamicLDLT: multifrontal, delaying happened, residual", dM.residual));

        // And the Hermitian one, on a genuinely Hermitian matrix: conjugate off-diagonals, real
        // diagonal, half of it zeroed so the pivot search has work to do. This is the extension
        // rather than a port, 0.9's complex LDL being symmetric only, so nothing here was checked
        // against a reference; the oracle is the residual and, in test_numfactor, reconstruction.
        const SparseMatrix<C> H = toSparseComplex(bandComplexHermitian(32, 3, 0.5, 7));

        const Outcome hL = run<C, FDC>(H, Ordering::Natural, Factorization::DynamicLDLH,
                                       Traversal::LeftLooking);
        const Outcome hR = run<C, FDC>(H, Ordering::Natural, Factorization::DynamicLDLH,
                                       Traversal::RightLooking);
        ck(hL.solved && hR.solved && hL.residual < tol && hR.residual < tol,
           with("complex DynamicLDLH: residual, both traversals", std::max(hL.residual, hR.residual)));
        ck(hL.ran && hL.delayed >= 1 && hL.pivots2x2 >= 1
                  && hR.delayed == hL.delayed && hR.pivots2x2 == hL.pivots2x2,
           counts("complex DynamicLDLH: delaying happened, traversals agree", hL));

        // The last cell of the whole matrix: complex Hermitian dynamic through the stack. This is the
        // trickiest path, the conjugating pivot kernel reached by the multifrontal driver, and the
        // one with no reference behind it. The residual is the oracle.
        const Outcome hM = run<C, FDC>(H, Ordering::Natural, Factorization::DynamicLDLH,
                                       Traversal::Multifrontal);
        ck(hM.solved && hM.residual < tol && hM.delayed >= 1 && hM.pivots2x2 >= 1,
           with("complex DynamicLDLH: multifrontal, delaying happened, residual", hM.residual));
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
            for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking,
                                 Traversal::Multifrontal}) {
                DirectSolver<double> solver(Ordering::Natural, fz, tr);
                if (!solver.analyze(A) || !solver.factor(A) || !solver.solve(b, x))
                    continue;
                ++reached;
                worst = std::max(worst, solver.relativeResidual(A, b, x));
            }

        ck(reached == 15, "DirectSolver      : all five factorizations reached, all three traversals");
        ck(reached == 15 && worst < tol,
           with("DirectSolver      : worst residual over all fifteen", worst));

        // Every ordering method, on the same matrix and through the same facade. An ordering can
        // be a valid permutation and still be one the rest of the pipeline cannot use, which is
        // the failure test_order's validity check cannot see, so what is asserted here is that a
        // factorization comes out the other end at machine precision. Fill is not asserted, for
        // any method: no ordering in this suite is checked for quality, and MMD1 is not held to
        // the vendored routines' output, being a different ordering rather than a copy of one.
        double worstOrder   = 0.0;
        int    reachedOrder = 0;
        for (Ordering om : {Ordering::Natural, Ordering::MMD, Ordering::MMD1,
                               Ordering::MMD2, Ordering::MMD3,
                               Ordering::AMD, Ordering::AMD1,
                               Ordering::AMD2, Ordering::AMD3, Ordering::AMD1B,
                               Ordering::AMD2B}) {
            DirectSolver<double> solver(om, Factorization::Cholesky, Traversal::LeftLooking);
            if (!solver.analyze(A) || !solver.factor(A) || !solver.solve(b, x))
                continue;
            ++reachedOrder;
            worstOrder = std::max(worstOrder, solver.relativeResidual(A, b, x));
        }

        // Two of the ten are the vendored routines, which are optional: without private/ they
        // refuse and the sweep skips them, so what is expected is every ordering the build has.
        // The count is written out rather than taken from the list's size, deliberately: adding an
        // enumerator should make this fail until someone has decided the new ordering belongs in
        // the sweep, which is exactly what happened when MMD3 was added.
#ifdef OBLIO_VENDORED_ORDERINGS
        const int expectedOrderings = 11;
#else
        const int expectedOrderings = 9;
#endif
        ck(reachedOrder == expectedOrderings, "Ordering       : every ordering built was reached");
        ck(reachedOrder == expectedOrderings && worstOrder < tol,
           with("Ordering       : worst residual over all orderings", worstOrder));

        // The multifrontal child ordering is computed during analyze, so the traversal has to be
        // known by then. Switching between left- and right-looking reads the same forest and must
        // not throw the analysis away; switching into or out of multifrontal must, because the
        // forest itself differs, its children reordered and its supernodes relabeled.
        DirectSolver<double> ds(Ordering::Natural, Factorization::Cholesky,
                                Traversal::LeftLooking);
        ck(ds.analyze(A) && ds.analyzed(), "DirectSolver      : analyze succeeds");

        ds.setTraversal(Traversal::RightLooking);
        ck(ds.analyzed(), "DirectSolver      : analysis survives left- to right-looking");

        ds.setTraversal(Traversal::Multifrontal);
        ck(!ds.analyzed(), "DirectSolver      : analysis invalidated switching to multifrontal");

        ck(ds.analyze(A) && ds.factor(A) && ds.solve(b, x)
               && ds.relativeResidual(A, b, x) < tol,
           "DirectSolver      : re-analyzed multifrontal solves");

        ds.setTraversal(Traversal::LeftLooking);
        ck(!ds.analyzed(), "DirectSolver      : analysis invalidated switching out of multifrontal");

        // Supernodes and amalgamation are constructor-only, so there is nothing to invalidate.
        // Fundamental and no amalgamation are the defaults; nodal and an amalgamating solver must
        // both reach the same answer, since neither changes what is being computed, only the block
        // structure it is computed in.
        DirectSolver<double> dsDefault(Ordering::Natural, Factorization::Cholesky,
                                       Traversal::LeftLooking);
        ck(dsDefault.supernodes() == Supernodes::Fundamental && !dsDefault.amalgamation().has_value(),
           "DirectSolver      : fundamental supernodes and no amalgamation by default");

        DirectSolver<double> dsNodal(Ordering::Natural, Factorization::Cholesky,
                                     Traversal::LeftLooking, Supernodes::Nodal);
        ck(dsNodal.supernodes() == Supernodes::Nodal
               && dsNodal.analyze(A) && dsNodal.factor(A) && dsNodal.solve(b, x)
               && dsNodal.relativeResidual(A, b, x) < tol,
           "DirectSolver      : nodal supernodes reachable and solve correct");

        DirectSolver<double> dsAmal(Ordering::Natural, Factorization::Cholesky,
                                    Traversal::Multifrontal, Supernodes::Fundamental, 8);
        ck(dsAmal.amalgamation().has_value() && *dsAmal.amalgamation() == 8
               && dsAmal.analyze(A) && dsAmal.factor(A) && dsAmal.solve(b, x)
               && dsAmal.relativeResidual(A, b, x) < tol,
           "DirectSolver      : amalgamation reachable, with multifrontal, and solve correct");

        // Both are settable afterwards too, and both invalidate the analysis, since the forest is
        // what they change and the analysis is the forest.
        ck(dsAmal.analyzed(), "DirectSolver      : analyzed after the amalgamating run");
        dsAmal.setSupernodes(Supernodes::Nodal);
        ck(!dsAmal.analyzed() && dsAmal.supernodes() == Supernodes::Nodal,
           "DirectSolver      : setSupernodes takes and invalidates the analysis");

        dsAmal.analyze(A);
        dsAmal.setAmalgamation(std::nullopt);
        ck(!dsAmal.analyzed() && !dsAmal.amalgamation().has_value(),
           "DirectSolver      : setAmalgamation takes and invalidates the analysis");

        ck(dsAmal.analyze(A) && dsAmal.factor(A) && dsAmal.solve(b, x)
               && dsAmal.relativeResidual(A, b, x) < tol,
           "DirectSolver      : solves again after both forest settings changed");
    }

    // =============================================================================================
    // What the facade reports about pivoting: the delay and pivot counts, and the inertia.
    //
    // These four accessors are the facade's only window onto the numeric factor, which it does not
    // expose, so they are checked two ways. The counts are checked against the by-hand sweep above,
    // which reaches the factor directly and whose numbers tier 1 already pins; agreement therefore
    // means the forwarders read what the factor holds. The inertia is checked against the closed
    // form for the grid's eigenvalues, which shares no code with anything here.
    // =============================================================================================
    {
        // The tier 1 matrix, whose counts are pinned above at 5 delayed and 4 two-by-twos.
        const SparseMatrix<double> A = toSparse(bandIndefinite(40, 3, 0.50, 7));
        const std::size_t n = A.size();

        const Outcome byHand = run<double, FD>(A, Ordering::Natural, Factorization::DynamicLDLT,
                                               Traversal::LeftLooking);

        DirectSolver<double> ds(Ordering::Natural, Factorization::DynamicLDLT,
                                Traversal::LeftLooking);
        ck(ds.analyze(A) && ds.factor(A), "facade counts     : dynamic factorization ran");

        ck(ds.numDelayedColumns() == static_cast<std::size_t>(byHand.delayed)
               && ds.numPivots1x1() == static_cast<std::size_t>(byHand.pivots1x1)
               && ds.numPivots2x2() == static_cast<std::size_t>(byHand.pivots2x2),
           "facade counts     : delayed, 1x1 and 2x2 agree with the by-hand sweep");

        // The columns are partitioned by the pivot choice, which is the invariant that would catch
        // a 2x2 counted in columns rather than in blocks.
        ck(ds.numPivots1x1() + 2 * ds.numPivots2x2() == n,
           "facade counts     : 1x1 + 2 * 2x2 covers every column exactly once");

        // A static factorization makes no pivot choices, so all three are zero and that is the
        // accurate report rather than a placeholder.
        DirectSolver<double> dsStatic(Ordering::Natural, Factorization::StaticLDLT,
                                      Traversal::LeftLooking);
        ck(dsStatic.analyze(A) && dsStatic.factor(A)
               && dsStatic.numDelayedColumns() == 0
               && dsStatic.numPivots1x1() == 0 && dsStatic.numPivots2x2() == 0,
           "facade counts     : a statically pivoted factor reports no pivot choices");

        // The three sizes, on all four classes that answer them and on the facade that forwards.
        //
        // These check against each other rather than against a recomputation, which is what makes
        // them worth having: the four classes describe one factorization from four points in the
        // pipeline, so they are constrained to agree except where the algorithm says they must
        // not. A wrong formula in one place breaks the agreement; a wrong formula copied into all
        // four would not, which is why the nesting and the delay relation are checked too.
        {
            const SparseMatrix<double> S = toSparse(bandIndefinite(40, 3, 0.50, 7));

            OrderEngine     ordEng(Ordering::MMD2);
            ElmForestEngine efEng;
            SymFactorEngine sfEng;
            Permutation     P;
            ElmForest       ef;
            SymFactor       sf;
            const bool analyzed = ordEng.compute(S, P) && efEng.compute(S, P, ef)
                               && sfEng.compute(S, P, ef, sf);
            ck(analyzed, "sizes             : analysis ran");

            // Nesting, on every class: an index set is one entry per row, nnz adds the values
            // against those rows, and the allocation adds the front's upper triangle on top.
            auto nests = [](std::size_t numNodeIdx, std::size_t nnz, std::size_t numVal) {
                return numNodeIdx <= nnz && nnz <= numVal;
            };

            // Recomputed here from the per-supernode sizes each class publishes, which is the only
            // check that pins the formulae rather than the relations between them. An inequality
            // cannot: a dynamic factor's own frontSize differs from the symbolic one once columns
            // have moved, so even a formula that ignored delaySize entirely still comes out larger
            // than predicted. `rows` is what the supernode's index set holds, which is where the
            // delays show.
            auto recompute = [](std::size_t snodes, auto front, auto rows) {
                struct { std::size_t numNodeIdx = 0, nnz = 0, numVal = 0; } sum;
                for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodes); ++jj) {
                    const std::size_t f = front(jj), r = rows(jj);
                    sum.numNodeIdx += r;
                    sum.nnz        += f * (f + 1) / 2 + f * (r - f);
                    sum.numVal     += r * f;
                }
                return sum;
            };
            auto matches = [](const auto& sum, const auto& object) {
                return sum.numNodeIdx == object.numNodeIdx() && sum.nnz == object.nnz()
                    && sum.numVal == object.numVal();
            };
            ck(nests(ef.numNodeIdx(), ef.nnz(), ef.numVal())
                   && nests(sf.numNodeIdx(), sf.nnz(), sf.numVal()),
               "sizes             : numNodeIdx <= nnz <= numVal, forest and symbolic");

            // The forest predicts and the symbolic factor materializes, from the same sizes, so
            // all three must match exactly. This is the assertion that catches a formula that
            // drifted in one of the two.
            ck(ef.numNodeIdx() == sf.numNodeIdx() && ef.nnz() == sf.nnz()
                   && ef.numVal() == sf.numVal(),
               "sizes             : forest and symbolic agree on all three");

            NumFactorStatic<double>  nfs;
            NumFactorDynamic<double> nfd;
            NumFactorEngine          neS(Factorization::StaticLDLT,  Traversal::LeftLooking);
            NumFactorEngine          neD(Factorization::DynamicLDLT, Traversal::LeftLooking);
            const bool factored = neS.compute(S, P, sf, nfs) && neD.compute(S, P, sf, nfd);
            ck(factored, "sizes             : both storages factored");

            ck(nests(nfs.numNodeIdx(), nfs.nnz(), nfs.numVal())
                   && nests(nfd.numNodeIdx(), nfd.nnz(), nfd.numVal()),
               "sizes             : numNodeIdx <= nnz <= numVal, both numeric factors");

            // A static factorization moves nothing, so its sizes are the symbolic ones exactly.
            ck(nfs.numNodeIdx() == sf.numNodeIdx() && nfs.nnz() == sf.nnz()
                   && nfs.numVal() == sf.numVal(),
               "sizes             : the static factor matches what symbolic predicted");

            // A dynamic one delays columns into their parents, which can only grow a front, so
            // its sizes are at least the predicted ones and strictly larger where a delay landed.
            // This matrix is tier 1's, which the suite already pins at 5 delayed columns.
            std::size_t delayed = 0;
            for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(nfd.snodeSize()); ++jj)
                delayed += nfd.delaySize(jj);
            ck(delayed > 0, "sizes             : the dynamic factor delayed at least one column");

            // Strictly larger, not merely no smaller. A delayed column adds a row to its parent's
            // front, which adds an index, values against it, and allocation for both, so with
            // delays on the books every one of the three must have grown. The weaker `>=` is
            // satisfied by a formula that ignores delaySize entirely, which is exactly the bug
            // worth catching here.
            ck(nfd.numNodeIdx() > sf.numNodeIdx() && nfd.nnz() > sf.nnz()
                   && nfd.numVal() > sf.numVal(),
               "sizes             : the dynamic factor is strictly larger than predicted");
            ck(nfd.numNodeIdx() == sf.numNodeIdx() + delayed,
               "sizes             : its index sets grew by exactly the delayed columns");

            // Each class against a recomputation from its own published sizes.
            ck(matches(recompute(ef.snodeSize(),
                                 [&](std::int32_t jj) { return ef.frontSize(jj); },
                                 [&](std::int32_t jj) {
                                     return ef.frontSize(jj) + ef.updateSize(jj);
                                 }), ef),
               "sizes             : forest agrees with a recomputation from its sizes");
            ck(matches(recompute(sf.snodeSize(),
                                 [&](std::int32_t jj) { return sf.frontSize(jj); },
                                 [&](std::int32_t jj) {
                                     return sf.frontSize(jj) + sf.updateSize(jj);
                                 }), sf),
               "sizes             : symbolic agrees with a recomputation from its sizes");
            ck(matches(recompute(nfs.snodeSize(),
                                 [&](std::int32_t jj) { return nfs.frontSize(jj); },
                                 [&](std::int32_t jj) {
                                     return nfs.frontSize(jj) + nfs.updateSize(jj);
                                 }), nfs),
               "sizes             : static factor agrees with a recomputation");
            ck(matches(recompute(nfd.snodeSize(),
                                 [&](std::int32_t jj) { return nfd.frontSize(jj); },
                                 [&](std::int32_t jj) { return nfd.frontSize(jj) + nfd.delaySize(jj)
                                                             + nfd.updateSize(jj); }), nfd),
               "sizes             : dynamic factor agrees with a recomputation, delays included");

            // The guard the facade puts in front of all three. A Cholesky handed this indefinite
            // matrix analyzes and then refuses, leaving whatever the attempt wrote behind; the
            // facade must report zero rather than that debris.
            DirectSolver<double> dsFailed(Ordering::MMD2, Factorization::Cholesky,
                                          Traversal::LeftLooking);
            const bool refused = dsFailed.analyze(S) && !dsFailed.factor(S);
            ck(refused && dsFailed.numNodeIdx() == 0 && dsFailed.nnz() == 0
                       && dsFailed.numVal() == 0,
               "sizes             : zero after a factorization that refused, not partial counts");

            // The facade forwards to whichever storage is live, and reports zero before there is
            // one. Zero is the count for a factor that does not exist, not a sentinel.
            DirectSolver<double> dsBefore(Ordering::MMD2, Factorization::DynamicLDLT,
                                          Traversal::LeftLooking);
            ck(dsBefore.numNodeIdx() == 0 && dsBefore.nnz() == 0 && dsBefore.numVal() == 0,
               "sizes             : the facade reports zero before factor()");
            ck(dsBefore.analyze(S) && dsBefore.numNodeIdx() == 0 && dsBefore.nnz() == 0,
               "sizes             : and still zero after analyze() alone");

            DirectSolver<double> dsS(Ordering::MMD2, Factorization::StaticLDLT,
                                     Traversal::LeftLooking);
            DirectSolver<double> dsD(Ordering::MMD2, Factorization::DynamicLDLT,
                                     Traversal::LeftLooking);
            const bool forwarded =
                   dsS.analyze(S) && dsS.factor(S) && dsD.analyze(S) && dsD.factor(S)
                && dsS.numNodeIdx() == nfs.numNodeIdx() && dsS.nnz() == nfs.nnz()
                && dsS.numVal()     == nfs.numVal()
                && dsD.numNodeIdx() == nfd.numNodeIdx() && dsD.nnz() == nfd.nnz()
                && dsD.numVal()     == nfd.numVal();
            ck(forwarded,
               "sizes             : the facade forwards the live storage's three exactly");
        }

        // Inertia, against the closed form. One pattern, three shifts, one analysis: the definite
        // matrix, a mildly indefinite one and a strongly indefinite one, the last of which takes
        // 2x2 pivots and so exercises the determinant branch rather than only the diagonal one.
        const std::size_t g = 8;
        const double      shift[] = { 0.0, 1.0, 3.0 };

        DirectSolver<double> dsIn(Ordering::MMD2, Factorization::DynamicLDLT);
        ck(dsIn.analyze(toSparse(shiftedGridLaplacian(g, 0.0))),
           "inertia           : one analysis serves every shift, the pattern being shared");

        bool inertiaMatches = true, inertiaSums = true, tookA2x2 = false;
        for (double sigma : shift) {
            const SparseMatrix<double> S = toSparse(shiftedGridLaplacian(g, sigma));
            Inertia in;
            if (!dsIn.factor(S) || !dsIn.inertia(in)) { inertiaMatches = false; break; }
            if (in.negative != negativeEigenvalues(g, sigma) || in.zero != 0) inertiaMatches = false;
            if (in.positive + in.negative + in.zero != g * g)                 inertiaSums = false;
            if (dsIn.numPivots2x2() > 0)                                      tookA2x2 = true;
        }
        ck(inertiaMatches, "inertia           : negative count matches the closed form at every shift");
        ck(inertiaSums,    "inertia           : positive + negative + zero is the matrix order");
        ck(tookA2x2,       "inertia           : at least one shift took a 2x2, so that branch ran");

        // **One branch of the 2x2 case is not exercised here, and probably cannot be.** A block
        // contributes one of each sign when its determinant is negative and two of the trace's sign
        // when positive, and only the first has ever been seen: 707 accepted blocks over 400 random
        // indefinite matrices were all negative, and mutating the positive case away leaves this
        // suite green. At a root that is a guarantee, since bounded Bunch-Kaufman reaches a 2x2 only
        // after both diagonals have failed their own tests, which forces det < 0. At a non-root it
        // is not: Figure 3.3's test reads |det| and says nothing about the sign, so the branch is
        // defensive rather than dead. Recorded in TESTING_SPECIFICATION under the gaps.

        // The same three matrices under Cholesky and static LDL, where they answer at all: a
        // congruence preserves inertia, so three different factors of one matrix must agree.
        DirectSolver<double> dsCh(Ordering::MMD2, Factorization::Cholesky);
        DirectSolver<double> dsSt(Ordering::MMD2, Factorization::StaticLDLT);
        const SparseMatrix<double> definite = toSparse(shiftedGridLaplacian(g, 0.0));
        Inertia inCh, inSt, inDy;
        const bool agree =
               dsCh.analyze(definite) && dsCh.factor(definite) && dsCh.inertia(inCh)
            && dsSt.analyze(definite) && dsSt.factor(definite) && dsSt.inertia(inSt)
            && dsIn.factor(definite)  && dsIn.inertia(inDy)
            && inCh.positive == g * g && inCh.negative == 0
            && inSt.positive == inCh.positive && inSt.negative == inCh.negative
            && inDy.positive == inCh.positive && inDy.negative == inCh.negative;
        ck(agree, "inertia           : Cholesky, static and dynamic agree on the definite matrix");

        // The two cases it declines rather than guessing.
        DirectSolver<double> dsUnfactored(Ordering::MMD2, Factorization::DynamicLDLT);
        Inertia unused;
        ck(!dsUnfactored.inertia(unused),
           "inertia           : refused before a factorization exists");

        DirectSolver<std::complex<double>> dsSym(Ordering::MMD2, Factorization::DynamicLDLT);
        ck(!dsSym.inertia(unused),
           "inertia           : refused for complex-symmetric LDLT, whose eigenvalues are complex");
    }

    std::cout << "\nPipeline tests: " << pass << "/" << (pass + fail) << " passed\n";
    return fail == 0 ? 0 : 1;
}
