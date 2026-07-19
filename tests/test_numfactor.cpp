// test_numfactor.cpp - the numeric factorization, against a dense Cholesky.
//
// The oracle is a textbook dense Cholesky, written independently of the sparse code. Every entry
// the sparse factor holds must match it. That checks *values*, not merely structure, which is the
// thing the symbolic tests could not reach.
//
// Four combinations are exercised, and the complex ones are the point: 0.9's complex Cholesky
// calls SYRK and TRSM('T'), which is the complex-*symmetric* pattern, on a factor its own POTRF
// treats as Hermitian. That is wrong whenever L has a complex entry, and it is what these tests
// would catch. See the factorization-space entry in DESIGN_DECISIONS.

#include "oblio/ElmForestEngine.h"
#include "oblio/NumFactorDynamic.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/NumFactorStatic.h"
#include "oblio/OrderEngine.h"
#include "oblio/SymFactorEngine.h"

#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace Oblio;
using Cplx = std::complex<double>;

namespace {

int pass = 0, fail = 0;
void ck(bool ok, const std::string& what) {
    if (ok) { ++pass; std::cout << "  PASS  " << what << "\n"; }
    else    { ++fail; std::cout << "  FAIL  " << what << "\n"; }
}

double conj_(double v) { return v; }
Cplx   conj_(Cplx v)   { return std::conj(v); }

// The oracle: dense Cholesky, A = L L^H (which is L L^T for real). Independent of everything the
// solver does, deliberately: it shares no code, no storage, no traversal.
template<class Val>
std::vector<std::vector<Val>> denseCholesky(const std::vector<std::vector<Val>>& A) {
    const std::size_t n = A.size();
    std::vector<std::vector<Val>> L(n, std::vector<Val>(n, Val(0)));
    for (std::size_t j = 0; j < n; ++j) {
        Val d = A[j][j];
        for (std::size_t k = 0; k < j; ++k)
            d -= L[j][k] * conj_(L[j][k]);
        L[j][j] = Val(std::sqrt(std::real(d)));   // Hermitian => the diagonal is real
        for (std::size_t i = j + 1; i < n; ++i) {
            Val s = A[i][j];
            for (std::size_t k = 0; k < j; ++k)
                s -= L[i][k] * conj_(L[j][k]);
            L[i][j] = s / L[j][j];
        }
    }
    return L;
}

// A random Hermitian positive definite matrix: a band, some scattered off-band entries, and a
// diagonal made dominant, which guarantees positive definiteness and keeps the diagonal real (as
// Hermitian requires).
template<class Val>
std::vector<std::vector<Val>> randomHpd(std::size_t n, std::mt19937& rng, int density) {
    std::vector<std::vector<Val>> A(n, std::vector<Val>(n, Val(0)));
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    auto entry = [&]() -> Val {
        if constexpr (std::is_same_v<Val, double>) return Val(u(rng));
        else                                       return Val(u(rng), u(rng));
    };

    for (std::size_t i = 1; i < n; ++i) {
        const Val v = entry();
        A[i][i - 1] = v;
        A[i - 1][i] = conj_(v);
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 2; j < n; ++j)
            if (rng() % 100 < static_cast<unsigned>(density)) {
                const Val v = entry();
                A[j][i] = v;
                A[i][j] = conj_(v);
            }
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j) off += std::abs(A[i][j]);
        A[i][i] = Val(off + 1.0);
    }
    return A;
}

template<class Val>
SparseMatrix<Val> toSparse(const std::vector<std::vector<Val>>& A) {
    const std::size_t n = A.size();
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<Val>          val;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i)
            if (std::abs(A[i][j]) != 0.0) {
                rowIdx.push_back(static_cast<std::int32_t>(i));
                val.push_back(A[i][j]);
            }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<Val>(n, colPtr, rowIdx, val);
}

// The worst discrepancy between the sparse factor and the dense one, over every entry the sparse
// factor holds. Returns -1 if the factorization failed outright.
//
// **The oracle must factor the permuted matrix.** Cholesky is not permutation-equivariant:
//
//     chol(P A P^T)  !=  P chol(A) P^T
//
// The solver factors P A P^T, and its indices are in the factor's (permuted) ordering. So the
// dense oracle has to be handed P A P^T too, and then the two are directly comparable. Factoring
// A and comparing against permuted indices compares the factors of two different matrices, which
// is invisible whenever P is the identity and wrong the moment it is not.
template<class Val>
double compare(const SparseMatrix<Val>& A, const Permutation& p,
               const std::vector<std::vector<Val>>& dense, Traversal traversal) {
    ElmForest f;
    ElmForestEngine fe;
    if (!fe.compute(A, p, f)) return -1;

    SymFactor s;
    SymFactorEngine se;
    if (!se.compute(A, p, f, s)) return -1;

    NumFactorStatic<Val> nf;
    NumFactorEngine ne(Factorization::Cholesky, traversal);
    if (!ne.compute(A, p, s, nf)) return -1;

    // P A P^T, in the factor's ordering: row li of the permuted matrix is row newToOld[li] of A.
    const std::size_t n = dense.size();
    const std::vector<std::int32_t>& newToOld = p.newToOld();
    std::vector<std::vector<Val>> permuted(n, std::vector<Val>(n, Val(0)));
    for (std::size_t li = 0; li < n; ++li)
        for (std::size_t lj = 0; lj < n; ++lj)
            permuted[li][lj] = dense[newToOld[li]][newToOld[lj]];

    const std::vector<std::vector<Val>> L = denseCholesky(permuted);

    double worst = 0;
    for (std::size_t kk = 0; kk < nf.snodeSize(); ++kk) {
        const std::size_t   frontSize = nf.frontSize(kk);
        const std::size_t   numNodeIdx = frontSize + nf.updateSize(kk);
        const std::int32_t* nodeIdx    = nf.nodeIdx().data() + nf.snodeNodeIdxPtr()[kk];
        const Val*          block     = nf.val().data() + nf.snodeValPtr()[kk];

        // The block is column-major, numNodeIdx by frontSize. Only its lower part is meaningful: the
        // strictly upper triangle of the front is allocated and left zero, so BLAS can take the
        // whole rectangle.
        for (std::size_t lclCol = 0; lclCol < frontSize; ++lclCol) {
            const std::size_t col = static_cast<std::size_t>(nodeIdx[lclCol]);
            for (std::size_t lclRow = lclCol; lclRow < numNodeIdx; ++lclRow) {
                const std::size_t row = static_cast<std::size_t>(nodeIdx[lclRow]);
                worst = std::max(worst,
                                 std::abs(block[lclCol * numNodeIdx + lclRow] - L[row][col]));
            }
        }
    }
    return worst;
}

// Every combination of scalar type and traversal, over many random matrices, natural and AMD
// ordered. One number: the worst error anywhere.
template<class Val>
double sweep(int trials, Traversal traversal, std::mt19937& rng, int& failures) {
    OrderEngine ord(OrderMethod::AMD);
    double worst = 0;
    for (int t = 0; t < trials; ++t) {
        const std::size_t n = 4 + rng() % 12;
        const auto dense = randomHpd<Val>(n, rng, 25);
        const SparseMatrix<Val> A = toSparse(dense);

        Permutation pNat(n), pAmd;
        if (!ord.compute(A, pAmd)) { ++failures; continue; }

        for (const Permutation& p : {pNat, pAmd}) {
            const double d = compare(A, p, dense, traversal);
            if (d < 0) ++failures;
            else       worst = std::max(worst, d);
        }
    }
    return worst;
}

// A 2D Laplacian on a grid. Genuinely sparse (five nonzeros a row), which the random matrices
// above are not: they are small and 25 percent dense off-band, so their fill is large, their
// forests are shallow, and the ordering has little to do. A grid is the opposite, and it is what
// exercises the structural machinery:
//
//   16x16 grid, natural order:  256 columns, 240 supernodes, forest height 240, fill x3.3
//
// A forest that deep means a supernode's updates travel a long way up the tree, so left-looking's
// `owed` lists get long and its `pos` bookkeeping is genuinely under test. Under AMD the fill
// drops and the supernodes grow, which is the other half of what the structural phases are for.
//
// Complex version: the same graph with complex Hermitian couplings, still diagonally dominant.
template<class Val>
std::vector<std::vector<Val>> gridLaplacian(std::size_t g, std::mt19937& rng) {
    const std::size_t n = g * g;
    std::vector<std::vector<Val>> A(n, std::vector<Val>(n, Val(0)));
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    auto coupling = [&]() -> Val {
        if constexpr (std::is_same_v<Val, double>) return Val(-1.0);
        else                                       return Val(-1.0, u(rng) * 0.5);
    };

    for (std::size_t y = 0; y < g; ++y)
        for (std::size_t x = 0; x < g; ++x) {
            const std::size_t j = y * g + x;
            A[j][j] = Val(4.5);   // dominant over the four couplings, so HPD
            if (x + 1 < g) { const Val v = coupling(); A[j + 1][j] = v; A[j][j + 1] = conj_(v); }
            if (y + 1 < g) { const Val v = coupling(); A[j + g][j] = v; A[j][j + g] = conj_(v); }
        }
    return A;
}

// Reconstruct the matrix from an LDL factor: L D L^T (or L D L^H) must equal P A P^T.
//
// A different kind of oracle from Cholesky's, and a better one where it applies. Cholesky is
// compared against an independently written dense Cholesky, which is a *second implementation*.
// LDL is checked by **reconstruction**, which needs no second implementation at all: multiply the
// factor back out and see whether the matrix comes back. That checks D, L, the storage layout and
// the supernodal assembly in one statement, and it works unchanged for all three symmetries.
//
// The factor's layout, which the reconstruction has to know:
//
//   the diagonal        holds D
//   below it            holds L, which is UNIT lower triangular (the 1s are implicit)
//   above it            holds U = D L^T, which we ignore: L and D are enough to rebuild A
template<class Val>
double reconstructLdl(const NumFactorStatic<Val>& f,
                      const std::vector<std::vector<Val>>& permuted, bool hermitian) {
    const std::size_t n = permuted.size();
    std::vector<std::vector<Val>> L(n, std::vector<Val>(n, Val(0)));
    std::vector<Val> D(n, Val(0));

    for (std::size_t kk = 0; kk < f.snodeSize(); ++kk) {
        const std::size_t   frontSize = f.frontSize(kk);
        const std::size_t   numNodeIdx = frontSize + f.updateSize(kk);
        const std::int32_t* nodeIdx    = f.nodeIdx().data() + f.snodeNodeIdxPtr()[kk];
        const Val*          block     = f.val().data() + f.snodeValPtr()[kk];

        for (std::size_t lclCol = 0; lclCol < frontSize; ++lclCol) {
            const std::size_t col = static_cast<std::size_t>(nodeIdx[lclCol]);
            D[col]      = block[lclCol * numNodeIdx + lclCol];
            L[col][col] = Val(1);
            for (std::size_t lclRow = lclCol + 1; lclRow < numNodeIdx; ++lclRow)
                L[static_cast<std::size_t>(nodeIdx[lclRow])][col] =
                    block[lclCol * numNodeIdx + lclRow];
        }
    }

    double worst = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            Val s(0);
            for (std::size_t k = 0; k <= j; ++k) {
                const Val ljk = hermitian ? conj_(L[j][k]) : L[j][k];
                s += L[i][k] * D[k] * ljk;
            }
            worst = std::max(worst, std::abs(s - permuted[i][j]));
        }
    return worst;
}

// Factor with LDL, then reconstruct. Returns -1 if the factorization failed.
template<class Val>
double compareLdl(const SparseMatrix<Val>& A, const Permutation& p,
                  const std::vector<std::vector<Val>>& dense,
                  Factorization factorization, Traversal traversal, std::size_t& perturbations) {
    ElmForest f;
    ElmForestEngine fe;
    if (!fe.compute(A, p, f)) return -1;

    SymFactor s;
    SymFactorEngine se;
    if (!se.compute(A, p, f, s)) return -1;

    NumFactorStatic<Val> nf;
    NumFactorEngine ne(factorization, traversal);
    if (!ne.compute(A, p, s, nf)) return -1;
    perturbations += std::as_const(nf).numPerturbations();

    const std::size_t n = dense.size();
    const std::vector<std::int32_t>& newToOld = p.newToOld();
    std::vector<std::vector<Val>> permuted(n, std::vector<Val>(n, Val(0)));
    for (std::size_t li = 0; li < n; ++li)
        for (std::size_t lj = 0; lj < n; ++lj)
            permuted[li][lj] = dense[newToOld[li]][newToOld[lj]];

    const bool hermitian = (factorization == Factorization::StaticLDLH
                            || factorization == Factorization::DynamicLDLH);
    return reconstructLdl(nf, permuted, hermitian);
}

// Factor the same matrix into the flat factor and the per-supernode factor, same factorization and
// same traversal, and return the largest block difference. Both run identical arithmetic through
// identical kernels; only where a block lives differs, so the two factors must come out bit for
// bit the same, and the difference must be exactly zero. Returns -1 if either factorization is
// refused, or if the two disagree on structure or perturbation count.
template<class Val>
double staticVsDynamic(const SparseMatrix<Val>& A, const Permutation& p,
                       Factorization factorization, Traversal traversal) {
    ElmForest f; ElmForestEngine fe;
    if (!fe.compute(A, p, f)) return -1;
    SymFactor s; SymFactorEngine se;
    if (!se.compute(A, p, f, s)) return -1;

    NumFactorEngine ne(factorization, traversal);
    NumFactorStatic<Val>  sfac;
    NumFactorDynamic<Val> dfac;
    if (!ne.compute(A, p, s, sfac)) return -1;
    if (!ne.compute(A, p, s, dfac)) return -1;

    const NumFactorStatic<Val>&  cs = sfac;
    const NumFactorDynamic<Val>& cd = dfac;
    if (cs.snodeSize() != cd.snodeSize()) return -1;
    if (cs.numPerturbations() != cd.numPerturbations()) return -1;

    double worst = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(cs.snodeSize()); ++kk) {
        const std::size_t len = (cs.frontSize(kk) + cs.updateSize(kk)) * cs.frontSize(kk);
        const Val* a = cs.val(kk);
        const Val* b = cd.val(kk);
        for (std::size_t t = 0; t < len; ++t)
            worst = std::max(worst, std::abs(a[t] - b[t]));
    }
    return worst;
}

// A random complex SYMMETRIC matrix (A = A^T, not Hermitian), which is what LDLT is for and what
// Cholesky may never be handed. Diagonally dominant, so the pivots stay away from zero.
std::vector<std::vector<Cplx>> randomComplexSymmetric(std::size_t n, std::mt19937& rng) {
    std::vector<std::vector<Cplx>> A(n, std::vector<Cplx>(n, Cplx(0)));
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (std::size_t i = 1; i < n; ++i) {
        const Cplx v(u(rng), u(rng));
        A[i][i - 1] = v;
        A[i - 1][i] = v;          // plain transpose, NOT conjugate
    }
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 2; j < n; ++j)
            if (rng() % 100 < 25) {
                const Cplx v(u(rng), u(rng));
                A[j][i] = v;
                A[i][j] = v;
            }
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0;
        for (std::size_t j = 0; j < n; ++j)
            if (i != j) off += std::abs(A[i][j]);
        A[i][i] = Cplx(off + 1.0, 0.0);
    }
    return A;
}

// Slice 1 check for dynamic LDL: factor a dense symmetric matrix with DynamicLDLT (left-looking),
// then reconstruct L D L^T from the block and pivotType and compare to the pivoted, permuted
// matrix. Returns the worst entry difference, -1 if the factorization was refused (not one front,
// or a column delayed with no ancestor to take it). twoByTwo counts the 2x2 pivots chosen.
double dynamicLdlWorst(const std::vector<std::vector<double>>& dense, int& twoByTwo) {
    const std::int32_t n = static_cast<std::int32_t>(dense.size());
    const SparseMatrix<double> A = toSparse(dense);

    OrderEngine ord(OrderMethod::AMD);
    Permutation p;
    if (!ord.compute(A, p)) return -1;

    ElmForest f; ElmForestEngine fe;
    if (!fe.compute(A, p, f)) return -1;
    SymFactor s; SymFactorEngine se;
    if (!se.compute(A, p, f, s)) return -1;

    NumFactorDynamic<double> nf;
    NumFactorEngine ne(Factorization::DynamicLDLT, Traversal::LeftLooking);
    if (!ne.compute(A, p, s, nf)) return -1;
    if (nf.snodeSize() != 1) return -1;

    const std::vector<std::int32_t>& newToOld = p.newToOld();
    std::vector<std::vector<double>> permuted(n, std::vector<double>(n));
    for (std::int32_t li = 0; li < n; ++li)
        for (std::int32_t lj = 0; lj < n; ++lj)
            permuted[li][lj] = dense[newToOld[li]][newToOld[lj]];

    const NumFactorDynamic<double>& cnf = nf;
    const double*       blk = cnf.val(0);
    const std::int32_t* idx = cnf.nodeIdx(0);
    const std::vector<std::int32_t>& pt = cnf.pivotType();
    const auto at = [n](std::int32_t r, std::int32_t c) {
        return static_cast<std::size_t>(c) * static_cast<std::size_t>(n) + static_cast<std::size_t>(r);
    };

    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0));
    for (std::int32_t j = 0; j < n; ) {
        if (pt[idx[j]] == 1) {
            L[j][j] = 1;
            D[j][j] = blk[at(j, j)];
            for (std::int32_t i = j + 1; i < n; ++i) L[i][j] = blk[at(i, j)];
            ++j;
        } else {                                    // 2x2 pivot: columns j, j+1
            ++twoByTwo;
            L[j][j] = 1; L[j + 1][j + 1] = 1;       // L[j+1][j] stays 0
            D[j][j]         = blk[at(j, j)];
            D[j + 1][j + 1] = blk[at(j + 1, j + 1)];
            D[j][j + 1] = D[j + 1][j] = blk[at(j, j + 1)];
            for (std::int32_t i = j + 2; i < n; ++i) { L[i][j] = blk[at(i, j)]; L[i][j + 1] = blk[at(i, j + 1)]; }
            j += 2;
        }
    }

    std::vector<std::vector<double>> LD(n, std::vector<double>(n, 0));
    for (std::int32_t i = 0; i < n; ++i)
        for (std::int32_t k = 0; k < n; ++k) {
            double acc = 0;
            for (std::int32_t m = 0; m < n; ++m) acc += L[i][m] * D[m][k];
            LD[i][k] = acc;
        }
    double worst = 0;
    for (std::int32_t i = 0; i < n; ++i)
        for (std::int32_t jc = 0; jc < n; ++jc) {
            double acc = 0;
            for (std::int32_t k = 0; k < n; ++k) acc += LD[i][k] * L[jc][k];
            worst = std::max(worst, std::abs(acc - permuted[idx[i]][idx[jc]]));
        }
    return worst;
}

} // namespace

int main() {
    std::mt19937 rng(20260713);
    const double tol = 1e-10;

    // A worked case, small enough to reason about. The 3x3 with A-edges 1-2 and 1-3 from the
    // notes: one fundamental supernode of all three columns.
    {
        std::vector<std::vector<double>> A = {{4, 2, 2}, {2, 5, 0}, {2, 0, 6}};
        // Hermitian, and diagonally dominant enough to be positive definite.
        const SparseMatrix<double> S = toSparse(A);
        Permutation p(3);
        const double d = compare(S, p, A, Traversal::LeftLooking);
        ck(d >= 0 && d < tol, "3x3 dense-ish       : matches dense Cholesky");
    }

    // A matrix that is NOT positive definite must be refused, not silently mangled.
    {
        std::vector<std::vector<double>> A = {{1, 2}, {2, 1}};   // eigenvalues 3 and -1
        const SparseMatrix<double> S = toSparse(A);
        Permutation p(2);
        ElmForest f; ElmForestEngine fe; fe.compute(S, p, f);
        SymFactor s; SymFactorEngine se; se.compute(S, p, f, s);
        NumFactorStatic<double> nf;
        NumFactorEngine ne(Factorization::Cholesky, Traversal::LeftLooking);
        ck(!ne.compute(S, p, s, nf), "indefinite 2x2      : refused, not factored");
    }

    // The sweeps. Both traversals must agree with the oracle, and therefore with each other.
    int failures = 0;
    const double rLeft  = sweep<double>(40, Traversal::LeftLooking,  rng, failures);
    const double rRight = sweep<double>(40, Traversal::RightLooking, rng, failures);
    const double cLeft  = sweep<Cplx>  (40, Traversal::LeftLooking,  rng, failures);
    const double cRight = sweep<Cplx>  (40, Traversal::RightLooking, rng, failures);

    ck(failures == 0, "random x40 x4       : every matrix factored");
    ck(rLeft  < tol, "real    left-looking : matches dense Cholesky");
    ck(rRight < tol, "real    right-looking: matches dense Cholesky");

    // These two are the ones 0.9 would fail. Its complex Cholesky uses SYRK and TRSM('T'), the
    // complex-symmetric pattern, while its POTRF is Hermitian. With genuinely complex data the
    // factor comes out wrong.
    ck(cLeft  < tol, "complex left-looking : matches dense Cholesky (Hermitian)");
    ck(cRight < tol, "complex right-looking: matches dense Cholesky (Hermitian)");

    // The grid. This is the sparse case: deep forest, real fill, an ordering that matters.
    {
        OrderEngine ord(OrderMethod::AMD);
        double worstR = 0, worstC = 0;
        int gridFail = 0;
        std::size_t supNat = 0, supAmd = 0, idxNat = 0, idxAmd = 0, heightNat = 0;

        for (std::size_t g : {6u, 10u}) {
            const std::size_t n = g * g;

            const auto denseR = gridLaplacian<double>(g, rng);
            const auto denseC = gridLaplacian<Cplx>(g, rng);
            const SparseMatrix<double> AR = toSparse(denseR);
            const SparseMatrix<Cplx>   AC = toSparse(denseC);

            Permutation pNat(n), pAmd;
            if (!ord.compute(AR, pAmd)) { ++gridFail; continue; }

            // Report the structure once, so the test shows it is doing real work.
            if (g == 10) {
                ElmForest f; ElmForestEngine fe; fe.compute(AR, pNat, f);
                SymFactor s; SymFactorEngine se; se.compute(AR, pNat, f, s);
                supNat = s.snodeSize(); idxNat = s.numNodeIdx(); heightNat = f.height();

                ElmForest fa; fe.compute(AR, pAmd, fa);
                SymFactor sa; se.compute(AR, pAmd, fa, sa);
                supAmd = sa.snodeSize(); idxAmd = sa.numNodeIdx();
            }

            for (const Permutation& p : {pNat, pAmd})
                for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
                    const double dr = compare(AR, p, denseR, tr);
                    const double dc = compare(AC, p, denseC, tr);
                    if (dr < 0 || dc < 0) ++gridFail;
                    else { worstR = std::max(worstR, dr); worstC = std::max(worstC, dc); }
                }
        }

        ck(gridFail == 0 && worstR < tol,
           "10x10 grid, real    : matches dense Cholesky, both traversals, natural and AMD");
        ck(gridFail == 0 && worstC < tol,
           "10x10 grid, complex : matches dense Cholesky, both traversals, natural and AMD");
        ck(heightNat > 50 && supAmd < supNat && idxAmd < idxNat,
           "10x10 grid          : structure exercised (forest height "
           + std::to_string(heightNat) + "; AMD cuts supernodes "
           + std::to_string(supNat) + "->" + std::to_string(supAmd)
           + ", indices " + std::to_string(idxNat) + "->" + std::to_string(idxAmd) + ")");
    }

    // Static LDL. Three symmetries, two traversals, checked by reconstruction.
    {
        OrderEngine ord(OrderMethod::AMD);
        int ldlFail = 0;
        std::size_t perturbations = 0;
        double wR = 0, wCs = 0, wCh = 0;

        for (int trial = 0; trial < 30; ++trial) {
            const std::size_t n = 5 + rng() % 12;

            const auto denseR  = randomHpd<double>(n, rng, 25);          // real symmetric
            const auto denseCs = randomComplexSymmetric(n, rng);         // complex symmetric
            const auto denseCh = randomHpd<Cplx>(n, rng, 25);            // complex Hermitian

            const SparseMatrix<double> AR  = toSparse(denseR);
            const SparseMatrix<Cplx>   ACs = toSparse(denseCs);
            const SparseMatrix<Cplx>   ACh = toSparse(denseCh);

            Permutation pNat(n), pAmd;
            if (!ord.compute(AR, pAmd)) { ++ldlFail; continue; }

            for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
                const double dr  = compareLdl(AR,  pAmd, denseR,  Factorization::StaticLDLT, tr, perturbations);
                const double dcs = compareLdl(ACs, pAmd, denseCs, Factorization::StaticLDLT, tr, perturbations);
                const double dch = compareLdl(ACh, pAmd, denseCh, Factorization::StaticLDLH, tr, perturbations);
                if (dr < 0 || dcs < 0 || dch < 0) ++ldlFail;
                else {
                    wR  = std::max(wR,  dr);
                    wCs = std::max(wCs, dcs);
                    wCh = std::max(wCh, dch);
                }
            }
        }

        ck(ldlFail == 0 && wR < tol,
           "LDLT real           : L D L^T reconstructs A, both traversals");
        ck(ldlFail == 0 && wCs < tol,
           "LDLT complex        : L D L^T reconstructs A (complex SYMMETRIC, D complex)");

        // 0.9 does not have this one at all: its complex LDL is symmetric only.
        ck(ldlFail == 0 && wCh < tol,
           "LDLH complex        : L D L^H reconstructs A (complex HERMITIAN, D real)");

        ck(perturbations == 0,
           "LDL                 : no perturbations needed (diagonally dominant input)");
    }

    // The static factorizations run into the dynamic factor. Same traversals, same kernels, only
    // per-supernode storage in place of one flat buffer, so the factor must come out identical to
    // the flat one, block for block. This exercises NumFactorDynamic's whole read/produce API
    // (setSymFactor, the accessors, the friend write path) without any of the growth verbs dynamic
    // LDL will add.
    {
        OrderEngine ord(OrderMethod::AMD);
        int dynFail = 0;
        double worst = 0;

        for (int trial = 0; trial < 30; ++trial) {
            const std::size_t n = 5 + rng() % 12;

            const auto denseR  = randomHpd<double>(n, rng, 25);          // real symmetric / Hermitian
            const auto denseCs = randomComplexSymmetric(n, rng);         // complex symmetric
            const auto denseCh = randomHpd<Cplx>(n, rng, 25);            // complex Hermitian

            const SparseMatrix<double> AR  = toSparse(denseR);
            const SparseMatrix<Cplx>   ACs = toSparse(denseCs);
            const SparseMatrix<Cplx>   ACh = toSparse(denseCh);

            Permutation pAmd;
            if (!ord.compute(AR, pAmd)) { ++dynFail; continue; }

            for (Traversal tr : {Traversal::LeftLooking, Traversal::RightLooking}) {
                const double dChR  = staticVsDynamic(AR,  pAmd, Factorization::Cholesky,   tr);
                const double dChC  = staticVsDynamic(ACh, pAmd, Factorization::Cholesky,   tr);
                const double dLdR  = staticVsDynamic(AR,  pAmd, Factorization::StaticLDLT, tr);
                const double dLdCs = staticVsDynamic(ACs, pAmd, Factorization::StaticLDLT, tr);
                const double dLdCh = staticVsDynamic(ACh, pAmd, Factorization::StaticLDLH, tr);
                for (double d : {dChR, dChC, dLdR, dLdCs, dLdCh}) {
                    if (d < 0) ++dynFail;
                    else worst = std::max(worst, d);
                }
            }
        }

        ck(dynFail == 0 && worst == 0.0,
           "static into dynamic : identical factor, flat vs per-supernode, all symmetries, both traversals");
    }

    // Dynamic LDL, slice 1: a dense front, real, left-looking. Reconstruct L D L^T and check it
    // matches the pivoted matrix. Dense PD inputs exercise 1x1 pivots; dense indefinite inputs
    // exercise the 2x2 pivots and the swap machinery.
    {
        int checked = 0, twoByTwo = 0;
        double worst = 0;
        std::uniform_real_distribution<double> u(-1.0, 1.0);

        for (int trial = 0; trial < 80; ++trial) {
            const std::size_t n = 3 + rng() % 8;
            std::vector<std::vector<double>> dense;

            if (trial % 2 == 0) {
                dense = randomHpd<double>(n, rng, 100);              // dense PD: 1x1 pivots
            } else {
                dense.assign(n, std::vector<double>(n, 0.0));         // dense symmetric indefinite
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t jc = i; jc < n; ++jc) {
                        const double v = u(rng);
                        dense[i][jc] = v;
                        dense[jc][i] = v;
                    }
            }

            const double w = dynamicLdlWorst(dense, twoByTwo);
            if (w < 0) continue;                                      // refused (multi-front or delayed): skip
            worst = std::max(worst, w);
            ++checked;
        }

        ck(checked > 20 && twoByTwo > 0 && worst < 1e-9,
           "dynamic LDL slice 1 : L D L^T reconstructs the dense front (1x1 and 2x2 pivots)");
    }

    std::cout << "\nNumFactor tests: " << pass << "/" << (pass + fail) << " passed\n";
    return fail == 0 ? 0 : 1;
}
