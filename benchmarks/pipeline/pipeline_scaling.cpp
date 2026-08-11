// pipeline_scaling.cpp - how a solve's cost grows with problem size, over every factorization
// Oblio has and all three value types.
//
//   make scale2d      the square ladder
//   make scale3d      the cubic ladder
//   ./pipeline_scaling_cpp 2d 100 200 400      any sides
//
// `pipeline_timing.cpp` beside this asks where a solve's time goes at a few fixed sizes.
// This asks how each phase GROWS, which needs a controlled family rather than a heterogeneous
// set: `../matrices` measures real matrices and cannot answer a scaling question, because two
// matrices of the same n there differ in fill by orders of magnitude. A grid ladder holds the
// structure fixed and moves one parameter.
//
// THE TWO FAMILIES ARE NOT ONE FAMILY WITH DIFFERENT NUMBERS, which is why both ladders exist.
// In 2D ten times the columns costs about 200 times the fill; in 3D sixty-four times the columns
// costs 660 times the fill. A 64 cube at n = 262144 fills five times more than a 1000 square at
// n = 1000000. A growth exponent measured on one says nothing about the other.
//
// EVERY FACTORIZATION, WHICH THE OTHER TWO BENCHMARK FOLDERS DO NOT DO. `../matrices` runs
// Cholesky only on the performance side, because a real matrix that is positive definite is the
// only kind Cholesky takes and complex definite matrices are effectively absent from the
// collection: of its 23 complex square pattern-symmetric matrices exactly one is marked positive
// definite, and that one is complex SYMMETRIC rather than Hermitian. A generated grid has no such
// problem, so this is the only place complex Cholesky can be measured at all.
//
//   real                Cholesky, static LDL^T, dynamic LDL^T
//   complex Hermitian   Cholesky, static LDL^H, dynamic LDL^H
//   complex symmetric   static LDL^T, dynamic LDL^T
//
// The two omissions are not gaps. Over the reals LDL^H IS LDL^T, the conjugate being the
// identity, so running both would repeat a column. And complex Cholesky needs Hermitian positive
// definiteness, so it does not apply to a complex symmetric matrix; worse, `zpotrf` reads one
// triangle and assumes the other is its conjugate, so handing it a complex symmetric matrix
// SUCCEEDS and returns a plausible wrong answer. `examples/example_pipeline_complex.cpp` says the
// same at greater length.
//
// NOTHING HERE SHOULD PERTURB OR PIVOT, and the run checks it. Every matrix is diagonally
// dominant by 25 percent: the diagonal is 1.25 times the sum of the off-diagonal magnitudes in
// its row, in all three arms, so no pivot is ever too small to divide by and no column is ever
// unusable. A static LDL that perturbs or a dynamic LDL that delays here means the construction
// is wrong, not that the matrix is hard, and the run says so loudly. The counts are printed
// beside the residual for exactly that reason.
//
// THE MARGIN IS DELIBERATE AND THE REAL 2D GRID DOES NOT HAVE IT. `pipeline_timing.cpp` uses a
// diagonal of 4 against four neighbors, so its interior rows have a dominance margin of ZERO;
// that matrix is positive definite by structure rather than by dominance, which is fine for
// Cholesky and is not enough for the complex symmetric arm, where dominance is the only guarantee
// available. So the ladder's matrices carry a margin and are NOT the same matrices the short
// ladder uses. Rows here and there are not comparable, deliberately.

#ifdef __APPLE__
#include <pthread.h>
#endif

#include "oblio/DirectSolver.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

namespace {

// ------------------------------------------------------------------------------------------
// The grids
// ------------------------------------------------------------------------------------------

// The pattern, shared by all three arms so that a difference between them is arithmetic and never
// structure. Each column's rows are SORTED with the diagonal, which is not tidiness: SparseMatrix
// requires ascending rows, and the order within a column is also a tie-break input, since it
// decides the content order of a clique. A 3D builder written the natural way is not ascending,
// and getting that wrong reads as an ordering divergence rather than as a harness fault.
struct Pattern {
    std::size_t               size = 0;
    std::size_t               degree = 0;   // interior degree: 4 in 2D, 6 in 3D
    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
};

Pattern gridPattern(int m, bool cubic) {
    const int n = cubic ? m * m * m : m * m;
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));

    if (!cubic) {
        auto id = [&](int r, int c) { return r * m + c; };
        for (int r = 0; r < m; ++r)
            for (int c = 0; c < m; ++c) {
                const int u = id(r, c);
                if (r > 0)     adj[u].push_back(id(r - 1, c));
                if (c > 0)     adj[u].push_back(id(r, c - 1));
                if (r + 1 < m) adj[u].push_back(id(r + 1, c));
                if (c + 1 < m) adj[u].push_back(id(r, c + 1));
            }
    } else {
        auto id = [&](int x, int y, int z) { return (z * m + y) * m + x; };
        for (int z = 0; z < m; ++z)
            for (int y = 0; y < m; ++y)
                for (int x = 0; x < m; ++x) {
                    const int u = id(x, y, z);
                    if (z > 0)     adj[u].push_back(id(x, y, z - 1));
                    if (y > 0)     adj[u].push_back(id(x, y - 1, z));
                    if (x > 0)     adj[u].push_back(id(x - 1, y, z));
                    if (x + 1 < m) adj[u].push_back(id(x + 1, y, z));
                    if (y + 1 < m) adj[u].push_back(id(x, y + 1, z));
                    if (z + 1 < m) adj[u].push_back(id(x, y, z + 1));
                }
    }

    Pattern p;
    p.size = static_cast<std::size_t>(n);
    p.degree = cubic ? 6 : 4;
    p.colPtr.assign(p.size + 1, 0);
    for (int j = 0; j < n; ++j) {
        std::vector<int> col = adj[static_cast<std::size_t>(j)];
        col.push_back(j);
        std::sort(col.begin(), col.end());
        for (int i : col)
            p.rowIdx.push_back(static_cast<std::int32_t>(i));
        p.colPtr[static_cast<std::size_t>(j) + 1] = p.rowIdx.size();
    }
    return p;
}

// THE DOMINANCE MARGIN, and the one number the whole construction turns on. Every arm gets a
// diagonal of `kMargin * degree * |offDiagonal|`, so an interior row's diagonal exceeds the sum of
// its off-diagonal magnitudes by 25 percent and a boundary row by more. That is what guarantees no
// pivot is ever too small to divide by, in any of the three arms, and it is the only guarantee the
// complex symmetric arm has: `x^T A x` is complex there, so there is no definiteness to appeal to.
const double kMargin = 1.25;

// The off-diagonal, magnitude sqrt(1.25) in the complex arms and 1 in the real one. The complex
// value carries a GENUINE imaginary part rather than a zero one, so the complex kernels do complex
// arithmetic rather than multiplying zeros in a complex type.
const double kOffReal = 1.0;
const std::complex<double> kOffComplex(-1.0, 0.5);

double offMagnitude(bool complexValues) {
    return complexValues ? std::abs(kOffComplex) : kOffReal;
}

SparseMatrix<double> realGrid(const Pattern& p) {
    const double diagonal = kMargin * double(p.degree) * offMagnitude(false);
    std::vector<double> val(p.rowIdx.size());
    for (std::size_t j = 0; j < p.size; ++j)
        for (std::size_t cp = p.colPtr[j]; cp < p.colPtr[j + 1]; ++cp)
            val[cp] = (p.rowIdx[cp] == static_cast<std::int32_t>(j)) ? diagonal : -kOffReal;
    return SparseMatrix<double>(p.size, p.colPtr, p.rowIdx, val);
}

// Hermitian: A = A^H, so the diagonal is REAL and the entry above the diagonal is the conjugate of
// the one below. Hermitian and diagonally dominant is Hermitian positive definite, which is what
// makes it a valid input for complex Cholesky.
SparseMatrix<std::complex<double>> hermitianGrid(const Pattern& p) {
    const double diagonal = kMargin * double(p.degree) * offMagnitude(true);
    std::vector<std::complex<double>> val(p.rowIdx.size());
    for (std::size_t j = 0; j < p.size; ++j)
        for (std::size_t cp = p.colPtr[j]; cp < p.colPtr[j + 1]; ++cp) {
            const std::int32_t i = p.rowIdx[cp];
            if (i == static_cast<std::int32_t>(j))
                val[cp] = std::complex<double>(diagonal, 0.0);
            else
                val[cp] = (i > static_cast<std::int32_t>(j)) ? kOffComplex : std::conj(kOffComplex);
        }
    return SparseMatrix<std::complex<double>>(p.size, p.colPtr, p.rowIdx, val);
}

// Complex symmetric: A = A^T, the SAME entry on both sides of the diagonal and no conjugation.
// Not Hermitian and not positive definite in any sense, which is why Cholesky is never handed it.
// The magnitudes match the Hermitian arm exactly, so a difference in time between the two is the
// conjugation and the kernel rather than a different matrix.
SparseMatrix<std::complex<double>> complexSymmetricGrid(const Pattern& p) {
    const double diagonal = kMargin * double(p.degree) * offMagnitude(true);
    std::vector<std::complex<double>> val(p.rowIdx.size());
    for (std::size_t j = 0; j < p.size; ++j)
        for (std::size_t cp = p.colPtr[j]; cp < p.colPtr[j + 1]; ++cp)
            val[cp] = (p.rowIdx[cp] == static_cast<std::int32_t>(j))
                          ? std::complex<double>(diagonal, 0.0)
                          : kOffComplex;
    return SparseMatrix<std::complex<double>>(p.size, p.colPtr, p.rowIdx, val);
}

// ------------------------------------------------------------------------------------------
// Measurement
// ------------------------------------------------------------------------------------------

// Best of three after a warm-up, as in ../ordering and ../pipeline and for the same reason.
template <class Work>
double bestOfThree(Work work) {
    work();
    double best = 1e30;
    for (int trial = 0; trial < 3; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        work();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

std::size_t fill(const SymFactor& sf) {
    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk);
        const std::size_t u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

template <class Val>
double infNorm(const SparseMatrix<Val>& A) {
    std::vector<double> rowSum(A.size(), 0.0);
    for (std::size_t j = 0; j < A.size(); ++j)
        for (std::size_t cp = A.colPtr()[j]; cp < A.colPtr()[j + 1]; ++cp)
            rowSum[static_cast<std::size_t>(A.rowIdx()[cp])] += std::abs(A.val()[cp]);
    double norm = 0.0;
    for (double sum : rowSum)
        norm = std::max(norm, sum);
    return norm;
}

template <class Val>
double infNorm(const Vector<Val>& v) {
    double norm = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i)
        norm = std::max(norm, std::abs(v[i]));
    return norm;
}

struct Result {
    bool        ok = false;
    double      factor = 0;
    double      solve = 0;
    double      residual = 0;
    double      backward = 0;
    std::size_t perturbed = 0;
    std::size_t delayed = 0;
};

// One factorization, one traversal, with the analysis outside the timed region. The residual and
// the backward error are computed from a single residual vector, as in ../matrices: the backward
// error is the verdict, since it does not inherit the conditioning, and the residual is reported
// beside it.
template <class Val>
Result measure(const SparseMatrix<Val>& A, Ordering ordering, Factorization factorization,
               Traversal traversal, double aNorm) {
    Result r;
    DirectSolver<Val> solver(ordering, factorization, traversal);
    if (!solver.analyze(A) || !solver.factor(A))
        return r;

    r.factor = bestOfThree([&] { solver.factor(A); });

    Vector<Val> b(A.size()), x(A.size()), residual(A.size());
    for (std::size_t i = 0; i < A.size(); ++i)
        b[i] = Val(1);

    if (!solver.solve(b, x))
        return r;
    r.solve = bestOfThree([&] { solver.solve(b, x); });

    const MultiplyEngine multiply;
    multiply.residual(A, x, b, residual);

    const double rNorm = infNorm(residual);
    const double xNorm = infNorm(x);
    const double bNorm = infNorm(b);

    r.residual = (bNorm > 0) ? rNorm / bNorm : 0.0;
    r.backward = rNorm / (aNorm * xNorm + bNorm);
    r.perturbed = solver.numPerturbations();
    r.delayed = solver.numDelayedColumns();
    r.ok = true;
    return r;
}

// The four orderings, and the three traversals, in the order the tables print them.
struct Named { Ordering method; const char* name; };
const Named kOrderings[] = {
    {Ordering::MMD,  "MMD"},
    {Ordering::MMD3, "MMD3"},
    {Ordering::AMD,  "AMD"},
    {Ordering::AMD3, "AMD3"},
};
const int kNumOrderings = 4;

struct NamedTraversal { Traversal traversal; const char* name; };
const NamedTraversal kTraversals[] = {
    {Traversal::LeftLooking,  "LL"},
    {Traversal::RightLooking, "RL"},
    {Traversal::Multifrontal, "MF"},
};
const int kNumTraversals = 3;

// The eight arms: a value type, a factorization, and the name the table prints.
enum class Arm { Real, Hermitian, ComplexSymmetric };

struct NamedArm {
    Arm           arm;
    Factorization factorization;
    const char*   name;
};

const NamedArm kArms[] = {
    {Arm::Real,             Factorization::Cholesky,    "real     Cholesky"},
    {Arm::Real,             Factorization::StaticLDLT,  "real     static LDLT"},
    {Arm::Real,             Factorization::DynamicLDLT, "real     dynamic LDLT"},
    {Arm::Hermitian,        Factorization::Cholesky,    "cplx herm Cholesky"},
    {Arm::Hermitian,        Factorization::StaticLDLH,  "cplx herm static LDLH"},
    {Arm::Hermitian,        Factorization::DynamicLDLH, "cplx herm dynamic LDLH"},
    {Arm::ComplexSymmetric, Factorization::StaticLDLT,  "cplx sym static LDLT"},
    {Arm::ComplexSymmetric, Factorization::DynamicLDLT, "cplx sym dynamic LDLT"},
};
const int kNumArms = 8;

} // namespace

int main(int argc, char** argv) {
#ifdef __APPLE__
    // Before anything is timed. A command-line process runs at the default quality-of-service
    // class, which prefers a performance core but permits the scheduler to park the thread on an
    // efficiency one, and that placement is sticky over long stretches rather than jittering per
    // sample, so a minimum over repeats does not filter it.
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    bool cubic = false;
    int  firstSide = 1;
    if (argc > 1) {
        const std::string word = argv[1];
        if      (word == "3d") { cubic = true;  firstSide = 2; }
        else if (word == "2d") { cubic = false; firstSide = 2; }
    }

    // THE LADDERS, six sizes each, chosen so the top of each is about 12 million entries of fill
    // in 2D and 51 million in 3D, which is 0.18 and 0.76 GB of complex values. Extending is one
    // number here, and the next rungs are 800 and 1000 in 2D, 56 and 64 in 3D, the last of which
    // needs 2.9 GB for the complex factor alone.
    std::vector<int> sides = cubic ? std::vector<int>{16, 20, 26, 32, 40, 48}
                                   : std::vector<int>{100, 140, 200, 300, 400, 600};
    if (argc > firstSide) {
        sides.clear();
        for (int k = firstSide; k < argc; ++k)
            sides.push_back(std::atoi(argv[k]));
    }

    std::printf("%s grids, best of three after a warm-up, milliseconds\n",
                cubic ? "cubic" : "square");
    std::printf("every matrix is diagonally dominant by %.0f percent, so NOTHING should perturb\n"
                "or delay; the two rightmost columns are the check\n\n", (kMargin - 1.0) * 100.0);

    std::size_t totalPerturbed = 0;
    std::size_t totalDelayed = 0;

    for (int m : sides) {
        const Pattern pattern = gridPattern(m, cubic);
        const SparseMatrix<double> A = realGrid(pattern);
        const SparseMatrix<std::complex<double>> AH = hermitianGrid(pattern);
        const SparseMatrix<std::complex<double>> AS = complexSymmetricGrid(pattern);

        std::printf("=== %s %d, n %zu, nnz(A) %zu ===\n",
                    cubic ? "cube" : "square", m, pattern.size, A.nnz());

        // The ordering and the analysis read the PATTERN and never a value, so they are timed
        // once per ordering and serve all eight arms.
        std::printf("  %-6s %12s %10s %10s\n", "order", "nnz(L)", "order", "analyze");
        for (int o = 0; o < kNumOrderings; ++o) {
            const OrderEngine oe(kOrderings[o].method);
            Permutation       P;
            if (!oe.compute(A, P)) {
                std::printf("  %-6s %12s\n", kOrderings[o].name, "refused");
                continue;
            }
            ElmForest ef;
            ElmForestEngine().compute(A, P, ef);
            SymFactor sf;
            SymFactorEngine().compute(A, P, ef, sf);

            const double order = bestOfThree([&] { Permutation Q; oe.compute(A, Q); });
            const double analyze = bestOfThree([&] {
                DirectSolver<double> s(kOrderings[o].method, Factorization::Cholesky,
                                       Traversal::LeftLooking);
                s.analyze(A);
            });
            std::printf("  %-6s %12zu %10.2f %10.2f\n", kOrderings[o].name, fill(sf), order,
                        analyze);
        }
        std::fflush(stdout);

        // The factorizations, at the first ordering, since the arm and the traversal are what
        // vary here and the ordering is priced above.
        const double aNorm = infNorm(A);
        const double ahNorm = infNorm(AH);
        const double asNorm = infNorm(AS);

        std::printf("  at %s: %-24s %8s %8s %8s %9s %9s %6s %6s\n", kOrderings[0].name,
                    "arm", "factLL", "factRL", "factMF", "bwd", "res", "pert", "delay");

        for (int a = 0; a < kNumArms; ++a) {
            Result byTraversal[kNumTraversals];
            bool   any = false;

            for (int t = 0; t < kNumTraversals; ++t) {
                switch (kArms[a].arm) {
                    case Arm::Real:
                        byTraversal[t] = measure(A, kOrderings[0].method, kArms[a].factorization,
                                                 kTraversals[t].traversal, aNorm);
                        break;
                    case Arm::Hermitian:
                        byTraversal[t] = measure(AH, kOrderings[0].method, kArms[a].factorization,
                                                 kTraversals[t].traversal, ahNorm);
                        break;
                    case Arm::ComplexSymmetric:
                        byTraversal[t] = measure(AS, kOrderings[0].method, kArms[a].factorization,
                                                 kTraversals[t].traversal, asNorm);
                        break;
                }
                if (byTraversal[t].ok)
                    any = true;
            }

            std::printf("         %-24s", kArms[a].name);
            for (int t = 0; t < kNumTraversals; ++t) {
                if (byTraversal[t].ok) std::printf(" %8.2f", byTraversal[t].factor);
                else                   std::printf(" %8s", "refused");
            }

            if (!any) {
                std::printf("\n");
                continue;
            }

            // The verdict columns come from the left-looking run, all three traversals computing
            // the same factorization of the same matrix.
            const Result& v = byTraversal[0];
            std::printf(" %9.1e %9.1e %6zu %6zu\n", v.backward, v.residual, v.perturbed, v.delayed);

            totalPerturbed += v.perturbed;
            totalDelayed += v.delayed;
        }
        std::fflush(stdout);
        std::printf("\n");
    }

    if (totalPerturbed == 0 && totalDelayed == 0)
        std::printf("no pivot was perturbed and no column was delayed, on any matrix, which is\n"
                    "what the dominance margin is for\n");
    else
        std::printf("PERTURBED %zu AND DELAYED %zu. The construction is wrong: these matrices are\n"
                    "meant to need neither.\n", totalPerturbed, totalDelayed);

    return 0;
}
