// pipeline_timing.cpp - where a solve's time actually goes, per ordering and per traversal.
//
// `benchmarks/ordering` measures one phase against itself: how long each ordering takes and how
// much it fills. This measures the phases against each other, which is a different question and
// the one that says whether the first question matters. An ordering that saves a millisecond of
// analysis and costs two per factorization is a loss for anyone who factors twice.
//
//   make run2d                      the default square-grid sides
//   make run3d                      the same on CUBIC grids
//   ./pipeline_timing_cpp 200       any square sides
//   ./pipeline_timing_cpp 3d 26     any cubic sides
//
// The target names and the ladders match `../ordering` exactly, and neither folder has a bare
// `make run`: the family is always named, because a claim measured on one of them is a claim
// about that one.
//
// -O3 -DNDEBUG, best of three after a warm-up, as in the ordering folder and for the same reason.
//
// **Grid Laplacians only, and that is a real limitation rather than a starting point.** A grid is
// where the MMD fill gap is largest, because nearly every live vertex has the same degree and the
// tie-break decides almost every pick; `experiments/ordering` measured a 16 percent spread in fill
// across four filing orders of one algorithm. So a break-even computed here may be a property of
// grids rather than of the orderings, and widening the matrices is the first thing this folder
// needs. Cholesky, real, left-looking against right-looking against multifrontal.
//
// **BOTH FAMILIES, added 2026-08-10, and the second is not an afterthought.** Every table this
// folder produced before that date was square, and on square grids MMD beats AMD outright. The
// ordering benchmark has run cubic grids since 2026-08-09 and they invalidated several claims that
// had been stated generally, among them which branch fills less. A break-even is a ratio of two
// phases and both move with the family, so one measured on squares says nothing about cubes. The
// 3D ladder stops shorter because fill grows far faster: a 26 cube is n = 17576 against a 140
// square's 19600, and its factor is several times larger. Both ladders match `../ordering`'s so
// that rows can be read across the two folders without converting sizes.

// ASK FOR A PERFORMANCE CORE, on the one platform where cores differ. Apple Silicon runs a
// command-line process at QOS_CLASS_DEFAULT, which prefers a performance core but permits the
// scheduler to park the thread on an efficiency one, and that placement is STICKY over long
// stretches rather than jittering per iteration. A minimum over repeats filters per-sample noise
// completely and filters a whole run placed on the wrong core not at all.
//
// **AND THIS BENCHMARK IS MORE EXPOSED THAN THE ORDERING ONE**, which sizes its repeat count from
// a timed probe and takes the best of fifteen to thirty. `bestOfThree` below takes the best of
// three, so a placement lasting a few hundred milliseconds covers every sample of a row. In the
// ordering benchmark this call turned a scattered null into a clean result on 2026-08-10 and
// explained a 4 percent disagreement between identical binaries that had been recorded as
// irreducible since 2026-08-08.
//
// QOS_CLASS_USER_INTERACTIVE is the class the scheduler will not put on an efficiency core. The
// alternative is `taskpolicy -c` at the shell, which has to be remembered every run and silently
// does nothing when it is not.
#ifdef __APPLE__
#include <pthread.h>
#endif

#include "oblio/DirectSolver.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"
#include "oblio/Types.h"
#include "oblio/Vector.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

static SparseMatrix<double> grid2D(int m) {
    const int n = m * m;
    std::vector<std::vector<int>> adj(n);
    auto id = [&](int r, int c) { return r * m + c; };
    for (int r = 0; r < m; ++r)
        for (int c = 0; c < m; ++c) {
            const int u = id(r, c);
            if (r > 0)     adj[u].push_back(id(r - 1, c));
            if (c > 0)     adj[u].push_back(id(r, c - 1));
            if (r + 1 < m) adj[u].push_back(id(r + 1, c));
            if (c + 1 < m) adj[u].push_back(id(r, c + 1));
        }
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (int j = 0; j < n; ++j) {
        std::vector<int> col = adj[j];
        col.push_back(j);
        std::sort(col.begin(), col.end());
        for (int i : col) {
            rowIdx.push_back(static_cast<std::int32_t>(i));
            val.push_back(i == j ? 4.0 : -1.0);
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<double>(static_cast<std::size_t>(n), colPtr, rowIdx, val);
}

// The cubic Laplacian, the same shape and the same assembly. Its neighbor list is built in a
// deliberate order and then SORTED with the diagonal, which is not tidiness: SparseMatrix requires
// a column's rows ascending, and the order within a column is also a TIE-BREAK INPUT, since it
// decides the content order of a clique. A 3D builder written the natural way is not ascending,
// and getting that wrong reads as an ordering divergence rather than as a harness fault. It cost a
// day in experiments/ordering on 2026-08-09.
static SparseMatrix<double> grid3D(int m) {
    const int n = m * m * m;
    std::vector<std::vector<int>> adj(n);
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
    std::vector<std::size_t>  colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double>       val;
    for (int j = 0; j < n; ++j) {
        std::vector<int> col = adj[j];
        col.push_back(j);
        std::sort(col.begin(), col.end());
        for (int i : col) {
            rowIdx.push_back(static_cast<std::int32_t>(i));
            val.push_back(i == j ? 7.0 : -1.0);        // diagonally dominant at degree six
        }
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<double>(static_cast<std::size_t>(n), colPtr, rowIdx, val);
}

// Best of three after a warm-up. The work is a callable so every phase is timed the same way.
template <class Work>
static double bestOfThree(Work work) {
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

// nnz(L), from the symbolic factor rather than from any ordering's estimate: a supernode
// contributes its own triangle plus its update rows. Identical to the ordering folder's, and
// deliberately duplicated rather than shared, since a benchmark should stand alone.
static std::size_t fill(const SparseMatrix<double>& A, Ordering method) {
    const OrderEngine oe(method);
    Permutation       P;
    if (!oe.compute(A, P)) return 0;
    const ElmForestEngine fe;
    ElmForest             ef;
    if (!fe.compute(A, P, ef)) return 0;
    const SymFactorEngine se;
    SymFactor             sf;
    if (!se.compute(A, P, ef, sf)) return 0;

    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk), u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

struct Row {
    std::string name;
    double      order = 0;      // OrderEngine::compute alone
    double      analyze = 0;    // the whole facade call, left-looking
    double      analyzeMf = 0;  // the same, multifrontal: it sorts children and relabels
    double      factLl = 0;
    double      factRl = 0;
    double      factMf = 0;
    double      solve = 0;
    std::size_t nnzL = 0;
};

// One traversal's factor time, with the analysis already done. Reported separately because
// analyze depends on the traversal (multifrontal orders the children and relabels the supernodes)
// while factor depends on it far more.
static double factorTime(const SparseMatrix<double>& A, Ordering method, Traversal traversal,
                         bool& ok) {
    DirectSolver<double> solver(method, Factorization::Cholesky, traversal);
    ok = solver.analyze(A);
    if (!ok) return 0;
    return bestOfThree([&] { solver.factor(A); });
}

static Row measure(const SparseMatrix<double>& A, const std::string& name, Ordering method) {
    Row row;
    row.name = name;
    row.nnzL = fill(A, method);

    const OrderEngine oe(method);
    row.order = bestOfThree([&] { Permutation P; oe.compute(A, P); });

    row.analyze = bestOfThree([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::LeftLooking);
        s.analyze(A);
    });
    row.analyzeMf = bestOfThree([&] {
        DirectSolver<double> s(method, Factorization::Cholesky, Traversal::Multifrontal);
        s.analyze(A);
    });

    bool ok = false;
    row.factLl = factorTime(A, method, Traversal::LeftLooking, ok);
    row.factRl = factorTime(A, method, Traversal::RightLooking, ok);
    row.factMf = factorTime(A, method, Traversal::Multifrontal, ok);

    DirectSolver<double> solver(method, Factorization::Cholesky, Traversal::LeftLooking);
    if (solver.analyze(A) && solver.factor(A)) {
        Vector<double> b(A.size()), x(A.size());
        for (std::size_t i = 0; i < A.size(); ++i) b[i] = 1.0;
        row.solve = bestOfThree([&] { solver.solve(b, x); });
    }
    return row;
}

int main(int argc, char** argv) {
#ifdef __APPLE__
    // Before anything is timed. See the note beside the include.
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    // A leading mode word rather than a side, as in the ordering benchmark, so the two folders are
    // driven the same way.
    bool cubic     = false;
    int  firstSide = 1;
    if (argc > 1) {
        const std::string word = argv[1];
        if      (word == "3d") { cubic = true;  firstSide = 2; }
        else if (word == "2d") { cubic = false; firstSide = 2; }
    }

    // THE SAME LADDERS THE ORDERING BENCHMARK USES, deliberately: 32, 64, 100, 140 on squares and
    // 6, 12, 16, 20, 26 on cubes, so a row here can be read beside a row there without converting
    // sizes. The 3D ladder is shorter for the reason that folder gives, fill growing far faster,
    // and the 6 cube is in it because it is the largest size at which the vendored AMD and its raw
    // order disagree on fill.
    //
    // 26 CUBED IS THE EXPENSIVE ROW and it is Natural's doing: its factorization there runs into
    // seconds where every ordered one is well under a hundred milliseconds, and each is timed
    // three times in three traversals. Pass explicit sides to skip it.
    std::vector<int> sides = cubic ? std::vector<int>{6, 12, 16, 20, 26}
                                   : std::vector<int>{32, 64, 100, 140};
    if (argc > firstSide) {
        sides.clear();
        for (int k = firstSide; k < argc; ++k) sides.push_back(std::atoi(argv[k]));
    }

    const std::vector<std::pair<std::string, Ordering>> methods = {
        {"Natural", Ordering::Natural}, {"MMD", Ordering::MMD},
        {"MMD1", Ordering::MMD1},       {"MMD2", Ordering::MMD2},
        {"MMD3", Ordering::MMD3},
        {"AMD", Ordering::AMD},         {"AMD1", Ordering::AMD1},
        {"AMD1B", Ordering::AMD1B},     {"AMD2", Ordering::AMD2},
        {"AMD3", Ordering::AMD3},
        {"AMD2B", Ordering::AMD2B},
    };

    for (int side : sides) {
        const SparseMatrix<double> A = cubic ? grid3D(side) : grid2D(side);
        if (cubic)
            std::printf("\n=== grid %d^3, n = %d, Cholesky, real ===\n\n", side,
                        side * side * side);
        else
            std::printf("\n=== grid %dx%d, n = %d, Cholesky, real ===\n\n", side, side,
                        side * side);
        std::printf("%-8s %8s %9s %10s %10s %9s %9s %9s %8s\n", "ordering", "order", "analyze",
                    "analyzeMF", "nnz(L)", "factLL", "factRL", "factMF", "solve");

        std::vector<Row> rows;
        for (const auto& m : methods) {
            const Row r = measure(A, m.first, m.second);
            rows.push_back(r);
            std::printf("%-8s %8.2f %9.2f %10.2f %10zu %9.2f %9.2f %9.2f %8.2f\n", r.name.c_str(),
                        r.order, r.analyze, r.analyzeMf, r.nnzL, r.factLl, r.factRl, r.factMf,
                        r.solve);
        }

        // Break-even against the vendored AMD, which is the default and the one ours would have to
        // displace. An ordering that analyzes faster and factors slower pays only while the factor
        // count stays below the ratio, so this is the number a caller actually faces: analyze runs
        // once per pattern, factor runs per Newton step or per time step.
        const Row* ref = nullptr;
        for (const Row& r : rows)
            if (r.name == "AMD") ref = &r;
        if (ref == nullptr) continue;

        std::printf("\n%-8s %12s %12s %14s   %s\n", "against", "analyze", "factor LL",
                    "break-even", "reading");
        std::printf("%-8s %12s %12s %14s\n", "AMD", "saved ms", "cost ms", "factorizations");
        for (const Row& r : rows) {
            if (&r == ref) continue;
            const double saved = ref->analyze - r.analyze;   // positive: this ordering analyzes faster
            const double cost  = r.factLl - ref->factLl;     // positive: and factors slower
            std::printf("%-8s %12.2f %12.2f", r.name.c_str(), saved, cost);
            if (saved > 0 && cost > 0)
                std::printf(" %14.1f   pays below that many factorizations\n", saved / cost);
            else if (saved < 0 && cost < 0)
                std::printf(" %14.1f   pays above that many factorizations\n", saved / cost);
            else if (saved >= 0 && cost <= 0)
                std::printf(" %14s   wins on both, at every count\n", "always");
            else
                std::printf(" %14s   loses on both, at every count\n", "never");
        }
    }
    return 0;
}
