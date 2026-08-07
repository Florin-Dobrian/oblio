// order_timing.cpp - ordering time and fill, per method, on grid Laplacians.
//
// It measures rather than asserts, and nothing in the tree depends on it.
//
//   make run                          the default grid sides, every method
//   make scale                        a wider ladder, the MMD branch only, with gap columns
//   ./order_timing_cpp 200 280        any sides, every method
//   ./order_timing_cpp mmd 200 280    any sides, the MMD branch, with gap columns
//
// The mmd mode exists because the question it answers is not the same one. `run` asks what each
// method costs; `scale` asks how MMD1's and MMD2's distance from the vendored routine MOVES with
// n, which is a ratio and is tedious to take off the wide table by hand. Same measurements either
// way, and the ratios are against MMD in the same row, so nothing is compared across machines.
//
// -O3 -DNDEBUG matters, and the Makefile carries it: the claim is about what an ordering costs,
// and at -O0 with asserts in it would measure something else. Every method is run three times at
// each size and the shortest wall time is kept, one warm-up run having gone first, since a single
// cold reading is a reading and not a result.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/SymFactor.h"
#include "oblio/SymFactorEngine.h"

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
    std::vector<std::size_t> colPtr(n + 1, 0);
    std::vector<std::int32_t> rowIdx;
    std::vector<double> val;
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

// Ordering time alone, in milliseconds: the best of three, after one warm-up.
static double orderTime(const SparseMatrix<double>& A, Ordering method) {
    const OrderEngine oe(method);
    Permutation warm;
    oe.compute(A, warm);

    double best = 1e30;
    for (int trial = 0; trial < 3; ++trial) {
        Permutation P;
        const auto t0 = std::chrono::steady_clock::now();
        oe.compute(A, P);
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

// nnz(L) under that ordering, from the symbolic factor: exact, and the honest way to compare two
// orderings, since a supernode contributes its own triangle plus its update rows.
static std::size_t fill(const SparseMatrix<double>& A, Ordering method) {
    const OrderEngine oe(method);
    Permutation P;
    if (!oe.compute(A, P)) return 0;
    const ElmForestEngine fe;
    ElmForest ef;
    if (!fe.compute(A, P, ef)) return 0;
    const SymFactorEngine se;
    SymFactor sf;
    if (!se.compute(A, P, ef, sf)) return 0;

    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk), u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

int main(int argc, char** argv) {
    // First argument may be a mode word rather than a side. Only one exists, and anything else is
    // read as a side, so the old spelling keeps working unchanged.
    int firstSide = 1;
    bool mmdOnly = false;
    if (argc > 1 && std::string(argv[1]) == "mmd") { mmdOnly = true; firstSide = 2; }

    std::vector<int> sides = mmdOnly ? std::vector<int>{32, 64, 100, 140, 200, 280, 400}
                                     : std::vector<int>{32, 64, 100, 140};
    if (argc > firstSide) {
        sides.clear();
        for (int k = firstSide; k < argc; ++k) sides.push_back(std::atoi(argv[k]));
    }

    const std::vector<std::pair<std::string, Ordering>> allMethods = {
        {"MMD",  Ordering::MMD},  {"MMD1", Ordering::MMD1},
        {"MMD2", Ordering::MMD2}, {"MMD3", Ordering::MMD3},
        {"AMD",  Ordering::AMD},  {"AMD1", Ordering::AMD1},
        {"AMD2", Ordering::AMD2},
        {"AMD1B", Ordering::AMD1B},
        {"AMD2B", Ordering::AMD2B},
    };
    const std::vector<std::pair<std::string, Ordering>> mmdMethods = {
        {"MMD",  Ordering::MMD},  {"MMD1", Ordering::MMD1},
        {"MMD2", Ordering::MMD2}, {"MMD3", Ordering::MMD3},
    };
    const auto& methods = mmdOnly ? mmdMethods : allMethods;

    std::printf("%-12s %8s", "grid", "n");
    for (const auto& m : methods) std::printf(" %10s", (m.first + " ms").c_str());
    for (const auto& m : methods) std::printf(" %10s", (m.first + " nnzL").c_str());
    if (mmdOnly)
        std::printf(" %10s %10s %10s %10s %10s %10s", "MMD1 time", "MMD1 fill",
                    "MMD2 time", "MMD2 fill", "MMD3 time", "MMD3 fill");
    std::printf("\n");

    for (int side : sides) {
        const SparseMatrix<double> A = grid2D(side);
        std::printf("%4dx%-7d %8zu", side, side, A.size());
        // Held rather than printed straight through, since the gap columns need them again. The
        // baseline is MMD in THIS row: a ratio against another machine's number would be a
        // different measurement wearing the same name.
        std::vector<double> times;
        std::vector<std::size_t> fills;
        for (const auto& m : methods) times.push_back(orderTime(A, m.second));
        for (const auto& m : methods) fills.push_back(fill(A, m.second));
        for (double t : times)       std::printf(" %10.2f", t);
        for (std::size_t f : fills)  std::printf(" %10zu", f);
        if (mmdOnly) {
            const double baseTime = times[0];
            const double baseFill = static_cast<double>(fills[0]);
            for (std::size_t k = 1; k < methods.size(); ++k) {
                const double f = static_cast<double>(fills[k]);
                std::printf(" %9.2fx", times[k] / baseTime);
                std::printf(" %9.1f%%", 100.0 * (f - baseFill) / baseFill);
            }
        }
        std::printf("\n");
        std::fflush(stdout);
    }
    return 0;
}
