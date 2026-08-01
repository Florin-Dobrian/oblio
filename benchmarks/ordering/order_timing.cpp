// order_timing.cpp - ordering time and fill, per method, on grid Laplacians.
//
// It measures rather than asserts, and nothing in the tree depends on it.
//
//   make run                     the default grid sides
//   ./order_timing_cpp 200 280   any sides
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
static double orderTime(const SparseMatrix<double>& A, OrderMethod method) {
    const OrderEngine oe(method);
    Permutation warm;
    oe.compute(A, warm);

    double best = 1e30;
    for (int trial = 0; trial < 3; ++trial) {
        Permutation p;
        const auto t0 = std::chrono::steady_clock::now();
        oe.compute(A, p);
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

// nnz(L) under that ordering, from the symbolic factor: exact, and the honest way to compare two
// orderings, since a supernode contributes its own triangle plus its update rows.
static std::size_t fill(const SparseMatrix<double>& A, OrderMethod method) {
    const OrderEngine oe(method);
    Permutation p;
    if (!oe.compute(A, p)) return 0;
    const ElmForestEngine fe;
    ElmForest ef;
    if (!fe.compute(A, p, ef)) return 0;
    const SymFactorEngine se;
    SymFactor sf;
    if (!se.compute(A, p, ef, sf)) return 0;

    std::size_t nnz = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(sf.snodeSize()); ++kk) {
        const std::size_t f = sf.frontSize(kk), u = sf.updateSize(kk);
        nnz += f * (f + 1) / 2 + f * u;
    }
    return nnz;
}

int main(int argc, char** argv) {
    std::vector<int> sides = {32, 64, 100, 140};
    if (argc > 1) {
        sides.clear();
        for (int k = 1; k < argc; ++k) sides.push_back(std::atoi(argv[k]));
    }

    const std::vector<std::pair<std::string, OrderMethod>> methods = {
        {"MMD",  OrderMethod::MMD},  {"MMD1", OrderMethod::MMD1},
        {"MMD2", OrderMethod::MMD2},
        {"AMD",  OrderMethod::AMD},  {"AMD1", OrderMethod::AMD1},
        {"AMD2", OrderMethod::AMD2},
        {"AMD1B", OrderMethod::AMD1B},
        {"AMD2B", OrderMethod::AMD2B},
    };

    std::printf("%-12s %8s", "grid", "n");
    for (const auto& m : methods) std::printf(" %10s", (m.first + " ms").c_str());
    for (const auto& m : methods) std::printf(" %10s", (m.first + " nnzL").c_str());
    std::printf("\n");

    for (int side : sides) {
        const SparseMatrix<double> A = grid2D(side);
        std::printf("%4dx%-7d %8zu", side, side, A.size());
        for (const auto& m : methods) std::printf(" %10.2f", orderTime(A, m.second));
        for (const auto& m : methods) std::printf(" %10zu", fill(A, m.second));
        std::printf("\n");
        std::fflush(stdout);
    }
    return 0;
}
