// order_profile.cpp - one ordering method, run in a loop, for a sampling profiler.
//
// The benchmark next door answers "how long"; this exists to answer "where". An ordering runs in
// a few milliseconds, which is far too short for a sampling profiler to see, so this repeats it
// until there is something to sample.
//
// Build with symbols, and with inlining left on: -O2 matches what the library actually runs, and
// a profile of -O0 code would attribute time to functions the optimizer deletes. -g adds the
// symbols without changing the code.
//
//   make profile
//   ./order_profile_cpp amd1 140 200      method, grid side, repeats
//
// The method may be one of ours (mmd1, amd1, amd2) or a vendored one (mmd, amd). Profiling both
// sides through the same driver is the point: the difference between them is what is being
// investigated, and a comparison across two programs measures their differences too.
//
// Then either:
//
//   sample $(pgrep order_profile_cpp) 10 -f /tmp/order.sample   # while it runs, no setup
//   open /tmp/order.sample
//
// or, for the real thing including cache behavior:
//
//   xcrun xctrace record --template 'Time Profiler' --launch -- ./order_profile_cpp amd1 140 200
//   open *.trace
//
// The method argument matters: mmd1 and amd1 spend their time in completely different places, so
// profile the one whose gap is being investigated rather than an average of both.

#include "oblio/Amd1.h"
#include "oblio/Amd2.h"
#include "oblio/Mmd1.h"
#include "oblio/Mmd2.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

// A square grid Laplacian's pattern, full-symmetric with the diagonal present, which is what
// SparseMatrix holds and what the orderings take.
static void gridPattern(int side,
                        std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = side * side;
    std::vector<std::vector<int>> adjacency(n);
    for (int r = 0; r < side; ++r)
        for (int c = 0; c < side; ++c) {
            const int u = r * side + c;
            if (r > 0)        adjacency[u].push_back(u - side);
            if (c > 0)        adjacency[u].push_back(u - 1);
            if (c + 1 < side) adjacency[u].push_back(u + 1);
            if (r + 1 < side) adjacency[u].push_back(u + side);
        }
    colPtr.assign(n + 1, 0);
    rowIdx.clear();
    for (int j = 0; j < n; ++j) {
        std::vector<int> column = adjacency[j];
        column.push_back(j);
        std::sort(column.begin(), column.end());
        for (int i : column) rowIdx.push_back(static_cast<std::int32_t>(i));
        colPtr[j + 1] = rowIdx.size();
    }
}

int main(int argc, char** argv) {
    const std::string method  = (argc > 1) ? argv[1] : "amd1";
    const int         side    = (argc > 2) ? std::atoi(argv[2]) : 140;
    const int         repeats = (argc > 3) ? std::atoi(argv[3]) : 200;

    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    gridPattern(side, colPtr, rowIdx);

    // The vendored pair is reached through the engine, which needs a matrix rather than a
    // pattern. The values are never read by an ordering; they are here because SparseMatrix holds
    // them.
    const bool vendored = (method == "mmd" || method == "amd");
    std::vector<double> val(rowIdx.size(), 1.0);
    for (std::size_t aj = 0; aj + 1 < colPtr.size(); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] == static_cast<std::int32_t>(aj)) val[cp] = 100.0;
    const SparseMatrix<double> A(colPtr.size() - 1, colPtr, rowIdx, val);
    const OrderEngine engine(method == "mmd" ? OrderMethod::MMD : OrderMethod::AMD);

    // The sum is only there to stop the optimizer deleting the calls.
    std::size_t sum = 0;
    for (int k = 0; k < repeats; ++k) {
        if      (method == "mmd1") sum += orderMmd1(colPtr, rowIdx).size();
        else if (method == "amd1") sum += orderAmd1(colPtr, rowIdx).size();
        else if (method == "amd2") sum += orderAmd2(colPtr, rowIdx).size();
        else if (method == "mmd2") sum += orderMmd2(colPtr, rowIdx).size();
        else if (vendored) { Permutation p; engine.compute(A, p); sum += p.size(); }
    }
    std::printf("%s, grid %dx%d, %d repeats, %zu\n", method.c_str(), side, side, repeats, sum);
    return 0;
}
