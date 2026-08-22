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
//   ./order_profile_cpp amd3 140 200         method, grid side, repeats
//   ./order_profile_cpp amd3 3d 26 200       the same on a CUBIC grid, a different family
//   ./order_profile_cpp mmd3c 140 200        one of the three non-enum layers
//
// BOTH FAMILIES, and the cubic mode is not an afterthought. Every profile taken through this
// driver before 2026-08-09 was square, because square was all it could build, and the amd hash
// pass profiled as diffuse there while doing eight times the work per pivot on a cube. Profile the
// family the question is about.
//
// The method may be one of ours (mmd3, amd3), one of the three layers
// reached as free functions (mmd3b, mmd3c, amd3b), or a vendored
// one (mmd, amd). An unrecognized name is REFUSED rather than ignored: this driver used to fall
// through such a name silently, order nothing, and produce a profile of process startup that looks
// like a real trace with the ordering missing from it. That cost a profiling session.
// Profiling both
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
//   xcrun xctrace record --template 'Time Profiler' --launch -- ./order_profile_cpp amd3 140 200
//   open *.trace
//
// The method argument matters: mmd3 and amd3 spend their time in completely different places, so
// profile the one whose gap is being investigated rather than an average of both.

// The same performance-core request the timing benchmark makes; see the note there. A profile
// taken on an efficiency core attributes time correctly and reports the wrong amount of it.
#ifdef __APPLE__
#include <pthread.h>
#endif

#include "oblio/Amd3.h"
#include "oblio/Amd3B.h"
#include "oblio/Mmd3.h"
#include "oblio/Mmd3B.h"
#include "oblio/Mmd3C.h"
#include "oblio/OrderEngine.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

// Adjacency lists to full-symmetric CSC with the diagonal present, which is what SparseMatrix
// holds and what the orderings take. Shared by both builders below, which differ only in which
// neighbors they list, exactly as order_timing.cpp's pair does.
//
// A column's rows come out ASCENDING, and that is a precondition rather than a tidiness: it is
// what SparseMatrix states, and it is also a tie-break input, since the order within a column
// decides the content order of C[pivot]. A 3D builder written the natural way is not ascending,
// which cost a day on 2026-08-09.
static void fromAdjacency(const std::vector<std::vector<int>>& adjacency,
                          std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = static_cast<int>(adjacency.size());
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
    fromAdjacency(adjacency, colPtr, rowIdx);
}

// The seven-point cubic grid, and it is a DIFFERENT PROBLEM FAMILY rather than a larger size: it
// fills faster, makes larger cliques, and mass-eliminates far more often.
//
// EVERY PROFILE THIS FOLDER HAS EVER TAKEN WAS SQUARE, and that is why this is here. `REPORT.md`
// finding 1 says 2D flattered us and all our published numbers are 2D; the warning was carried
// across to the fill columns and not to the instrument. It cost something concrete: the amd hash
// pass profiled as diffuse at 140 a side, with no line above 378 ms, while it was testing eight
// times as many pairs per pivot on a cube. A profile of the flattering family is a profile of the
// flattering family.
static void grid3dPattern(int side,
                          std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = side * side * side;
    std::vector<std::vector<int>> adjacency(n);
    auto id = [&](int x, int y, int z) { return (z * side + y) * side + x; };
    for (int z = 0; z < side; ++z)
        for (int y = 0; y < side; ++y)
            for (int x = 0; x < side; ++x) {
                const int u = id(x, y, z);
                if (z > 0)        adjacency[u].push_back(id(x, y, z - 1));
                if (y > 0)        adjacency[u].push_back(id(x, y - 1, z));
                if (x > 0)        adjacency[u].push_back(id(x - 1, y, z));
                if (x + 1 < side) adjacency[u].push_back(id(x + 1, y, z));
                if (y + 1 < side) adjacency[u].push_back(id(x, y + 1, z));
                if (z + 1 < side) adjacency[u].push_back(id(x, y, z + 1));
            }
    fromAdjacency(adjacency, colPtr, rowIdx);
}

int main(int argc, char** argv) {
#ifdef __APPLE__
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    // An optional shape word between the method and the side, as order_timing.cpp takes:
    //     ./order_profile_cpp amd3 140 200        square, the default
    //     ./order_profile_cpp amd3 3d 26 200      cubic
    // The shape is named rather than implied, for the reason that folder's targets are: a default
    // invisible in the invocation is the one nobody revisits.
    std::vector<std::string> words;
    for (int k = 1; k < argc; ++k) words.emplace_back(argv[k]);

    bool cubic = false;
    for (auto it = words.begin(); it != words.end(); )
        if (*it == "3d")      { cubic = true;  it = words.erase(it); }
        else if (*it == "2d") { cubic = false; it = words.erase(it); }
        else                  ++it;

    const std::string method  = words.size() > 0 ? words[0] : std::string("amd3");
    const int         side    = words.size() > 1 ? std::atoi(words[1].c_str()) : (cubic ? 26 : 140);
    const int         repeats = words.size() > 2 ? std::atoi(words[2].c_str()) : 200;

    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    if (cubic) grid3dPattern(side, colPtr, rowIdx); else gridPattern(side, colPtr, rowIdx);

    // The vendored pair is reached through the engine, which needs a matrix rather than a
    // pattern. The values are never read by an ordering; they are here because SparseMatrix holds
    // them.
    const bool vendored = (method == "mmd" || method == "amd");
    std::vector<double> val(rowIdx.size(), 1.0);
    for (std::size_t aj = 0; aj + 1 < colPtr.size(); ++aj)
        for (std::size_t cp = colPtr[aj]; cp < colPtr[aj + 1]; ++cp)
            if (rowIdx[cp] == static_cast<std::int32_t>(aj)) val[cp] = 100.0;
    const SparseMatrix<double> A(colPtr.size() - 1, colPtr, rowIdx, val);
    const OrderEngine engine(method == "mmd" ? Ordering::MMD : Ordering::AMD);

    // Refuse an unknown method rather than falling through it. A profile of a method that never
    // ran is not empty, it is a trace of dyld and process startup, and it reads as a real one.
    if (!vendored && method != "mmd3" &&
        method != "amd3" &&
        method != "mmd3b" && method != "mmd3c" && method != "amd3b") {
        std::fprintf(stderr,
                     "order_profile: unknown method \"%s\"\n"
                     "  ours:     mmd3 amd3\n"
                     "  layers:   mmd3b mmd3c amd3b\n"
                     "  vendored: mmd amd\n", method.c_str());
        return 2;
    }

    // The sum is only there to stop the optimizer deleting the calls.
    std::size_t sum = 0;
    for (int k = 0; k < repeats; ++k) {
        if      (method == "mmd3") sum += orderMmd3(colPtr, rowIdx).size();
        else if (method == "amd3") sum += orderAmd3(colPtr, rowIdx).size();
        // THE THREE NON-ENUM LAYERS, added 2026-08-17. Each is its original computed differently,
        // so each is a thing whose cost is a question in its own right: Mmd3B and Amd3B carry the
        // vendored clique storage schemes and Mmd3C carries the port of the amd array folds onto
        // the mmd side. Profiling them through THIS driver rather than another is the same
        // argument the header already makes for the vendored pair: a comparison across two
        // programs measures their differences too.
        else if (method == "mmd3b") sum += orderMmd3B(colPtr, rowIdx).size();
        else if (method == "mmd3c") sum += orderMmd3C(colPtr, rowIdx).size();
        else if (method == "amd3b") sum += orderAmd3B(colPtr, rowIdx).size();
        else if (vendored) { Permutation P; engine.compute(A, P); sum += P.size(); }
    }
    if (cubic)
        std::printf("%s, grid3d %d^3 (n = %d), %d repeats, %zu\n",
                    method.c_str(), side, side * side * side, repeats, sum);
    else
        std::printf("%s, grid %dx%d (n = %d), %d repeats, %zu\n",
                    method.c_str(), side, side, side * side, repeats, sum);
    return 0;
}
