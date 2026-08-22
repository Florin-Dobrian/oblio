// tools/sweep.cpp - a FINE SWEEP over grid side, looking for SELF-ALIASING POINTS in the drivers
// that ship: a size at which the ordering costs noticeably more per vertex than at the sizes either
// side of it, because arrays indexed by vertex have landed in the same cache sets.
//
// WHY. On 2026-08-17 `MmdCompacted` was found to read 1.28x `MmdFlat` at exactly 200 a side and
// about 0.99x at 199 and 201, and padding every size-n allocation by one page removed it: size-n arrays landing
// in the same cache sets at that n. `MmdCompacted` is transitional and its addresses are about to
// change, so the question that outlives it is whether `QuotientGraph` has such points, since `MmdFlat` and
// `AmdFlat` allocate the same shaped set of vectors and would collide by the same mechanism at
// whatever n aligns. The scaling ladder cannot see it: it visits twelve sides, and this one was
// found only because 200 happens to be a rung.
//
// HOW TO READ IT. A self-aliasing point is a SPIKE AT ONE SIDE with its neighbors flat, which is
// what makes it distinguishable from the smooth growth of the ordering itself. The sweep prints
// time per vertex, which removes that growth over a range this narrow, and flags any side more
// than a threshold above the median of its neighbors. Sides are consecutive, so nothing here
// depends on a ladder.
//
// A FLAGGED SIDE IS A CANDIDATE AND NOT A RESULT. Confirming one means padding the allocations by a
// page and seeing it go, exactly as `MmdCompacted`'s `OBLIO_MMD3C_PAD` did; nothing here proves a
// cause.
//
// REPEATS MATTER MORE THAN USUAL, and the default is 25 for a measured reason. At best-of-3 in the
// sandbox the series scattered by about 10 percent, which is the same size as the effect, and one
// side was flagged; at best-of-25 the scatter fell to about 5 percent and the flag went away. The
// threshold below is 10 percent, so anything under about 15 repeats produces candidates that are
// noise. If a real one is found it should hold across two invocations before anyone believes it.
//
//   c++ -std=c++17 -O3 -DNDEBUG -I../include sweep.cpp ../src/*.cpp -framework Accelerate -o sweep
//   ./sweep mmd3 190 210          one driver, a range of sides
//   ./sweep amd3 500 520 40       and more repeats where the series is noisy
//
// It links the whole library because the drivers do, and nothing in an ordering calls BLAS; on a
// machine without Accelerate, any BLAS will do, or the stub in tmp/blasstub.cpp.
//
// RUN IT WHERE THE QUESTION IS. This measures whether one side stands out from its neighbors ON THE
// MACHINE IT RUNS ON, and cache geometry differs between machines: the effect that prompted it does
// not reproduce in the Linux sandbox at all. A clean sweep elsewhere says nothing about alpamayo.
//
// WHAT IT FOUND, 2026-08-17, the run it was written for: NOTHING, which was the wanted answer.
// `MmdFlat` and `AmdFlat` are flat across sides 190 to 210 on alpamayo, and 200 in particular is
// unremarkable in both, `AmdFlat` reading BELOW both its neighbors there. So the point that
// prompted this belongs to `MmdCompacted`'s allocations and not to the shared class's. One band of 21 sides, so a
// point at some other n is unobserved rather than excluded.

#include "oblio/AmdFlat.h"
#include "oblio/MmdFlat.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

using namespace Oblio;

namespace {

// The 5-point Laplacian on an m by m grid, built without a container per vertex: the first version
// of tmp/one.cpp allocated twice per vertex and put 17 percent of a cachegrind profile in malloc.
void grid2d(int m, std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = m * m;
    colPtr.assign(1, 0);
    rowIdx.clear();
    rowIdx.reserve(static_cast<std::size_t>(n) * 4);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            const int v = i * m + j;
            if (i > 0)     rowIdx.push_back(v - m);
            if (j > 0)     rowIdx.push_back(v - 1);
            if (j < m - 1) rowIdx.push_back(v + 1);
            if (i < m - 1) rowIdx.push_back(v + m);
            colPtr.push_back(rowIdx.size());
            (void) v;
        }
    }
}

double bestMs(const std::string& which, const std::vector<std::size_t>& colPtr,
              const std::vector<std::int32_t>& rowIdx, int repeats) {
    double best = 1e300;
    for (int r = 0; r < repeats; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        const std::vector<std::int32_t> p = which == "amd3" ? orderAmdFlat(colPtr, rowIdx)
                                                            : orderMmdFlat(colPtr, rowIdx);
        const auto t1 = std::chrono::steady_clock::now();
        if (p.empty()) std::abort();                       // keeps the call from being elided
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

} // namespace

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "mmd3";
    const int lo      = argc > 2 ? std::atoi(argv[2]) : 190;
    const int hi      = argc > 3 ? std::atoi(argv[3]) : 210;
    const int repeats = argc > 4 ? std::atoi(argv[4]) : 25;

    if (which != "mmd3" && which != "amd3") {
        std::fprintf(stderr, "sweep: unknown driver \"%s\"; mmd3 or amd3\n", which.c_str());
        return 2;
    }
    if (lo < 2 || hi < lo) { std::fprintf(stderr, "sweep: bad range\n"); return 2; }

    std::printf("%s, sides %d to %d, best of %d\n\n", which.c_str(), lo, hi, repeats);
    std::printf("%6s %10s %10s %12s\n", "side", "n", "ms", "ns/vertex");

    std::vector<double> perVertex;
    std::vector<int>    sides;
    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    for (int m = lo; m <= hi; ++m) {
        grid2d(m, colPtr, rowIdx);
        const int    n  = m * m;
        const double ms = bestMs(which, colPtr, rowIdx, repeats);
        const double pv = ms * 1e6 / n;
        std::printf("%6d %10d %10.3f %12.2f\n", m, n, ms, pv);
        perVertex.push_back(pv);
        sides.push_back(m);
    }

    // A candidate is a side standing above BOTH neighbors by more than the threshold. Against the
    // neighbors rather than against a fit, because consecutive sides differ in n by under two
    // percent here and the ordering's own growth over that is far below what is being looked for.
    const double threshold = 1.10;
    std::printf("\ncandidates, more than %.0f percent above both neighbors:\n",
                (threshold - 1.0) * 100.0);
    int found = 0;
    for (std::size_t k = 1; k + 1 < perVertex.size(); ++k) {
        if (perVertex[k] > threshold * perVertex[k - 1] &&
            perVertex[k] > threshold * perVertex[k + 1]) {
            std::printf("  side %d: %.2f ns/vertex against %.2f and %.2f\n",
                        sides[k], perVertex[k], perVertex[k - 1], perVertex[k + 1]);
            ++found;
        }
    }
    if (found == 0) std::printf("  none\n");
    std::printf("\nA candidate is not a result. Confirm by padding the size-n allocations by one\n"
                "page and re-running: the spike must go and the permutation must not move.\n");
    return 0;
}
