// width.cpp - what the integer width costs, measured on the VENDORED routines rather than on ours.
//
// THE QUESTION. Our quotient graph holds about half its arrays as std::size_t where an
// std::int32_t would do: mSourcePtr, mAdjacencySize, mIncidenceSize, mCliquePtr, mCliqueSize,
// mWeight, and the drivers' degrees, outside and cliqueDegree. Every one is a count, one
// dimensional and bounded by n, not a position into an n x n object. Widening them was a
// deliberate decision recorded in docs/DESIGN_DECISIONS.md, taken for uniformity and for the
// reason that a count has no type of its own. The question this file answers is what that costs,
// and it cannot be answered by changing our own code: retyping ours would change the code AND the
// width at once, so a slowdown could be either.
//
// THE METHOD, and it is the useful part. Widen the ORACLE instead. The vendored routines are
// driven by a single `using Int = int32_t`, so retyping the kernel to int64_t gives a program
// doing byte-for-byte identical work with every array twice as wide and nothing else changed. The
// difference is the width and can be nothing else. Fill is asserted equal on every row, which is
// what says the work really was identical.
//
// The conversion at the API boundary is O(n + nnz) against a kernel that is not, so it does not
// move the measurement; it exists so the kernel runs on int64_t throughout rather than on a mix.
//
// WHY IT IS AN UPPER BOUND AND NOT AN ESTIMATE. This doubles every array. Ours are mixed: the
// pools and the mark array are already int32_t. So the number here bounds what a count sweep could
// recover rather than predicting it. It bounds it loosely from one side and tightly from the
// other: the arrays we hold wide are the SIZE arrays, read once per element rather than once per
// vertex, which are the ones a width penalty hits hardest.
//
// BOTH BRANCHES, IN ONE RUN. genmmd is measured here too, and that is the point rather than a
// convenience. The 17 to 26 percent recorded for genmmd in REPORT.md came from a different harness
// at a different time, so comparing the two branches across those two numbers would be measuring
// the harnesses as well. Same machine, same method, same run, or the comparison means nothing.
//
// genmmd has no `using Int` line to change: it is plain `int` throughout, all 145 of them the
// arithmetic type. So its transform is a blanket rename to a typedef instead, which is mechanical
// for that file and would not be for a file mixing widths.
//
// Run:  make width && ./width_cpp                 the default ladder
//       ./width_cpp 32 100 400                    chosen sides
//
// Needs private/. The Makefile target reports and skips when it is absent.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" int amd_order(std::int32_t n, const std::int32_t Ap[], const std::int32_t Ai[],
                         std::int32_t P[], double Control[], double Info[]);
extern "C" int amd_order64(std::int32_t n, const std::int32_t Ap[], const std::int32_t Ai[],
                           std::int32_t P[], double Control[], double Info[]);
void mmd_order(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);
extern "C" void mmd_order64(std::int32_t n, const std::int32_t colPtr[], const std::int32_t rowIdx[],
                            std::int32_t perm[], std::int32_t invp[]);

// A 2D grid Laplacian's pattern, off-diagonal only, which is what both routines want.
static void grid(int side, std::vector<std::int32_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = side * side;
    colPtr.assign(n + 1, 0);
    rowIdx.clear();
    for (int j = 0; j < n; ++j) {
        const int r = j / side, c = j % side;
        if (r > 0)        rowIdx.push_back(j - side);
        if (c > 0)        rowIdx.push_back(j - 1);
        if (c + 1 < side) rowIdx.push_back(j + 1);
        if (r + 1 < side) rowIdx.push_back(j + side);
        colPtr[j + 1] = static_cast<std::int32_t>(rowIdx.size());
    }
}

// The BEST of several runs rather than the mean. A mean measures the machine's other work as well
// as ours, and the question here is what the code costs when nothing is in its way.
//
// The untimed warm-up is not decoration. Without it the FIRST routine timed at the smallest size
// carried the cost of the core ramping up its clock, and reported 0.19 ms against the second
// routine's 0.07 on the same work: a 2.7x difference that is entirely the machine waking up. Taking
// the minimum does not remove it, because the ramp lasts longer than the whole loop at n = 1024.
// One discarded call per routine per row costs nothing and removes it.
//
// The repeat count is CHOSEN FROM A TIMED PROBE rather than fixed, so every row runs for about the
// same wall time whatever its size. A fixed count gave 25 repeats of a 0.17 ms ordering at
// n = 1024, which is four milliseconds of measurement and is not enough to see through the
// machine: two runs of that row disagreed by a factor of two and both reported the WIDER build as
// faster. The same fixed count gave seven repeats of a 29 ms ordering at n = 160000, which is the
// opposite problem and a slow benchmark. Targeting a duration fixes both.
template <class F>
static double best(F run, int repeats) {
    run();                             // warm-up, discarded: see below
    double lowest = 1e18;
    for (int k = 0; k < repeats; ++k) {
        const auto t0 = std::chrono::steady_clock::now();
        run();
        const auto t1 = std::chrono::steady_clock::now();
        lowest = std::min(lowest, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return lowest;
}

// How many repeats to reach roughly `targetMs` of measurement, from one timed call. The warm-up
// inside best() means the probe's own call is also the ramp, so it is discarded here too.
template <class F>
static int repeatsFor(F run, double targetMs) {
    run();
    const auto t0 = std::chrono::steady_clock::now();
    run();
    const auto t1 = std::chrono::steady_clock::now();
    const double one = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const int wanted = static_cast<int>(targetMs / std::max(one, 1e-3));
    return std::min(std::max(wanted, 20), 20000);
}

int main(int argc, char** argv) {
    std::vector<int> sides;
    for (int a = 1; a < argc; ++a) sides.push_back(std::atoi(argv[a]));
    if (sides.empty()) sides = {32, 64, 100, 140, 200, 280, 400};

    std::printf("%-11s %9s %11s %11s %10s %11s %11s %10s %7s\n", "grid", "n",
                "AMD int32", "AMD int64", "penalty", "MMD int32", "MMD int64", "penalty", "reps");
    for (int side : sides) {
        std::vector<std::int32_t> colPtr, rowIdx;
        grid(side, colPtr, rowIdx);
        const std::int32_t n = side * side;
        std::vector<std::int32_t> perm(n);

        // AMD_DENSE far above sqrt(n) drives the dense threshold into its MIN(n, dense) clamp, so
        // dense removal is unreachable and both builds do the same work. AMD_AGGRESSIVE at its
        // default. See the README's method section on why the alpha knob is the wrong one.
        double control[5] = {1e30, 1, 0, 0, 0};
        double info32[20], info64[20];

        // One count for the row, taken from the int32 AMD build, so both routines and both widths
        // are measured over the same number of orderings.
        const int repeats = repeatsFor([&] {
            amd_order(n, colPtr.data(), rowIdx.data(), perm.data(), control, info32);
        }, 400.0);
        const double t32 = best([&] {
            amd_order(n, colPtr.data(), rowIdx.data(), perm.data(), control, info32);
        }, repeats);
        const double t64 = best([&] {
            amd_order64(n, colPtr.data(), rowIdx.data(), perm.data(), control, info64);
        }, repeats);

        // The check that makes the row mean anything: identical fill is what says the two builds
        // did the same work and differ only in width.
        if (info32[9] != info64[9]) {
            std::printf("%4dx%-6d %9d   FILL DIFFERS, %.0f against %.0f: the two builds are not\n"
                        "                        doing the same work and the row means nothing\n",
                        side, side, n, info32[9], info64[9]);
            continue;
        }
        // genmmd, the same transform on the other branch. Its permutations are compared rather
        // than its fill, that being what this routine returns.
        std::vector<std::int32_t> permM(n), invpM(n), permM64(n), invpM64(n);
        const double m32 = best([&] {
            mmd_order(n, colPtr.data(), rowIdx.data(), permM.data(), invpM.data());
        }, repeats);
        const double m64 = best([&] {
            mmd_order64(n, colPtr.data(), rowIdx.data(), permM64.data(), invpM64.data());
        }, repeats);
        const bool mmdSame = (permM == permM64);

        std::printf("%4dx%-6d %9d %11.2f %11.2f %9.1f%% %11.2f %11.2f %9.1f%%  %6d%s\n",
                    side, side, n, t32, t64, 100.0 * (t64 - t32) / t32,
                    m32, m64, 100.0 * (m64 - m32) / m32, repeats,
                    mmdSame ? "" : "   MMD PERMUTATIONS DIFFER, the row means nothing");
        std::fflush(stdout);
    }
    return 0;
}
