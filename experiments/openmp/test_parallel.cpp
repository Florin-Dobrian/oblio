// How much parallelism is already there, and how much is left to take.
//
// Two questions on one kernel, in that order of importance.
//
// First: at what front size does the BLAS stop threading a call on its own? One
// product, timed, with the library's thread cap on and off. Where the two agree, the
// cores are idle and only we can reach them; where they diverge, the BLAS is already
// using the machine and OpenMP would be redundant. That crossover is a direct design
// input for any future scheduler and nothing in the tree knows it today.
//
// Second, and only interesting below that crossover: what do two whole products run
// side by side buy over the same two run back to back? That is task parallelism, whole
// independent units rather than slices of one, and in the sparse-direct vocabulary it is
// tree-level parallelism: what factoring two disjoint forest branches looks like. See the parallelism section of
// docs/ARCHITECTURE.md and docs/DESIGN_DECISIONS.md (2026-07-22).
//
// Two kernels throughout, because a speedup alone cannot be read. The hand loop uses
// only per-core resources, so if two of them do not scale, OpenMP never gave us two
// threads and no other row on the page means anything. It is a control, not a
// competitor. Where the hand loop scales and the BLAS does not, under identical
// scheduling, the contention is in hardware the BLAS reaches for and we do not.

#include "Gemm.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

using Oblio::gemmBlas;
using Oblio::gemmHand;

using Kernel = void (*)(std::size_t, const double*, const double*, double*);

// Every configuration is warmed up once and then timed this many times, at every
// order. An earlier version tapered these with order to keep the sweep short, and the
// large rows came out incoherent: a single unwarmed call absorbed cold-page first
// touch and reported one product as slower than two. Cheap measurements are worse
// than no measurements, because they look like results.
const int kTrials = 3;

int gChecks = 0;
int gFailures = 0;

void check(bool ok, const char* what)
{
    ++gChecks;
    if (!ok) ++gFailures;
    std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what);
}

// The OpenMP runtime API does not degrade the way the pragmas do. A pragma is ignored
// without -fopenmp; a call to omp_get_max_threads() fails to link. So every runtime
// call sits behind _OPENMP, and that guard is what keeps the serial build buildable.
// Demonstrated rather than asserted: this file compiles both ways.
int maxThreads()
{
#if defined(_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
}

bool openmpEnabled()
{
#if defined(_OPENMP)
    return true;
#else
    return false;
#endif
}

void fillRandom(std::vector<double>& v, unsigned seed)
{
    std::mt19937 gen(seed);
    for (double& x : v) x = static_cast<double>(gen() % 2048) / 2048.0 - 0.5;
}

// One task's worth of work: reps products of order n. At small n this is the leaf
// swarm, many small fronts; at large n it is a single root front.
void runUnit(Kernel kernel, std::size_t n, std::size_t reps, const double* a,
             const double* b, double* c)
{
    for (std::size_t r = 0; r < reps; ++r) kernel(n, a, b, c);
}

using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point t0, Clock::time_point t1)
{
    return std::chrono::duration<double>(t1 - t0).count();
}

// wall is what the caller waited. taskA and taskB are what each unit took on its own,
// measured inside the parallel region. Two very different task times mean the pair
// landed on unlike cores, which on Apple Silicon means one performance and one
// efficiency core: the wall time is then set by the slower of the two and says
// nothing about contention. There is no way to pin threads on Darwin, so this is
// measured rather than prevented.
struct Sample {
    double wall = 0.0;
    double taskA = 0.0;
    double taskB = 0.0;

    double skew() const
    {
        const double lo = taskA < taskB ? taskA : taskB;
        const double hi = taskA < taskB ? taskB : taskA;
        return lo > 0.0 ? hi / lo : 1.0;
    }
};

struct Workspace {
    std::vector<double> a1, b1, c1;
    std::vector<double> a2, b2, c2;

    explicit Workspace(std::size_t n)
        : a1(n * n), b1(n * n), c1(n * n, 0.0), a2(n * n), b2(n * n), c2(n * n, 0.0)
    {
        fillRandom(a1, 11);
        fillRandom(b1, 22);
        fillRandom(a2, 33);
        fillRandom(b2, 44);
    }

    void zeroResults()
    {
        std::fill(c1.begin(), c1.end(), 0.0);
        std::fill(c2.begin(), c2.end(), 0.0);
    }
};

// One unit alone. This is the data-parallelism probe: whatever the BLAS does inside a
// single call, with no help and no competition from us.
Sample timeOne(Kernel kernel, std::size_t n, std::size_t reps, Workspace& w)
{
    w.zeroResults();
    Sample s;
    const auto t0 = Clock::now();
    runUnit(kernel, n, reps, w.a1.data(), w.b1.data(), w.c1.data());
    s.wall = seconds(t0, Clock::now());
    s.taskA = s.wall;
    s.taskB = s.wall;
    return s;
}

// The two units back to back.
Sample timeSerial(Kernel kernel, std::size_t n, std::size_t reps, Workspace& w)
{
    w.zeroResults();
    Sample s;
    const auto t0 = Clock::now();

    const auto a0 = Clock::now();
    runUnit(kernel, n, reps, w.a1.data(), w.b1.data(), w.c1.data());
    s.taskA = seconds(a0, Clock::now());

    const auto b0 = Clock::now();
    runUnit(kernel, n, reps, w.a2.data(), w.b2.data(), w.c2.data());
    s.taskB = seconds(b0, Clock::now());

    s.wall = seconds(t0, Clock::now());
    return s;
}

// The two units as two tasks. The parallel/single/task skeleton is the one a postorder
// walk of the elimination forest would use, which is why it is written with tasks
// rather than the simpler `sections`: with no dependency to express here, `sections`
// would read more cleanly and generalize to nothing.
//
// Without -fopenmp all four pragmas are ignored and this runs the two units back to
// back, giving a correct answer and a speedup of one. That is the property the
// recommendation in DESIGN_DECISIONS rests on.
Sample timeParallel(Kernel kernel, std::size_t n, std::size_t reps, Workspace& w)
{
    w.zeroResults();
    Sample s;
    const auto t0 = Clock::now();
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task shared(s, w)
            {
                const auto a0 = Clock::now();
                runUnit(kernel, n, reps, w.a1.data(), w.b1.data(), w.c1.data());
                s.taskA = seconds(a0, Clock::now());
            }
#pragma omp task shared(s, w)
            {
                const auto b0 = Clock::now();
                runUnit(kernel, n, reps, w.a2.data(), w.b2.data(), w.c2.data());
                s.taskB = seconds(b0, Clock::now());
            }
        }
    }
    s.wall = seconds(t0, Clock::now());
    return s;
}

using Timer = Sample (*)(Kernel, std::size_t, std::size_t, Workspace&);

// Warm up once, discarded, then keep the trial with the shortest wall time. The
// warm-up matters most at the large orders, where the buffers are big enough that
// first touch dominates a single call.
Sample best(Timer timer, Kernel kernel, std::size_t n, std::size_t reps, Workspace& w)
{
    timer(kernel, n, reps, w);
    Sample s = timer(kernel, n, reps, w);
    for (int i = 1; i < kTrials; ++i) {
        const Sample u = timer(kernel, n, reps, w);
        if (u.wall < s.wall) s = u;
    }
    return s;
}

double maxAbsDiff(const std::vector<double>& x, const std::vector<double>& y)
{
    double d = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double e = std::fabs(x[i] - y[i]);
        if (e > d) d = e;
    }
    return d;
}

// The arithmetic in one unit: reps products of order n, at 2 n^3 flops each.
double unitFlops(std::size_t n, std::size_t reps)
{
    const double d = static_cast<double>(n);
    return 2.0 * d * d * d * static_cast<double>(reps);
}

// The number of products per task, chosen so each task does roughly the same amount of
// arithmetic whatever the order.
std::size_t repsFor(std::size_t n)
{
    const double target = 2.0e8;
    const double per = 2.0 * static_cast<double>(n) * static_cast<double>(n) *
                       static_cast<double>(n);
    const double r = target / per;
    return r < 1.0 ? std::size_t{1} : static_cast<std::size_t>(r);
}

void reportEnvironment()
{
    std::printf("OpenMP compiled in    %s\n", openmpEnabled() ? "yes" : "no");
    std::printf("omp_get_max_threads   %d\n", maxThreads());

    // The BLAS thread cap is set from outside the process, so nothing here has to
    // guess at an Accelerate API. The Makefile runs this binary twice, once with the
    // cap and once without, and the GF/s column of the two runs is the comparison.
    const char* names[] = {"VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "OMP_NUM_THREADS"};
    for (const char* name : names) {
        const char* value = std::getenv(name);
        std::printf("%-21s %s\n", name, value ? value : "(unset)");
    }
    std::printf("\n");
}

void checkKernelsAgree()
{
    const std::size_t n = 96;
    Workspace w(n);

    std::vector<double> hand(n * n, 0.0);
    std::vector<double> blas(n * n, 0.0);
    gemmHand(n, w.a1.data(), w.b1.data(), hand.data());
    gemmBlas(n, w.a1.data(), w.b1.data(), blas.data());

    const double d = maxAbsDiff(hand, blas);
    std::printf("     hand against BLAS, max abs difference %.3e\n", d);
    check(d < 1e-12, "the two kernels compute the same product");
}

// The parallel run must reproduce the serial one exactly. Each task performs the same
// arithmetic in the same order as it would serially, and the tasks touch disjoint
// memory, so there is no reduction and no summation-order freedom: anything short of
// bit-identical is a data race, not rounding.
void checkParallelMatchesSerial(Kernel kernel, const char* name)
{
    const std::size_t n = 96;
    const std::size_t reps = 3;
    Workspace w(n);

    timeSerial(kernel, n, reps, w);
    const std::vector<double> serial1 = w.c1;
    const std::vector<double> serial2 = w.c2;

    timeParallel(kernel, n, reps, w);

    const bool same =
        maxAbsDiff(serial1, w.c1) == 0.0 && maxAbsDiff(serial2, w.c2) == 0.0;
    char what[128];
    std::snprintf(what, sizeof(what), "%s, parallel is bit-identical to serial", name);
    check(same, what);
}

// Both kernels sweep the same orders so the two tables can be read side by side, which
// is the whole point of having a control. The sweep stops at 512: the hand loop costs
// 130 ms a call at 1024, and a row present in one table and absent from the other is
// worse than a row missing from both.
void sweep(Kernel kernel, const char* name)
{
    const std::size_t sizes[] = {32, 64, 128, 256, 512};

    std::printf("%s\n", name);
    std::printf(
        "    n   reps    one call     GF/s   two serial  two parallel   tree   skew\n");
    for (std::size_t n : sizes) {
        const std::size_t reps = repsFor(n);
        Workspace w(n);

        const Sample one = best(timeOne, kernel, n, reps, w);
        const Sample ser = best(timeSerial, kernel, n, reps, w);
        const Sample par = best(timeParallel, kernel, n, reps, w);

        std::printf("%5zu %6zu %9.2f ms %8.2f %9.2f ms %10.2f ms %6.2fx %6.2f\n", n,
                    reps, one.wall * 1e3, unitFlops(n, reps) / one.wall / 1e9,
                    ser.wall * 1e3, par.wall * 1e3, ser.wall / par.wall, par.skew());
    }
    std::printf("\n");
}

}  // namespace

int main()
{
    std::printf("Dense products: what the BLAS already threads, and what is left\n\n");
    reportEnvironment();

    checkKernelsAgree();
    checkParallelMatchesSerial(gemmHand, "hand loop");
    checkParallelMatchesSerial(gemmBlas, "BLAS");
    std::printf("\n");

    sweep(gemmHand, "gemmHand, one core's own resources (the control)");
    sweep(gemmBlas, "gemmBlas, whatever the BLAS reaches for");

    std::printf("GF/s is the single-call rate. Compare that column between the capped\n");
    std::printf("and uncapped runs: where they agree, the BLAS is not threading and\n");
    std::printf("the idle cores are ours to take. tree is two serial / two parallel.\n");
    std::printf("skew is the slower task over the faster one within the parallel run:\n");
    std::printf("well above 1 means the pair landed on unlike cores, not contention.\n\n");

    std::printf("%d checks, %d failures\n", gChecks, gFailures);
    return gFailures == 0 ? 0 : 1;
}
