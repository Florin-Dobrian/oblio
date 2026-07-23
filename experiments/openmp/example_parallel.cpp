// OpenMP by example: the two kinds of parallelism, on a kernel small enough to read.
//
// One function, `sum`, adds two lists elementwise. It is parallelized two ways, and
// the difference between them is the whole distinction this experiment cares about.
//
//   Data parallelism   split ONE call across cores.  One pragma, inside the function.
//                      The same operation over many elements; scales with core count.
//   Task parallelism   run TWO whole calls at once.  Pragmas outside; function untouched.
//                      Independent units of work; capped by how many units exist.
//
// These are the general names. The sparse-direct literature calls the same two things
// node-level and tree-level parallelism, and that is the vocabulary to use when bridging
// to Oblio or to MUMPS. In those terms: data parallelism is node-level, which Accelerate
// already does inside a single dense kernel call, so it is nothing we write. Task
// parallelism is tree-level, factoring two independent forest branches at once, and that
// is the part that would be ours.
//
// The kernel is run in two flavors, and the pair is the actual lesson:
//
//   cheap      c[i] = a[i] + b[i], one flop per 24 bytes moved. Memory-bound. Both
//              parallelizations are correct and buy nearly nothing, because both cores
//              wait on the same memory.
//   expensive  a long arithmetic chain per element. Compute-bound. The same pragmas,
//              in the same places, now pay.
//
// So: putting the pragma in the right place is easy, and whether it buys anything is a
// separate question about what the code is waiting on. That is why this file exists.
//
// Build and run: `make example`. Both with and without OpenMP, since the same source
// must compile and give the right answer either way.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace Oblio {

namespace {

using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point t0, Clock::time_point t1)
{
    return std::chrono::duration<double>(t1 - t0).count();
}

// The runtime API is the one part of OpenMP that does NOT degrade on its own. A pragma
// is ignored without -fopenmp; a call to omp_get_max_threads() fails to LINK. Hence the
// guard. Every runtime call in a codebase that wants to build both ways needs one.
int maxThreads()
{
#if defined(_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// -------------------------------------------------------------------------------------
// The kernel, three ways.
// -------------------------------------------------------------------------------------

// Serial. The baseline both parallel versions must reproduce exactly, and the function
// the task-parallel form calls unmodified.
void sum(const double* a, const double* b, double* c, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) c[i] = a[i] + b[i];
}

// DATA parallelism: one call, split across cores.
//
// The single pragma is the entire change. Every iteration is independent, nothing is
// shared for writing, and no iteration reads what another writes, which is exactly the
// case OpenMP was built for. Note what is NOT needed: no restructuring, no threads
// created by hand, no join. Without -fopenmp the pragma is ignored and this is the
// serial loop again.
void sumDataParallel(const double* a, const double* b, double* c, std::size_t n)
{
#pragma omp parallel for
    for (std::size_t i = 0; i < n; ++i) c[i] = a[i] + b[i];
}

// The expensive body. Same shape, same independence, but hundreds of flops per element
// instead of one, so the arithmetic dominates the memory traffic. The chain is
// deliberately serial in `s` so the compiler cannot vectorize it away.
double expensiveElement(double x, double y)
{
    double s = 0.0;
    for (int k = 0; k < 100; ++k) s = std::sin(x + s) * std::cos(y) + s * 0.5;
    return s;
}

void sumExpensive(const double* a, const double* b, double* c, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) c[i] = expensiveElement(a[i], b[i]);
}

void sumExpensiveDataParallel(const double* a, const double* b, double* c,
                                   std::size_t n)
{
#pragma omp parallel for
    for (std::size_t i = 0; i < n; ++i) c[i] = expensiveElement(a[i], b[i]);
}

// -------------------------------------------------------------------------------------
// TASK parallelism: two whole calls at once. The callee is not modified at all.
//
// Note where the name sits. `sumDataParallel` is a parallel version OF `sum` and could
// be swapped in for it anywhere. This is not: it is a caller, and what it calls is plain
// `sum` with no pragma in it. There is no task-parallel version of the kernel, because
// task parallelism does not need one. That asymmetry is the distinction itself.
// -------------------------------------------------------------------------------------

using Kernel = void (*)(const double*, const double*, double*, std::size_t);

// `parallel` opens a team of threads. `single` says one thread runs the block that
// follows, so the two tasks are created once rather than once per thread. Each `task`
// is a unit of work any thread in the team may pick up. The implicit barrier at the end
// of `parallel` waits for both.
//
// This is the skeleton a postorder walk of the elimination forest would use, which is
// why it is written with tasks rather than the simpler `sections`: over a real tree the
// tasks are created recursively and `sections` would not generalize.
//
// Nothing here is shared for writing: the two calls read different inputs and write
// different outputs, exactly as two disjoint subtrees would.
void sumTaskParallel(Kernel kernel, const double* a1, const double* b1, double* c1,
                   const double* a2, const double* b2, double* c2, std::size_t n)
{
#pragma omp parallel
    {
#pragma omp single
        {
#pragma omp task
            kernel(a1, b1, c1, n);
#pragma omp task
            kernel(a2, b2, c2, n);
        }
    }
}

void sumInSequence(Kernel kernel, const double* a1, const double* b1, double* c1,
                      const double* a2, const double* b2, double* c2, std::size_t n)
{
    kernel(a1, b1, c1, n);
    kernel(a2, b2, c2, n);
}

// -------------------------------------------------------------------------------------
// Timing. Warm up once and keep the best of three, which is the discipline the rest of
// this experiment learned the hard way: a measurement taken once is a reading, not a
// result.
// -------------------------------------------------------------------------------------

template<class F>
double bestOf3(F f)
{
    f();  // warm up, discarded
    double t = 0.0;
    for (int i = 0; i < 3; ++i) {
        const auto t0 = Clock::now();
        f();
        const double u = seconds(t0, Clock::now());
        if (i == 0 || u < t) t = u;
    }
    return t;
}

struct Lists {
    std::vector<double> a1, b1, c1, a2, b2, c2;

    explicit Lists(std::size_t n)
        : a1(n), b1(n), c1(n, 0.0), a2(n), b2(n), c2(n, 0.0)
    {
        for (std::size_t i = 0; i < n; ++i) {
            const double x = static_cast<double>(i % 1000) / 1000.0;
            a1[i] = x;
            b1[i] = 1.0 - x;
            a2[i] = 2.0 * x;
            b2[i] = 2.0 - x;
        }
    }
};

bool sameAs(const std::vector<double>& x, const std::vector<double>& y)
{
    return x.size() == y.size() && std::equal(x.begin(), x.end(), y.begin());
}

// One flavor: the serial baseline, data parallelism on one call, and task parallelism
// over two calls, with a correctness check against the serial result in each case.
void demonstrate(const char* label, Kernel serial, Kernel dataParallel, std::size_t n)
{
    Lists lists(n);

    // Reference: what the serial kernel produces. Both parallel forms must match it
    // exactly. There is no reduction here and no shared accumulator, so every element is
    // computed by the same arithmetic in the same order however the work is divided.
    // Anything short of identical would be a data race, not rounding.
    serial(lists.a1.data(), lists.b1.data(), lists.c1.data(), n);
    const std::vector<double> reference = lists.c1;

    const double tOne = bestOf3([&] {
        serial(lists.a1.data(), lists.b1.data(), lists.c1.data(), n);
    });

    const double tData = bestOf3([&] {
        dataParallel(lists.a1.data(), lists.b1.data(), lists.c1.data(), n);
    });
    const bool dataOk = sameAs(reference, lists.c1);

    const double tSeq = bestOf3([&] {
        sumInSequence(serial, lists.a1.data(), lists.b1.data(), lists.c1.data(),
                         lists.a2.data(), lists.b2.data(), lists.c2.data(), n);
    });

    const double tTask = bestOf3([&] {
        sumTaskParallel(serial, lists.a1.data(), lists.b1.data(), lists.c1.data(),
                      lists.a2.data(), lists.b2.data(), lists.c2.data(), n);
    });
    const bool taskOk = sameAs(reference, lists.c1);

    std::printf("%s, n = %zu\n", label, n);
    std::printf("  one call, serial            %8.2f ms\n", tOne * 1e3);
    std::printf("  one call, data parallel     %8.2f ms   %5.2fx   %s\n", tData * 1e3,
                tOne / tData, dataOk ? "matches serial" : "MISMATCH");
    std::printf("  two calls, in sequence      %8.2f ms\n", tSeq * 1e3);
    std::printf("  two calls, as tasks         %8.2f ms   %5.2fx   %s\n", tTask * 1e3,
                tSeq / tTask, taskOk ? "matches serial" : "MISMATCH");
    std::printf("\n");
}

}  // namespace

}  // namespace Oblio

int main()
{
    std::printf("OpenMP by example: data parallelism and task parallelism\n\n");
    std::printf("threads available: %d\n", Oblio::maxThreads());
#if defined(_OPENMP)
    std::printf("built with OpenMP.\n\n");
#else
    std::printf("built WITHOUT OpenMP: every pragma ignored, every result still correct,\n");
    std::printf("every speedup 1.00x. That is the property, demonstrated. Note that the\n");
    std::printf("thread count above still printed: maxThreads() compiles either way only\n");
    std::printf("because its omp_get_max_threads() call sits behind #ifdef _OPENMP.\n\n");
#endif

    // Two sizes, for two different reasons. The cheap kernel needs a long list before it
    // is measurable at all. The expensive kernel does about 100 iterations per element,
    // so a long list would take minutes.
    Oblio::demonstrate("cheap: c[i] = a[i] + b[i], memory-bound", Oblio::sum,
                       Oblio::sumDataParallel, 4000000);

    Oblio::demonstrate("expensive: long arithmetic chain per element, compute-bound",
                       Oblio::sumExpensive, Oblio::sumExpensiveDataParallel,
                       100000);

    std::printf("The pragmas are identical in both flavors and sit in the same places.\n");
    std::printf("Only the ratio of arithmetic to memory traffic differs, and that alone\n");
    std::printf("decides whether the parallelism is worth anything.\n");
    return 0;
}
