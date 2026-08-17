// matrix_ordering.cpp -- what the ordering step alone costs on real matrices, genmmd against ours.
//
// ONE PHASE, TWO ROUTINES, AND NOTHING ELSE TIMED. matrix_performance.cpp already reports an
// ordering column, and it is a different measurement for a different purpose: it prices ordering
// against the rest of a solve, over four orderings and three traversals, on the 107 positive
// definite matrices a factorization can actually be run on. This file asks the narrower question
// the ordering work of 2026-08-15 raised, whether `Mmd3` reproduces genmmd's SPEED on real
// structure as `mmdmatrices` showed it reproduces genmmd's PERMUTATION, and that question wants a
// different set: every matrix whose pattern can be read, positive definite or not, values or not.
// 246 files rather than 107.
//
// WHY IT IS A BENCHMARK AND `mmdmatrices` IS NOT. docs/CODING_RULES.md draws the line at whether
// something can FAIL. An alignment check is a verdict and lives in experiments/ordering beside its
// twin; a timing table is a table to read, and the three benchmark directories deliberately carry
// no `test` target for exactly that reason. This file prints numbers and returns 0 whatever they
// say.
//
// WHAT THE TWO COLUMNS MEAN, and the pairing matters. On every matrix `mmdmatrices` checks, the two
// routines return the SAME permutation, so they do the same work and produce the same factor and
// the only thing that separates them is how long they take. That makes this the cleanest timing
// comparison in the tree: no fill difference to trade against, no tie-break to argue about, one
// number each. The `nnz(L)` column is therefore printed once and is both of them.
//
// It is measured through the same path for both, which is a mistake this tree has already paid
// for: a free-function column timed by one helper against a standing method timed by another
// differed by up to 2.4 percent on 2026-08-10, and the difference was read as a result. Here both
// go through `bestOf` with the same warm-up, and neither builds a Permutation the other does not.
//
// `pattern` files are read, since an ordering needs no values; see MatrixMarket.h's banner note.
// That is most of why this set is twice the size of the performance one.

// ASK FOR A PERFORMANCE CORE, on the one platform where cores differ. Apple Silicon runs a
// command-line process at QOS_CLASS_DEFAULT, which prefers a performance core but permits the
// scheduler to park the thread on an efficiency one, and that placement is STICKY over long
// stretches rather than jittering per iteration, so a minimum over repeats does not filter it.
// The same call in ../ordering turned a scattered null into a clean result on 2026-08-10.
#ifdef __APPLE__
#include <pthread.h>
#endif

#include "MatrixMarket.h"

#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/Mmd3.h"
#include "oblio/Permutation.h"
#include "oblio/SparseMatrix.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace Oblio;

// The vendored routine, from ../../private/Mmd.cpp. Absent, this file still builds and the MMD
// column reports a refusal, which is how every driver in this folder treats private/.
#ifdef OBLIO_VENDORED_ORDERINGS
void mmd_order(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);
#endif

namespace {

// Best of N after a warm-up, the work a callable so both routines are timed the same way. Taken
// from matrix_performance.cpp, including the warm-up, which matters here for its reason: a matrix
// read from disk has just been touched by the reader and its pages are cold in a way a generated
// grid's are not.
template <class Work>
double bestOf(Work work, int repeats) {
    work();
    double best = 1e30;
    for (int trial = 0; trial < repeats; ++trial) {
        const auto t0 = std::chrono::steady_clock::now();
        work();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

// How many times to repeat, from one timed pass. A grid benchmark can fix this because its
// matrices span a known range; here the set spans n = 1000 to n = 1062400 and an ordering from
// under a millisecond to minutes, so a fixed count either wastes an hour on the large end or
// measures noise on the small one. `--repeats=` overrides it outright.
int repeatsFor(double oneMs, double targetMs) {
    if (oneMs <= 0.0) return 3;
    const int wanted = static_cast<int>(targetMs / oneMs);
    return std::min(std::max(wanted, 1), 200);
}

// nnz(L) with the diagonal. C is NOT taken from here: the ordering reports its own arena size, and
// this used to approximate it by summing update parts over supernodes, which overcounted by 12 to
// 18 percent in 2D and up to 40 on cubes. The two groupings differ, supernodes being found from L's
// structure and supervariables during the ordering, and a clique's entries are supervariable
// representatives where `updateSize` counts rows. The measured number needs no such argument.
std::size_t fillOf(const SparseMatrix<double>& matrix, const std::vector<std::int32_t>& order) {
    Permutation P(matrix.size());
    if (!P.setNewToOld(order)) return 0;
    const ElmForestEngine engine;
    ElmForest forest;
    if (!engine.compute(matrix, P, forest)) return 0;
    return forest.nnz();
}

// n and m, from the pattern rather than from a formula: a Matrix Market file may or may not carry
// every diagonal entry, so `m` is counted off the off-diagonal entries and halved, each edge
// appearing in both endpoints' lists. `tril(A) = n + m` is A in its stored form and `A+I = 2m` is
// what mSource holds.
std::size_t edgesOf(const SparseMatrix<double>& matrix) {
    std::size_t offDiagonal = 0;
    const std::vector<std::size_t>&  colPtr = matrix.colPtr();
    const std::vector<std::int32_t>& rowIdx = matrix.rowIdx();
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(matrix.size()); ++j)
        for (std::size_t p = colPtr[j]; p < colPtr[j + 1]; ++p)
            if (rowIdx[p] != j) ++offDiagonal;
    return offDiagonal / 2;
}

// The pattern without the diagonal, which is what genmmd takes. Ours takes the matrix's own CSC.
// Guarded because it has no other caller: without private/ there is no vendored routine to feed
// and the function is dead, which -Wunused-function says out loud.
#ifdef OBLIO_VENDORED_ORDERINGS
void toVendored(const SparseMatrix<double>& matrix,
                std::vector<int>& colPtr, std::vector<int>& rowIdx) {
    const std::vector<std::size_t>&  ap = matrix.colPtr();
    const std::vector<std::int32_t>& ai = matrix.rowIdx();
    const std::size_t size = matrix.size();

    colPtr.assign(size + 1, 0);
    rowIdx.clear();
    rowIdx.reserve(ai.size());
    for (std::size_t j = 0; j < size; ++j) {
        for (std::size_t p = ap[j]; p < ap[j + 1]; ++p)
            if (ai[p] != static_cast<std::int32_t>(j))
                rowIdx.push_back(ai[p]);
        colPtr[j + 1] = static_cast<int>(rowIdx.size());
    }
}
#endif

struct Options {
    int         repeats  = 0;          // 0 means choose from a timed pass
    double      targetMs = 200.0;      // how long one matrix's repeats should take
    std::size_t maxN     = 2000000;
    std::size_t maxNnz   = 50000000;
};

void run(const std::string& path, const Options& options) {
    const std::size_t slash = path.find_last_of('/');
    const std::size_t prior = (slash == std::string::npos)
                                  ? std::string::npos
                                  : path.find_last_of('/', slash - 1);
    const std::string name = (prior == std::string::npos) ? path : path.substr(prior + 1);

    const MatrixBenchmark::ReadResult read = MatrixBenchmark::readMatrixMarket(path, true);
    if (!read.ok) {
        std::printf("  %-38s %8s %11s %13s %9s %9s %8s  %s\n",
                    name.c_str(), "-", "-", "-", "-", "-", "-", read.reason.c_str());
        return;
    }

    const SparseMatrix<double>& A = read.matrix;
    const std::size_t size = A.size();
    const std::size_t nnz  = A.rowIdx().size();
    const std::size_t m    = edgesOf(A);
    const std::size_t tril = size + m;
    if (size == 0 || size > options.maxN || nnz > options.maxNnz) {
        std::printf("  %-38s %8zu %10zu %11zu %11zu %11s %13s %8s %8s %9s %9s %8s  past the cap\n",
                    name.c_str(), size, m, tril, 2 * m, "-", "-", "-", "-", "-", "-", "-");
        return;
    }

    const std::vector<std::size_t>&  colPtr = A.colPtr();
    const std::vector<std::int32_t>& rowIdx = A.rowIdx();

    // One pass to choose the repeat count, then the measurement. The pass is ours, and the same
    // count is given to both routines so neither gets a longer or shorter minimum than the other.
    const double one = bestOf([&] { const auto p = orderMmd3(colPtr, rowIdx); (void) p; }, 1);
    const int repeats = options.repeats > 0 ? options.repeats
                                            : repeatsFor(one, options.targetMs);

    const double ours = bestOf([&] { const auto p = orderMmd3(colPtr, rowIdx); (void) p; }, repeats);
    // C AS THE ORDERING ACTUALLY FILLED IT, not as the forest implies. `orderMmd3`'s four-argument
    // overload reports the arena's entry count, which is a size rather than a capacity and, this
    // arena never shrinking, also its peak.
    std::size_t arena = 0;
    const std::vector<std::int32_t> order = orderMmd3(colPtr, rowIdx, 0, arena);
    const std::size_t nnzL = fillOf(A, order);
    // Held as strings so a zero denominator, which a forest that failed to build would give, prints
    // as a dash rather than as a number nobody can read.
    char cTril[16] = "-", cFill[16] = "-";
    if (tril != 0) std::snprintf(cTril, sizeof cTril, "%.2f", (double) arena / (double) tril);
    if (nnzL != 0) std::snprintf(cFill, sizeof cFill, "%.3f", (double) arena / (double) nnzL);

#ifdef OBLIO_VENDORED_ORDERINGS
    std::vector<int> ap, ai;
    toVendored(A, ap, ai);
    const int n = static_cast<int>(size);
    std::vector<int> perm(n), invp(n);
    const double vendored = bestOf(
        [&] { mmd_order(n, ap.data(), ai.data(), perm.data(), invp.data()); }, repeats);

    // NO RATIO ON A READING TOO SMALL TO CARRY ONE. The set spans four orders of magnitude in
    // ordering time and the smallest matrices finish inside the clock's own resolution, where a
    // quotient of two such readings is noise printed to two decimals. Below a hundredth of a
    // millisecond the column is a dash, which says "not measured here" where `1.00x` would have
    // said "the same".
    if (vendored >= 0.01 && ours >= 0.01)
        std::printf("  %-38s %8zu %10zu %11zu %11zu %11zu %13zu %8s %8s %9.3f %9.3f %7.2fx\n",
                    name.c_str(), size, m, tril, 2 * m, arena, nnzL, cTril, cFill,
                    vendored, ours, ours / vendored);
    else
        std::printf("  %-38s %8zu %10zu %11zu %11zu %11zu %13zu %8s %8s %9.3f %9.3f %8s\n",
                    name.c_str(), size, m, tril, 2 * m, arena, nnzL, cTril, cFill,
                    vendored, ours, "-");
#else
    std::printf("  %-38s %8zu %10zu %11zu %11zu %11zu %13zu %8s %8s %9s %9.3f %8s  MMD needs private/\n",
                name.c_str(), size, m, tril, 2 * m, arena, nnzL, cTril, cFill,
                "-", ours, "-");
#endif
}

void usage() {
    std::printf("the ordering step alone, genmmd against Mmd3, on real matrices\n\n");
    std::printf("  ./matrix_ordering_cpp [--repeats=N] [--target-ms=X]\n");
    std::printf("                        [--max-n=N] [--max-nnz=N] <file.mtx> ...\n\n");
    std::printf("  The two return the same permutation on every matrix experiments/ordering's\n");
    std::printf("  `make mmdmatrices` checks, so the fill column is both of them and the only\n");
    std::printf("  difference between the routines is time.\n");
}

} // namespace

int main(int argc, char** argv) {
#ifdef __APPLE__
    // Before anything is timed. See the note beside the include.
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    Options options;
    std::vector<std::string> paths;

    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        if (arg.rfind("--repeats=", 0) == 0)
            options.repeats = std::atoi(arg.c_str() + 10);
        else if (arg.rfind("--target-ms=", 0) == 0)
            options.targetMs = std::atof(arg.c_str() + 12);
        else if (arg.rfind("--max-n=", 0) == 0)
            options.maxN = std::strtoull(arg.c_str() + 8, nullptr, 10);
        else if (arg.rfind("--max-nnz=", 0) == 0)
            options.maxNnz = std::strtoull(arg.c_str() + 10, nullptr, 10);
        else if (arg == "--help" || arg == "-h") { usage(); return 0; }
        else
            paths.push_back(arg);
    }

    if (paths.empty()) {
        usage();
        std::printf("\nNo matrices given, so nothing was measured.\n");
        return 0;
    }

    std::sort(paths.begin(), paths.end());

    std::printf("the ordering step alone, genmmd against Mmd3, on real matrices\n");
    std::printf("  (best of N after a warm-up, both through the same helper; nnz(L) is both\n");
    std::printf("   routines', the two returning the same permutation. MMD3/MMD below 1 is ours\n");
    std::printf("   faster)\n\n");
    // THE SPACE COLUMNS COME FIRST because they are a property of the matrix and the ordering, not
    // of a run: `tril(A) = n + m` is A in its stored form, `A+I = 2m` is what mSource holds, `C` is
    // the clique arena and `nnz(L)` includes the diagonal. C/tril(A) and C/nnz(L) are the two
    // ratios worth reading, and they say whether our arena tracked the input or the factor on this
    // matrix. On grids it is about 2x tril(A) in 2D and up to 4.5x on cubes; the compression is
    // bought by mass elimination and so is a property of the MATRIX, which is exactly why a real
    // set is worth measuring. See docs/DESIGN_DECISIONS.md (2026-08-16).
    //
    // The timing half is still genmmd against Mmd3 alone. The amd ladder belongs here too and is
    // not in yet.
    std::printf("  %-38s %8s %10s %11s %11s %11s %13s %8s %8s %9s %9s %8s\n",
                "matrix", "n", "m", "tril(A)", "A+I", "C", "nnz(L)", "C/tril", "C/nnzL",
                "MMD ms", "MMD3 ms", "MMD3/MMD");

    for (const std::string& path : paths) run(path, options);

    std::printf("\n%zu files\n", paths.size());
    return 0;                          // a table, not a verdict
}
