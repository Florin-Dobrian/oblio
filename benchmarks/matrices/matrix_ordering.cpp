// matrix_ordering.cpp -- what the ordering step alone costs on real matrices, genmmd against ours.
//
// ONE PHASE, TWO ROUTINES, AND NOTHING ELSE TIMED. matrix_performance.cpp already reports an
// ordering column, and it is a different measurement for a different purpose: it prices ordering
// against the rest of a solve, over four orderings and three traversals, on the 107 positive
// definite matrices a factorization can actually be run on. This file asks the narrower question
// the ordering work of 2026-08-15 raised, whether `MmdFlat` reproduces genmmd's SPEED on real
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

#include "oblio/AmdFlat.h"
#include "oblio/Amd3B.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/MmdFlat.h"
#include "oblio/Mmd3B.h"
#include "oblio/Mmd3C.h"
#include "oblio/Permutation.h"
#include "oblio/QuotientGraph.h"
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

// TWO ENTRY POINTS ON THE AMD SIDE, AND THEY MEASURE DIFFERENT THINGS.
//
//   amd_order      SuiteSparse as it ships. TIMED, because it is what a caller selecting AMD
//                  actually pays: `AMD_aat`, `AMD_1`, `AMD_2` and `AMD_postorder`, of which the
//                  first and last have no counterpart here.
//   amd_order_raw  the hooked copy, ../../tools/hook_amd.py, reporting the elimination order
//                  BEFORE the postorder. NOT timed: it carries a vector per vertex, so it is an
//                  oracle rather than a stopwatch.
//
// Everything but the clock therefore comes from the raw copy, because `AmdFlat` does no postorder
// and so returns the pre-postorder order. Comparing against `amd_order`'s output would compare two
// different permutations and report a fill difference that is the postorder's, not ours.
extern "C" int amd_order(int n, const int* colPtr, const int* rowIdx, int* perm,
                         double* control, double* info);
extern "C" int amd_order_raw(int n, const int* colPtr, const int* rowIdx, int* perm,
                             double* control, double* info);
// The raw copy hands its order back through this callback rather than through a parameter, the
// generator having no way to widen the signature it copies. Same arrangement as
// ../ordering/order_timing.cpp.
std::vector<int> gRawOrder;
extern "C" void pbRawOrder(const int* raw, int n) { gRawOrder.assign(raw, raw + n); }
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

enum class Branch { Mmd, Amd };

struct Options {
    Branch      branch   = Branch::Mmd;
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
    const bool   amd = options.branch == Branch::Amd;
    const double one = amd
        ? bestOf([&] { const auto p = orderAmdFlat(colPtr, rowIdx); (void) p; }, 1)
        : bestOf([&] { const auto p = orderMmdFlat(colPtr, rowIdx); (void) p; }, 1);
    const int repeats = options.repeats > 0 ? options.repeats
                                            : repeatsFor(one, options.targetMs);

    const double ours = amd
        ? bestOf([&] { const auto p = orderAmdFlat(colPtr, rowIdx); (void) p; }, repeats)
        : bestOf([&] { const auto p = orderMmdFlat(colPtr, rowIdx); (void) p; }, repeats);
    // TWO CLIQUE FIGURES, AND THEY ANSWER DIFFERENT QUESTIONS.
    //
    //   cC   CUMULATIVE, from the ordering's own arena-reporting overload: every entry the arena
    //        ever handed
    //        out. Nothing here is reclaimed, so it holds dead cliques AND the members that
    //        contractions dropped out of live ones. It is what the flat layout actually costs.
    //   pC   PEAK LIVE MEMBERS, read from `gPeakCliqueMembers`: the most that was alive at any one
    //        instant. It is what a CHUNKED clique store would ask the allocator for at its worst
    //        moment, and PAYLOAD ONLY, with no per-clique header, no allocator rounding and no
    //        capacity above size.
    //
    // pC IS A PROPERTY OF THE ALGORITHM, NOT OF THE LAYOUT. `MmdFlat`, `Mmd3B` and `Mmd3C` return
    // the same permutation, so they form the same cliques and lose the same members at the same
    // moments; the figure is identical for all three, and the same holds for `AmdFlat` against
    // `Amd3B`. Neither vendored routine can report one. Only cC is ours, and a chunked version
    // would show this same pC with a cC of nothing.
    std::size_t cC = 0;
    gPeakCliqueMembers = 0;
    const std::vector<std::int32_t> order = amd ? orderAmdFlat(colPtr, rowIdx, cC)
                                                : orderMmdFlat(colPtr, rowIdx, 0, cC);
    const std::size_t pC   = gPeakCliqueMembers;
    const std::size_t nnzL = fillOf(A, order);

    // THE B SIBLING RIDES ALONG ON BOTH BRANCHES, and only its time is reported: `Amd3B` on the amd
    // side, `Mmd3B` on the mmd side. Each is its branch's driver on a DIFFERENT CLIQUE LAYOUT,
    // `AMD_2`'s compacting pool and genmmd's chained segments respectively, so in neither case is
    // there anything for a `cC` column to compare; see experiments/ordering/README.md on the three
    // layouts. Everything else about a sibling MUST match its driver, the two being the same
    // algorithm with the same encodings, so the order, the fill and the peak live members are
    // CHECKED rather than printed and a mismatch is flagged at the end of the row.
    //
    // THE `pC` CHECK IS THE SHARP ONE, and it is why this exists. Peak live clique members is a
    // property of the ALGORITHM and not of the layout, so two drivers agreeing on the permutation
    // can still be caught doing different work. It found two defects in `Amd3B`'s mid-walk
    // collector on 2026-08-19, neither of which moved a permutation or a fill figure and neither
    // of which the digest, the vendored acceptance checks, `test_order` or the sanitizers could
    // see. Order and fill compare the ANSWER; a B layer exists to differ in MECHANISM while
    // agreeing on the answer, so until this check there was nothing watching the thing it varies.
    // TWO SIBLINGS ON THE MMD BRANCH, ONE ON THE AMD BRANCH, and every one of them gets all three
    // checks. `Mmd3C` is the mmd algorithm on `AMD_2`'s compacting pool, the second layout this
    // branch has, so it belongs beside `Mmd3B` rather than instead of it. There is no `Amd3C`.
    double      sibMs[2] = {0.0, 0.0};
    std::size_t nCmp[2]  = {0, 0};
    const int   sibCount = amd ? 1 : 2;
    std::string extra;
    for (int si = 0; si < sibCount; ++si) {
        const char* tag = amd ? "AMD3B" : (si == 0 ? "MMD3B" : "MMD3C");
        auto sibling = [&] {
            if (amd) return orderAmd3B(colPtr, rowIdx);
            return si == 0 ? orderMmd3B(colPtr, rowIdx) : orderMmd3C(colPtr, rowIdx);
        };
        sibMs[si] = bestOf([&] { const auto p = sibling(); (void) p; }, repeats);
        gPeakCliqueMembers = 0;
        const std::vector<std::int32_t> sibOrder = sibling();
        const std::size_t sibPC = gPeakCliqueMembers;
        // HOW OFTEN THE POOL HAD TO BE COMPACTED, for the two siblings that have one. Both size it
        // to `nzaat + nzaat/5 + n`, which is `AMD_2`'s `iwlen`, so the figure is directly against
        // that routine's Info[AMD_NCMPA] below. `Mmd3B` chains and never compacts, so it has none.
        if (amd)          nCmp[si] = gAmd3BCompactions;
        else if (si == 1) nCmp[si] = gMmd3CCompactions;
        if (sibOrder != order)                extra += std::string("  ") + tag + " order differs";
        else if (fillOf(A, sibOrder) != nnzL) extra += std::string("  ") + tag + " fill differs";
        if (sibPC != pC)                      extra += std::string("  ") + tag + " pC differs";
    }
    // Held as strings so a zero denominator, which a forest that failed to build would give, prints
    // as a dash rather than as a number nobody can read.
    char cTril[16] = "-", cFill[16] = "-", pOverC[16] = "-";
    if (tril != 0) std::snprintf(cTril, sizeof cTril, "%.2f", (double) cC / (double) tril);
    if (nnzL != 0) std::snprintf(cFill, sizeof cFill, "%.3f", (double) cC / (double) nnzL);
    if (cC   != 0) std::snprintf(pOverC, sizeof pOverC, "%.2f", (double) pC / (double) cC);

#ifdef OBLIO_VENDORED_ORDERINGS
    std::vector<int> ap, ai;
    toVendored(A, ap, ai);
    const int n = static_cast<int>(size);
    std::vector<int> perm(n), invp(n);
    double      vendored    = 0.0;
    std::size_t vendoredCmp = 0;
    if (amd) {
        // THE CLOCK IS ON `amd_order` WHOLE and everything else comes from the raw copy; see the
        // note beside the declarations. The raw copy is called once, outside the timed loop.
        vendored = bestOf(
            [&] { amd_order(n, ap.data(), ai.data(), perm.data(), nullptr, nullptr); }, repeats);
        // AND THE REFERENCE'S OWN COUNT, one untimed call with an Info array. Info[AMD_NCMPA] is
        // the same quantity `Amd3B` reports, over a workspace of the same size, so the two are
        // comparable entry for entry.
        {
            double info[20] = {0.0};
            amd_order(n, ap.data(), ai.data(), perm.data(), nullptr, info);
            vendoredCmp = static_cast<std::size_t>(info[8]);          // AMD_NCMPA
        }
        gRawOrder.clear();
        amd_order_raw(n, ap.data(), ai.data(), perm.data(), nullptr, nullptr);
    } else {
        vendored = bestOf(
            [&] { mmd_order(n, ap.data(), ai.data(), perm.data(), invp.data()); }, repeats);
    }

    // AND THE FILL COLUMN'S PREMISE IS CHECKED RATHER THAN ASSUMED, ON BOTH BRANCHES. It is printed
    // once and called both routines', which is only true while the two orders agree. Nothing on
    // this set was checking that until 2026-08-18, and when the amd half got the check it lit up
    // about sixty rows and cost four fixes, one of them a correctness bug in supervariable
    // detection that had been quietly costing fill. The mmd half is now checked the same way.
    //
    // WHAT IS COMPARED DIFFERS BY BRANCH, and the asymmetry is genmmd's rather than ours. `AMD_2`
    // relabels its output through `AMD_postorder`, so the comparable object has to be
    // reconstructed upstream of it by the hooked copy. genmmd does no postorder at all, so
    // `mmd_order`'s own `perm` IS the elimination order, 0-based new-to-old exactly as
    // `orderMmdFlat` returns it. Same reasoning as experiments/ordering's two checkers; see the
    // header of mmdorder.cpp there.
    const std::vector<int>& theirs = amd ? gRawOrder : perm;
    std::string note = extra;
    if (theirs.size() != order.size() ||
        !std::equal(theirs.begin(), theirs.end(), order.begin()))
        note += amd ? "  raw order differs" : "  order differs";

    // NO RATIO ON A READING TOO SMALL TO CARRY ONE. The set spans four orders of magnitude in
    // ordering time and the smallest matrices finish inside the clock's own resolution, where a
    // quotient of two such readings is noise printed to two decimals. Below a hundredth of a
    // millisecond the column is a dash, which says "not measured here" where `1.00x` would have
    // said "the same".
    // ONE ROW, PRINTED IN THREE PIECES so the branch-specific columns need no second format
    // string: the shared columns, then the amd branch's two, then whatever the checks had to say.
    char ratio[16] = "-";
    if (vendored >= 0.01 && ours >= 0.01)
        std::snprintf(ratio, sizeof ratio, "%.2fx", ours / vendored);
    std::printf("  %-38s %8zu %10zu %11zu %11zu %11zu %11zu %13zu %8s %8s %7s %9.3f %9.3f %8s",
                name.c_str(), size, m, tril, 2 * m, cC, pC, nnzL, cTril, cFill, pOverC,
                vendored, ours, ratio);
    for (int si = 0; si < sibCount; ++si) {
        char sibRatio[16] = "-";
        if (vendored >= 0.01 && sibMs[si] >= 0.01)
            std::snprintf(sibRatio, sizeof sibRatio, "%.2fx", sibMs[si] / vendored);
        std::printf(" %9.3f %9s", sibMs[si], sibRatio);
    }
    // The compaction counts last: the reference's, then the compacting sibling's. On the mmd side
    // the reference is genmmd, which chains and has no such figure, so only ours is printed.
    if (amd) std::printf(" %8zu %8zu", vendoredCmp, nCmp[0]);
    else     std::printf(" %8s %8zu", "-", nCmp[1]);
    std::printf("%s\n", note.c_str());
#else
    std::printf("  %-38s %8zu %10zu %11zu %11zu %11zu %11zu %13zu %8s %8s %7s %9s %9.3f %8s",
                name.c_str(), size, m, tril, 2 * m, cC, pC, nnzL, cTril, cFill, pOverC,
                "-", ours, "-");
    for (int si = 0; si < sibCount; ++si) std::printf(" %9.3f %9s", sibMs[si], "-");
    std::printf(" %8s %8zu", "-", amd ? nCmp[0] : nCmp[1]);
    std::printf("%s  vendored needs private/\n", extra.c_str());
#endif
}

void usage() {
    std::printf("the ordering step alone, a vendored routine against ours, on real matrices\n\n");
    std::printf("  ./matrix_ordering_cpp mmd|amd [--repeats=N] [--target-ms=X]\n");
    std::printf("                        [--max-n=N] [--max-nnz=N] <file.mtx> ...\n\n");
    std::printf("  mmd  genmmd against MmdFlat. They return the same permutation on every "
                "matrix\n");
    std::printf("       here, checked per matrix below, so the fill column is both of theirs\n");
    std::printf("       and the only difference between them is time.\n");
    std::printf("  amd  amd_order against AmdFlat. The CLOCK is on amd_order whole, what a "
                "caller\n");
    std::printf("       pays; EVERYTHING ELSE comes from the hooked pre-postorder copy, because\n");
    std::printf("       AmdFlat does no postorder and returns that order.\n\n");
    std::printf("  Each branch also times its B sibling, Mmd3B or Amd3B, and CHECKS its order,\n");
    std::printf("  fill and pC against the driver's rather than printing them. A mismatch is\n");
    std::printf("  flagged at the end of the row.\n");
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
        else if (arg == "mmd") options.branch = Branch::Mmd;
        else if (arg == "amd") options.branch = Branch::Amd;
        else
            paths.push_back(arg);
    }

    if (paths.empty()) {
        usage();
        std::printf("\nNo matrices given, so nothing was measured.\n");
        return 0;
    }

    std::sort(paths.begin(), paths.end());

    const bool amd = options.branch == Branch::Amd;
    if (amd) {
        std::printf("the ordering step alone, amd_order against AmdFlat, on real matrices\n");
        std::printf("  (best of N after a warm-up, through the same helper. THE CLOCK IS ON\n");
        std::printf("   amd_order WHOLE, which is what a caller pays and includes AMD_aat and\n");
        std::printf("   AMD_postorder, neither of which we do; nnz(L) is both routines', from\n");
        std::printf("   the hooked PRE-POSTORDER copy, AmdFlat returning it. AmdFlat/AMD "
                    "below 1 is\n");
        std::printf("   ours faster)\n\n");
    } else {
        std::printf("the ordering step alone, genmmd against MmdFlat, on real matrices\n");
        std::printf("  (best of N after a warm-up, both through the same helper; nnz(L) is both\n");
        std::printf("   routines', the two returning the same permutation. MmdFlat/MMD "
                    "below 1 is\n");
        std::printf("   ours faster)\n\n");
    }
    // THE SPACE COLUMNS COME FIRST because they are a property of the matrix and the ordering, not
    // of a run: `tril(A) = n + m` is A in its stored form, `A+I = 2m` is what mSource holds, `cC`
    // and `pC` are the clique figures described where they are computed, and `nnz(L)` includes the
    // diagonal. cC/tril(A) and cC/nnz(L) say whether our arena tracked the input or the factor on
    // this matrix. On grids it is about 2x tril(A) in 2D and up to 4.5x on cubes; the compression
    // is bought by mass elimination and so is a property of the MATRIX, which is exactly why a
    // real set is worth measuring. See docs/DESIGN_DECISIONS.md (2026-08-16).
    //
    // pC/cC IS THE ONE THAT DECIDES ANYTHING: the fraction of the flat arena a CHUNKED clique store
    // would hold at its worst moment, and so the saving that scheme would buy. On grids it is
    // about 0.43 in 2D and 0.28 on cubes.
    //
    // READ IT CONDITIONED ON cC/nnzL, or it will mislead. Where the arena is one per cent of the
    // factor, and this set has such matrices, clique storage is noise beside what the
    // factorization will need and no clique scheme is worth building for it whatever pC/cC says.
    // Where the arena is half the factor, and this set has those too, it decides.
    //
    // Each branch's siblings follow its own columns: `Mmd3B` and `Mmd3C` on the mmd side, `Amd3B`
    // on the amd side. Only their times are printed; their order, fill and pC are checked against
    // the branch driver's and a mismatch is flagged at the end of the row.
    std::printf("  %-38s %8s %10s %11s %11s %11s %11s %13s %8s %8s %7s %9s %9s %8s",
                "matrix", "n", "m", "tril(A)", "A+I", "cC", "pC", "nnz(L)", "cC/tril", "cC/nnzL",
                "pC/cC", amd ? "AMD ms" : "MMD ms", amd ? "AmdFlat ms" : "MmdFlat ms",
                amd ? "AmdFlat/AMD" : "MmdFlat/MMD");
    std::printf(" %9s %9s", amd ? "AMD3B ms" : "MMD3B ms",
                amd ? "AMD3B/AMD" : "MMD3B/MMD");
    if (!amd) std::printf(" %9s %9s", "MMD3C ms", "MMD3C/MMD");
    std::printf(" %8s %8s", amd ? "AMD cmp" : "-", amd ? "AMD3B cmp" : "MMD3C cmp");
    std::printf("\n");

    for (const std::string& path : paths) run(path, options);

    std::printf("\n%zu files\n", paths.size());
    return 0;                          // a table, not a verdict
}
