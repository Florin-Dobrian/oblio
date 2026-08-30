// mmdmatrices.cpp -- production MmdFlat against MmdCorrected's elimination order, on REAL
// MATRICES.
//
// The same assertion mmdorder.cpp makes, on different input, and the difference in input is the
// whole point. mmdorder runs four shapes: seven hand-built examples, square grids, cubic grids and
// random patterns. Every one of them is generated here, and a generated matrix has the structure
// somebody chose to give it. Widening a grid exercises SCALE and never MECHANISM, which is a thing
// this tree has already been caught by once: the 2D-only version of the amd check was green while
// production AmdFlat carried a stale clique degree that a 3D grid at 16 a side finds in a second.
//
// Real matrices bring structures nothing here generates. Disconnected components. Isolated
// vertices. Rows adjacent to nearly everything, which is the dense-row pathology
// benchmarks/matrices/README.md records: one column of degree 10000 among 9992 of degree 5 takes
// MMD from 0.83 ms to 70.7. Supernodal structure that grids have only in the trivial form. And
// duplicate entries, which our orderings are on record as not tolerating and which the vendored
// AMD cleans in `amd_preprocess`.
//
// SO A DIVERGENCE HERE IS A FINDING, NOT A FAILURE OF THE EXERCISE. `docs/NEXT.md` item 6 records
// that AmdFlat and the vendored AMD already differ on fill on a minority of the 107-matrix
// performance set, once by 4 percent, and calls it "a divergence the acceptance tests cannot see".
// This is the check that would see it. The mmd branch is done first because it is the aligned one:
// MmdFlat reproduces the reference exactly on all 38 generated cases, so anything this finds is
// new.
//
// WHY IT LIVES HERE AND NOT IN benchmarks/matrices, WHERE THE MATRICES ARE. An alignment check is
// a VERDICT. docs/CODING_RULES.md states that `test` exists exactly where something can fail,
// which is why the three benchmark directories have none: a benchmark prints a table to read. A
// PASS or FAIL over a hundred matrices is not a table. And this file asserts the same property
// against the same oracle as mmdorder.cpp, so the two belong where a change to one obliges the
// other; split across directories, nothing would detect them drifting apart.
//
// WHAT CROSSES THE BOUNDARY IS THE DATA, AND IT CROSSES AS AN ARGUMENT. docs/WRITING_RULES.md
// warns about a Makefile in one directory naming another directory's files: nothing binds them,
// either may move, and the pointer dies without a sound. So no path is baked into a rule here.
// The matrices arrive on the command line exactly as they do for matrix_accuracy_cpp, the shell
// choosing the subset, and with none given this prints where to get them and exits clean. The one
// thing that does cross is an #include of benchmarks/matrices/MatrixMarket.h, which is a
// compile-time dependency: if that file moves, this stops building, loudly, which is the safe half
// of the same rule.
//
// THE READER CONDITIONS THE INPUT, and that is load bearing rather than incidental. Its
// `fromTriplets` deduplicates, sorts, mirrors the stored triangle and inserts a diagonal where
// nothing landed. Our orderings need all four and do none of them. What is compared here is
// therefore two orderings of the same CONDITIONED pattern, which is the honest comparison: it
// asks whether we reproduce the reference on real structure, not whether either of us survives a
// malformed file.
//
// `pattern` files are what this wants and what nothing else could use, so the reader takes them
// on request; see its banner note. They carry structure and no values, which is exactly an
// ordering's input.

#include "oblio/MmdFlat.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/Permutation.h"

#include "../../benchmarks/matrices/MatrixMarket.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// The reference routine, from private/MmdCorrected.cpp, linked rather than hooked: it returns the
// order it eliminates in and does no postorder, so its output vector IS the object to compare.
// mmdorder.cpp's header explains why the amd side needs a generated copy and this side does not.
void mmd_order_corrected(int n, const int colPtr[], const int rowIdx[], int perm[], int invp[]);

namespace {

// Two caps, and they are about PATIENCE rather than capability. Nothing here cannot order a large
// matrix: the reference's `int` arrays would refuse a pattern past 2^31 nonzeros, and the largest
// file in
// the collection as fetched is 28.7 million, seventy times under that. What a cap buys is that one
// matrix cannot quietly become the whole run, minimum degree paying for fill as it goes and
// PARSEC/Si87H76 predicting 5.68 billion entries of it under MmdFlat.
//
// THE DEFAULTS COVER THE COLLECTION AS FETCHED, which is the point: an alignment check should see
// everything on disk, and the three matrices these once excluded, PARSEC/Ga41As41H72,
// PARSEC/Si87H76 and Schenk/nlpkkt80, are exactly the ones with structure furthest from a grid.
// benchmarks/matrices caps them out of ITS runs for a different reason, that a factorization of
// them does not fit in memory, and that reason does not reach an ordering.
//
// Lower them to bound a run: `--max-n=` and `--max-nnz=` take any value, and a matrix past either
// is reported as skipped with its size rather than silently dropped.
constexpr std::size_t defaultMaxN   = 2000000;
constexpr std::size_t defaultMaxNnz = 50000000;

struct Options {
    std::size_t maxN   = defaultMaxN;
    std::size_t maxNnz = defaultMaxNnz;
};

// A conditioned SparseMatrix, two ways. Ours takes the pattern with the diagonal present, the
// reference
// takes it without. These are mmdorder.cpp's two conversions with a matrix in place of a Graph,
// and they stay in the driver for the same reason they do there: what a driver feeds its routine
// is the driver's own business.
void toOurs(const Oblio::SparseMatrix<double>& matrix,
            std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    colPtr = matrix.colPtr();
    rowIdx = matrix.rowIdx();
}

void toVendored(const Oblio::SparseMatrix<double>& matrix,
                std::vector<int>& colPtr, std::vector<int>& rowIdx) {
    const std::vector<std::size_t>&  ap = matrix.colPtr();
    const std::vector<std::int32_t>& ai = matrix.rowIdx();
    const std::size_t size = matrix.size();

    colPtr.assign(size + 1, 0);
    rowIdx.clear();
    rowIdx.reserve(ai.size());
    for (std::size_t j = 0; j < size; ++j) {
        for (std::size_t p = ap[j]; p < ap[j + 1]; ++p)
            if (ai[p] != static_cast<std::int32_t>(j))     // the reference wants no diagonal
                rowIdx.push_back(ai[p]);
        colPtr[j + 1] = static_cast<int>(rowIdx.size());
    }
}

// nnz(L) FROM THE FOREST ALONE, and that choice is load bearing rather than a shortcut.
// `ElmForest::nnz()` sums `f*(f+1)/2 + f*u` over the supernodes, which is the same formula
// benchmarks/ordering computes from the symbolic factor, and `test_pipeline` asserts the two agree
// on it. The forest reaches it from the front and update SIZES, so its own storage is O(n); the
// symbolic factor would materialize the index sets, which for PARSEC/Si87H76 at a predicted 5.68
// billion entries is about 22 GB. A fill figure must not cost more than the ordering it describes.
//
// Zero on refusal, which no matrix here should produce and which reads as "no forest" rather than
// "no fill".
std::size_t fillOf(const Oblio::SparseMatrix<double>& matrix,
                   const std::vector<std::int32_t>& order) {
    Oblio::Permutation P(matrix.size());
    if (!P.setNewToOld(order)) return 0;
    const Oblio::ElmForestEngine engine;
    Oblio::ElmForest forest;
    if (!engine.compute(matrix, P, forest)) return 0;
    return forest.nnz();
}

enum class Outcome { Matched, Differed, Skipped };

Outcome check(const std::string& path, const Options& options) {
    // The file name without its directory, since a full SuiteSparse path is long and the group
    // is the part worth seeing.
    const std::size_t slash = path.find_last_of('/');
    const std::size_t prior = (slash == std::string::npos)
                                  ? std::string::npos
                                  : path.find_last_of('/', slash - 1);
    const std::string name =
        (prior == std::string::npos) ? path : path.substr(prior + 1);

    const MatrixBenchmark::ReadResult read = MatrixBenchmark::readMatrixMarket(path, true);
    if (!read.ok) {
        std::printf("  %-38s %8s %11s %13s  skipped: %s\n",
                    name.c_str(), "-", "-", "-", read.reason.c_str());
        return Outcome::Skipped;
    }

    const std::size_t size = read.matrix.size();
    const std::size_t nnz  = read.matrix.rowIdx().size();
    if (size == 0) {
        std::printf("  %-38s %8s %11s %13s  skipped: empty\n", name.c_str(), "-", "-", "-");
        return Outcome::Skipped;
    }
    if (size > options.maxN || nnz > options.maxNnz) {
        std::printf("  %-38s %8zu %11zu %13s  skipped: past the cap\n",
                    name.c_str(), size, nnz, "-");
        return Outcome::Skipped;
    }

    std::vector<std::size_t>  colPtr;
    std::vector<std::int32_t> rowIdx;
    toOurs(read.matrix, colPtr, rowIdx);
    const std::vector<std::int32_t> ours = Oblio::orderMmdFlat(colPtr, rowIdx).order();

    std::vector<int> ap, ai;
    toVendored(read.matrix, ap, ai);
    const int n = static_cast<int>(size);
    std::vector<int> perm(n), invp(n);
    mmd_order_corrected(n, ap.data(), ai.data(), perm.data(), invp.data());

    if (ours.size() != perm.size()) {
        std::printf("  %-38s %8zu %11zu %13s  SIZE MISMATCH: ours %zu, reference %zu\n",
                    name.c_str(), size, nnz, "-", ours.size(), perm.size());
        return Outcome::Differed;
    }

    // BOTH FILLS COMPUTED, ONE PRINTED WHEN THEY MATCH. Matching permutations give equal fill
    // necessarily, so a second identical column on every row would be noise down the page; what
    // the pair is for is the divergence, where it says whether the disagreement COSTS anything.
    // A minimum degree ordering is a heuristic and two of them may legitimately differ while
    // filling the same, which is the distinction docs/NEXT.md item 6 turns on for the amd branch,
    // where ours fills LESS on the minority of matrices where the two disagree.
    //
    // Both are computed either way, and cheaply: the forest is O(n) and the fill comes off it, so
    // the second one also serves as a check that identical permutations really do produce
    // identical fill, which would be a defect below this file if they ever did not.
    const std::size_t fillOurs   = fillOf(read.matrix, ours);
    const std::size_t fillVendor = fillOf(read.matrix, std::vector<std::int32_t>(perm.begin(),
                                                                                perm.end()));

    std::size_t firstDiffer = ours.size();
    for (std::size_t k = 0; k < ours.size(); ++k)
        if (ours[k] != perm[k]) { firstDiffer = k; break; }

    if (firstDiffer != ours.size()) {
        std::printf("  %-38s %8zu %11zu   nnz(L) ours %zu reference %zu"
                    "  DIFFER at pivot %zu of %zu: ours %d, reference %d\n",
                    name.c_str(), size, nnz, fillOurs, fillVendor,
                    firstDiffer, ours.size(), ours[firstDiffer], perm[firstDiffer]);
        return Outcome::Differed;
    }
    std::printf("  %-38s %8zu %11zu %13zu  matches\n", name.c_str(), size, nnz, fillOurs);
    return Outcome::Matched;
}

void usage() {
    std::printf("production MmdFlat against MmdCorrected's elimination order, "
                "on real matrices\n\n");
    std::printf("  ./mmdmatrices_cpp [--max-n=N] [--max-nnz=N] <file.mtx> ...\n\n");
    std::printf("  The matrices are not in the repository. From benchmarks/matrices:\n");
    std::printf("      ./ssget.py list --max-nnz 500000 > candidates.txt\n");
    std::printf("      ./ssget.py fetch candidates.txt\n");
    std::printf("  They land in data/, which is gitignored. `make mmdmatrices` runs this over\n");
    std::printf("  all of them; any subset can be given directly, the shell choosing it.\n");
}

} // namespace

int main(int argc, char** argv) {
    Options options;
    std::vector<std::string> paths;

    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        if (arg.rfind("--max-n=", 0) == 0)
            options.maxN = std::strtoull(arg.c_str() + 8, nullptr, 10);
        else if (arg.rfind("--max-nnz=", 0) == 0)
            options.maxNnz = std::strtoull(arg.c_str() + 10, nullptr, 10);
        else if (arg == "--help" || arg == "-h") {
            usage();
            return 0;
        } else
            paths.push_back(arg);
    }

    if (paths.empty()) {
        usage();
        std::printf("\nNo matrices given, so nothing was checked.\n");
        return 0;                      // not a failure: an empty data/ is the ordinary state
    }

    std::sort(paths.begin(), paths.end());

    std::printf("production MmdFlat against MmdCorrected's elimination order, on real matrices\n");
    std::printf("  (the permutation itself, entry for entry; see mmdorder.cpp for why the mmd\n");
    std::printf("   side needs no hook, and this file's header for what real structure adds)\n\n");
    // A matching row carries one fill, the two being equal by construction; a differing row
    // carries both, which is the case the second one exists for.
    std::printf("  %-38s %8s %11s %13s\n", "matrix", "n", "nnz(A)", "nnz(L)");

    int matched = 0, differed = 0, skipped = 0;
    for (const std::string& path : paths) {
        switch (check(path, options)) {
            case Outcome::Matched:  ++matched;  break;
            case Outcome::Differed: ++differed; break;
            case Outcome::Skipped:  ++skipped;  break;
        }
    }

    std::printf("\n%zu files: %d matched, %d differed, %d skipped\n",
                paths.size(), matched, differed, skipped);
    return differed == 0 ? 0 : 1;
}
