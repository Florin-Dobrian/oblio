// order_digest.cpp - one hash per driver per grid, compared against a recorded baseline.
//
// WHAT IT ANSWERS. "Did anything move?", asked of every driver at once. Almost every change to the
// ordering code is a RE-ENCODING rather than an algorithm change: an array folded into another, a
// fact read off a sign instead of a stamp, a pass reordered. Those must leave every permutation
// exactly as it was, and when the change is in `QuotientGraph` that means every driver's, not the
// one being worked on. This asks all of them in about two seconds.
//
// WHAT IT IS NOT, and this matters more than what it is. It is SELF-REFERENTIAL: it says "same as
// last time", not "correct". `make amdorder` and `make mmdorder` are stronger IN KIND, comparing
// our elimination order entry for entry against the vendored routines, which is an external
// oracle. This has breadth and speed instead: every driver rather than one, 73 grids rather than
// 38, no `private/` needed, two seconds rather than a build. Use it BETWEEN those runs, not
// instead of them.
//
// It checks the permutation and nothing else. Not fill, not memory, not the numeric path.
//
// THE BASELINE IS DELIBERATELY NOT COMMITTED, and the reason is a failure mode rather than a
// preference. A checked-in baseline catches drift over months, which is real value; but when it
// fails, the cheapest response is to regenerate it, and a baseline someone regenerates to silence
// a failure is worse than no baseline at all. This tree has already recorded that shape for other
// oracles. So: record at the start of a session, compare through it, and let the durable external
// check stay with `make amdorder` and `make mmdorder`, which cannot be quietly re-blessed.
//
//     make digest-record     writes .digest-baseline, gitignored
//     make digest            compares against it
//
// A driver added to the list below simply has no baseline entry until the next record, and the run
// says so rather than passing silently.

#include "oblio/AmdFlat.h"
#include "oblio/Amd3B.h"
#include "oblio/MmdFlat.h"
#include "oblio/Mmd3B.h"
#include "oblio/Mmd3C.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace Oblio;

namespace {

// The pattern both families are built from, five point or seven point, sorted within a column with
// the diagonal in place, which is what a SparseMatrix holds. Built without a per-vertex container:
// this runs 73 times per driver and the allocation would be most of it.
void grid(int side, bool cubic, std::vector<std::size_t>& colPtr, std::vector<std::int32_t>& rowIdx) {
    const int n = cubic ? side * side * side : side * side;
    colPtr.clear(); colPtr.reserve(n + 1); colPtr.push_back(0);
    rowIdx.clear(); rowIdx.reserve(7 * n);
    const auto id = [&](int x, int y, int z) { return (z * side + y) * side + x; };
    for (int j = 0; j < n; ++j) {
        if (cubic) {
            const int x = j % side, y = (j / side) % side, z = j / (side * side);
            if (z > 0)        rowIdx.push_back(id(x, y, z - 1));
            if (y > 0)        rowIdx.push_back(id(x, y - 1, z));
            if (x > 0)        rowIdx.push_back(id(x - 1, y, z));
            rowIdx.push_back(j);
            if (x + 1 < side) rowIdx.push_back(id(x + 1, y, z));
            if (y + 1 < side) rowIdx.push_back(id(x, y + 1, z));
            if (z + 1 < side) rowIdx.push_back(id(x, y, z + 1));
        } else {
            const int r = j / side, c = j % side;
            if (r > 0)        rowIdx.push_back(j - side);
            if (c > 0)        rowIdx.push_back(j - 1);
            rowIdx.push_back(j);
            if (c + 1 < side) rowIdx.push_back(j + 1);
            if (r + 1 < side) rowIdx.push_back(j + side);
        }
        colPtr.push_back(rowIdx.size());
    }
}

// FNV-1a over the permutation. Any hash would do; what matters is that a single moved entry
// changes it, and that the value is stable across machines and compilers so a baseline recorded
// on one is readable on another.
unsigned long long digest(const std::vector<std::int32_t>& perm) {
    unsigned long long h = 1469598103934665603ull;
    for (const std::int32_t v : perm) {
        h ^= static_cast<unsigned long long>(static_cast<std::uint32_t>(v));
        h *= 1099511628211ull;
    }
    return h;
}

using OrderFn = std::vector<std::int32_t> (*)(const std::vector<std::size_t>&,
                                              const std::vector<std::int32_t>&);

// EVERY MMD DRIVER TAKES A `delta` WITH A DEFAULT, so none of them has type OrderFn: a default
// argument is not part of a function's type and cannot be bound through a pointer. Forwarded
// rather than widening OrderFn, because delta is the mmd layers' business and the amd ones have no
// counterpart to it. Zero is what OrderEngine passes.
#define OBLIO_DIGEST_FORWARD(NAME, CALL)                                          \
    std::vector<std::int32_t> NAME(const std::vector<std::size_t>&  colPtr,       \
                                   const std::vector<std::int32_t>& rowIdx) {     \
        return CALL(colPtr, rowIdx);                                              \
    }
OBLIO_DIGEST_FORWARD(mmd3Default,  orderMmdFlat)
OBLIO_DIGEST_FORWARD(mmd3bDefault, orderMmd3B)
OBLIO_DIGEST_FORWARD(mmd3cDefault, orderMmd3C)
#undef OBLIO_DIGEST_FORWARD

struct Driver { const char* name; OrderFn fn; };

// Every driver reached as a free function. The vendored routines are NOT here: they need
// `private/`, and what this instrument is for is our own code not moving under a re-encoding.
// `make amdorder` and `make mmdorder` are where the vendored comparison lives.
const std::vector<Driver>& drivers() {
    static const std::vector<Driver> d = {
        {"MmdFlat",  mmd3Default}, {"MMD3B", mmd3bDefault}, {"MMD3C", mmd3cDefault},
        {"AmdFlat",  orderAmdFlat},   {"AMD3B", orderAmd3B},
    };
    return d;
}

// SMALL GRIDS, AND MANY OF THEM, which is the opposite of what the timing benchmark wants. A
// re-encoding that breaks anything breaks it on a small graph too, and usually on the smallest
// where the shape first appears; what catches it is COVERAGE of shapes, not size. 2 a side is
// worth having: a 2x2 grid is a clique, and an ordering that mishandles a fully connected
// component fails there and nowhere else.
const char* const kBaseline = ".digest-baseline";

} // namespace

int main(int argc, char** argv) {
    const bool record = argc > 1 && std::strcmp(argv[1], "record") == 0;

    std::map<std::string, unsigned long long> baseline;
    if (!record) {
        std::FILE* f = std::fopen(kBaseline, "r");
        if (!f) {
            std::printf("order_digest: no %s. Run `make digest-record` first.\n", kBaseline);
            return 2;
        }
        char key[64];
        unsigned long long value;
        while (std::fscanf(f, "%63s %llu", key, &value) == 2) baseline[key] = value;
        std::fclose(f);
    }

    std::FILE* out = record ? std::fopen(kBaseline, "w") : nullptr;
    if (record && !out) { std::printf("order_digest: cannot write %s\n", kBaseline); return 2; }

    std::vector<std::size_t> colPtr;
    std::vector<std::int32_t> rowIdx;
    int checked = 0, moved = 0, missing = 0;

    for (int cubic = 0; cubic < 2; ++cubic) {
        const int last = cubic ? 20 : 55;
        for (int side = 2; side <= last; ++side) {
            grid(side, cubic != 0, colPtr, rowIdx);
            for (const Driver& d : drivers()) {
                const unsigned long long got = digest(d.fn(colPtr, rowIdx));
                char key[64];
                std::snprintf(key, sizeof key, "%s/%s%d", d.name, cubic ? "cube" : "square", side);
                if (record) {
                    std::fprintf(out, "%s %llu\n", key, got);
                    ++checked;
                    continue;
                }
                const auto it = baseline.find(key);
                if (it == baseline.end()) {
                    std::printf("  NO BASELINE  %s\n", key);
                    ++missing;
                    continue;
                }
                ++checked;
                if (it->second != got) {
                    std::printf("  MOVED        %s\n", key);
                    ++moved;
                }
            }
        }
    }

    if (record) {
        std::fclose(out);
        std::printf("recorded %d digests over %zu drivers into %s\n",
                    checked, drivers().size(), kBaseline);
        return 0;
    }

    if (missing != 0)
        std::printf("\n%d entries had no baseline, which is what a NEW DRIVER looks like.\n"
                    "Re-record if that is expected.\n", missing);
    std::printf("%s: %d digests over %zu drivers\n",
                moved ? "PERMUTATIONS MOVED" : "all identical", checked, drivers().size());
    return moved ? 1 : 0;
}
