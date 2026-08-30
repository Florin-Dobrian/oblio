#include "oblio/AmdCompacted.h"

#include "oblio/AmdEngine.h"
#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

// AmdCompacted.cpp - the free-function form of amd over a COMPACTED CLIQUE STORE: one pool with a
// free cursor and a compaction, where the flat form keeps an arena that only grows.
//
// AN ADAPTER, NOT AN IMPLEMENTATION. The ordering is `AmdEngine`, instantiated here on
// `QuotientGraphCompacted`; this file exists so that a caller with a pattern can ask for the order
// in one call, and so that `Ordering::AmdCompacted` has a function to dispatch to.
//
// IT IS THE PREDICTABLE-SPACE VERSION. Given a machine you know whether A fits, but you cannot know
// whether L fits, nnz(L) depending on the ordering being computed. A method that stays within
// `O(n + m)` carries a guarantee no amount of speed substitutes for: IF THE INPUT FITS, THE ANSWER
// IS REACHABLE. Its obligation is to return AmdFlat's permutation exactly, the two being one engine
// over two stores, and it is the only one of the pair that can report compactions.

namespace Oblio {

ElmOrder orderAmdCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx) {
    ElmOrder eo;
    AmdEngine<QuotientGraphCompacted>().compute(colPtr, rowIdx, eo);
    return eo;
}

} // namespace Oblio
