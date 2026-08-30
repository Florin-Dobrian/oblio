#include "oblio/MmdCompacted.h"

#include "oblio/ElmOrder.h"
#include "oblio/MmdEngine.h"

#include <cstddef>
#include <cstdint>
#include <vector>

// MmdCompacted.cpp - the free-function form of mmd over a COMPACTED CLIQUE STORE: one pooled
// workspace with a free cursor and a compaction, where the flat form keeps an append-only arena.
//
// AN ADAPTER, NOT AN IMPLEMENTATION. The ordering is `MmdEngine`, instantiated here on
// `QuotientGraphCompacted`; this file exists so that a caller with a pattern and no engine to
// configure can ask for the order in one call, and so that `Ordering::MmdCompacted` has a function
// to dispatch to.
//
// Its obligation is to return MmdFlat's permutation exactly, the two being one engine over two
// stores. It is the only one of the pair that can report compactions.

namespace Oblio {

ElmOrder orderMmdCompacted(const std::vector<std::size_t>&  colPtr,
                           const std::vector<std::int32_t>& rowIdx,
                           std::int32_t delta) {
    ElmOrder eo;
    MmdEngine<QuotientGraphCompacted>(delta).compute(colPtr, rowIdx, eo);
    return eo;
}

} // namespace Oblio
