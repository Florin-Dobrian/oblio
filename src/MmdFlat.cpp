#include "oblio/MmdFlat.h"

#include "oblio/ElmOrder.h"
#include "oblio/MmdEngine.h"

#include <cstddef>
#include <cstdint>
#include <vector>

// MmdFlat.cpp - the free-function form of mmd over OUR OWN CLIQUE ARENA, a store that only appends,
// in elimination order, and never reclaims.
//
// AN ADAPTER, NOT AN IMPLEMENTATION. The ordering is `MmdEngine`, instantiated here on
// `QuotientGraphFlat`; this file exists so that a caller with a pattern and no engine to configure can
// ask for the order in one call, and so that `Ordering::MmdFlat` has a function to dispatch to.
//
// The store is what the name selects. This one is the unbounded member of the pair, and it is what
// makes the bounded one checkable: `MmdCompacted` runs the same engine over a pooled workspace and
// must return this file's permutation entry for entry.

namespace Oblio {

ElmOrder orderMmdFlat(const std::vector<std::size_t>&  colPtr,
                      const std::vector<std::int32_t>& rowIdx,
                      std::int32_t delta) {
    ElmOrder eo;
    MmdEngine<QuotientGraphFlat>(delta).compute(colPtr, rowIdx, eo);
    return eo;
}

} // namespace Oblio
