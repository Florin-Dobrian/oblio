#include "oblio/AmdFlat.h"

#include "oblio/AmdEngine.h"
#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

// AmdFlat.cpp - the free-function form of amd over OUR OWN CLIQUE ARENA, a store that only
// appends, in elimination order, and never reclaims.
//
// AN ADAPTER, NOT AN IMPLEMENTATION. The ordering is `AmdEngine`, instantiated here on
// `QuotientGraphFlat`; this file exists so that a caller with a pattern can ask for the order in
// one call, and so that `Ordering::AmdFlat` has a function to dispatch to.
//
// The store is what the name selects. This one is the unbounded member of the pair, and it is what
// makes the bounded one checkable: `AmdCompacted` runs the same engine over a pooled workspace and
// must return this file's permutation entry for entry. It is also the amd branch's reference
// against the vendored `AMD_2`, which `make amdorder` asserts over 38 cases.

namespace Oblio {

ElmOrder orderAmdFlat(const std::vector<std::size_t>&  colPtr,
                      const std::vector<std::int32_t>& rowIdx) {
    ElmOrder eo;
    AmdEngine<QuotientGraphFlat>().compute(colPtr, rowIdx, eo);
    return eo;
}

} // namespace Oblio
