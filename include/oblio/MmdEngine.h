#pragma once

// MmdEngine.h - multiple minimum degree, over any of our clique stores.
//
// ONE ENGINE, TWO STORES. `MmdFlatEngine` runs over the append-only arena and `MmdCompactedEngine`
// over `AMD_2`'s pooled workspace with its free cursor and its compaction. They are the same
// algorithm and MUST return the same permutation, entry for entry, which is what makes the pair
// worth having: neither can be checked against the other unless both are reachable at once.
//
// THE STORE IS THE TEMPLATE PARAMETER, and the only one. Nothing above the store varies, so the
// body is written once and instantiated twice, both instantiations sitting in MmdEngine.cpp
// alongside the two graph classes, which is the same arrangement NumFactorEngine uses for its two
// factor storages.
//
// `delta` IS CONFIGURATION rather than an argument, as an engine's settings are: zero keeps a batch
// to true minima, negative takes one pivot per round. Callers that never change it never mention
// it.
//
// The ordering itself, the prepass, the filing convention, the clique-by-clique refresh and
// pairwise merging with outmatched marking, is described in MmdFlat.h.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// DECLARED, NOT DEFINED. Naming `MmdEngine<QuotientGraphFlat>` needs no more than this, the engine
// holding no member of the graph's type, so a caller that only wants an ordering does not pull in
// either store's bodies. The definitions meet the template in MmdEngine.cpp.
class QuotientGraphFlat;
class QuotientGraphCompacted;

template<class QuotientGraph>
class MmdEngine {
public:
    MmdEngine() = default;
    explicit MmdEngine(std::int32_t delta) : mDelta(delta) {}

    void         setDelta(std::int32_t delta) { mDelta = delta; }
    std::int32_t delta() const                { return mDelta; }

    // Order the pattern into `eo`.
    //
    // NO RETURN VALUE, unlike the other engines. Theirs is false where the job could not be done,
    // and this one has no such case: an ordering reads a pattern and always produces one. An empty
    // pattern gives an empty order, which is an answer rather than a failure.
    //
    // It takes a PATTERN rather than a matrix, which is all an ordering reads, and builds its own
    // quotient graph from it.
    void compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 ElmOrder& eo) const;

private:
    std::int32_t mDelta = 0;   // see the class comment
};

} // namespace Oblio
