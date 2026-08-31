#pragma once

// AmdEngine.h - approximate minimum degree, over any of our clique stores.
//
// ONE ENGINE, TWO STORES. Instantiated on `QuotientGraphFlat` it runs over the append-only arena,
// and on `QuotientGraphCompacted` over a pooled workspace with a free cursor and a compaction.
// They are the same algorithm and MUST return the same permutation, entry for entry,
// which is what makes the pair worth having: neither can be checked against the other unless both
// are reachable at once.
//
// THE STORE IS THE TEMPLATE PARAMETER, and the only one. Nothing above the store varies, so the
// body is written once and instantiated twice, both instantiations sitting in AmdEngine.cpp
// alongside the two graph classes, which is the same arrangement NumFactorEngine uses for its two
// factor storages.
//
// NO SETTINGS, which is why there is no constructor taking one. The mmd engine carries `delta`
// because its batching is a choice; this branch's one tunable, the dense-row threshold, is fixed
// rather than exposed, so the class holds nothing. `SymFactorEngine` is the same shape.
//
// The ordering itself, the prepass and the dense-row rule, the bound in place of an exact refresh,
// aggressive absorption, mass elimination and supervariable detection by hash, is described in
// AmdFlat.h.

#include "oblio/ElmOrder.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

// DECLARED, NOT DEFINED. Naming `AmdEngine<QuotientGraphFlat>` needs no more than this, the engine
// holding no member of the graph's type, so a caller that only wants an ordering does not pull in
// either store's bodies. The definitions meet the template in AmdEngine.cpp.
class QuotientGraphFlat;
class QuotientGraphCompacted;

template<class QuotientGraph>
class AmdEngine {
public:
    AmdEngine() = default;

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
};

} // namespace Oblio
