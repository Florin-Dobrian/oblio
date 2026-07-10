#pragma once

// ElmForest.h — the elimination forest of a factored sparse matrix.
//
// "Forest" (not tree) because a reducible/disconnected matrix yields several
// trees. Nodes are factor columns in the permuted (factor) order. This is the
// minimal core: the parent links (the elimination tree/etree). Child/sibling
// links, supernode structure, and statistics are layered on later.
//
//   parent[i] = parent node of i in the forest, or NIL if i is a root.
// By construction parent[i] > i (a node's parent is always a later factor column).
// Nodes are IDs -> std::int32_t; loop counters over them are std::size_t.

#include "oblio/Types.h"

#include <vector>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class ElmForestEngine;

class ElmForest {
public:
    ElmForest() = default;

    std::size_t size() const { return mParent.size(); }   // number of nodes

    const std::vector<std::int32_t>& parent() const { return mParent; }

private:
    std::vector<std::int32_t> mParent;   // parent[i] = parent of i, or NIL

    friend class ElmForestEngine;   // fills the forest; add engines as needed
};

} // namespace Oblio
