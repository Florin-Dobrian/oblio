#pragma once

// ElmForest.h — the elimination forest of a factored sparse matrix.
//
// "Forest" (not tree) because a reducible/disconnected matrix yields several
// trees. Nodes are factor columns in the permuted (factor) order. This is the
// minimal core: the parent links (the elimination tree/etree). Child/sibling
// links, supernode structure, and statistics are layered on later.
//
//   parent[i] = parent node of i in the forest, or cNoParent if i is a root.
// By construction parent[i] > i (a node's parent is always a later factor column).

#include <vector>
#include <cstddef>

namespace Oblio {

class ElmForestEngine;

class ElmForest {
public:
    // Sentinel for "no parent" (a tree root). (A permutation index is always < n,
    // so the max value is safe as an out-of-band marker.)
    static constexpr std::size_t cNoParent = static_cast<std::size_t>(-1);

    ElmForest() = default;

    std::size_t size() const { return mParent.size(); }   // number of nodes

    const std::vector<std::size_t>& parent() const { return mParent; }

private:
    std::vector<std::size_t> mParent;   // parent[i] = parent of i, or cNoParent

    friend class ElmForestEngine;   // fills the forest; add engines as needed
};

} // namespace Oblio
