#pragma once

// Permutation.h — holds a symmetric reordering of a sparse matrix as two index maps
// that are inverses of each other:
//   oldToNew[i] = new index of old index i
//   newToOld[i] = old index of new index i
// An ordering engine fills these (via friend access); other engines read them.
//
// A permutation has structure (a bijection of 0..n-1), so unlike a matrix it can hold
// well-defined but invalid data. The constructor therefore builds a *valid* state —
// the identity — never a sized-but-garbage one. (An ordering engine could still
// produce an invalid permutation; validate() checks for that.)

#include <vector>
#include <cstddef>

namespace Oblio {

class Permutation {
public:
    Permutation() = default;                 // empty (size 0) — trivially valid
    explicit Permutation(std::size_t size);  // identity of the given size

    std::size_t size() const { return mOldToNew.size(); }

    const std::vector<std::size_t>& oldToNew() const { return mOldToNew; }
    const std::vector<std::size_t>& newToOld() const { return mNewToOld; }

    // Reset to the identity permutation of the given size.
    void setIdentity(std::size_t size);

    // Rebuild newToOld as the inverse of oldToNew. An ordering engine fills oldToNew,
    // then calls this to complete the other direction consistently.
    void rebuildInverse();

    // Check that the maps form a genuine bijection of 0..size()-1 and are consistent
    // inverses. Returns true if valid. (Used to check ordering-engine output.)
    bool validate() const;

private:
    std::vector<std::size_t> mOldToNew;
    std::vector<std::size_t> mNewToOld;

    friend class OrderEngine;   // fills the maps; add other engine friends as needed
};

} // namespace Oblio
