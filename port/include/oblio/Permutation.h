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
#include <cstdint>

namespace Oblio {

class Permutation {
public:
    Permutation() = default;                 // empty (size 0) — trivially valid
    explicit Permutation(std::size_t size);  // identity of the given size

    std::size_t size() const { return mOldToNew.size(); }

    const std::vector<std::int32_t>& oldToNew() const { return mOldToNew; }
    const std::vector<std::int32_t>& newToOld() const { return mNewToOld; }

    // Reset to the identity permutation of the given size.
    void setIdentity();   // fill the existing maps (size fixed at construction) to identity

    // Adopt a chosen permutation, given in either direction, and complete the other
    // direction from it. The size is taken from the map. Returns false, leaving the
    // permutation unchanged (so still valid), if the map is not a bijection of
    // 0..size-1. This is 0.9's set(), whose direction flag becomes two named methods.
    //
    // Ordering is not always ours to compute: an ordering may come from a file, from
    // an external tool, or from the problem's own numbering. This is also the only way
    // to test the permutation maps against a known answer rather than against whatever
    // an ordering engine happened to produce.
    bool setOldToNew(const std::vector<std::int32_t>& map);
    bool setNewToOld(const std::vector<std::int32_t>& map);

    // Compose with p, applying this permutation first and p second: p reorders the
    // indices this one has already produced, so afterwards
    //     oldToNew[i] == p.oldToNew[oldToNew[i]]   (with the old value on the right)
    // The usual case is refining an ordering: AMD, then a post-order of the resulting
    // elimination forest. Returns false, leaving this permutation unchanged, if the
    // sizes differ. No revalidation is needed, the composition of two bijections is a
    // bijection.
    bool compose(const Permutation& p);

    // Rebuild newToOld as the inverse of oldToNew. An ordering engine fills oldToNew,
    // then calls this to complete the other direction consistently.
    void rebuildInverse();

    // Check that the maps form a genuine bijection of 0..size()-1 and are consistent
    // inverses. Returns true if valid. (Used to check ordering-engine output.)
    bool validate() const;

private:
    // Is the map a bijection of 0..map.size()-1? (Range and duplicate check, as in
    // 0.9's validate. Does not look at the opposite direction.)
    static bool isBijection(const std::vector<std::int32_t>& map);

    // Fill inv with the inverse of map. Assumes map is a bijection.
    static void invert(const std::vector<std::int32_t>& map,
                       std::vector<std::int32_t>& inv);

    std::vector<std::int32_t> mOldToNew;
    std::vector<std::int32_t> mNewToOld;

    friend class OrderEngine;   // fills the maps; add other engine friends as needed
};

} // namespace Oblio
