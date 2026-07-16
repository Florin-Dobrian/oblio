#include "oblio/Permutation.h"

#include "oblio/Types.h"   // checkIndexRange

namespace Oblio {

// checkIndexRange guards the size in the first member's initializer, so it throws before either map
// allocates; the maps are then sized in the init list and setIdentity fills them. Indices are
// std::int32_t, so the bounded size keeps the casts in setIdentity from wrapping.
Permutation::Permutation(std::size_t size)
    : mOldToNew(checkIndexRange(size, "Permutation size")),
      mNewToOld(size) {
    setIdentity();
}

void Permutation::setIdentity() {
    // Fill the existing maps to the identity. The size was bounded by MAX_IDX at
    // construction, so the counter cannot overflow an index.
    const std::size_t size = mOldToNew.size();
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i) {
        mOldToNew[i] = i;
        mNewToOld[i] = i;
    }
}

bool Permutation::isBijection(const std::vector<std::int32_t>& map) {
    const std::size_t size = map.size();
    if (size > MAX_IDX)
        return false;

    // Every entry in range, and each value hit exactly once. The bounds are checked as
    // indices, explicitly: a widening cast would turn a negative entry into a huge position
    // and the range test would catch it by accident, which is a trick a reader has to know.
    // Two comparisons say it outright.
    std::vector<bool> seen(size, false);
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i) {
        const std::int32_t j = map[i];
        if (j < 0 || j >= static_cast<std::int32_t>(size) || seen[j])
            return false;
        seen[j] = true;
    }
    return true;
}

void Permutation::invert(const std::vector<std::int32_t>& map,
                         std::vector<std::int32_t>& inv) {
    const std::size_t size = map.size();
    inv.assign(size, 0);
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i)
        inv[map[i]] = i;
}

bool Permutation::setOldToNew(const std::vector<std::int32_t>& map) {
    if (!isBijection(map))
        return false;   // rejected before any state changes, so the permutation stays valid
    mOldToNew = map;
    invert(mOldToNew, mNewToOld);
    return true;
}

bool Permutation::setNewToOld(const std::vector<std::int32_t>& map) {
    if (!isBijection(map))
        return false;
    mNewToOld = map;
    invert(mNewToOld, mOldToNew);
    return true;
}

bool Permutation::compose(const Permutation& p) {
    const std::size_t size = mOldToNew.size();
    if (p.size() != size)
        return false;

    // This permutation runs first, p second, so p maps the indices this one produces.
    // The update is in place: slot i is read before it is written, and no other slot is
    // touched, so there is no aliasing. The new value then indexes newToOld, which
    // rebuilds the inverse as we go (every slot written exactly once, the composition
    // being a bijection), so no separate inversion pass is needed.
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i) {
        mOldToNew[i] = p.mOldToNew[mOldToNew[i]];
        mNewToOld[mOldToNew[i]] = i;
    }
    return true;
}

void Permutation::rebuildInverse() {
    invert(mOldToNew, mNewToOld);
}

bool Permutation::validate() const {
    const std::size_t size = mOldToNew.size();
    if (mNewToOld.size() != size)
        return false;

    // Every old index maps into range, each new index hit exactly once, and the two
    // maps are consistent inverses.
    std::vector<bool> seen(size, false);
    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i) {
        const std::int32_t j = mOldToNew[i];
        if (j < 0 || j >= static_cast<std::int32_t>(size) || seen[j])
            return false;
        seen[j] = true;
        if (mNewToOld[j] != i)
            return false;
    }
    return true;
}

} // namespace Oblio
