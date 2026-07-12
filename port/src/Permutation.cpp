#include "oblio/Permutation.h"

#include "oblio/Types.h"   // MAX_IDX

#include <stdexcept>   // std::length_error

namespace Oblio {

Permutation::Permutation(std::size_t size) {
    // Size is fixed here (a permutation of `size` elements); guard it once so the
    // index casts in setIdentity()/fills can never wrap. Indices are std::int32_t.
    if (size > MAX_IDX)
        throw std::length_error(
            "Permutation: size exceeds the std::int32_t index range");
    mOldToNew.resize(size);
    mNewToOld.resize(size);
    setIdentity();
}

void Permutation::setIdentity() {
    // Fill the existing maps to the identity. Size is whatever was set at
    // construction (already bounded by MAX_IDX), so the casts here cannot wrap.
    const std::size_t size = mOldToNew.size();
    for (std::size_t i = 0; i < size; ++i) {
        mOldToNew[i] = static_cast<std::int32_t>(i);
        mNewToOld[i] = static_cast<std::int32_t>(i);
    }
}

bool Permutation::isBijection(const std::vector<std::int32_t>& map) {
    const std::size_t size = map.size();
    if (size > MAX_IDX)
        return false;

    // Every entry in range, and each value hit exactly once. A negative entry casts to
    // a huge std::size_t, so the range test catches it.
    std::vector<bool> seen(size, false);
    for (std::size_t i = 0; i < size; ++i) {
        const std::size_t j = static_cast<std::size_t>(map[i]);
        if (j >= size || seen[j])
            return false;
        seen[j] = true;
    }
    return true;
}

void Permutation::invert(const std::vector<std::int32_t>& map,
                         std::vector<std::int32_t>& inv) {
    const std::size_t size = map.size();
    inv.assign(size, 0);
    for (std::size_t i = 0; i < size; ++i)
        inv[static_cast<std::size_t>(map[i])] = static_cast<std::int32_t>(i);
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
    for (std::size_t i = 0; i < size; ++i) {
        mOldToNew[i] = p.mOldToNew[static_cast<std::size_t>(mOldToNew[i])];
        mNewToOld[static_cast<std::size_t>(mOldToNew[i])] = static_cast<std::int32_t>(i);
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
    for (std::size_t i = 0; i < size; ++i) {
        const std::size_t j = static_cast<std::size_t>(mOldToNew[i]);
        if (j >= size || seen[j])
            return false;
        seen[j] = true;
        if (static_cast<std::size_t>(mNewToOld[j]) != i)
            return false;
    }
    return true;
}

} // namespace Oblio
