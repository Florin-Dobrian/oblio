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
    const std::size_t n = mOldToNew.size();
    for (std::size_t i = 0; i < n; ++i) {
        mOldToNew[i] = static_cast<std::int32_t>(i);
        mNewToOld[i] = static_cast<std::int32_t>(i);
    }
}

void Permutation::rebuildInverse() {
    const std::size_t n = mOldToNew.size();
    mNewToOld.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        mNewToOld[static_cast<std::size_t>(mOldToNew[i])] = static_cast<std::int32_t>(i);
}

bool Permutation::validate() const {
    const std::size_t n = mOldToNew.size();
    if (mNewToOld.size() != n)
        return false;

    // Every old index maps into range, each new index hit exactly once, and the two
    // maps are consistent inverses.
    std::vector<bool> seen(n, false);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = static_cast<std::size_t>(mOldToNew[i]);
        if (j >= n || seen[j])
            return false;
        seen[j] = true;
        if (static_cast<std::size_t>(mNewToOld[j]) != i)
            return false;
    }
    return true;
}

} // namespace Oblio
