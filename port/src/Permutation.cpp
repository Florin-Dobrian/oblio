#include "oblio/Permutation.h"

namespace Oblio {

Permutation::Permutation(std::size_t size) {
    setIdentity(size);
}

void Permutation::setIdentity(std::size_t size) {
    mOldToNew.resize(size);
    mNewToOld.resize(size);
    for (std::size_t i = 0; i < size; ++i) {
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
