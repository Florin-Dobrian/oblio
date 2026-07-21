#include "oblio/UpdateMatrix.h"

namespace Oblio {

template<class Val>
void UpdateMatrix<Val>::allocate(std::size_t size) {
    mSize = size;
    mNodeIdx.assign(size, 0);          // sized; the engine fills the global indices
    mVal.assign(size * size, Val(0));  // zeroed, because assembly adds into it
}

template<class Val>
void UpdateMatrix<Val>::discard() {
    mSize = 0;
    // Free rather than clear: clear keeps the capacity, which would let the stack's peak grow to
    // the sum of every block ever allocated. Swapping with an empty vector releases the storage.
    std::vector<std::int32_t>().swap(mNodeIdx);
    std::vector<Val>().swap(mVal);
}

// Explicit instantiation. See the note in NumFactorStatic.cpp.
template class UpdateMatrix<double>;
template class UpdateMatrix<std::complex<double>>;

} // namespace Oblio
