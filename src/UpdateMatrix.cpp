#include "oblio/UpdateMatrix.h"

namespace Oblio {

template<class Val>
void UpdateMatrix<Val>::allocate(std::size_t size) {
    mSize = size;
    mRowIdx.assign(size, 0);           // sized; the engine fills the global row indices
    mVal.assign(size * size, Val(0));  // zeroed, because assembly adds into it
}

template<class Val>
void UpdateMatrix<Val>::discard() {
    // Back to the default-constructed state, which is where the empty values are stated: mSize's
    // own initializer and two empty vectors. Move-assigning from a fresh object deallocates both
    // buffers, where clear() would keep the capacity and let the stack's peak grow to the sum of
    // every block ever allocated. Freeing is what keeps the peak bounded.
    //
    // `std::vector<T>().swap(v)` did the same job here until 2026-08-13 and is the pre-C++11 idiom
    // for it, from when clear() did not free and shrink_to_fit did not exist. shrink_to_fit is not
    // the replacement: it is a non-binding request, so it cannot be relied on to free.
    *this = UpdateMatrix<Val>();
}

// Explicit instantiation. See the note in NumFactorStatic.cpp.
template class UpdateMatrix<double>;
template class UpdateMatrix<std::complex<double>>;

} // namespace Oblio
