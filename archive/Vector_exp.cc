// Vector_exp.cc — explicit instantiation style

#include "Vector_exp.h"
#include <cassert>

namespace Oblio {

template<class Val>
Vector<Val>::Vector()
    : mSize(0) {}

template<class Val>
Vector<Val>::Vector(std::size_t size)
    : mSize(size), mVals(size, Val{0}) {}

template<class Val>
Val& Vector<Val>::operator[](std::size_t i) {
    assert(i < mSize);
    return mVals[i];
}

template<class Val>
const Val& Vector<Val>::operator[](std::size_t i) const {
    assert(i < mSize);
    return mVals[i];
}

template<class Val>
std::size_t Vector<Val>::size() const { return mSize; }

// ── Explicit instantiations ───────────────────────────────────────────────

template class Vector<double>;
template class Vector<std::complex<double>>;

} // namespace Oblio
