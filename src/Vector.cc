#include "oblio/Vector.h"
#include <algorithm>
#include <complex>

namespace Oblio {

template<class Val>
Vector<Val>::Vector() {}

template<class Val>
Vector<Val>::Vector(Size n)
    : mData(n, Val{0}) {}

template<class Val>
Vector<Val>::Vector(Size n, Val fill)
    : mData(n, fill) {}

template<class Val>
Vector<Val>::Vector(Size n, const Val* src)
    : mData(src, src + n) {}

template<class Val>
Vector<Val>::Vector(std::initializer_list<Val> il)
    : mData(il) {}

template<class Val>
Vector<Val>::Vector(const std::vector<Val>& v)
    : mData(v) {}

template<class Val>
Vector<Val>::Vector(const Vector& o)
    : mData(o.mData) {}

template<class Val>
Vector<Val>::Vector(Vector&& o) noexcept
    : mData(std::move(o.mData)) {}

template<class Val>
Vector<Val>& Vector<Val>::operator=(const Vector& o) {
    if (this != &o) mData = o.mData;
    return *this;
}

template<class Val>
Vector<Val>& Vector<Val>::operator=(Vector&& o) noexcept {
    if (this != &o) mData = std::move(o.mData);
    return *this;
}

template<class Val>
Size Vector<Val>::size() const { return mData.size(); }

template<class Val>
bool Vector<Val>::empty() const { return mData.empty(); }

template<class Val>
Val* Vector<Val>::data() { return mData.data(); }

template<class Val>
const Val* Vector<Val>::data() const { return mData.data(); }

template<class Val>
Val& Vector<Val>::operator[](Size i) { return mData[i]; }

template<class Val>
const Val& Vector<Val>::operator[](Size i) const { return mData[i]; }

template<class Val>
void Vector<Val>::resize(Size n) {
    mData.resize(n, Val{0});
}

template<class Val>
void Vector<Val>::resize(Size n, Val fill) {
    Size old = mData.size();
    mData.resize(n);
    if (n > old) std::fill(mData.begin() + old, mData.end(), fill);
}

template<class Val>
void Vector<Val>::setZero() {
    std::fill(mData.begin(), mData.end(), Val{0});
}

template<class Val>
void Vector<Val>::fill(Val v) {
    std::fill(mData.begin(), mData.end(), v);
}

template<class Val>
std::vector<Val> Vector<Val>::toStdVector() const {
    return mData;
}

// ---- explicit instantiations ------------------------------------------------
template class Vector<double>;
template class Vector<std::complex<double>>;

} // namespace Oblio
