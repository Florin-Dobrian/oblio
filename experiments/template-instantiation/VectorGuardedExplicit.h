#pragma once

// VectorGuardedExplicit.h — Guarded explicit: plain explicit + extern template guard
//
// Declaration only. Implementation in VectorGuardedExplicit.cpp.
// extern template lines suppress implicit instantiation everywhere else.

#include <vector>
#include <complex>
#include <cstddef>

namespace Oblio {

template<class Val>
class Vector {
public:
    Vector();
    explicit Vector(std::size_t size);

    Val&       operator[](std::size_t i);
    const Val& operator[](std::size_t i) const;

    std::size_t size() const;

private:
    std::size_t      mSize;
    std::vector<Val> mVals;
};

extern template class Vector<double>;
extern template class Vector<std::complex<double>>;

} // namespace Oblio
