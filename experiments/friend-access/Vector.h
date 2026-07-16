#pragma once

// Vector.h - dense vector (guarded-explicit instantiation).
//
// Public API: operator[] for convenient, bounds-checked access. MultiplyEngine is a
// `friend` for direct access to the contiguous storage (mVals) on the fast path.

#include <vector>
#include <complex>
#include <cstddef>
#include <cassert>

namespace Oblio {

class MultiplyEngine;  // befriended below

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

    friend class MultiplyEngine;
};

extern template class Vector<double>;
extern template class Vector<std::complex<double>>;

} // namespace Oblio
