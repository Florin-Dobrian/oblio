#pragma once

// VectorImplicit.h — Implicit: template inclusion (body in header)
//
// Same pattern as MatrixImplicit.h. Full implementation in the header.
// Implicitly instantiated for every Val type in every .cpp that includes this.

#include <vector>
#include <complex>
#include <cassert>
#include <cstddef>

namespace Oblio {

template<class Val>
class Vector {
public:
    // Constructs an empty vector.
    Vector()
        : mSize(0) {}

    // Constructs a zero-initialised vector of the given size.
    explicit Vector(std::size_t size)
        : mSize(size), mVals(size, Val{0}) {}

    Val& operator[](std::size_t i) {
        assert(i < mSize);
        return mVals[i];
    }

    const Val& operator[](std::size_t i) const {
        assert(i < mSize);
        return mVals[i];
    }

    std::size_t size() const { return mSize; }

private:
    std::size_t      mSize;
    std::vector<Val> mVals;
};

} // namespace Oblio
