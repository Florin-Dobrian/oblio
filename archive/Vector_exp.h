#pragma once

// Vector_exp.h — Case 2: explicit instantiation (forcing only)
//
// Declaration only; implementation in Vector_exp.cc. No member bodies are visible
// here, so other translation units cannot implicitly instantiate Vector<Val> —
// they link the explicit instantiations forced in the .cc. No extern template.

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

} // namespace Oblio
