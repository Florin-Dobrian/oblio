#pragma once
#include "oblio/Types.h"
#include <vector>
#include <complex>
#include <initializer_list>

namespace Oblio {

// Vector<Val> — dense vector of length n, owns its storage.
// Thin wrapper over a contiguous Val array; provides the interface
// consumed by SolveEngine and OblioEngine.

template<class Val>
class Vector {
public:
    // ---- construction ----
    Vector();
    explicit Vector(Size n);                      // zero-initialised
    Vector(Size n, Val fill);                     // filled with fill
    Vector(Size n, const Val* data);              // copy from raw array
    Vector(std::initializer_list<Val> il);
    explicit Vector(const std::vector<Val>& v);   // copy from std::vector
    Vector(const Vector&);
    Vector(Vector&&) noexcept;
    Vector& operator=(const Vector&);
    Vector& operator=(Vector&&) noexcept;
    ~Vector() = default;

    // ---- size / data access ----
    Size        size()                  const;
    bool        empty()                 const;
    Val*        data();
    const Val*  data()                  const;
    Val&        operator[](Size i);
    const Val&  operator[](Size i)      const;

    // ---- mutators ----
    void        resize(Size n);                   // preserves existing values up to min(old,new)
    void        resize(Size n, Val fill);         // new slots filled with fill
    void        setZero();
    void        fill(Val v);

    // ---- conversion ----
    std::vector<Val> toStdVector()      const;

private:
    std::vector<Val> mData;
};

// Explicit instantiation declarations — definitions in Vector.cc
extern template class Vector<double>;
extern template class Vector<std::complex<double>>;

} // namespace Oblio
