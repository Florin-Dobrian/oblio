#include "oblio/Vector.h"

#include "oblio/Types.h"   // checkIndexRange

#include <utility>   // std::move

namespace Oblio {

// A vector's length is the system dimension, which no valid matrix can push past MAX_IDX (indices
// are std::int32_t, and the dimension is walked by int32 loop counters), so a longer vector is
// meaningless and is rejected at construction, as SparseMatrix and Permutation reject theirs.
// checkIndexRange guards mSize before mVal is constructed (members initialize in declaration order),
// so the first constructor checks before it allocates. Both constructors stay here rather than in
// the header so the exception path stays out of the units that compile the hot multiply and solve
// kernels.

template<class Val>
Vector<Val>::Vector(std::size_t size)
    : mSize(checkIndexRange(size, "Vector size")),
      mVal(size, Val(0)) {}

template<class Val>
Vector<Val>::Vector(std::size_t size, std::vector<Val> val)
    : mSize(checkIndexRange(size, "Vector size")),
      mVal(std::move(val)) {}

template class Vector<double>;
template class Vector<std::complex<double>>;

} // namespace Oblio
