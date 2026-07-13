#include "oblio/NumFactorDynamic.h"

namespace Oblio {

// Explicit instantiation. See the note in NumFactorStatic.cpp.
template class NumFactorDynamic<double>;
template class NumFactorDynamic<std::complex<double>>;

} // namespace Oblio
