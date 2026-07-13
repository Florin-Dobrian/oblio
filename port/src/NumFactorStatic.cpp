#include "oblio/NumFactorStatic.h"

namespace Oblio {

// Explicit instantiation. The class is data only, so there is nothing to define here, but the
// instantiations must live in one translation unit and not be emitted per-TU: that is what the
// `extern template` lines in the header suppress, and it is the whole point of the
// explicit-instantiation decision. Adding a scalar type is one line.
template class NumFactorStatic<double>;
template class NumFactorStatic<std::complex<double>>;

} // namespace Oblio
