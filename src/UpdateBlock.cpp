#include "oblio/UpdateBlock.h"

namespace Oblio {

// Explicit instantiation. See the note in NumFactorStatic.cpp.
template class UpdateBlock<double>;
template class UpdateBlock<std::complex<double>>;

} // namespace Oblio
