#include "oblio/Types.h"

#include <stdexcept>   // std::length_error
#include <string>      // std::string, to build the message

namespace Oblio {

// The single home for the index-range guard and its message. Kept out of Types.h so the throw does
// not sit in a header that the numeric kernels include (see the codegen note in DESIGN_DECISIONS).
std::size_t checkIndexRange(std::size_t size, const char* what) {
    if (size > MAX_IDX)
        throw std::length_error(std::string(what) + " exceeds the std::int32_t index range");
    return size;
}

} // namespace Oblio
