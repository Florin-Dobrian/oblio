#pragma once
#include <cstddef>
#include <complex>
#include <limits>

namespace Oblio {

using Size  = std::size_t;
using RVal  = double;

inline constexpr Size cNullIdx = std::numeric_limits<Size>::max();

enum class FactorType  { eCholesky, eStaticLDL, eDynamicLDL };
enum class FactorAlg   { eLeftLooking, eRightLooking, eMultifrontal };
enum class OrderAlg    { eNatural, eMMD, eAMD };
enum class SolveScope  { eFull, eLower, eDiagonal, eUpper };
enum class Err         { eNone=0, eInvArg, eInvMat, eInvPerm, eInvPivot, eDimMismatch, eNotImpl };

} // namespace Oblio
