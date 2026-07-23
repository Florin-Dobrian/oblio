#ifndef OBLIO_EXPERIMENTS_OPENMP_GEMM_H
#define OBLIO_EXPERIMENTS_OPENMP_GEMM_H

#include <cstddef>

namespace Oblio {

// C := C + A B, on n by n column-major blocks with leading dimension n.
//
// Two kernels computing the same thing by different means, which is the whole
// point of the experiment: the hand loop uses only per-core resources (registers,
// cache, NEON), while the BLAS call may reach the AMX coprocessor and Accelerate's
// own threading. Running the pair of them across two cores is what separates
// "OpenMP gave us two cores" from "the two cores are contending underneath".
//
// The hand loop is deliberately untuned. It is a control, not a competitor.

void gemmHand(std::size_t n, const double* a, const double* b, double* c);
void gemmBlas(std::size_t n, const double* a, const double* b, double* c);

}  // namespace Oblio

#endif  // OBLIO_EXPERIMENTS_OPENMP_GEMM_H
