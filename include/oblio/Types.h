#pragma once

// Types.h - shared index-type conventions (see the "index types" design decision),
// plus the enums that select what the numeric factorization computes and how.
//
// IDs (vertices/rows/columns/supernodes, which may carry a "none" marker) are
// std::int32_t; offsets/counts/sizes are std::size_t. NIL is the shared "none"
// sentinel for ID arrays (parent/child/sibling links, supernode maps, etc.).
// MAX_IDX is the largest representable ID; construction guards check inputs against
// it so an over-range value fails loudly rather than wrapping at a narrowing cast.

#include <cstddef>
#include <complex>
#include <cstdint>
#include <limits>

namespace Oblio {

// "No such ID" sentinel for std::int32_t link/ID arrays (e.g. a tree root's parent).
inline constexpr std::int32_t NIL = -1;

// Largest valid ID, as a std::size_t for comparison against sizes/offsets. Equal to
// INT32_MAX. Bounds both the matrix dimension (indices are std::int32_t) and, for now,
// nnz (the entry count is narrowed to int at the vendored AMD/MMD boundary); those two
// caps happen to coincide at this value today but have different origins.
inline constexpr std::size_t MAX_IDX =
    static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

// Guard a size or count against the index range. Throws std::length_error if `size` exceeds
// MAX_IDX, naming `what` (e.g. "Vector size", "SparseMatrix nnz") in the message; otherwise returns
// `size` unchanged, so it can guard an init-list member directly:
//
//   Vector(std::size_t size) : mSize(checkIndexRange(size, "Vector size")), mVal(size, Val(0)) {}
//
// Members initialize in declaration order, so a guard in the first member's initializer runs before
// any later member allocates. That is how a class that allocates from its size (Vector) checks
// before allocating; a class that moves or resizes instead (SparseMatrix, Permutation) can call it
// as a plain statement in the body. Defined in Types.cpp, not here: the throw is kept out of this
// widely-included header so the exception path stays out of the units that compile the numeric
// kernels (an in-header throw was measured to perturb a hot loop's codegen).
std::size_t checkIndexRange(std::size_t size, const char* what);

// What the numeric factorization computes. **The symmetry is part of this, not a separate
// setting.**
//
// For Cholesky there is nothing to choose: `A = CC^H` always, and over the reals `C^H` is `C^T`,
// so real symmetric and real Hermitian are the same case. But for LDL, a complex matrix may be
// symmetric (`A = A^T`, `D` complex) or Hermitian (`A = A^H`, `D` real), and those are genuinely
// different factorizations. So the transpose is named:
//
//   ...LDLT     A = LDL^T    plain transpose.        Complex: D is complex.
//   ...LDLH     A = LDL^H    conjugate transpose.    Complex: D is real.
//
// Over the reals the two coincide, and either name computes the same thing.
//
// Folding symmetry into this enum, rather than adding a `Symmetry` setting alongside it, is
// deliberate and it is the same principle as BlasLapack's operation-named wrappers: **make the
// wrong thing unwriteable.** A separate flag would let a caller ask for a complex *symmetric*
// Cholesky, which we would then have to reject at runtime. Here there is no such value to name.
// The combination is not forbidden; it does not exist.
//
// (Why it could not exist: positive definiteness means `x* A x > 0`, which requires that quantity
// to be real, which holds for all `x` exactly when `A` is Hermitian. For a complex symmetric `A`,
// `x^T A x` is a complex number and the inequality does not typecheck. Cholesky has no pivoting to
// recover with, because not needing pivoting is the entire point of Cholesky, and that guarantee
// comes from positive definiteness. LAPACK has no complex-symmetric Cholesky for this reason.)
enum class Factorization {
    Cholesky,      // A = CC^H  (CC^T for real). Requires positive definiteness, hence Hermitian.
    StaticLDLT,    // A = LDL^T, pivots fixed by the symbolic structure. Complex: D complex.
    StaticLDLH,    // A = LDL^H, ditto.                                  Complex: D real.
    DynamicLDLT,   // A = LDL^T, with 2x2 pivoting and delayed columns.
    DynamicLDLH    // A = LDL^H, ditto.
};

// Two facts a factorization implies, answerable from the enum alone. `hermitian` is true when the
// factorization conjugates: Cholesky always (A = CC^H), LDLH yes, LDLT no. Over the reals the
// conjugate is the identity, so the flag is harmless there. `separateDiagonal` is true when the
// diagonal is a pass of its own: LDL keeps L unit and holds D apart, while Cholesky folds its own
// diagonal into C.
// Conjugate a value, or leave it. `hermitian` decides, and for `double` both overloads are the
// identity, which is what lets one templated kernel serve the symmetric and Hermitian
// factorizations alike rather than being written twice.
//
// These sit beside the predicate that drives them because three separate parts of the library need
// them: the dense kernels in BlasLapack, the dynamic pivot and swap code, and the solve.
inline double               maybeConjugate(double v, bool)                    { return v; }
inline std::complex<double> maybeConjugate(std::complex<double> v, bool herm) {
    return herm ? std::conj(v) : v;
}

// **A correction for floating point, not for the algorithm.** A Hermitian matrix has a real
// diagonal, so the imaginary part of such an entry is zero mathematically and there is nothing to
// decide. It is not zero in storage. So this restores a value that is already determined rather
// than choosing one: nothing is tuned, there is no threshold, and the true imaginary part is known
// exactly.
//
// Where the nonzero comes from is worth knowing, because it is not generic rounding but a
// cancellation that exact arithmetic guarantees and floating point does not. The kernels form a
// diagonal update as `l * u` with `u = d * conj(l)` materialized first, `d` real. Writing
// `l = a + bi`, the imaginary part of that product is
//
//     b * fl(d*a) - a * fl(d*b) = abd(1 + e1) - abd(1 + e2) = abd(e1 - e2)
//
// which is zero only if the two roundings agree. Reassociating to `l * conj(l) * d` would cancel
// bitwise, the two cross terms then being the same product with opposite signs, but that is not
// the association the update uses: materializing `u` is what lets the multiply be a plain gemm.
// The residues are of order epsilon apiece and accumulate over a supernode's descendants, and the
// entry they land on is about to become a pivot that its whole column is divided by.
//
// Applied to pivots only, which is the whole of D, off-diagonal entries having no such constraint
// to restore. Four call sites, and they divide into two kinds:
//
//   ldl, factor1x1, factor2x2's diagonals    **the last moment.** These correct a value that is
//                                            about to be divided by, and then store the corrected
//                                            one where the solve will divide by it again. Nothing
//                                            downstream could put it right, having no way to know
//                                            it should be real.
//   readPivotBlock2x2 from acceptPivot2x2    not the last moment, but it matters all the same: the
//                                            acceptance test compares magnitudes against the 2x2
//                                            determinant, and that determinant is real only if the
//                                            diagonals are. A residue there perturbs the test that
//                                            decides whether the pivot is used at all.
inline double               forceReal(double v, bool)                    { return v; }
inline std::complex<double> forceReal(std::complex<double> v, bool herm) {
    return herm ? std::complex<double>(v.real(), 0.0) : v;
}

inline bool hermitian(Factorization factorization) {
    return factorization == Factorization::Cholesky
        || factorization == Factorization::StaticLDLH
        || factorization == Factorization::DynamicLDLH;
}
inline bool separateDiagonal(Factorization factorization) {
    return factorization != Factorization::Cholesky;
}

// Whether pivots are chosen while the arithmetic runs. This is the other axis the enum names, and
// it is about *pivoting*, not storage: Cholesky is static by nature (a positive definite matrix
// needs no pivot search), while LDL can go either way, statically for matrices that tolerate the
// pivots symbolic assigned, dynamically for the hard ones, where an unstable pivot is delayed to an
// ancestor.
//
// **Dynamic pivoting implies dynamic storage**, because delaying a column expands a front and a
// flat
// buffer cannot expand; the reverse does not hold, since static pivoting runs in either storage
// (and
// prefers the flat one, for locality). This predicate is where that rule is stated once: the engine
// asks it to pick a traversal, DirectSolver asks it to pick a storage.
inline bool dynamicPivoting(Factorization factorization) {
    return factorization == Factorization::DynamicLDLT
        || factorization == Factorization::DynamicLDLH;
}

// How the supernodes are traversed. All three compute the same factor; they differ in when the
// updates are formed and where they are kept.
//
//   LeftLooking   for each supernode: pull every update owed to it, then factor it.
//   RightLooking  for each supernode: factor it, then push its updates to every ancestor.
//   Multifrontal  for each supernode: assemble its children's update matrices from a stack,
//                 factor, and push its own Schur complement back on the stack.
enum class Traversal {
    LeftLooking,
    RightLooking,
    Multifrontal
};

} // namespace Oblio
