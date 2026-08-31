#pragma once

// OrderEngine.h, computes a fill-reducing permutation of a SparseMatrix.
//
//   Natural,       identity (no reordering)
//   MmdVendored,   Multiple Minimum Degree (Liu/Sparspak)
//   MmdCorrected,  the same code with its degree scale repaired         (private/MmdCorrected.cpp)
//   MmdFlat,       ours, matching MmdCorrected's permutation      (src/MmdFlat.cpp)
//   MmdChained,    the same, on the mmd oracle's chained store   (src/MmdChained.cpp)
//   MmdCompacted,  the same, on the amd oracle's store           (src/MmdCompacted.cpp)
//   AmdVendored,   Approximate Minimum Degree (SuiteSparse 3.3.4, BSD-3)
//   AmdFlat,       ours, matching AmdVendored's permutation      (src/AmdFlat.cpp)
//   AmdCompacted,  the same, on the amd oracle's store           (src/AmdCompacted.cpp)  DEFAULT
//
// NINE ENUMERATORS, SIX OF THEM OURS, and TWO AXES rather than one. The BRANCH decides the
// permutation: every mmd enumerator returns `MmdCorrected`'s and every amd one `AmdVendored`'s,
// exactly, on every matrix the benchmarks cover. The STORE decides only what the computation costs.
// So a caller choosing among the three mmd enumerators, or between the two amd ones, is choosing an
// implementation and not an ordering, and the fill is identical by construction.
//
// THE EARLIER LADDER LAYERS WERE RETIRED ON 2026-08-21 and are in retired/. `MMD1`, `MMD2`, `AMD1`
// and `AMD2` were ours too, each carrying the base algorithm without its reference's later
// refinements, so they ordered differently and were not drop-in replacements. They were the
// ladder that got us to MmdFlat and AmdFlat, and the same ladder lives as working C++ and Python
// twins in experiments/ordering/. See retired/README.md.
//
// THE THREE ALTERNATIVE-STORE LAYERS JOINED THIS ENUM ON 2026-08-21, having been free functions
// before it. The entry this replaces argued they were measuring instruments rather than orderings
// and that an enumerator would put a benchmark's oracle into the public enum. The premise was
// sound and the conclusion did not follow: `MmdFlat` against `MmdCorrected` is already the same
// kind of choice, one ordering computed two ways, and it has been offered from the start. What
// the store changes is cost, which is exactly the sort of thing a caller may want to choose.
//
// EACH IS ITS ORIGINAL ON A DIFFERENT CLIQUE STORE, the chained one the mmd oracle uses and the
// compacted one the amd oracle uses, both of which hold the ordering inside `O(n + m)` where our
// flat arena does
// not. Each must return exactly its branch's permutation, and a difference is a defect in one of
// them rather than a variation. That obligation is what makes them useful as instruments, and it
// is also what makes them safe to offer.
//
// `AmdCompacted` IS THE DEFAULT, changed from `MmdFlat` on 2026-08-21. Two decisions, and only the
// second is free. AMD over MMD is the branch decision: on the 246 real matrices the amd branch
// orders in 4.25 s against the mmd branch's 50.7 s, and it is the decision that changes what a
// caller sees, since the two branches return different permutations and different fill. Compacted
// over flat is the store decision and changes nothing but time, the two returning one permutation.
// See benchmarks/matrices/ORDERING.md and docs/DESIGN_DECISIONS.md.
//
// Two lineages sit behind these names. `MmdVendored`, `MmdCorrected` and `AmdVendored` operate on
// raw int CSC arrays and live in private/. The first and third are vendored and untouched; the
// second is `MmdVendored` with its degree scale repaired, so it is inherited code we have
// modified rather than code we wrote, and it is what our three mmd drivers are checked against.
// It files at the degree in its setup and at the degree PLUS ONE in its refresh, two scales in one
// bucket array, so a refreshed vertex is penalised by one against one no pivot has reached and the
// minimum is not always the minimum; see private/MmdCorrected.cpp. The other six are ours, built
// over the quotient graph classes in include/oblio/. Either way this engine is the seam that reads
// the matrix structure and fills the Permutation. Returns true on success.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {


enum class Ordering { Natural, MmdVendored, MmdCorrected, MmdFlat, MmdChained, MmdCompacted,
                      AmdVendored, AmdFlat, AmdCompacted };

class OrderEngine {
public:
    OrderEngine() = default;
    explicit OrderEngine(Ordering ordering) : mOrdering(ordering) {}

    void     setOrdering(Ordering ordering) { mOrdering = ordering; }
    Ordering ordering() const               { return mOrdering; }

    // Order A into P.
    //
    // An ordering is a pure graph operation: AMD and MMD read the sparsity pattern and
    // would not know what to do with a value. So the implementation is the non-templated
    // overload, taking colPtr and rowIdx, and it is compiled once. The templated overload
    // is an adapter over it, for the common case of holding a matrix. The pattern overload
    // is public: a caller holding a graph with no numbers attached can order it without
    // inventing a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, Permutation& P) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 Permutation& P) const;

private:
    // Ours, not a vendored one, because the vendored orderings are optional: a default has to be an
    // ordering that is always there. That constraint has held since 2026-08-07 and still does.
    //
    // `AmdCompacted`, replacing `MmdFlat`. The old entry's argument was that a
    // default is a bet on the cases nobody has run, so reproducing a reference beats a tie-break
    // of our own. That argument is unchanged and `AmdCompacted` satisfies it: it returns the amd
    // oracle's
    // permutation EXACTLY on every example, every square grid, and all 246 matrices in
    // benchmarks/matrices, compaction for compaction. What changed is that the amd branch now has
    // a layer that meets the bar, where in 2026-08-07 only the mmd one did.
    //
    // Two decisions, and only the second is free. AMD over MMD is the branch, and it is the one a
    // caller sees: different permutation, different fill, and on the 246 real matrices 4.25 s of
    // ordering against 50.7 s. Compacted over flat is the store, and it changes nothing but time,
    // the two returning one permutation. See benchmarks/matrices/ORDERING.md.
    Ordering mOrdering = Ordering::AmdCompacted;

    bool orderNatural(std::size_t size, Permutation& P) const;
    // The two vendored orderings, declared only when private/ supplies their sources. Everything
    // else here is ours and is always present.
#ifdef OBLIO_VENDORED_ORDERINGS
    bool orderMmdVendored(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
    bool orderMmdCorrected(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
    bool orderAmdVendored(std::size_t size,
                  const std::vector<std::size_t>&  colPtr,
                  const std::vector<std::int32_t>& rowIdx,
                  Permutation& P) const;
#endif
    bool orderMmdFlat(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAmdFlat(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderMmdChained(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderMmdCompacted(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
    bool orderAmdCompacted(std::size_t size,
                   const std::vector<std::size_t>&  colPtr,
                   const std::vector<std::int32_t>& rowIdx,
                   Permutation& P) const;
};

extern template bool OrderEngine::compute(const SparseMatrix<double>&, Permutation&) const;
extern template bool OrderEngine::compute(const SparseMatrix<std::complex<double>>&, Permutation&) const;

} // namespace Oblio
