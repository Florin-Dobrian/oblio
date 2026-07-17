#pragma once

// NumFactorDynamic.h - the numeric factorization of a sparse matrix, dynamically stored.
//
// **A placeholder. No engine writes it yet.** It is here because the storage split is a settled
// decision and the class that embodies it should exist next to its sibling, not be summoned into
// being when dynamic LDL is written. It also fixes what the two classes do and do not share, at
// the moment that is cheap to get right.
//
// Dynamic means the structure changes while the arithmetic runs. Dynamic LDL delays an unstable
// pivot by passing its column up to an ancestor, which makes that ancestor's front grow by an
// amount symbolic factorization never predicted. The growth is *local*: one front grows, its
// siblings do not.
//
// Which is why the storage differs, and only the storage differs:
//
//   NumFactorStatic    one flat buffer,  snodeValPtr[kk] an offset into it  nothing grows
//   NumFactorDynamic   one buffer per supernode, val[kk] its own vector     one may grow
//
// Everything above the values is identical, and identically copied from SymFactor. There is
// deliberately **no common base class**. A base exists so that one algorithm can serve two
// storages, and experiments/storage-options measured that a plain array of pointers already does
// that, for a single compiled function, at about a one percent cost. A base would buy nothing the
// pointer array does not, while costing a vtable and, worse, forcing accessors where the engines
// use friendship and direct field access. The engine's kernels take (Val* block, rows, cols, ld)
// and never see either class, so both storages reach them unchanged.
//
// Note what a growable front does to the flat layout, and why it cannot simply be kept: a
// supernode whose block grows must either be reallocated (moving every later supernode in the
// buffer) or preallocated to a worst case that symbolic cannot bound tightly. One vector per
// supernode makes the growth local, which is what the algorithm already is.

#include "oblio/Types.h"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Oblio {

class NumFactorEngine;

template<class Val>
class NumFactorDynamic {
public:
    NumFactorDynamic() = default;

    std::size_t   size()          const { return mSize; }
    std::size_t   snodeSize()     const { return mSnodeSize; }
    Factorization factorization() const { return mFactorization; }

    std::size_t numNodeIdx() const { return mNumNodeIdx; }

    // How many pivots the factorization had to replace. Meaningful when a *static* factorization
    // (Cholesky, static LDL) runs into this storage, exactly as on the static factor. Dynamic LDL,
    // when it lands, delays an unstable pivot instead of replacing it, and leaves this zero.
    std::size_t numPerturbations() const { return mNumPerturbations; }

    const std::vector<std::int32_t>& nodeToSnode() const { return mNodeToSnode; }

    const std::vector<std::size_t>& frontSize()  const { return mFrontSize; }
    const std::vector<std::size_t>& updateSize() const { return mUpdateSize; }
    std::size_t frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

    const std::vector<std::size_t>&  snodeNodeIdxPtr() const { return mSnodeNodeIdxPtr; }
    const std::vector<std::int32_t>& nodeIdx()         const { return mNodeIdx; }

    // The values, one dense block per supernode. Supernode kk's block is val()[kk], column-major,
    // with the same shape as the static case: indexSize rows by frontSize columns. It differs in
    // that it may be resized while the factorization runs.
    const std::vector<std::vector<Val>>& val() const { return mVal; }

    // Where supernode kk's node indices start. Flat, exactly as in the static factor: the index
    // sets are stored identically between the two, only the values differ.
    const std::int32_t* nodeIdxPtr(std::int32_t kk) const { return mNodeIdx.data() + mSnodeNodeIdxPtr[kk]; }

    // Where supernode kk's dense block lives. The static factor's counterpart computes an offset
    // into one flat buffer; this one hands over the column's own vector. Same question, same
    // signature, different storage, and the engines cannot tell which they are talking to.
    //
    // Reading is public, as in the static factor: a consumer that only reads (SolveEngine) reaches
    // these const overloads and needs no friendship; the mutable overloads stay private.
    //
    // **Call it at the moment of use, never hoist it.** Here the warning is not theoretical: a
    // delayed pivot grows an ancestor's front, which resizes its vector, which dangles every
    // pointer previously taken into it.
    const Val*          valPtr(std::int32_t kk)     const { return mVal[kk].data(); }

private:
    // The write path, reached only through friendship: NumFactorEngine fills each block.
    std::int32_t*       nodeIdxPtr(std::int32_t kk) { return mNodeIdx.data() + mSnodeNodeIdxPtr[kk]; }
    Val*                valPtr(std::int32_t kk)     { return mVal[kk].data(); }

    // Also the write path: the engine accumulates the perturbation count through this reference
    // (factorSupernode increments it). The const read overload above is public.
    std::size_t& numPerturbations() { return mNumPerturbations; }

    std::size_t   mSize      = 0;
    std::size_t   mSnodeSize = 0;
    Factorization mFactorization = Factorization::DynamicLDLT;

    // Copied from SymFactor, exactly as in NumFactorStatic. Under delayed pivoting the index sets
    // themselves grow, which is the second reason the factor owns a copy rather than referring
    // back: SymFactor's sets are the *predicted* ones and must not be disturbed.
    std::vector<std::int32_t> mNodeToSnode;
    std::vector<std::size_t>  mFrontSize;
    std::vector<std::size_t>  mUpdateSize;
    std::size_t               mNumNodeIdx = 0;
    std::vector<std::size_t>  mSnodeNodeIdxPtr;
    std::vector<std::int32_t> mNodeIdx;

    // One buffer per supernode, so a front can grow without moving its neighbours.
    std::vector<std::vector<Val>> mVal;

    std::size_t mNumPerturbations = 0;

    friend class NumFactorEngine;
};

extern template class NumFactorDynamic<double>;
extern template class NumFactorDynamic<std::complex<double>>;

} // namespace Oblio
