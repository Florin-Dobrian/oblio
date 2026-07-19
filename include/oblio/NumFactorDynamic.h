#pragma once

// NumFactorDynamic.h - the numeric factorization of a sparse matrix, dynamically stored.
//
// **Engines write it now.** The static factorizations (Cholesky, static LDL) run into it unchanged
// through the templated traversals, and dynamic LDL, the reason it exists, writes it through the
// growth verbs below. It began as a placeholder, built next to its sibling before anything filled
// it, because the storage split was a settled decision and because that was the cheap moment to fix
// what the two classes do and do not share.
//
// Dynamic means the structure changes while the arithmetic runs. Dynamic LDL delays an unstable
// pivot by passing its column up to an ancestor, which makes that ancestor's front grow by an
// amount symbolic factorization never predicted. The growth is *local*: one front grows, its
// siblings do not.
//
// Which is why the storage differs:
//
//   NumFactorStatic    flat index and value buffers, an offset per supernode   nothing grows
//   NumFactorDynamic   one index vector and one value vector per supernode      either may grow
//
// The sizes and the node-to-supernode map are identical and identically copied from SymFactor. The
// index sets and value blocks are held one vector per supernode, because a delayed column grows
// both a front's index set and its block, and the growth must stay local. There is deliberately
// **no common base class**. A base exists so that one algorithm can serve two
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
#include <utility>
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

    // How many pivots the factorization had to replace. Meaningful when a *static* factorization
    // (Cholesky, static LDL) runs into this storage, exactly as on the static factor. Dynamic LDL,
    // when it lands, delays an unstable pivot instead of replacing it, and leaves this zero.
    std::size_t numPerturbations() const { return mNumPerturbations; }

    const std::vector<std::int32_t>& nodeToSnode() const { return mNodeToSnode; }

    const std::vector<std::size_t>& frontSize()  const { return mFrontSize; }
    const std::vector<std::size_t>& updateSize() const { return mUpdateSize; }
    std::size_t frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

    // Per supernode: how many of its columns were delayed up to its parent, unpivotable where they
    // stood. Zero until dynamic LDL runs. (No flat nodeIdx()/snodeNodeIdxPtr() here: unlike the
    // static factor, the index sets are stored one vector per supernode, since a delayed column
    // grows a front's index set, so there is no single flat buffer to point into.)
    const std::vector<std::int32_t>& numberOfDelayedColumns() const { return mNumberOfDelayedColumns; }
    std::int32_t numberOfDelayedColumns(std::int32_t kk) const { return mNumberOfDelayedColumns[kk]; }

    // Per column (global): 0 = not yet pivoted, 1 = 1x1 pivot, 2 and 3 = the two halves of a 2x2.
    // The solve reads this to apply 1x1 divisions and 2x2 block solves in the diagonal pass.
    const std::vector<std::int32_t>& pivotType() const { return mPivotType; }

    // The values, one dense block per supernode. Supernode kk's block is val()[kk], column-major,
    // with the same shape as the static case: indexSize rows by frontSize columns. It differs in
    // that it may be resized while the factorization runs.
    const std::vector<std::vector<Val>>& val() const { return mVal; }

    // Where supernode kk's node indices start. One vector per supernode (not a flat buffer as in
    // the static factor), so a delayed column can grow the set; the signature is the same either
    // way, and the engines cannot tell which storage they are reading.
    const std::int32_t* nodeIdx(std::int32_t kk) const { return mNodeIdx[kk].data(); }

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
    const Val*          val(std::int32_t kk)     const { return mVal[kk].data(); }

private:
    // The write path, reached only through friendship: NumFactorEngine fills each block.
    std::int32_t*       nodeIdx(std::int32_t kk) { return mNodeIdx[kk].data(); }
    Val*                val(std::int32_t kk)     { return mVal[kk].data(); }

    // Also the write path: the engine accumulates the perturbation count through this reference
    // (factorSupernode increments it). The const read overload above is public.
    std::size_t& numPerturbations() { return mNumPerturbations; }

    // The growth verbs, the engine's alone, exercised as dynamic LDL delays and pivots. Like the
    // other engine-internal steps (factorSupernode, assembleFromA) they are validated through the
    // factorization's residual rather than in isolation.

    // Grow supernode jj's index set by n slots, for n columns delayed into it. The existing indices
    // stay at the front; the new slots (zero for now) are for the caller to fill with the delayed
    // columns' global indices. Ported from 0.9 extendIndex_.
    void extendIndex(std::int32_t jj, std::int32_t n) {
        mNodeIdx[jj].resize(mNodeIdx[jj].size() + static_cast<std::size_t>(n));
    }

    // Swap columns j_ and k_ of supernode jj, symmetrically. The block is a symmetric matrix stored
    // column-major with leading dimension the index size, so exchanging two pivot columns exchanges
    // the matching rows too. Also swaps the two node indices and repairs the global-to-local map.
    // Ported from 0.9 swap_.
    void swap(std::int32_t jj, std::int32_t j_, std::int32_t k_, std::vector<std::int32_t>& gblToLcl) {
        if (k_ < j_)
            std::swap(j_, k_);

        Val*               block = mVal[jj].data();
        std::int32_t*      idx   = mNodeIdx[jj].data();
        const std::int32_t ld    = static_cast<std::int32_t>(mFrontSize[jj] + mUpdateSize[jj]);

        // Column-major position of (row r, column c) in jj's block, in size_t to avoid overflow.
        const auto at = [ld](std::int32_t r, std::int32_t c) {
            return static_cast<std::size_t>(c) * static_cast<std::size_t>(ld) + static_cast<std::size_t>(r);
        };

        for (std::int32_t i_ = 0; i_ < j_; ++i_)            // rows j_ and k_, in the columns left of j_
            std::swap(block[at(j_, i_)], block[at(k_, i_)]);

        for (std::int32_t i_ = j_ + 1; i_ < k_; ++i_)       // column j_ against row k_, between j_ and k_
            std::swap(block[at(i_, j_)], block[at(k_, i_)]);

        for (std::int32_t i_ = k_ + 1; i_ < ld; ++i_)       // columns j_ and k_, in the rows below k_
            std::swap(block[at(i_, j_)], block[at(i_, k_)]);

        std::swap(block[at(j_, j_)], block[at(k_, k_)]);    // the two diagonal entries

        std::swap(idx[j_], idx[k_]);                        // the node indices and the global-to-local map
        gblToLcl[idx[j_]] = j_;
        gblToLcl[idx[k_]] = k_;
    }

    // Drop the n delayed columns from supernode jj's block, reclaiming their column storage. Called
    // once factorDynamicLDL_ has reduced frontSize[jj] by n and the delayed columns' values have
    // been assembled into the parent: the block loses its n trailing front columns and keeps all
    // its rows (frontSize + n + updateSize), so a column-major truncation does it. The delayed
    // columns live on only as rows, counted now in numberOfDelayedColumns. Ported from 0.9
    // shrinkEntry_.
    void shrinkEntry(std::int32_t jj, std::int32_t n) {
        const std::size_t rows = mFrontSize[jj] + static_cast<std::size_t>(n) + mUpdateSize[jj];
        mVal[jj].resize(mFrontSize[jj] * rows);
    }

    std::size_t   mSize      = 0;
    std::size_t   mSnodeSize = 0;
    Factorization mFactorization = Factorization::DynamicLDLT;

    // Copied from SymFactor, exactly as in NumFactorStatic. Under delayed pivoting the index sets
    // themselves grow, which is the second reason the factor owns a copy rather than referring
    // back: SymFactor's sets are the *predicted* ones and must not be disturbed.
    std::vector<std::int32_t> mNodeToSnode;
    std::vector<std::size_t>  mFrontSize;
    std::vector<std::size_t>  mUpdateSize;

    // Per supernode, filled by dynamic LDL: how many columns it delayed up to its parent.
    std::vector<std::int32_t> mNumberOfDelayedColumns;

    // Per column (global): the pivot kind, 0 / 1x1 / the two halves of a 2x2.
    std::vector<std::int32_t> mPivotType;

    // The index sets and value blocks, one vector per supernode so a front can grow without moving
    // its neighbors. In the static factor both are flat buffers with offsets; here each supernode
    // owns its own, because a delayed column extends both.
    std::vector<std::vector<std::int32_t>> mNodeIdx;
    std::vector<std::vector<Val>>          mVal;

    std::size_t mNumPerturbations = 0;

    friend class NumFactorEngine;
};

extern template class NumFactorDynamic<double>;
extern template class NumFactorDynamic<std::complex<double>>;

} // namespace Oblio
