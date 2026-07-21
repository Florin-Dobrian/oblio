#pragma once

// NumFactorDynamic.h - the numeric factorization of a sparse matrix, dynamically stored.
//
// **Engines write it now.** The static factorizations (Cholesky, static LDL) run into it unchanged
// through the templated traversals, and dynamic LDL, the reason it exists, writes it through the
// expansion and contraction verbs below. It began as a placeholder, built next to its sibling
// before anything filled
// it, because the storage split was a settled decision and because that was the cheap moment to fix
// what the two classes do and do not share.
//
// Dynamic means the structure changes while the arithmetic runs. Dynamic LDL delays an unstable
// pivot by passing its column up to an ancestor, which makes that ancestor's front expand by an
// amount symbolic factorization never predicted. The expansion is *local*: one front expands, its
// siblings do not.
//
// Which is why the storage differs:
//
//   NumFactorStatic    flat index and value buffers, an offset per supernode   nothing expands
//   NumFactorDynamic   one index vector and one value vector per supernode      either may expand
//
// The sizes and the node-to-supernode map are identical and identically copied from SymFactor. The
// index sets and value blocks are held one vector per supernode, because a delayed column expands
// both a front's index set and its block, and the expansion must stay local. There is deliberately
// **no common base class**. A base exists so that one algorithm can serve two
// storages, and experiments/storage-options measured that a plain array of pointers already does
// that, for a single compiled function, at about a one percent cost. A base would buy nothing the
// pointer array does not, while costing a vtable and, worse, forcing accessors where the engines
// use friendship and direct field access. The engine's kernels take (Val* block, rows, cols, ld)
// and never see either class, so both storages reach them unchanged.
//
// Note what an expandable front does to the flat layout, and why it cannot simply be kept: a
// supernode whose block expands must either be reallocated (moving every later supernode in the
// buffer) or preallocated to a worst case that symbolic cannot bound tightly. One vector per
// supernode makes the expansion local, which is what the algorithm already is.

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

    // Node to supernode: the whole map, and the supernode a single node belongs to. Indexed by a
    // *node*, unlike frontSize and its kin, which are indexed by a supernode.
    const std::vector<std::int32_t>& nodeToSnode()                  const { return mNodeToSnode; }
    std::int32_t                     nodeToSnode(std::int32_t node) const { return mNodeToSnode[node]; }

    // **The three regions of a supernode's block, in the order they occupy it.** The height is
    // frontSize + delaySize + updateSize and is conserved throughout: expanding moves rows into the
    // front, factoring moves columns from the front into the delay region, and updateSize is never
    // rewritten. Every leading dimension in the dynamic path is that sum.
    //
    // delaySize is how many of this supernode's columns were delayed up to its *parent*, unpivotable
    // where they stood, never how many were delayed into it, which is a sum over its children and
    // is never stored. Zero until dynamic LDL runs, and never cleared afterwards: the delayed
    // columns leave, but their rows stay, so the count remains part of this block's geometry.
    const std::vector<std::size_t>&  frontSize()  const { return mFrontSize; }
    const std::vector<std::size_t>&  delaySize()  const { return mDelaySize; }
    const std::vector<std::size_t>&  updateSize() const { return mUpdateSize; }

    std::size_t  frontSize(std::int32_t kk)  const { return mFrontSize[kk]; }
    std::size_t  delaySize(std::int32_t kk)  const { return mDelaySize[kk]; }
    std::size_t  updateSize(std::int32_t kk) const { return mUpdateSize[kk]; }

    // Per column (global): 0 = not yet pivoted, 1 = 1x1 pivot, 2 and 3 = the two halves of a 2x2.
    // The solve reads this to apply 1x1 divisions and 2x2 block solves in the diagonal pass.
    const std::vector<std::int32_t>& pivotType() const { return mPivotType; }

    // The values, one dense block per supernode. Supernode kk's block is val()[kk], column-major,
    // with the same shape as the static case: indexSize rows by frontSize columns. It differs in
    // that it may be resized while the factorization runs.
    const std::vector<std::vector<Val>>& val() const { return mVal; }

    // Where supernode kk's node indices start. One vector per supernode (not a flat buffer as in
    // the static factor), so a delayed column can expand the set; the signature is the same either
    // way, and the engines cannot tell which storage they are reading.
    //
    // There is deliberately no `snodeNodeIdxPtr()` counterpart. The static factor offers one
    // because its index sets live end to end in a single buffer and an offset per supernode is
    // meaningful; here there is no such buffer to point into, and inventing an accessor that
    // returned something offset-shaped would promise a layout this class does not have.
    const std::int32_t* nodeIdx(std::int32_t kk) const { return mNodeIdx[kk].data(); }

    // Where supernode kk's dense block lives. The static factor's counterpart computes an offset
    // into one flat buffer; this one hands over the column's own vector. Same question, same
    // signature, different storage, and the engines cannot tell which they are talking to.
    //
    // Reading is public, as in the static factor: a consumer that only reads (SolveEngine) reaches
    // these const overloads and needs no friendship; the mutable overloads stay private.
    //
    // **Call it at the moment of use, never hoist it.** Here the warning is not theoretical: a
    // delayed pivot expands an ancestor's front, which resizes its vector, which dangles every
    // pointer previously taken into it.
    const Val*          val(std::int32_t kk)     const { return mVal[kk].data(); }

private:
    // The write path, reached only through friendship: NumFactorEngine fills each block.
    std::int32_t*       nodeIdx(std::int32_t kk) { return mNodeIdx[kk].data(); }
    Val*                val(std::int32_t kk)     { return mVal[kk].data(); }

    // Also the write path: the engine accumulates the perturbation count through this reference
    // (factorStaticSupernode increments it). The const read overload above is public.
    std::size_t& numPerturbations() { return mNumPerturbations; }

    // The expansion and contraction verbs, the engine's alone, exercised as dynamic LDL delays
    // and pivots. Like the other engine-internal steps (factorStaticSupernode, assembleFromA)
    // they are validated through the factorization's residual rather than in isolation.

    // Grow supernode jj's index set by n slots, for n columns delayed into it. The existing indices
    // stay at the front; the new slots (zero for now) are for the caller to fill with the delayed
    // columns' global indices. Ported from 0.9 extendIndex_.
    void expandNodeIdx(std::int32_t jj, std::size_t n) {
        mNodeIdx[jj].resize(mNodeIdx[jj].size() + n);
    }

    // Size supernode jj's block to its *current* shape and zero it: frontSize columns by
    // frontSize + updateSize rows. The other half of expanding a front, called once frontSize has
    // absorbed the columns delayed into it, and only there, since the old contents are discarded.
    //
    // 0.9 spells this as three calls, discardEntry_ then allocateEntry_ then zeroEntry_, because it
    // manages the block by hand. One vector per supernode collapses all three into an assign, which
    // is the whole of what the storage choice buys here. Left-looking never expands a block in
    // place,
    // which is why it can discard; expandVal below is right-looking's counterpart, and 0.9's
    // extendEntry_ is what that ports.
    void resetVal(std::int32_t jj) {
        mVal[jj].assign(mFrontSize[jj] * (mFrontSize[jj] + mUpdateSize[jj]), Val(0));
    }

    // The other way to expand a block, and the one right-looking needs: size jj's block to its
    // current shape while **keeping what is already in it**, the n new columns arriving empty at
    // the left. Ported from 0.9 extendEntry_.
    //
    // The two traversals differ in exactly this. Left-looking assembles A into a front immediately
    // before factoring it, so when the front expands there is nothing in it worth keeping and
    // resetVal discards. Right-looking assembles A into every front at the start and then pushes
    // each supernode's update into its ancestors, so by the time a front expands it already holds
    // A's
    // values and every update from every descendant already factored. Discarding there would throw
    // the factorization away.
    //
    // The delayed columns are prepended, so the old contents move down and right by n: old (i, j)
    // becomes new (i + n, j + n). Only the lower triangle is carried, which is all that is occupied
    // before a front is factored.
    //
    // **Called after mFrontSize has been widened**, like resetVal, so both verbs mean the same
    // thing: make the block match the shape the fields already describe. 0.9 calls its version
    // before widening and passes the old width implicitly; the arithmetic is identical either way,
    // and one convention across the two verbs is worth more than matching that call order.
    void expandVal(std::int32_t jj, std::size_t n) {
        const std::int32_t nInt     = static_cast<std::int32_t>(n);
        const std::int32_t newFront = static_cast<std::int32_t>(mFrontSize[jj]);
        const std::int32_t update   = static_cast<std::int32_t>(mUpdateSize[jj]);
        const std::int32_t oldFront = newFront - nInt;
        const std::int32_t oldRows  = oldFront + update;
        const std::int32_t newRows  = newFront + update;

        std::vector<Val> expanded(static_cast<std::size_t>(newFront) * static_cast<std::size_t>(newRows),
                               Val(0));
        const std::vector<Val>& old = mVal[jj];

        for (std::int32_t j_ = 0; j_ < oldFront; ++j_)
            for (std::int32_t i_ = j_; i_ < oldRows; ++i_)
                expanded[static_cast<std::size_t>(j_ + nInt) * static_cast<std::size_t>(newRows)
                      + static_cast<std::size_t>(i_ + nInt)]
                    = old[static_cast<std::size_t>(j_) * static_cast<std::size_t>(oldRows)
                          + static_cast<std::size_t>(i_)];

        mVal[jj] = std::move(expanded);
    }

    // Swap columns j_ and k_ of supernode jj, symmetrically. The block is a symmetric matrix stored
    // column-major with leading dimension the index size, so exchanging two pivot columns exchanges
    // the matching rows too. Also swaps the two node indices and repairs the global-to-local map.
    // Ported from 0.9 swap_.
    // **The middle loop conjugates for a Hermitian factorization, and that is not decoration.** It
    // exchanges a column entry with a row entry, and those two are reflections across the diagonal:
    // for a symmetric factor A(i,j) == A(j,i) and a raw swap is right, but for a Hermitian one
    // A(k,i) == conj(A(i,k)), so a value crossing the diagonal has to be conjugated on the way. The
    // other two loops move entries between columns or between rows without crossing, and the
    // diagonal entries are real, so this is the only place it arises. For `double` conj is the
    // identity and all three loops are plain swaps again.
    void swap(std::int32_t jj, std::int32_t j_, std::int32_t k_, std::vector<std::int32_t>& gblToLcl) {
        if (k_ < j_)
            std::swap(j_, k_);

        const bool withHermitian = hermitian(mFactorization);

        Val*               block = mVal[jj].data();
        std::int32_t*      idx   = mNodeIdx[jj].data();
        const std::int32_t ld    = static_cast<std::int32_t>(mFrontSize[jj] + mUpdateSize[jj]);

        // Column-major position of (row r, column c) in jj's block, in size_t to avoid overflow.
        const auto at = [ld](std::int32_t r, std::int32_t c) {
            return static_cast<std::size_t>(c) * static_cast<std::size_t>(ld) + static_cast<std::size_t>(r);
        };

        for (std::int32_t i_ = 0; i_ < j_; ++i_)            // rows j_ and k_, in the columns left of j_
            std::swap(block[at(j_, i_)], block[at(k_, i_)]);

        for (std::int32_t i_ = j_ + 1; i_ < k_; ++i_) {     // column j_ against row k_, between j_ and k_
            const Val a = block[at(i_, j_)];                // A(i, j)
            const Val b = block[at(k_, i_)];                // A(k, i)
            block[at(i_, j_)] = maybeConjugate(b, withHermitian);   // becomes A(i, k)
            block[at(k_, i_)] = maybeConjugate(a, withHermitian);   // becomes A(j, i)
        }

        for (std::int32_t i_ = k_ + 1; i_ < ld; ++i_)       // columns j_ and k_, in the rows below k_
            std::swap(block[at(i_, j_)], block[at(i_, k_)]);

        std::swap(block[at(j_, j_)], block[at(k_, k_)]);    // the two diagonal entries

        // **The entry between the two swapped positions is its own reflection**, and so is touched
        // by none of the loops above: under the permutation A(k, j) becomes A(j, k), which is the
        // same stored position. For a symmetric factor those are equal and there is nothing to do,
        // which is why 0.9 leaves it alone. For a Hermitian one they are conjugates, so it has to
        // be conjugated in place. Leaving it out produces a factor that reconstructs the conjugate
        // of the matrix in the affected rows, with no other symptom.
        block[at(k_, j_)] = maybeConjugate(block[at(k_, j_)], withHermitian);

        std::swap(idx[j_], idx[k_]);                        // the node indices and the global-to-local map
        gblToLcl[idx[j_]] = j_;
        gblToLcl[idx[k_]] = k_;
    }

    // Drop the n delayed columns from supernode jj's block, reclaiming their column storage. Called
    // once factorDynamicLDL_ has reduced frontSize[jj] by n and the delayed columns' values have
    // been assembled into the parent: the block loses its n trailing front columns and keeps all
    // its rows (frontSize + n + updateSize), so a column-major truncation does it. The delayed
    // columns live on only as rows, counted now in delaySize. Ported from 0.9
    // shrinkEntry_.
    void contractVal(std::int32_t jj, std::size_t n) {
        const std::size_t rows = mFrontSize[jj] + n + mUpdateSize[jj];
        mVal[jj].resize(mFrontSize[jj] * rows);
    }

    std::size_t   mSize      = 0;
    std::size_t   mSnodeSize = 0;
    Factorization mFactorization = Factorization::DynamicLDLT;

    // Copied from SymFactor, exactly as in NumFactorStatic. Under delayed pivoting the index sets
    // themselves expand, which is the second reason the factor owns a copy rather than referring
    // back: SymFactor's sets are the *predicted* ones and must not be disturbed.
    std::vector<std::int32_t> mNodeToSnode;
    std::vector<std::size_t>  mFrontSize;
    std::vector<std::size_t>  mUpdateSize;

    // Per supernode, filled by dynamic LDL: how many columns it delayed up to its parent.
    std::vector<std::size_t>  mDelaySize;

    // Per column (global): the pivot kind, 0 / 1x1 / the two halves of a 2x2.
    std::vector<std::int32_t> mPivotType;

    // The index sets and value blocks, one vector per supernode so a front can expand without
    // moving
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
