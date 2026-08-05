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

    // Three sizes of the factor, in the same order on all four classes that can answer. They count
    // different things and the difference grows with front width, so the names are not
    // interchangeable:
    //
    //   numNodeIdx  the index sets, one entry per row a supernode touches
    //   numVal      values *allocated*: the full rectangle, whose front block is stored whole with
    //               its strict upper triangle left zero so BLAS can take it in one call
    //   nnz         entries of L: each supernode's own lower triangle plus its update rectangle
    //
    // In magnitude they nest, numNodeIdx <= nnz <= numVal, with equality only when every
    // supernode is one column. 0.9 spelled the last two numberOfEntries and
    // numberOfAllocatedEntries.
    //
    // **nnz counts entries, not nonzeros.** Amalgamation pads index sets with structurally zero
    // rows and those are counted here; the two agree exactly when the patterns are exact.
    //
    // Summed over the per-supernode vectors, which is the only way here and also the most direct:
    // these are the lengths actually allocated rather than a formula reproducing them. They are
    // therefore the *delayed* sizes, larger than the symbolic factorization predicted wherever a
    // column was delayed into a front, which is what delaySize records. A supernode's block is
    // frontSize columns by frontSize + delaySize + updateSize rows.
    std::size_t numNodeIdx() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj) sum += mNodeIdx[jj].size();
        return sum;
    }
    std::size_t numVal() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj) sum += mVal[jj].size();
        return sum;
    }
    std::size_t nnz() const {
        std::size_t sum = 0;
        for (std::size_t jj = 0; jj < mSnodeSize; ++jj)
            sum += mFrontSize[jj] * (mFrontSize[jj] + 1) / 2
                 + mFrontSize[jj] * (mDelaySize[jj] + mUpdateSize[jj]);
        return sum;
    }

    // How many pivots the factorization had to replace. Meaningful when a *static* factorization
    // (Cholesky, static LDL) runs into this storage, exactly as on the static factor. Dynamic LDL,
    // when it lands, delays an unstable pivot instead of replacing it, and leaves this zero.
    std::size_t numPerturbations() const { return mNumPerturbations; }

    // Numerical rank, tracked as 0.9's rank_ is: it starts at the factor's size and each 1x1 pivot
    // accepted with an exactly zero diagonal drops it by one, so after factorization it is the
    // matrix dimension minus the count of zero pivots. Only dynamic LDL touches it (only dynamic
    // pivoting can accept a zero 1x1); the value is meaningful there and stays at full size on a
    // factorization that never hits one.
    std::size_t rank() const { return mRank; }

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

    // The write path for rank, mirroring numPerturbations: the engine decrements this reference at a
    // zero-diagonal 1x1. The const read overload above is public.
    std::size_t& rank() { return mRank; }

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

        for (std::int32_t j = 0; j < oldFront; ++j)
            for (std::int32_t i = j; i < oldRows; ++i)
                expanded[static_cast<std::size_t>(j + nInt) * static_cast<std::size_t>(newRows)
                      + static_cast<std::size_t>(i + nInt)]
                    = old[static_cast<std::size_t>(j) * static_cast<std::size_t>(oldRows)
                          + static_cast<std::size_t>(i)];

        mVal[jj] = std::move(expanded);
    }

    // Swap columns j and k of supernode jj, symmetrically. The block is a symmetric matrix stored
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
    void swap(std::int32_t jj, std::int32_t j, std::int32_t k, std::vector<std::int32_t>& gblToLcl) {
        if (k < j)
            std::swap(j, k);

        const bool withHermitian = hermitian(mFactorization);

        Val*               block = mVal[jj].data();
        std::int32_t*      idx   = mNodeIdx[jj].data();
        const std::int32_t ld    = static_cast<std::int32_t>(mFrontSize[jj] + mUpdateSize[jj]);

        // Column-major position of (row r, column c) in jj's block, in size_t to avoid overflow.
        const auto at = [ld](std::int32_t r, std::int32_t c) {
            return static_cast<std::size_t>(c) * static_cast<std::size_t>(ld) + static_cast<std::size_t>(r);
        };

        for (std::int32_t i = 0; i < j; ++i)                // rows j and k, in the columns left of j
            std::swap(block[at(j, i)], block[at(k, i)]);

        for (std::int32_t i = j + 1; i < k; ++i) {          // column j against row k, between j and k
            const Val a = block[at(i, j)];                  // A(i, j)
            const Val b = block[at(k, i)];                  // A(k, i)
            block[at(i, j)] = maybeConjugate(b, withHermitian);     // becomes A(i, k)
            block[at(k, i)] = maybeConjugate(a, withHermitian);     // becomes A(j, i)
        }

        for (std::int32_t i = k + 1; i < ld; ++i)           // columns j and k, in the rows below k
            std::swap(block[at(i, j)], block[at(i, k)]);

        std::swap(block[at(j, j)], block[at(k, k)]);        // the two diagonal entries

        // **The entry between the two swapped positions is its own reflection**, and so is touched
        // by none of the loops above: under the permutation A(k, j) becomes A(j, k), which is the
        // same stored position. For a symmetric factor those are equal and there is nothing to do,
        // which is why 0.9 leaves it alone. For a Hermitian one they are conjugates, so it has to
        // be conjugated in place. Leaving it out produces a factor that reconstructs the conjugate
        // of the matrix in the affected rows, with no other symptom.
        block[at(k, j)] = maybeConjugate(block[at(k, j)], withHermitian);

        std::swap(idx[j], idx[k]);                          // the node indices and the global-to-local map
        gblToLcl[idx[j]] = j;
        gblToLcl[idx[k]] = k;
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
    std::size_t mRank             = 0;   // set to size at setup, decremented per zero 1x1 pivot

    friend class NumFactorEngine;
};

extern template class NumFactorDynamic<double>;
extern template class NumFactorDynamic<std::complex<double>>;

} // namespace Oblio
