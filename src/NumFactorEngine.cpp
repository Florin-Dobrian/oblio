#include "oblio/NumFactorEngine.h"

#include "oblio/BlasLapack.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <type_traits>
#include <utility>

namespace Oblio {

namespace {

// The 2x2 pivot val D, read out of a front at local columns k1_ and k2_.
//
// **This is the one place the symmetry of D is decided, and that is the point of it existing.**
// `d12 = d21` is the symmetric statement, and it is what `LDL^T` means over the reals and over the
// complex field alike. A Hermitian factorization wants the conjugate here instead, and wants `d11`
// and `d22` known real. Everything downstream, the acceptance test and the elimination, works from
// what this returns, so that change is made once rather than in four places that must not drift.
template<class Val>
struct PivotBlock2x2 {
    Val d11, d22, d21, d12;
    Val det;
};

template<class Val>
PivotBlock2x2<Val> readPivotBlock2x2(const Val* val, std::ptrdiff_t ld,
                                     std::int32_t k1_, std::int32_t k2_, bool withHermitian) {
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    PivotBlock2x2<Val> d;
    d.d11 = forceReal(val[at(k1_, k1_)], withHermitian);
    d.d22 = forceReal(val[at(k2_, k2_)], withHermitian);

    // Only the lower triangle is occupied before the front is factored, so exactly one of the two
    // off-diagonal positions is stored, and **which one it is depends on the order of k1 and k2**.
    // That distinction is invisible in the symmetric case, where the two are equal, and matters
    // here: the stored entry is d21 when k2 sits below k1 and d12 when it sits above, and the other
    // is its conjugate.
    if (k2_ > k1_) {
        d.d21 = val[at(k2_, k1_)];
        d.d12 = maybeConjugate(d.d21, withHermitian);
    } else {
        d.d12 = val[at(k1_, k2_)];
        d.d21 = maybeConjugate(d.d12, withHermitian);
    }

    // Hermitian: d11 * d22 is real and d12 * d21 is |d21|^2, so the determinant is real, as the
    // Bunch-Kaufman test downstream assumes when it compares magnitudes against it.
    d.det = d.d11 * d.d22 - d.d12 * d.d21;
    return d;
}

} // namespace

// =================================================================================================
// Naming, as elsewhere. Supernodes are doubled letters: jj the supernode being factored, kk an
// ancestor it updates, and jj < kk since a descendant's label is below its ancestor's. Single
// letters are columns and rows: lk a column, li a row, in L's ordering; ak, ai the same in A's.
// Positions carry the initials of the pointer array they walk: cp into A's colPtr, sp into a
// supernode's node indices, vp into the factor's val. `nextUpdateSp[jj]` is one such sp, held
// across visits rather than recomputed: the position in jj's index set where its next update
// begins. (A supernode's indices live under snodeNodeIdxPtr on the factor and snodePtr on
// SymFactor; sp is the settled name for a position into them whichever class holds the array, and
// snip/svp are reserved for the day a position into the paired index and value arrays needs
// distinguishing.)
//
// Local coordinates are the exception worth naming. Inside a supernode's val, a row is
// identified by its *position in that supernode's index set*, not by its global row index. Those
// are the two things gblToLcl converts between, and the code says which it means: `li` is global,
// `lcl` is local.
// =================================================================================================

// The global-to-local map is a scratchpad the size of the matrix, allocated once per factorization
// with every entry NIL, and touched only at the positions a supernode names. Both directions cost
// O(|Idx(jj)|), never O(n), so the map is proportional to the supernode rather than to the matrix
// however sparse it is.
//
// **The clear is required, not hygiene.** assembleFromA reads NIL as an input check: every
// lower-triangle row of A's column must appear in the supernode's index set, because symbolic made
// the same cut, and a row that maps to NIL means A carries a nonzero the symbolic structure did not
// predict. Leave a stale local index behind and that check passes on a value from a previous
// supernode, putting the entry at a plausible offset in the wrong place, silently, on exactly the
// malformed input the check exists to catch.
//
// The escape, should the map ever show up in a profile, is a generation counter: store a stamp
// beside each local index and treat a stale stamp as NIL, trading the clear pass for a wider array
// and a comparison per lookup. Not worth it at this cost level.
void NumFactorEngine::setGlobalToLocal(std::size_t numNodeIdx, const std::int32_t* nodeIdx,
                                       std::vector<std::int32_t>& gblToLcl) const {
    for (std::size_t sp = 0; sp < numNodeIdx; ++sp)
        gblToLcl[nodeIdx[sp]] = static_cast<std::int32_t>(sp);
}

void NumFactorEngine::clearGlobalToLocal(std::size_t numNodeIdx, const std::int32_t* nodeIdx,
                                         std::vector<std::int32_t>& gblToLcl) const {
    for (std::size_t sp = 0; sp < numNodeIdx; ++sp)
        gblToLcl[nodeIdx[sp]] = NIL;
}

template<class Val>
void NumFactorEngine::setSymFactor(const SymFactor& sf, NumFactorStatic<Val>& nf) const {
    nf.mSize          = sf.size();
    nf.mSnodeSize     = sf.snodeSize();
    nf.mFactorization = mFactorization;

    // The structure, copied. The factor owns it, so SymFactor may be discarded afterwards, and so
    // dynamic LDL can expand its copy without disturbing the prediction.
    nf.mNodeToSnode     = sf.nodeToSnode();
    nf.mFrontSize       = sf.frontSize();
    nf.mUpdateSize      = sf.updateSize();
    nf.mNumNodeIdx      = sf.numNodeIdx();
    nf.mSnodeNodeIdxPtr = sf.snodePtr();
    nf.mNodeIdx         = sf.nodeIdx();

    // The value blocks. Supernode kk's is a dense column-major rectangle, indexSize rows by
    // frontSize columns, so it holds indexSize * frontSize values. Offsets accumulated the usual
    // way: an exclusive prefix sum, so snodeValPtr[kk] is where kk's val starts.
    nf.mSnodeValPtr.resize(nf.mSnodeSize + 1);
    nf.mSnodeValPtr[0] = 0;
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.mSnodeSize); ++kk) {
        const std::size_t numNodeIdx = nf.mFrontSize[kk] + nf.mUpdateSize[kk];
        nf.mSnodeValPtr[kk + 1] = nf.mSnodeValPtr[kk] + numNodeIdx * nf.mFrontSize[kk];
    }
    nf.mNumVal = nf.mSnodeValPtr[nf.mSnodeSize];

    // Zeroed, because assembly *adds* into it: A's original values first, then every descendant's
    // update.
    nf.mVal.assign(nf.mNumVal, Val(0));
}

template<class Val>
void NumFactorEngine::setSymFactor(const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    nf.mSize          = sf.size();
    nf.mSnodeSize     = sf.snodeSize();
    nf.mFactorization = mFactorization;

    // The structure, copied, exactly as for the static factor.
    nf.mNodeToSnode = sf.nodeToSnode();
    nf.mFrontSize   = sf.frontSize();
    nf.mUpdateSize  = sf.updateSize();

    // No columns delayed and no pivots chosen yet: dynamic LDL fills these as it runs, and a static
    // factorization into this storage leaves them untouched. pivotType is per column, the rest per
    // supernode.
    nf.mDelaySize.assign(nf.mSnodeSize, 0);
    nf.mPivotType.assign(nf.mSize, 0);

    // The index sets and value blocks, one vector per supernode so a front can later expand without
    // moving its neighbors. The index set is copied from SymFactor's flat buffer, sliced per
    // supernode; each val is indexSize rows by frontSize columns, zeroed because assembly adds
    // into it.
    const std::vector<std::int32_t>& sfNodeIdx  = sf.nodeIdx();
    const std::vector<std::size_t>&  sfSnodePtr = sf.snodePtr();
    nf.mNodeIdx.resize(nf.mSnodeSize);
    nf.mVal.resize(nf.mSnodeSize);
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.mSnodeSize); ++kk) {
        nf.mNodeIdx[kk].assign(sfNodeIdx.begin() + static_cast<std::ptrdiff_t>(sfSnodePtr[kk]),
                               sfNodeIdx.begin() + static_cast<std::ptrdiff_t>(sfSnodePtr[kk + 1]));

        const std::size_t numNodeIdx = nf.mFrontSize[kk] + nf.mUpdateSize[kk];
        nf.mVal[kk].assign(numNodeIdx * nf.mFrontSize[kk], Val(0));
    }
}

template<class Val>
bool NumFactorEngine::assembleFromA(const SparseMatrix<Val>& A, const Permutation& p,
                                    const std::vector<std::int32_t>& gblToLcl,
                                    std::size_t delaySize,
                                    std::size_t frontSize, std::size_t numNodeIdx,
                                    const std::int32_t* nodeIdx, Val* val) const {
    const std::vector<std::size_t>&  colPtr   = A.colPtr();
    const std::vector<std::int32_t>& aRowIdx  = A.rowIdx();
    const std::vector<Val>&          aVal     = A.val();
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // For each front column of the supernode that holds an entry of A. Its local column position
    // is lcl, and its val column starts at lcl * numNodeIdx (column-major). The run starts at
    // delaySize: zero for a static factorization, and under dynamic pivoting the
    // columns delayed into this front from its children, which A knows nothing about.
    for (std::size_t lcl = delaySize; lcl < frontSize; ++lcl) {
        const std::int32_t lk = nodeIdx[lcl];        // the global column, in L's ordering
        const std::int32_t ak = newToOld[lk];       // the same column, in A's

        for (std::size_t cp = colPtr[ak]; cp < colPtr[ak + 1]; ++cp) {
            const std::int32_t ai = aRowIdx[cp];
            const std::int32_t li = oldToNew[ai];   // the same row, in L's ordering

            // A is stored full-symmetric, so each entry appears twice. Take the lower triangle:
            // rows at or below the column. (The symbolic factorization made the same cut, which
            // is why the local map below is guaranteed to find these and only these.)
            if (li < lk)
                continue;

            const std::int32_t lclRow = gblToLcl[li];
            if (lclRow == NIL)
                return false;   // A has an entry the symbolic structure does not predict

            val[lcl * numNodeIdx + static_cast<std::size_t>(lclRow)] = aVal[cp];
        }
    }
    return true;
}

template<class Val>
void NumFactorEngine::assembleUpdate(const std::vector<std::int32_t>& gblToLcl,
                                     const UpdateBlock<Val>& updateBlock,
                                     std::size_t numNodeIdx, Val* val) const {
    // The update val's rows and columns carry global row indices; gblToLcl maps them into the
    // ancestor's local coordinates. Only the lower triangle of the val is meaningful (row at or
    // below column), which is exactly the part the two BLAS calls filled. Entry (si, sj) of the
    // update sits at updateBlock.mVal[sj * updateBlock.mHeight + si], the same column-major layout
    // the solve reads; (di, dj) is where that entry lands in the destination front. The global
    // (L-ordering) index in between is what gblToLcl converts.
    //
    // si, sj, di, dj are column/row indices (int32_t), as in the solve; the two positions here are
    // the flat column offsets into the source and destination vals, both size_t in the cp family.
    for (std::int32_t sj = 0; sj < static_cast<std::int32_t>(updateBlock.mWidth); ++sj) {
        const std::int32_t lj = updateBlock.mRowIdx[sj];   // this column's global (L-ordering) index
        const std::int32_t dj = gblToLcl[lj];              // its column in the destination front

        const std::size_t scp = static_cast<std::size_t>(sj) * updateBlock.mHeight;  // source column
        const std::size_t dcp = static_cast<std::size_t>(dj) * numNodeIdx;            // dest column

        for (std::int32_t si = sj; si < static_cast<std::int32_t>(updateBlock.mHeight); ++si) {
            const std::int32_t li = updateBlock.mRowIdx[si];
            const std::int32_t di = gblToLcl[li];

            val[dcp + di] += updateBlock.mVal[scp + si];
        }
    }
}

template<class Val>
bool NumFactorEngine::factorStaticSupernode(std::size_t frontSize, std::size_t numNodeIdx,
                                            Val* val,
                                      std::size_t& numPerturbations) const {
    const int f  = static_cast<int>(frontSize);
    const int u  = static_cast<int>(numNodeIdx - frontSize);
    const int ld = static_cast<int>(numNodeIdx);

    if (mFactorization == Factorization::Cholesky) {
        // The front, which is the diagonal val: A11 = L11 L11^H.
        int info = 0;
        potrf('L', f, val, ld, &info);
        if (info > 0)
            return false;   // not positive definite: the leading minor of order `info` failed

        // The update rows: L21 = A21 (L11^H)^-1, solved in place.
        //
        // Blas<Val>::conjTrans is 'T' for real and 'C' for complex, which is the whole of what the
        // scalar type decides here. 0.9 writes 'T' unconditionally, which is wrong for a complex
        // Hermitian factor, and there is nothing at its call site to reveal that.
        if (u > 0)
            trsm('R', 'L', Blas<Val>::conjTrans, 'N', u, f, Val(1), val, ld, val + f, ld);

        return true;
    }

    // LDL. The kernel is ours: LAPACK has no unpivoted LDL^T (its ?sytrf pivots, and pivoting is
    // what a *static* factorization refuses to do). It cannot fail, because there is no positive
    // definiteness to violate; a pivot too small to divide by is perturbed and counted.
    int numPert = 0;
    ldl(f, val, ld, mPerturbation, &numPert, hermitian(mFactorization));
    numPerturbations += static_cast<std::size_t>(numPert);

    // The update rows: L21 = A21 U11^-1, where U11 = D11 L11^H sits in the front's *upper*
    // triangle. So the solve is against the upper, untransposed, which is exactly what storing U
    // buys: Cholesky would have to transpose, and does.
    if (u > 0)
        trsm('R', 'U', 'N', 'N', u, f, Val(1), val, ld, val + f, ld);

    return true;
}

template<class Val>
void NumFactorEngine::updateStaticSupernode(std::size_t frontSize, std::size_t numNodeIdx,
                                            const Val* val,
                                      std::size_t offset, UpdateBlock<Val>& updateBlock) const {
    const int f      = static_cast<int>(frontSize);
    const int ld     = static_cast<int>(numNodeIdx);
    const int height = static_cast<int>(updateBlock.mHeight);
    const int width  = static_cast<int>(updateBlock.mWidth);
    const int tld    = height;

    Val*       tVal = updateBlock.mVal.data();
    const Val* L21  = val + offset;   // the update rows that reach this ancestor, and below

    if (mFactorization == Factorization::Cholesky) {
        // The square part: the block's (0..width, 0..width) -= L21' L21'^H, where L21' is the `width` rows
        // that land in the ancestor. Symmetric, so HERK, which touches only the lower triangle.
        //
        // `herk` means "A times A-conjugate-transpose": dsyrk_ for real, zherk_ for complex. The
        // engine never names either, which is what makes 0.9's bug (SYRK on a Hermitian factor)
        // impossible to write here.
        herk('L', 'N', width, f, -1.0, L21, ld, 1.0, tVal, tld);

        // The rectangle below: the block's (width.., 0..width) -= L21'' L21'^H, where L21'' is the rows of
        // the supernode's update val that lie *below* the ancestor's. Not symmetric, so GEMM.
        if (height > width)
            gemm('N', Blas<Val>::conjTrans, height - width, width, f,
                 Val(-1), L21 + width, ld,
                 L21, ld,
                 Val(1), tVal + width, tld);
        return;
    }

    // LDL. The update is `block -= L21 D L21^H`, and **no BLAS routine computes it**: the D in the
    // middle rules out a rank-k call, which is why Cholesky gets one and LDL does not.
    //
    // So form U := D L21'^H explicitly, into a scratch, and then two multiplies. The scratch is
    // f by width, and it is the whole price of the D.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(width), Val(0));
    formUpper(width, f, L21, ld, upper.data(), f, val, ld, hermitian(mFactorization));

    // The square part: symmetric, so only its lower triangle is filled. BLAS has nothing for this
    // either (syrk does A A^T, not A B with B known to make the product symmetric), so gemmLower
    // is ours as well.
    gemmLower(width, f, L21, ld, upper.data(), f, tVal, tld);

    // The rectangle below: not symmetric, so a plain GEMM. Note 'N','N': the transpose is already
    // baked into U.
    if (height > width)
        gemm('N', 'N', height - width, width, f,
             Val(-1), L21 + width, ld,
             upper.data(), f,
             Val(1), tVal + width, tld);
}

// =================================================================================================
// Left-looking. For each supernode: pull in every update to it from below, then factor.
//
// The bookkeeping is the interesting part. When jj is factored it must update every supernode that
// owns one of jj's update rows, and those ancestors are visited in increasing order. So jj is placed
// on the queue of the *next* ancestor it must update, and moves to the next queue each time it
// delivers one. nextUpdateSp[jj] is the sp -- the position in jj's index set -- of the first row of
// that next update; nodeToSnode of the row there names the ancestor.
//
// descendantUpdateQueue[kk] therefore holds, at the moment kk is reached, exactly the supernodes
// still queued to update kk. No search, no scan over descendants.
//
// **The relay is not an optimization over telling every ancestor at once**, and the reason it works
// is the elimination tree's absorption property: every update row of jj also appears in its parent's
// index set, so the ancestor jj is re-queued onto always has more work above it and the chain never
// breaks. jj does still update each of its ancestors individually; the relay only defers *finding*
// them, and it defers finding them without re-visiting anyone: each (descendant, ancestor) pair is
// popped exactly once, and a hop skips straight past ancestors jj has no rows in. What it buys is
// that one supernode sits on one queue at a time, so the record costs O(N) in the number of
// supernodes rather than O(pairs) -- peak storage, not operations, which are the same either way.
// "The life of an update" in docs/ARCHITECTURE.md has the full account, including why right-looking
// needs none of this.
// =================================================================================================

template<class Val, class Factor>
bool NumFactorEngine::factorStaticLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                        const SymFactor& sf, Factor& nf) const {
    setSymFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // Assemble A's original values into every supernode first. Cheaper than doing it inside the
    // main loop: the local map is set and cleared once per supernode either way, but this keeps
    // the traversal's own bookkeeping uncluttered.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        const bool ok = assembleFromA(A, p, gblToLcl, 0, kkFrontSize,
                                      kkNumNodeIdx, kkNodeIdx, kkVal);
        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        if (!ok)
            return false;
    }

    // For each kk, the descendants queued to update it; and how far each has got.
    std::vector<std::list<std::int32_t>> descendantUpdateQueue(snodeSize);
    std::vector<std::size_t>             nextUpdateSp(snodeSize, 0);   // sp of jj's next update row

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // Every descendant jj queued to update kk.
        while (!descendantUpdateQueue[kk].empty()) {
            const std::int32_t jj = descendantUpdateQueue[kk].front();
            descendantUpdateQueue[kk].pop_front();

            const std::size_t   jjFrontSize  = nf.frontSize(jj);
            const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
            const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
            const Val*          jjVal        = nf.val(jj);

            // How many of jj's remaining rows belong to kk. They are contiguous, because jj's
            // index set is sorted and the supernodes partition it in increasing order.
            const std::size_t jjUpdateSp = nextUpdateSp[jj];
            const std::size_t jjHeight       = jjNumNodeIdx - jjUpdateSp;
            std::size_t       jjWidth        = 0;
            while (jjUpdateSp + jjWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjUpdateSp + jjWidth]) == kk)
                ++jjWidth;

            UpdateBlock<Val> updateBlock(jjHeight, jjWidth);
            std::copy(jjNodeIdx + jjUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateStaticSupernode(jjFrontSize, jjNumNodeIdx, jjVal, jjUpdateSp, updateBlock);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);

            // jj has updated kk. Queue it for the next ancestor it must update.
            nextUpdateSp[jj] = jjUpdateSp + jjWidth;
            if (nextUpdateSp[jj] < jjNumNodeIdx)
                descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        }

        if (!factorStaticSupernode(kkFrontSize, kkNumNodeIdx, kkVal,
                                   nf.numPerturbations())) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;   // not positive definite (Cholesky only; LDL perturbs instead)
        }

        // kk is factored, so it now has updates of its own to deliver. Its front rows are its own
        // columns and update nobody; the first update row names the first ancestor it must update.
        nextUpdateSp[kk] = kkFrontSize;
        if (nextUpdateSp[kk] < kkNumNodeIdx)
            descendantUpdateQueue[nf.nodeToSnode(kkNodeIdx[nextUpdateSp[kk]])].push_back(kk);

        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
    }

    return true;
}

// =================================================================================================
// Right-looking. For each supernode: factor it, then push its update to every ancestor.
//
// The mirror of left-looking, and simpler: a supernode's ancestors are found by walking its own
// update rows, so no lists are needed and nothing has to be remembered between supernodes. The
// cost is that the local map is set and cleared once per (descendant, ancestor) pair rather than
// once per supernode.
// =================================================================================================

template<class Val, class Factor>
bool NumFactorEngine::factorStaticRightLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                               const SymFactor& sf, Factor& nf) const {
    setSymFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        const bool ok = assembleFromA(A, p, gblToLcl, 0, kkFrontSize,
                                      kkNumNodeIdx, kkNodeIdx, kkVal);
        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        if (!ok)
            return false;
    }

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        Val*                jjVal        = nf.val(jj);

        if (!factorStaticSupernode(jjFrontSize, jjNumNodeIdx, jjVal, nf.numPerturbations()))
            return false;   // not positive definite (Cholesky only; LDL perturbs instead)

        // Walk jj's update rows. Each run of them belonging to one ancestor is one update.
        std::size_t jjUpdateSp = jjFrontSize;
        while (jjUpdateSp < jjNumNodeIdx) {
            const std::int32_t kk = nf.nodeToSnode(jjNodeIdx[jjUpdateSp]);

            const std::size_t   kkFrontSize  = nf.frontSize(kk);
            const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
            const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
            Val*                kkVal        = nf.val(kk);

            const std::size_t jjHeight = jjNumNodeIdx - jjUpdateSp;
            std::size_t       jjWidth  = 0;
            while (jjUpdateSp + jjWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjUpdateSp + jjWidth]) == kk)
                ++jjWidth;

            UpdateBlock<Val> updateBlock(jjHeight, jjWidth);
            std::copy(jjNodeIdx + jjUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateStaticSupernode(jjFrontSize, jjNumNodeIdx, jjVal, jjUpdateSp, updateBlock);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjUpdateSp += jjWidth;
        }
    }

    return true;
}

// =================================================================================================
// The two pivot eliminations, each applied once a selection loop has accepted it.
//
// These were duplicated across the kernel's two passes until now, character for character: 0.9
// splits factorDynamicLDL_ on whether the supernode has update rows and writes both bodies out, and
// the port followed it. What the transcription showed is that the split is entirely in the
// *selection*: which candidates are eligible, which partner may be paired, and which acceptance
// test applies. Once a pivot is accepted the arithmetic is the same, so it lives here once and the
// two selection loops call it.
//
// The selection loops stay separate, deliberately. They are genuinely two algorithms rather than
// one with flags, and merging them behind parameters would save lines and cost the reader the
// ability to see what each pass actually does.
// =================================================================================================

template<class Val>
void NumFactorEngine::applyPivot1x1(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t j_,
                                    std::int32_t k1_, std::int32_t k1, std::int32_t jjFrontSize,
                                    std::int32_t rows,
                                    std::vector<std::int32_t>& gblToLcl) const {
    Val*                 val = nf.val(jj);
    const std::ptrdiff_t ld  = rows;

    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    const bool withHermitian = hermitian(nf.factorization());

    // Read before the swap, which is why it is not simply val[at(j_, j_)]: the swap is what puts
    // the pivot there. Forced real for a Hermitian factorization, and written back so the solve
    // divides by the same value the elimination used.
    const Val diagonal1 = forceReal(val[at(k1_, k1_)], withHermitian);

    if (j_ != k1_) nf.swap(jj, j_, k1_, gblToLcl);

    val[at(j_, j_)] = diagonal1;

    for (std::int32_t i_ = j_ + 1; i_ < rows; ++i_)               // L column: divide by the pivot
        val[at(i_, j_)] /= diagonal1;

    for (std::int32_t k_ = j_ + 1; k_ < jjFrontSize; ++k_)        // D L^H row, in the upper part
        val[at(j_, k_)] = diagonal1 * maybeConjugate(val[at(k_, j_)], withHermitian);

    for (std::int32_t k_ = j_ + 1; k_ < jjFrontSize; ++k_)        // rank-1 trailing update
        for (std::int32_t i_ = k_; i_ < rows; ++i_)
            val[at(i_, k_)] -= val[at(i_, j_)] * val[at(j_, k_)];

    nf.mPivotType[k1] = 1;
}

template<class Val>
void NumFactorEngine::applyPivot2x2(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t j_,
                                    std::int32_t k1_, std::int32_t k2_, std::int32_t k1,
                                    std::int32_t k2, std::int32_t jjFrontSize, std::int32_t rows,
                                    std::vector<std::int32_t>& gblToLcl) const {
    Val*                 val = nf.val(jj);
    const std::ptrdiff_t ld  = rows;

    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    const bool withHermitian = hermitian(nf.factorization());

    // Read before the swaps, as above. The caller read the same val to decide on it; nothing
    // touches the front in between.
    const PivotBlock2x2<Val> d = readPivotBlock2x2(val, ld, k1_, k2_, withHermitian);

    const std::int32_t j1_ = j_;
    const std::int32_t j2_ = j_ + 1;

    // Bring k1 and k2 to the front's next two columns. Which swaps are needed depends on where they
    // already are relative to each other, and doing them in the wrong order would undo one.
    if (k1_ < k2_) {
        if (!(j1_ == k1_ && j2_ == k2_)) {
            if (j1_ != k1_) nf.swap(jj, j1_, k1_, gblToLcl);
            nf.swap(jj, j2_, k2_, gblToLcl);
        }
    } else {
        if (j1_ == k2_ && j2_ == k1_) {
            nf.swap(jj, j1_, j2_, gblToLcl);
        } else {
            if (j2_ != k2_) nf.swap(jj, j2_, k2_, gblToLcl);
            nf.swap(jj, j1_, k1_, gblToLcl);
        }
    }

    // D's own four entries, written back where the solve expects them. The lower pair are already
    // in place, the swaps having carried them (conjugating on the way, for Hermitian); the diagonal
    // pair is rewritten because forceReal may have changed them, and the upper off-diagonal has no
    // other home.
    val[at(j1_, j1_)] = d.d11;
    val[at(j2_, j2_)] = d.d22;
    val[at(j1_, j2_)] = d.d12;

    for (std::int32_t i_ = j1_ + 2; i_ < rows; ++i_) {        // L columns: solve against D
        const Val t1 = val[at(i_, j1_)];
        const Val t2 = val[at(i_, j2_)];
        val[at(i_, j1_)] = (t1 * d.d22 - t2 * d.d21) / d.det;
        val[at(i_, j2_)] = (t2 * d.d11 - t1 * d.d12) / d.det;
    }

    for (std::int32_t k_ = j1_ + 2; k_ < jjFrontSize; ++k_) { // D L^H rows, in the upper part
        const Val l1 = maybeConjugate(val[at(k_, j1_)], withHermitian);
        const Val l2 = maybeConjugate(val[at(k_, j2_)], withHermitian);
        val[at(j1_, k_)] = d.d11 * l1 + d.d12 * l2;
        val[at(j2_, k_)] = d.d21 * l1 + d.d22 * l2;
    }

    for (std::int32_t k_ = j1_ + 2; k_ < jjFrontSize; ++k_)   // rank-2 trailing update
        for (std::int32_t i_ = k_; i_ < rows; ++i_)
            val[at(i_, k_)] -= val[at(i_, j1_)] * val[at(j1_, k_)];
    for (std::int32_t k_ = j2_ + 1; k_ < jjFrontSize; ++k_)
        for (std::int32_t i_ = k_; i_ < rows; ++i_)
            val[at(i_, k_)] -= val[at(i_, j2_)] * val[at(j2_, k_)];

    nf.mPivotType[k1] = 2;
    nf.mPivotType[k2] = 3;
}


template<class Val>
bool NumFactorEngine::factorDynamicSupernode(NumFactorDynamic<Val>& nf, std::int32_t jj,
                                       std::vector<std::int32_t>& gblToLcl) const {
    Val*          val = nf.val(jj);
    std::int32_t* idx = nf.mNodeIdx[jj].data();

    const std::int32_t   jjFrontSize = static_cast<std::int32_t>(nf.mFrontSize[jj]);
    const std::int32_t   rows        = jjFrontSize + static_cast<std::int32_t>(nf.mUpdateSize[jj]);
    const std::ptrdiff_t ld          = rows;
    const double         threshold   = mPivotThreshold;

    // Column-major position of (row r, column c), in ptrdiff_t to avoid overflow.
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    // The candidate pivot columns, by global index, in front order.
    std::list<std::int32_t> pivotList;
    for (std::int32_t j_ = 0; j_ < jjFrontSize; ++j_)
        pivotList.push_back(idx[j_]);

    std::int32_t j_ = 0;

    // Two passes, and 0.9 splits them on whether the supernode has update rows. Pass 1 is the
    // dense-front case, where a column that cannot pivot has nowhere to go; pass 2 is the general
    // one, where an ancestor is waiting. The bodies below say how they differ.
    //
    // **What remains below is selection only.** The eliminations themselves are applyPivot1x1 and
    // applyPivot2x2, shared by both passes, because once a pivot is accepted the arithmetic is
    // identical and 0.9's two copies of it were an artifact of writing the passes out separately.
    // The three real differences are the ones this loop and the next make visible:
    //
    //   No forced 1x1.  Pass 1 accepts the last remaining candidate whatever it looks like, since
    //                   a dense front has nowhere to delay it to. Pass 2 never does.
    //   Front partners. Pass 2's 2x2 partner scan stops at jjFrontSize, so a partner is always a
    //                   front column and never an update row.
    //   A real test.    Pass 1 accepts a 2x2 on max1 == max2, on the magnitudes alone, without ever
    //                   reading the 2x2 val. Pass 2 applies the Bunch-Kaufman test to its
    //                   determinant against the growth bound.
    if (nf.mUpdateSize[jj] == 0) {
        while (!pivotList.empty()) {
            bool         pivotFound = false;
            std::int32_t trials     = static_cast<std::int32_t>(pivotList.size());

            while (trials > 0) {
                const std::int32_t k1  = pivotList.front(); pivotList.pop_front();
                const std::int32_t k1_ = gblToLcl[k1];

                if (pivotList.empty()) {                    // only candidate left: a forced 1x1
                    pivotFound = true;
                    nf.mPivotType[k1] = 1;
                    ++j_;
                    break;
                }

                // max1 = k1's largest off-diagonal magnitude (its row on the left, its column
                // below);
                // k2_ the local row where it occurs.
                std::int32_t k2_  = -1;
                double       max1 = -1;
                for (std::int32_t i_ = j_; i_ < k1_; ++i_)
                    if (max1 < std::abs(val[at(k1_, i_)])) { k2_ = i_; max1 = std::abs(val[at(k1_, i_)]); }
                for (std::int32_t i_ = k1_ + 1; i_ < rows; ++i_)
                    if (max1 < std::abs(val[at(i_, k1_)])) { k2_ = i_; max1 = std::abs(val[at(i_, k1_)]); }

                const Val diagonal1 = val[at(k1_, k1_)];

                if (max1 == 0) {                            // isolated column: 1x1, nothing to eliminate
                    pivotFound = true;
                    if (j_ != k1_) nf.swap(jj, j_, k1_, gblToLcl);
                    nf.mPivotType[k1] = 1;
                    ++j_;
                    break;
                }

                if (std::abs(diagonal1) > 0 && std::abs(diagonal1) >= threshold * max1) {   // accept 1x1
                    pivotFound = true;
                    applyPivot1x1(nf, jj, j_, k1_, k1, jjFrontSize, rows, gblToLcl);
                    ++j_;
                    break;
                }
                else {                                      // try a 2x2 with k1 and its max partner k2
                    const std::int32_t k2 = idx[k2_];

                    double max2 = -1;
                    for (std::int32_t i_ = j_; i_ < k2_; ++i_)
                        if (max2 < std::abs(val[at(k2_, i_)])) max2 = std::abs(val[at(k2_, i_)]);
                    for (std::int32_t i_ = k2_ + 1; i_ < rows; ++i_)
                        if (max2 < std::abs(val[at(i_, k2_)])) max2 = std::abs(val[at(i_, k2_)]);

                    // Note what is *not* read here. Pass 1 decides on the magnitudes alone and
                    // never examines the 2x2 val itself; pass 2 tests its determinant. That
                    // difference was invisible while the two passes each carried their own copy of
                    // the elimination.
                    if (max1 == max2) {                     // accept 2x2
                        pivotFound = true;
                        pivotList.remove(k2);

                        applyPivot2x2(nf, jj, j_, k1_, k2_, k1, k2, jjFrontSize, rows, gblToLcl);
                        j_ += 2;
                        break;
                    }
                    else {                                  // neither k1 nor the 2x2 acceptable: delay k1
                        pivotList.push_back(k1);
                    }
                }

                --trials;
            }

            if (!pivotFound)
                break;
        }
    } else {
        // Pass 2: the supernode has update rows, so an ancestor exists to take a column this front
        // cannot pivot. Ported from 0.9 factorDynamicLDL_, the jjUpdateSize != 0 branch.
        //
        // **Not pass 1 with different bounds**, which is the guess worth naming so nobody makes it
        // twice. The two arithmetic bodies below, 1x1 and 2x2, are identical to pass 1's, character
        // for character. Everything that differs is in the *selection*:
        //
        //   No forced 1x1.  Pass 1 accepts the last remaining candidate whatever it looks like,
        //                   because there is nowhere to delay it to. Here there is, so the last
        //                   candidate falls through and is delayed like any other.
        //   Two scans.      max1 measures k1's largest off-diagonal over the whole column height,
        //                   update rows included. The partner scan `max` repeats it stopping at
        //                   jjFrontSize, so a 2x2 partner is always a front column and never an
        //                   update row. Hence max <= max1, and max == max1 says the largest entry
        //                   in k1's line is a front column after all.
        //   A real test.    Pass 1 accepts a 2x2 on max1 == max2. Here the test is the
        //                   Bunch-Kaufman one on the 2x2 determinant against the growth bound
        //                   maxmax, with the symmetric-maximum case kept as a separate disjunct.
        while (!pivotList.empty()) {
            bool         pivotFound = false;
            std::int32_t trials     = static_cast<std::int32_t>(pivotList.size());

            while (trials > 0) {
                const std::int32_t k1  = pivotList.front(); pivotList.pop_front();
                const std::int32_t k1_ = gblToLcl[k1];

                // k1's largest off-diagonal magnitude: its row on the left, its column below, all
                // the way down through the update rows. No argmax here, unlike pass 1; the partner
                // is chosen by the narrower scan further down.
                double max1 = -1;
                for (std::int32_t i_ = j_; i_ < k1_; ++i_)
                    if (max1 < std::abs(val[at(k1_, i_)])) max1 = std::abs(val[at(k1_, i_)]);
                for (std::int32_t i_ = k1_ + 1; i_ < rows; ++i_)
                    if (max1 < std::abs(val[at(i_, k1_)])) max1 = std::abs(val[at(i_, k1_)]);

                const Val diagonal1 = val[at(k1_, k1_)];

                if (max1 == 0) {                            // isolated column: 1x1, nothing to eliminate
                    pivotFound = true;
                    if (j_ != k1_) nf.swap(jj, j_, k1_, gblToLcl);
                    nf.mPivotType[k1] = 1;
                    ++j_;
                    break;
                }

                if (std::abs(diagonal1) > 0 && std::abs(diagonal1) >= threshold * max1) {   // accept 1x1
                    pivotFound = true;
                    applyPivot1x1(nf, jj, j_, k1_, k1, jjFrontSize, rows, gblToLcl);
                    ++j_;
                    break;
                }
                else if (!pivotList.empty()) {              // try a 2x2 with k1 and a front partner
                    // The same scan as max1's, stopped at the end of the front. k2_ is where the
                    // largest such entry sits, and it is a front column by construction.
                    std::int32_t k2_ = -1;
                    double       max = -1;
                    for (std::int32_t i_ = j_; i_ < k1_; ++i_)
                        if (max < std::abs(val[at(k1_, i_)])) { k2_ = i_; max = std::abs(val[at(k1_, i_)]); }
                    for (std::int32_t i_ = k1_ + 1; i_ < jjFrontSize; ++i_)
                        if (max < std::abs(val[at(i_, k1_)])) { k2_ = i_; max = std::abs(val[at(i_, k1_)]); }

                    const std::int32_t k2 = idx[k2_];

                    double max2 = -1;
                    for (std::int32_t i_ = j_; i_ < k2_; ++i_)
                        if (max2 < std::abs(val[at(k2_, i_)])) max2 = std::abs(val[at(k2_, i_)]);
                    for (std::int32_t i_ = k2_ + 1; i_ < rows; ++i_)
                        if (max2 < std::abs(val[at(i_, k2_)])) max2 = std::abs(val[at(i_, k2_)]);

                    const PivotBlock2x2<Val> d =
                        readPivotBlock2x2(val, ld, k1_, k2_, hermitian(nf.factorization()));

                    // The growth bound the 2x2 determinant is tested against: the larger of the two
                    // ways of pairing each diagonal with the other column's maximum.
                    const double maxmax = std::max(std::abs(d.d22) * max1 + max * max2,
                                                   std::abs(d.d11) * max2 + max * max1);

                    if ((max == max1 && max == max2 && max != 0)
                        || (std::abs(d.det) > 0 && std::abs(d.det) >= threshold * maxmax)) {   // accept 2x2
                        pivotFound = true;
                        pivotList.remove(k2);

                        applyPivot2x2(nf, jj, j_, k1_, k2_, k1, k2, jjFrontSize, rows, gblToLcl);
                        j_ += 2;
                        break;
                    }
                    else {                                  // neither k1 nor the 2x2 acceptable: delay k1
                        pivotList.push_back(k1);
                    }
                }
                else {                                      // no partner available: delay k1
                    pivotList.push_back(k1);
                }

                --trials;
            }

            if (!pivotFound)
                break;
        }
    }

    // Whatever is left could not be pivoted here; delay it to an ancestor. frontSize contracts by
    // that
    // count, and the val height (frontSize + delaySize + updateSize) is preserved.
    const std::int32_t delayed = static_cast<std::int32_t>(pivotList.size());
    nf.mDelaySize[jj] = static_cast<std::size_t>(delayed);
    nf.mFrontSize[jj] -= static_cast<std::size_t>(delayed);

    return true;
}

template<class Val>
void NumFactorEngine::updateDynamicSupernode(const NumFactorDynamic<Val>& nf, std::int32_t jj,
                                             std::size_t offset, UpdateBlock<Val>& updateBlock) const {
    const int f = static_cast<int>(nf.mFrontSize[jj]);
    if (f == 0)
        return;   // every column of jj was delayed: there is no pivot here to update anyone with

    // The height, and the one number that differs from the static twin, where it is just
    // frontSize + updateSize. A delayed column keeps its row, so the stride still counts it.
    const int ld     = f + static_cast<int>(nf.mDelaySize[jj])
                         + static_cast<int>(nf.mUpdateSize[jj]);
    const int height = static_cast<int>(updateBlock.mHeight);
    const int width  = static_cast<int>(updateBlock.mWidth);
    const int tld    = height;

    const bool          withHermitian = hermitian(nf.factorization());
    const Val*          val           = nf.val(jj);
    const std::int32_t* idx           = nf.mNodeIdx[jj].data();
    // jj's rows from `offset` down: the block this ancestor is about to receive.
    const Val*          L21           = val + offset;
    Val*                tVal          = updateBlock.mVal.data();

    // Column-major positions: `at` into jj's val (and equally into L21, which shares its leading
    // dimension), `atU` into the scratch.
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };
    const auto atU = [f](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * f + static_cast<std::ptrdiff_t>(r);
    };

    // U := D L21^H, f by width, the conjugate being the identity for a symmetric factorization. The
    // static twin gets this from formUpper in one call; here D is val-diagonal, so the front is
    // walked a pivot at a time.
    //
    // For a 2x2 the four entries of D sit where the factorization left them: the 2x2 elimination
    // starts at row j1 + 2 and never touches its own corner, so the lower pair are the original
    // matrix entries and the upper one was written back explicitly.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(width), Val(0));

    for (std::int32_t j_ = 0; j_ < f; ) {
        if (nf.mPivotType[idx[j_]] == 1) {
            const Val d = val[at(j_, j_)];
            for (std::int32_t tc = 0; tc < width; ++tc)
                upper[atU(j_, tc)] = d * maybeConjugate(L21[at(tc, j_)], withHermitian);
            ++j_;
        } else {                                   // a 2x2 pivot: two columns solved together
            const std::int32_t j1_ = j_;
            const std::int32_t j2_ = j_ + 1;

            const Val d11 = val[at(j1_, j1_)];
            const Val d12 = val[at(j1_, j2_)];   // the upper part, written back by the pivot
            const Val d21 = val[at(j2_, j1_)];   // the lower part, never overwritten
            const Val d22 = val[at(j2_, j2_)];

            for (std::int32_t tc = 0; tc < width; ++tc) {
                const Val l1 = maybeConjugate(L21[at(tc, j1_)], withHermitian);
                const Val l2 = maybeConjugate(L21[at(tc, j2_)], withHermitian);
                upper[atU(j1_, tc)] = d11 * l1 + d12 * l2;
                upper[atU(j2_, tc)] = d21 * l1 + d22 * l2;
            }
            j_ += 2;
        }
    }

    // From here it is the static twin exactly: the square part is symmetric and gets gemmLower,
    // the rectangle below it is not and gets a plain GEMM, with the transpose already baked into U.
    gemmLower(width, f, L21, ld, upper.data(), f, tVal, tld);

    if (height > width)
        gemm('N', 'N', height - width, width, f,
             Val(-1), L21 + width, ld,
             upper.data(), f,
             Val(1), tVal + width, tld);
}

template<class Val>
void NumFactorEngine::assembleDelayed(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t kk,
                                      const std::vector<std::int32_t>& gblToLcl) const {
    const std::int32_t jjFrontSize = static_cast<std::int32_t>(nf.mFrontSize[jj]);
    const std::int32_t jjDelayed   = static_cast<std::int32_t>(nf.mDelaySize[jj]);
    const std::int32_t jjRows      = jjFrontSize + jjDelayed
                                   + static_cast<std::int32_t>(nf.mUpdateSize[jj]);
    const std::int32_t kkRows      = static_cast<std::int32_t>(nf.mFrontSize[kk] + nf.mUpdateSize[kk]);

    const std::int32_t* jjNodeIdx = nf.mNodeIdx[jj].data();
    const Val*          jjVal     = nf.val(jj);
    Val*                kkVal     = nf.val(kk);

    const auto atJj = [jjRows](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * jjRows + static_cast<std::ptrdiff_t>(r);
    };
    const auto atKk = [kkRows](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * kkRows + static_cast<std::ptrdiff_t>(r);
    };

    // The delayed columns are the run just past jj's new front, and each carries its rows from the
    // diagonal down. Every one of those rows is in kk's index set: the delayed columns because kk
    // was expanded to hold exactly them, and the update rows below by the property that makes the
    // elimination forest work, a supernode's update indices lying inside its parent's index set.
    // So gblToLcl finds all of them, and none is NIL.
    for (std::int32_t sj_ = jjFrontSize; sj_ < jjFrontSize + jjDelayed; ++sj_) {
        const std::int32_t dj_ = gblToLcl[jjNodeIdx[sj_]];

        for (std::int32_t si_ = sj_; si_ < jjRows; ++si_) {
            const std::int32_t di_ = gblToLcl[jjNodeIdx[si_]];

            kkVal[atKk(di_, dj_)] = jjVal[atJj(si_, sj_)];
        }
    }
}

// =================================================================================================
// Dynamic LDL, left-looking. The same shape as the static traversal above, pull every update in
// then factor, with three additions that all follow from one fact: a column that cannot be pivoted
// where it stands is passed up to the parent, so a front's width is no longer what symbolic
// predicted.
//
// Per supernode kk, in order, and the order is load bearing twice:
//
//   Grow.       Sum what kk's children delayed. If nonzero, extend kk's index set, shift its own
//               indices right to make room, prepend the children's delayed globals, widen the
//               front, and discard-and-rezero the val.
//   Assemble A. Starting past the prepended columns, which hold no entry of A.
//   Update.     For each jj owing kk: if kk is jj's parent, fold jj's delayed columns in and only
//               then contract them away. Then the ordinary update.
//   Factor.     Which may itself delay, reducing frontSize[kk] and setting
//               delaySize[kk].
//
// The two orderings that matter: the delayed columns must be assembled into the parent *before*
// contractVal drops them, and kk must be expanded *before* A is assembled into it, since the offset
// assumes the wider front.
//
// **The height is conserved throughout.** frontSize + delaySize + updateSize is the
// val's row count and never changes: expanding moves rows from nowhere into the front, factoring
// moves them from the front into delayed, and updateSize is never rewritten. When a residual comes
// out wrong, that identity is the first thing to check.
// =================================================================================================

template<class Val>
bool NumFactorEngine::factorDynamicLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                               const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    setSymFactor(sf, nf);

    const std::int32_t snodeSize = static_cast<std::int32_t>(nf.snodeSize());

    const std::vector<std::int32_t>& parent      = sf.parent();
    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(nf.size(), NIL);

    // For each kk, the descendants queued to update it, and how far each has got, exactly as in
    // the static traversal. A is not
    // assembled up front here as it is there: kk's val does not reach its final width until its
    // children have been factored, so the assemble has to happen inside the loop.
    std::vector<std::list<std::int32_t>> descendantUpdateQueue(nf.snodeSize());
    std::vector<std::size_t>             nextUpdateSp(nf.snodeSize(), 0);   // sp of jj's next update row

    for (std::int32_t kk = 0; kk < snodeSize; ++kk) {
        std::int32_t delayedIntoKk = 0;
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj])
            delayedIntoKk += static_cast<std::int32_t>(nf.mDelaySize[jj]);

        if (delayedIntoKk > 0) {
            // The height after expanding, which is also the height before it: the new columns come
            // from rows the val already had. Only the front/update split moves.
            const std::int32_t expandedRows =
                static_cast<std::int32_t>(nf.mFrontSize[kk] + nf.mUpdateSize[kk]) + delayedIntoKk;

            // The index set. Extend first, then shift kk's own indices right by the delayed count,
            // descending so the copy does not overwrite its own source.
            nf.expandNodeIdx(kk, delayedIntoKk);
            std::vector<std::int32_t>& kkNodeIdxVec = nf.mNodeIdx[kk];

            for (std::int32_t sij_ = expandedRows - delayedIntoKk - 1, dij_ = expandedRows - 1;
                 sij_ >= 0; --sij_, --dij_)
                kkNodeIdxVec[dij_] = kkNodeIdxVec[sij_];

            // Then the vacated slots at the left, filled from the children in sibling order. Each
            // child's delayed columns are the run just past its (already reduced) front.
            std::int32_t dij_ = 0;
            for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
                const std::int32_t  jjFrontSize = static_cast<std::int32_t>(nf.mFrontSize[jj]);
                const std::int32_t* jjNodeIdx   = nf.nodeIdx(jj);

                for (std::int32_t sij_ = jjFrontSize;
                     sij_ < jjFrontSize + static_cast<std::int32_t>(nf.mDelaySize[jj]);
                     ++sij_, ++dij_)
                    kkNodeIdxVec[dij_] = jjNodeIdx[sij_];
            }

            // And the val. The front is wider, so the old contents are discarded rather than
            // moved: nothing has been written into kk yet, which is why left-looking never needs
            // 0.9's extendEntry_.
            nf.mFrontSize[kk] += static_cast<std::size_t>(delayedIntoKk);
            nf.resetVal(kk);
        }

        // The full height, captured before the factorization reclassifies part of the front as
        // delayed. Every use below wants this number and not the contracted front.
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        if (!assembleFromA(A, p, gblToLcl, static_cast<std::size_t>(delayedIntoKk),
                           kkFrontSize, kkNumNodeIdx, kkNodeIdx, nf.val(kk))) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Every descendant jj queued to update kk.
        while (!descendantUpdateQueue[kk].empty()) {
            const std::int32_t jj = descendantUpdateQueue[kk].front();
            descendantUpdateQueue[kk].pop_front();

            const std::int32_t  jjDelayed    = static_cast<std::int32_t>(nf.mDelaySize[jj]);
            const std::size_t   jjNumNodeIdx = nf.frontSize(jj)
                                             + static_cast<std::size_t>(jjDelayed)
                                             + nf.updateSize(jj);
            const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);

            // How many of jj's remaining rows belong to kk. Contiguous, as in the static case.
            const std::size_t jjUpdateSp = nextUpdateSp[jj];
            const std::size_t jjHeight       = jjNumNodeIdx - jjUpdateSp;
            std::size_t       jjWidth        = 0;
            while (jjUpdateSp + jjWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjUpdateSp + jjWidth]) == kk)
                ++jjWidth;

            // jj's delayed columns go to its parent and nowhere else, so this fires at most once
            // per jj across the whole traversal, and only if jj delayed something.
            //
            // 0.9 tests only the parent, letting both calls run as no-ops when jj delayed nothing:
            // assembleDelayed loops zero times and contractVal resizes a val to the size it
            // already has. Correct, but it makes the verb fire on every non-root supernode, which
            // is misleading under measurement and asymmetric with the right-looking driver, where
            // the same calls sit under a delay test already.
            if (jjDelayed > 0 && parent[jj] == kk) {
                assembleDelayed(nf, jj, kk, gblToLcl);
                nf.contractVal(jj, jjDelayed);
            }

            UpdateBlock<Val> updateBlock(jjHeight, jjWidth);
            std::copy(jjNodeIdx + jjUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateDynamicSupernode(nf, jj, jjUpdateSp, updateBlock);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, nf.val(kk));

            nextUpdateSp[jj] = jjUpdateSp + jjWidth;
            if (nextUpdateSp[jj] < jjNumNodeIdx)
                descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        }

        if (!factorDynamicSupernode(nf, kk, gblToLcl)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // kk now has updates of its own to deliver. **The advance is over the front and the delayed columns
        // together**, not the front alone: both are kk's own rows, neither updates an ancestor, and
        // the delayed ones are handed over by assembleDelayed rather than as an update. Getting
        // this wrong sends the delayed rows into a temporary and corrupts the parent quietly.
        nextUpdateSp[kk] = nf.frontSize(kk) + nf.delaySize(kk);
        if (nextUpdateSp[kk] < kkNumNodeIdx)
            descendantUpdateQueue[nf.nodeToSnode(kkNodeIdx[nextUpdateSp[kk]])].push_back(kk);

        // A root has no parent to take a delayed column, so a delay there is unrecoverable. 0.9
        // treats it as an error rather than a numeric failure, and so do we: it means the pivoting
        // strategy did not do its job, not that the matrix is singular.
        const bool delayedAtRoot = nf.delaySize(kk) > 0 && parent[kk] == NIL;

        // Clear the whole height, delayed columns included. 0.9 clears only frontSize + updateSize
        // here, and since frontSize has just contracted it leaves the delayed entries stale.
        // Harmless
        // there, because the map is only ever read at indices known to be in the current
        // supernode's set, but it costs the array its stated invariant (NIL everywhere outside the
        // supernode in hand), and assembleFromA's NIL test is exactly a check that relies on it.
        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        if (delayedAtRoot)
            return false;
    }

    return true;
}


// =================================================================================================
// Dynamic LDL, right-looking. Factor a supernode, then push its update to every ancestor: the
// mirror of the traversal above, and it uses the same two kernels unchanged, because 0.9's are
// byte-identical between its two engines.
//
// One thing genuinely differs, and it is the reason expandVal exists. A is assembled into every
// front here *before* the traversal starts, and ancestors accumulate updates as their descendants
// are factored, so a front that expands is never empty. It must keep what it holds while the
// delayed
// columns are inserted at its left, which is expandVal; left-looking, whose fronts are still
// empty when they expand, calls resetVal instead.
//
// The second difference is smaller and follows from the direction. Left-looking folds a child's
// delayed columns into the parent when it reaches the parent's update list; here the parent takes
// them from all its children at once, at the moment it expands, since by then every child is
// finished. The order that matters is unchanged: assemble the delayed columns, then contract them
// away.
//
// **The cost trade is the same one the static pair makes**, and neither side is free. Both do work
// proportional to the number of (descendant, ancestor) pairs; they differ only in what they spend
// per pair. This traversal sets and clears the global-to-local map, O(|Idx(kk)|) but a contiguous
// sweep over an array the assembly is about to touch anyway. Left-looking spends a list node,
// allocated and freed, plus a position write: O(1), but an allocation and a pointer chase.
//
// So it is O(1)-with-an-allocation against O(|Idx|)-with-locality, and which wins depends on how
// large the index sets are. Empirical and unmeasured; see docs/TODO.md, which also notes that 0.9
// pooled its queue nodes where this port uses std::list, so left-looking's side of the trade is
// probably worse here than in the reference.
// =================================================================================================

template<class Val>
bool NumFactorEngine::factorDynamicRightLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                                const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    setSymFactor(sf, nf);

    const std::int32_t snodeSize = static_cast<std::int32_t>(nf.snodeSize());

    const std::vector<std::int32_t>& parent      = sf.parent();
    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(nf.size(), NIL);

    // A first, into every front, while every front is still the width symbolic predicted. Nothing
    // has expanded yet, so the delayed-column offset is zero everywhere.
    for (std::int32_t kk = 0; kk < snodeSize; ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        const bool ok = assembleFromA(A, p, gblToLcl, 0, kkFrontSize,
                                      kkNumNodeIdx, kkNodeIdx, nf.val(kk));
        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
        if (!ok)
            return false;
    }

    for (std::int32_t jj = 0; jj < snodeSize; ++jj) {
        std::int32_t delayedIntoJj = 0;
        for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii])
            delayedIntoJj += static_cast<std::int32_t>(nf.mDelaySize[ii]);

        if (delayedIntoJj > 0) {
            const std::int32_t expandedRows =
                static_cast<std::int32_t>(nf.mFrontSize[jj] + nf.mUpdateSize[jj]) + delayedIntoJj;

            // The index set, exactly as in the left-looking driver: extend, shift right, prepend
            // the children's delayed globals in sibling order.
            nf.expandNodeIdx(jj, delayedIntoJj);
            std::vector<std::int32_t>& jjNodeIdxVec = nf.mNodeIdx[jj];

            for (std::int32_t sij_ = expandedRows - delayedIntoJj - 1, dij_ = expandedRows - 1;
                 sij_ >= 0; --sij_, --dij_)
                jjNodeIdxVec[dij_] = jjNodeIdxVec[sij_];

            std::int32_t dij_ = 0;
            for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii]) {
                const std::int32_t  iiFrontSize = static_cast<std::int32_t>(nf.mFrontSize[ii]);
                const std::int32_t* iiNodeIdx   = nf.nodeIdx(ii);

                for (std::int32_t sij_ = iiFrontSize;
                     sij_ < iiFrontSize + static_cast<std::int32_t>(nf.mDelaySize[ii]);
                     ++sij_, ++dij_)
                    jjNodeIdxVec[dij_] = iiNodeIdx[sij_];
            }

            // And the val, keeping what A and the descendants already put there.
            nf.mFrontSize[jj] += static_cast<std::size_t>(delayedIntoJj);
            nf.expandVal(jj, delayedIntoJj);
        }

        // The full height, captured before the factorization reclassifies part of the front as
        // delayed. The map is computed once and serves both the delayed assembly and the factor,
        // since jj's index set does not change between them.
        const std::size_t   jjNumNodeIdx = nf.frontSize(jj) + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);

        setGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        // Take every child's delayed columns, then let the child reclaim their storage. Tested per
        // child, not just on the total: one child delaying does not mean its siblings did, and the
        // calls are no-ops for those that did not.
        if (delayedIntoJj > 0)
            for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii])
                if (nf.mDelaySize[ii] > 0) {
                    assembleDelayed(nf, ii, jj, gblToLcl);
                    nf.contractVal(ii, static_cast<std::int32_t>(nf.mDelaySize[ii]));
                }

        if (!factorDynamicSupernode(nf, jj, gblToLcl)) {
            clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
            return false;
        }

        clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        // Push. The walk starts past the front *and* the delayed columns, for the same reason the
        // left-looking nextUpdateSp seed does: both are jj's own rows and neither updates an
        // ancestor, the delayed ones going up by assembleDelayed instead.
        std::size_t jjUpdateSp = nf.frontSize(jj) + nf.delaySize(jj);
        while (jjUpdateSp < jjNumNodeIdx) {
            const std::int32_t kk = nf.nodeToSnode(jjNodeIdx[jjUpdateSp]);

            const std::size_t jjHeight = jjNumNodeIdx - jjUpdateSp;
            std::size_t       jjWidth  = 0;
            while (jjUpdateSp + jjWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjUpdateSp + jjWidth]) == kk)
                ++jjWidth;

            UpdateBlock<Val> updateBlock(jjHeight, jjWidth);
            std::copy(jjNodeIdx + jjUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateDynamicSupernode(nf, jj, jjUpdateSp, updateBlock);

            // kk has not expanded yet, and need not have: jj's update rows are kk's own nodes, which
            // its index set already holds. When kk later expands, expandVal carries these values
            // along with the rest.
            const std::size_t   kkFrontSize  = nf.frontSize(kk);
            const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
            const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, nf.val(kk));
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjUpdateSp += jjWidth;
        }

        if (nf.delaySize(jj) > 0 && parent[jj] == NIL)
            return false;   // delayed at a root: nowhere left to put it
    }

    return true;
}

template<class Val>
bool NumFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                              NumFactorStatic<Val>& nf) const {
    if (A.size() != p.size() || A.size() != sf.size())
        return false;

    // **Dynamic pivoting cannot go into this storage, and never will.** A delayed column expands
    // its
    // parent's front, and this factor's value buffer is one flat array sized once from the symbolic
    // factorization; expanding a front in the middle of it would mean moving everything after it.
    // So
    // this is a design refusal rather than a missing feature, and the combination is asserted to be
    // refused in test_pipeline. Callers wanting a dynamic factorization pass NumFactorDynamic,
    // which
    // the overload below takes.
    switch (mFactorization) {
        case Factorization::Cholesky:
        case Factorization::StaticLDLT:
        case Factorization::StaticLDLH:
            break;
        case Factorization::DynamicLDLT:
        case Factorization::DynamicLDLH:
            return false;   // by design, see above; not a gap to be filled
    }

    switch (mTraversal) {
        case Traversal::LeftLooking:  return factorStaticLeftLooking(A, p, sf, nf);
        case Traversal::RightLooking: return factorStaticRightLooking(A, p, sf, nf);
        case Traversal::Multifrontal: return false;   // not implemented
    }
    return false;
}

template bool NumFactorEngine::compute(const SparseMatrix<double>&, const Permutation&,
                                       const SymFactor&, NumFactorStatic<double>&) const;
template bool NumFactorEngine::compute(const SparseMatrix<std::complex<double>>&,
                                       const Permutation&, const SymFactor&,
                                       NumFactorStatic<std::complex<double>>&) const;

template<class Val>
bool NumFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                              NumFactorDynamic<Val>& nf) const {
    if (A.size() != p.size() || A.size() != sf.size())
        return false;

    // Dynamic LDL is the reason this storage exists. Everything runs except complex `LDL^H`;
    // multifrontal is unported for all of them. The static factorizations run unchanged, below.
    //
    // **The kernels needed nothing to become complex.** 0.9's complex `factorDynamicLDL_` differs
    // from its real one in six lines, all the same edit: the pivot *magnitudes* (`max1`, `max2`,
    // `maxmax`) are declared real rather than scalar. This port declared them `double` from the
    // start, so it was already the complex form; `updateDynamicLDL_` is byte-identical between
    // 0.9's two engines to begin with. Everything else is `Val` arithmetic and `std::abs`, which
    // means modulus for complex and is the right comparison either way.
    //
    // **Complex `LDL^H` is the one cell still missing, and it is an extension, not a port.** 0.9's
    // complex LDL is symmetric only, so there is nothing to transcribe. `factorDynamicSupernode`
    // writes `diagonal12 = diagonal21`, which is exactly the complex-symmetric assumption; the
    // Hermitian form needs the conjugate there and a Bunch-Kaufman test that knows the 2x2's
    // diagonal is real. Refused here rather than silently computing the symmetric factorization
    // under a Hermitian name. See docs/TODO.md.
    //
    // Over the reals the question does not arise: the two transposes are the same computation, and
    // `test_pipeline` asserts they agree bit for bit.
    if (dynamicPivoting(mFactorization)) {
        switch (mTraversal) {
            case Traversal::LeftLooking:  return factorDynamicLeftLooking(A, p, sf, nf);
            case Traversal::RightLooking: return factorDynamicRightLooking(A, p, sf, nf);
            case Traversal::Multifrontal: return false;   // not ported yet
        }
        return false;
    }

    // The static factorizations, which this storage holds too. Every enumerator is named rather
    // than collected under a default: the two dynamic ones are unreachable here, the guard above
    // having returned for them, but the compiler cannot know that, and naming them is what keeps
    // -Wswitch pointing at this switch if a sixth factorization is ever added.
    switch (mFactorization) {
        case Factorization::Cholesky:
        case Factorization::StaticLDLT:
        case Factorization::StaticLDLH:
            break;
        case Factorization::DynamicLDLT:
        case Factorization::DynamicLDLH:
            return false;   // unreachable: handled by the guard above
    }

    switch (mTraversal) {
        case Traversal::LeftLooking:  return factorStaticLeftLooking(A, p, sf, nf);
        case Traversal::RightLooking: return factorStaticRightLooking(A, p, sf, nf);
        case Traversal::Multifrontal: return false;   // not ported yet
    }
    return false;
}

template bool NumFactorEngine::compute(const SparseMatrix<double>&, const Permutation&,
                                       const SymFactor&, NumFactorDynamic<double>&) const;
template bool NumFactorEngine::compute(const SparseMatrix<std::complex<double>>&,
                                       const Permutation&, const SymFactor&,
                                       NumFactorDynamic<std::complex<double>>&) const;

} // namespace Oblio
