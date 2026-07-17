#include "oblio/NumFactorEngine.h"

#include "oblio/BlasLapack.h"

#include <algorithm>
#include <cassert>
#include <list>
#include <utility>

namespace Oblio {

// =================================================================================================
// Naming, as elsewhere. Supernodes are doubled letters: jj the supernode being factored, kk an
// ancestor it updates, and jj < kk since a descendant's label is below its ancestor's. Single
// letters are columns and rows: lk a column, li a row, in L's ordering; ak, ai the same in A's.
// Positions carry the initials of the pointer array they walk: cp into A's colPtr, sp into
// SymFactor's supPtr, vp into the factor's valPtr.
//
// Local coordinates are the exception worth naming. Inside a supernode's block, a row is
// identified by its *position in that supernode's index set*, not by its global row index. Those
// are the two things gblToLcl converts between, and the code says which it means: `li` is global,
// `lcl` is local.
// =================================================================================================

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
    // dynamic LDL can grow its copy without disturbing the prediction.
    nf.mNodeToSnode     = sf.nodeToSnode();
    nf.mFrontSize       = sf.frontSize();
    nf.mUpdateSize      = sf.updateSize();
    nf.mNumNodeIdx      = sf.numNodeIdx();
    nf.mSnodeNodeIdxPtr = sf.snodePtr();
    nf.mNodeIdx         = sf.nodeIdx();

    // The value blocks. Supernode kk's is a dense column-major rectangle, indexSize rows by
    // frontSize columns, so it holds indexSize * frontSize values. Offsets accumulated the usual
    // way: an exclusive prefix sum, so snodeValPtr[kk] is where kk's block starts.
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
    nf.mNodeToSnode     = sf.nodeToSnode();
    nf.mFrontSize       = sf.frontSize();
    nf.mUpdateSize      = sf.updateSize();
    nf.mNumNodeIdx      = sf.numNodeIdx();
    nf.mSnodeNodeIdxPtr = sf.snodePtr();
    nf.mNodeIdx         = sf.nodeIdx();

    // The value blocks. Same shape as the static factor, indexSize rows by frontSize columns, but
    // one vector per supernode rather than offsets into one flat buffer, so a front can later grow
    // without moving its neighbours. Each zeroed, because assembly adds into it.
    nf.mVal.resize(nf.mSnodeSize);
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(nf.mSnodeSize); ++kk) {
        const std::size_t numNodeIdx = nf.mFrontSize[kk] + nf.mUpdateSize[kk];
        nf.mVal[kk].assign(numNodeIdx * nf.mFrontSize[kk], Val(0));
    }
}

template<class Val>
bool NumFactorEngine::assembleFromA(const SparseMatrix<Val>& A, const Permutation& p,
                                    const std::vector<std::int32_t>& gblToLcl,
                                    std::size_t frontSize, std::size_t numNodeIdx,
                                    const std::int32_t* nodeIdx, Val* block) const {
    const std::vector<std::size_t>&  colPtr   = A.colPtr();
    const std::vector<std::int32_t>& aRowIdx  = A.rowIdx();
    const std::vector<Val>&          aVal     = A.val();
    const std::vector<std::int32_t>& oldToNew = p.oldToNew();
    const std::vector<std::int32_t>& newToOld = p.newToOld();

    // For each front column of the supernode. Its local column position is lcl, and its block
    // column starts at lcl * numNodeIdx (column-major).
    for (std::size_t lcl = 0; lcl < frontSize; ++lcl) {
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

            block[lcl * numNodeIdx + static_cast<std::size_t>(lclRow)] = aVal[cp];
        }
    }
    return true;
}

template<class Val>
void NumFactorEngine::assembleUpdate(const std::vector<std::int32_t>& gblToLcl,
                                     const UpdateBlock<Val>& t,
                                     std::size_t numNodeIdx, Val* block) const {
    // The update block's rows and columns carry global row indices; gblToLcl maps them into the
    // ancestor's local coordinates. Only the lower triangle of the block is meaningful (row at or
    // below column), which is exactly the part the two BLAS calls filled.
    for (std::size_t tCol = 0; tCol < t.mWidth; ++tCol) {
        const std::int32_t li = t.mRowIdx[tCol];             // this column's global row index
        const std::size_t  lclCol = static_cast<std::size_t>(gblToLcl[li]);

        for (std::size_t tRow = tCol; tRow < t.mHeight; ++tRow) {
            const std::int32_t liRow  = t.mRowIdx[tRow];
            const std::size_t  lclRow = static_cast<std::size_t>(gblToLcl[liRow]);

            block[lclCol * numNodeIdx + lclRow] += t.mVal[tCol * t.mHeight + tRow];
        }
    }
}

template<class Val>
bool NumFactorEngine::factorSupernode(std::size_t frontSize, std::size_t numNodeIdx, Val* block,
                                      std::size_t& numPerturbations) const {
    const int f  = static_cast<int>(frontSize);
    const int u  = static_cast<int>(numNodeIdx - frontSize);
    const int ld = static_cast<int>(numNodeIdx);

    if (mFactorization == Factorization::Cholesky) {
        // The front, which is the diagonal block: A11 = L11 L11^H.
        int info = 0;
        potrf('L', f, block, ld, &info);
        if (info > 0)
            return false;   // not positive definite: the leading minor of order `info` failed

        // The update rows: L21 = A21 (L11^H)^-1, solved in place.
        //
        // Blas<Val>::conjTrans is 'T' for real and 'C' for complex, which is the whole of what the
        // scalar type decides here. 0.9 writes 'T' unconditionally, which is wrong for a complex
        // Hermitian factor, and there is nothing at its call site to reveal that.
        if (u > 0)
            trsm('R', 'L', Blas<Val>::conjTrans, 'N', u, f, Val(1), block, ld, block + f, ld);

        return true;
    }

    // LDL. The kernel is ours: LAPACK has no unpivoted LDL^T (its ?sytrf pivots, and pivoting is
    // what a *static* factorization refuses to do). It cannot fail, because there is no positive
    // definiteness to violate; a pivot too small to divide by is perturbed and counted.
    int numPert = 0;
    ldl(f, block, ld, mPerturbation, &numPert, hermitian(mFactorization));
    numPerturbations += static_cast<std::size_t>(numPert);

    // The update rows: L21 = A21 U11^-1, where U11 = D11 L11^H sits in the front's *upper*
    // triangle. So the solve is against the upper, untransposed, which is exactly what storing U
    // buys: Cholesky would have to transpose, and does.
    if (u > 0)
        trsm('R', 'U', 'N', 'N', u, f, Val(1), block, ld, block + f, ld);

    return true;
}

template<class Val>
void NumFactorEngine::updateSupernode(std::size_t frontSize, std::size_t numNodeIdx, const Val* block,
                                      std::size_t offset, UpdateBlock<Val>& t) const {
    const int f      = static_cast<int>(frontSize);
    const int ld     = static_cast<int>(numNodeIdx);
    const int height = static_cast<int>(t.mHeight);
    const int width  = static_cast<int>(t.mWidth);
    const int tld    = height;

    Val*       tVal = t.mVal.data();
    const Val* L21  = block + offset;   // the update rows that reach this ancestor, and below

    if (mFactorization == Factorization::Cholesky) {
        // The square part: t(0..width, 0..width) -= L21' L21'^H, where L21' is the `width` rows
        // that land in the ancestor. Symmetric, so HERK, which touches only the lower triangle.
        //
        // `herk` means "A times A-conjugate-transpose": dsyrk_ for real, zherk_ for complex. The
        // engine never names either, which is what makes 0.9's bug (SYRK on a Hermitian factor)
        // impossible to write here.
        herk('L', 'N', width, f, -1.0, L21, ld, 1.0, tVal, tld);

        // The rectangle below it: t(width.., 0..width) -= L21'' L21'^H, where L21'' is the rows of
        // the supernode's update block that lie *below* the ancestor's. Not symmetric, so GEMM.
        if (height > width)
            gemm('N', Blas<Val>::conjTrans, height - width, width, f,
                 Val(-1), L21 + width, ld,
                 L21, ld,
                 Val(1), tVal + width, tld);
        return;
    }

    // LDL. The update is t -= L21 D L21^H, and **no BLAS routine computes it**: the D in the
    // middle rules out a rank-k call, which is why Cholesky gets one and LDL does not.
    //
    // So form U := D L21'^H explicitly, into a scratch, and then two multiplies. The scratch is
    // f by width, and it is the whole price of the D.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(width), Val(0));
    formUpper(width, f, L21, ld, upper.data(), f, block, ld, hermitian(mFactorization));

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
// Left-looking. For each supernode: pull every update owed to it, then factor.
//
// The bookkeeping is the interesting part. When jj is factored it owes an update to every
// supernode that owns one of jj's update rows, and those ancestors are visited in increasing
// order. So jj is placed in a list belonging to the *next* ancestor it must update, and moves to
// the next list each time it discharges one. pos[jj] records how far through jj's index set that
// has got.
//
// The lists therefore hold, at the moment kk is reached, exactly the supernodes that still owe kk
// an update. No search, no scan over descendants.
// =================================================================================================

template<class Val, class Factor>
bool NumFactorEngine::factorLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                        const SymFactor& sf, Factor& nf) const {
    setSymFactor(sf, nf);

    const std::size_t size    = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // Assemble A's original values into every supernode first. Cheaper than doing it inside the
    // main loop: the local map is set and cleared once per supernode either way, but this keeps
    // the traversal's own bookkeeping uncluttered.
    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   numNodeIdx = nf.frontSize(kk) + nf.updateSize(kk);
        const std::int32_t* nodeIdx    = nf.nodeIdxPtr(kk);
        Val*                block  = nf.valPtr(kk);

        setGlobalToLocal(numNodeIdx, nodeIdx, gblToLcl);
        const bool ok = assembleFromA(A, p, gblToLcl, nf.frontSize(kk), numNodeIdx, nodeIdx, block);
        clearGlobalToLocal(numNodeIdx, nodeIdx, gblToLcl);
        if (!ok)
            return false;
    }

    // Who still owes whom, and how far each has got.
    std::vector<std::list<std::int32_t>> owed(snodeSize);
    std::vector<std::size_t>             pos(snodeSize, 0);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkNumNodeIdx = nf.frontSize(kk) + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdxPtr(kk);
        Val*                kkBlock  = nf.valPtr(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // Every supernode jj that owes kk an update.
        while (!owed[kk].empty()) {
            const std::int32_t jj = owed[kk].front();
            owed[kk].pop_front();

            const std::size_t   jjFrontSize = nf.frontSize(jj);
            const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
            const std::int32_t* jjNodeIdx    = nf.nodeIdxPtr(jj);
            const Val*          jjBlock     = std::as_const(nf).valPtr(jj);

            // How many of jj's remaining rows belong to kk. They are contiguous, because jj's
            // index set is sorted and the supernodes partition it in increasing order.
            const std::size_t from   = pos[jj];
            const std::size_t height = jjNumNodeIdx - from;
            std::size_t       width  = 0;
            while (from + width < jjNumNodeIdx
                   && nf.nodeToSnode()[jjNodeIdx[from + width]] == kk)
                ++width;

            UpdateBlock<Val> t(height, width);
            std::copy(jjNodeIdx + from, jjNodeIdx + jjNumNodeIdx, t.mRowIdx.begin());

            updateSupernode(jjFrontSize, jjNumNodeIdx, jjBlock, from, t);
            assembleUpdate(gblToLcl, t, kkNumNodeIdx, kkBlock);

            // jj has discharged kk. Queue it against the next ancestor it owes.
            pos[jj] = from + width;
            if (pos[jj] < jjNumNodeIdx)
                owed[nf.nodeToSnode()[jjNodeIdx[pos[jj]]]].push_back(jj);
        }

        if (!factorSupernode(nf.frontSize(kk), kkNumNodeIdx, kkBlock, nf.numPerturbations())) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;   // not positive definite (Cholesky only; LDL perturbs instead)
        }

        // kk is factored, so it now owes updates of its own. Its front rows are its own columns
        // and update nobody; the first update row names the first ancestor it owes.
        pos[kk] = nf.frontSize(kk);
        if (pos[kk] < kkNumNodeIdx)
            owed[nf.nodeToSnode()[kkNodeIdx[pos[kk]]]].push_back(kk);

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
bool NumFactorEngine::factorRightLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                         const SymFactor& sf, Factor& nf) const {
    setSymFactor(sf, nf);

    const std::size_t size    = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   numNodeIdx = nf.frontSize(kk) + nf.updateSize(kk);
        const std::int32_t* nodeIdx    = nf.nodeIdxPtr(kk);
        Val*                block  = nf.valPtr(kk);

        setGlobalToLocal(numNodeIdx, nodeIdx, gblToLcl);
        const bool ok = assembleFromA(A, p, gblToLcl, nf.frontSize(kk), numNodeIdx, nodeIdx, block);
        clearGlobalToLocal(numNodeIdx, nodeIdx, gblToLcl);
        if (!ok)
            return false;
    }

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        const std::size_t   jjFrontSize = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdxPtr(jj);
        Val*                jjBlock     = nf.valPtr(jj);

        if (!factorSupernode(jjFrontSize, jjNumNodeIdx, jjBlock, nf.numPerturbations()))
            return false;   // not positive definite (Cholesky only; LDL perturbs instead)

        // Walk jj's update rows. Each run of them belonging to one ancestor is one update.
        std::size_t from = jjFrontSize;
        while (from < jjNumNodeIdx) {
            const std::int32_t kk = nf.nodeToSnode()[jjNodeIdx[from]];

            const std::size_t   kkNumNodeIdx = nf.frontSize(kk) + nf.updateSize(kk);
            const std::int32_t* kkNodeIdx    = nf.nodeIdxPtr(kk);
            Val*                kkBlock  = nf.valPtr(kk);

            const std::size_t height = jjNumNodeIdx - from;
            std::size_t       width  = 0;
            while (from + width < jjNumNodeIdx
                   && nf.nodeToSnode()[jjNodeIdx[from + width]] == kk)
                ++width;

            UpdateBlock<Val> t(height, width);
            std::copy(jjNodeIdx + from, jjNodeIdx + jjNumNodeIdx, t.mRowIdx.begin());

            updateSupernode(jjFrontSize, jjNumNodeIdx, jjBlock, from, t);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdate(gblToLcl, t, kkNumNodeIdx, kkBlock);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            from += width;
        }
    }

    return true;
}

template<class Val>
bool NumFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                              NumFactorStatic<Val>& nf) const {
    if (A.size() != p.size() || A.size() != sf.size())
        return false;

    // Cholesky and static LDL, in both transposes. Dynamic LDL and multifrontal follow.
    switch (mFactorization) {
        case Factorization::Cholesky:
        case Factorization::StaticLDLT:
        case Factorization::StaticLDLH:
            break;
        case Factorization::DynamicLDLT:
        case Factorization::DynamicLDLH:
            return false;   // not implemented
    }

    switch (mTraversal) {
        case Traversal::LeftLooking:  return factorLeftLooking(A, p, sf, nf);
        case Traversal::RightLooking: return factorRightLooking(A, p, sf, nf);
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

    // The static factorizations run unchanged into per-supernode storage: same traversals, same
    // kernels, only the block addresses differ. Dynamic LDL, which needs this storage to grow a
    // front under a delayed pivot, is next; until then its cases are refused here as *not yet*,
    // unlike the static overload, where they are refused permanently (a flat buffer cannot grow).
    switch (mFactorization) {
        case Factorization::Cholesky:
        case Factorization::StaticLDLT:
        case Factorization::StaticLDLH:
            break;
        case Factorization::DynamicLDLT:
        case Factorization::DynamicLDLH:
            return false;   // not implemented yet
    }

    switch (mTraversal) {
        case Traversal::LeftLooking:  return factorLeftLooking(A, p, sf, nf);
        case Traversal::RightLooking: return factorRightLooking(A, p, sf, nf);
        case Traversal::Multifrontal: return false;   // not implemented
    }
    return false;
}

template bool NumFactorEngine::compute(const SparseMatrix<double>&, const Permutation&,
                                       const SymFactor&, NumFactorDynamic<double>&) const;
template bool NumFactorEngine::compute(const SparseMatrix<std::complex<double>>&,
                                       const Permutation&, const SymFactor&,
                                       NumFactorDynamic<std::complex<double>>&) const;

} // namespace Oblio
