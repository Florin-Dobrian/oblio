#include "oblio/NumFactorEngine.h"

#include "oblio/BlasLapack.h"

#include <algorithm>
#include <cmath>
#include <list>

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
void NumFactorEngine::initNumFactor(const SymFactor& sf, NumFactorStatic<Val>& nf) const {
    nf.mSize          = sf.size();
    nf.mSnodeSize     = sf.snodeSize();
    nf.mFactorization = mFactorization;

    // This run's perturbation count starts at zero; the static LDL factorization accumulates into
    // it (Cholesky never perturbs, dynamic LDL delays instead). Reset here, not at construction, so
    // a reused factor does not carry a previous run's count.
    nf.mNumPerturbations = 0;

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
void NumFactorEngine::initNumFactor(const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    nf.mSize          = sf.size();
    nf.mSnodeSize     = sf.snodeSize();
    nf.mFactorization = mFactorization;

    // This run's perturbation count starts at zero. Dynamic LDL never perturbs (it delays), but a
    // static factorization run into this storage does, and a reused factor must not carry a previous
    // run's count.
    nf.mNumPerturbations = 0;

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

template<class Val, class Factor>
void NumFactorEngine::assembleUpdateMatrix(const UpdateMatrix<Val>& childUpdate, std::int32_t kk,
                                           Factor& nf, UpdateMatrix<Val>& kkUpdate,
                                           const std::vector<std::int32_t>& gblToLcl) const {
    const std::int32_t  childSize = static_cast<std::int32_t>(childUpdate.size());
    const std::int32_t* childIdx  = childUpdate.nodeIdx();
    const Val*          childVal  = childUpdate.val();   // childSize by childSize, ld == childSize

    const std::int32_t kkFrontSize  = static_cast<std::int32_t>(nf.frontSize(kk));
    const std::int32_t kkUpdateSize = static_cast<std::int32_t>(nf.updateSize(kk));
    const std::size_t  kkIndexSize  = static_cast<std::size_t>(kkFrontSize) + kkUpdateSize;
    Val*               kkVal        = nf.val(kk);         // lu block, ld == kkIndexSize
    Val*               kkUpdateVal  = kkUpdate.val();     // kk's contribution block, ld == kkUpdateSize

    // Each column of the child's block carries a global index that lies somewhere in kk's index set.
    // Where it lands decides which block receives it: a pivot column of kk (local position below
    // kkFrontSize) goes into the lu block, an update row of kk into kk's contribution block. The
    // child's indices are sorted and gblToLcl preserves order, so within a column the rows run at or
    // below the diagonal, and only the lower triangle is written.
    for (std::int32_t c = 0; c < childSize; ++c) {
        const std::int32_t dj = gblToLcl[childIdx[c]];
        const std::size_t  sc = static_cast<std::size_t>(c) * childSize;   // child column offset

        if (dj < kkFrontSize) {
            const std::size_t dc = static_cast<std::size_t>(dj) * kkIndexSize;
            for (std::int32_t r = c; r < childSize; ++r) {
                const std::int32_t di = gblToLcl[childIdx[r]];
                kkVal[dc + di] += childVal[sc + r];
            }
        } else {
            const std::size_t dc = static_cast<std::size_t>(dj - kkFrontSize) * kkUpdateSize;
            for (std::int32_t r = c; r < childSize; ++r) {
                const std::int32_t di = gblToLcl[childIdx[r]];
                kkUpdateVal[dc + di - kkFrontSize] += childVal[sc + r];
            }
        }
    }
}

template<class Val, class Factor>
bool NumFactorEngine::factorStaticSupernode(Factor& nf, std::int32_t jj) const {
    const std::size_t frontSize   = nf.frontSize(jj);
    const std::size_t numNodeIdx  = frontSize + nf.updateSize(jj);
    Val*              val         = nf.val(jj);

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
    nf.numPerturbations() += static_cast<std::size_t>(numPert);

    // The update rows: L21 = A21 U11^-1, where U11 = D11 L11^H sits in the front's *upper*
    // triangle. So the solve is against the upper, untransposed, which is exactly what storing U
    // buys: Cholesky would have to transpose, and does.
    if (u > 0)
        trsm('R', 'U', 'N', 'N', u, f, Val(1), val, ld, val + f, ld);

    return true;
}

template<class Val, class Factor>
void NumFactorEngine::updateStaticSupernode(const Factor& nf, std::int32_t jj,
                                      std::size_t jjKkUpdateSp, UpdateBlock<Val>& updateBlock) const {
    // `jjKkUpdateSp` is a position into jj's index set, and it is what selects the ancestor kk this
    // update reaches: kk == nf.nodeToSnode(nf.nodeIdx(jj)[jjKkUpdateSp]). The kernel never names kk
    // (it forms jj's outer product from that row down, and the caller has already sized the block to
    // the rows kk owns), but that is the position's meaning and it is worth seeing first.
    const std::size_t jjFrontSize  = nf.frontSize(jj);
    const std::size_t jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
    const Val*        jjVal        = nf.val(jj);

    const int f          = static_cast<int>(jjFrontSize);
    const int sld        = static_cast<int>(jjNumNodeIdx);
    const int jjKkHeight = static_cast<int>(updateBlock.mHeight);
    const int jjKkWidth  = static_cast<int>(updateBlock.mWidth);
    const int dld        = jjKkHeight;

    const bool withHermitian = hermitian(mFactorization);
    const Val* l21Val        = jjVal + jjKkUpdateSp;   // the update rows that reach this ancestor, and below
    Val*       u22Val        = updateBlock.mVal.data();

    // `jjKkHeight` and `jjKkWidth` are dimensions of the (jj, kk) *edge*, not of either supernode
    // alone. `jjKkWidth` is how many of jj's rows land in kk, a two-supernode count: change jj or kk
    // and it moves. `jjKkHeight` is all of jj's rows from `jjKkUpdateSp` down, kk's run plus any that
    // reach higher ancestors, so it too is anchored at kk's start. Only `f`, jj's own front, is a
    // one-supernode dimension here. This kernel is the edge's work: it never forms all of jj's L, only
    // the `f`-by-`jjKkWidth` slice D L21_kk^H restricted to jj's kk-landing rows, which is why it runs
    // once per (jj, kk) the traversal visits rather than once per jj.
    //
    // The outer product U22 = L21 (D) L21^H splits into two pieces by symmetry, not by anything about
    // the ancestor. The top jjKkWidth-by-jjKkWidth square is a block times its own conjugate transpose,
    // symmetric, so only its lower triangle is computed. The (jjKkHeight - jjKkWidth)-by-jjKkWidth
    // rectangle below is jj's higher-reaching rows against kk's rows, not symmetric, a plain multiply.
    // The block is stored as a full rectangle for convenience but no arithmetic is spent on the
    // square's upper triangle.
    //
    // `jjKkHeight == jjKkWidth` means jj has no rows beyond kk, so there is no rectangle and only the
    // square runs. It does *not* mean kk is a root: it is a fact about how far jj reaches, not about
    // kk's place in the tree. A root kk forces it for every jj that reaches it, but it also happens at
    // interior ancestors jj simply does not extend past.

    if (mFactorization == Factorization::Cholesky) {
        // The square part: the block's (0..jjKkWidth, 0..jjKkWidth) -= L21' L21'^H, where L21' is the `jjKkWidth` rows
        // that land in the ancestor. Symmetric, so HERK, which touches only the lower triangle.
        //
        // `herk` means "A times A-conjugate-transpose": dsyrk_ for real, zherk_ for complex. The
        // engine never names either, which is what makes 0.9's bug (SYRK on a Hermitian factor)
        // impossible to write here.
        herk('L', 'N', jjKkWidth, f, -1.0, l21Val, sld, 1.0, u22Val, dld);

        // The rectangle below, present only when jj reaches past this ancestor: the block's
        // (jjKkWidth.., 0..jjKkWidth) -= L21'' L21'^H, where L21'' is the rows of jj's update val below the
        // ancestor's. Not symmetric, so GEMM.
        if (jjKkHeight > jjKkWidth)
            gemm('N', Blas<Val>::conjTrans, jjKkHeight - jjKkWidth, jjKkWidth, f,
                 Val(-1), l21Val + jjKkWidth, sld,
                 l21Val, sld,
                 Val(1), u22Val + jjKkWidth, dld);
        return;
    }

    // LDL. The update is `block -= L21 D L21^H`, and **no BLAS routine computes it**: the D in the
    // middle rules out a rank-k call, which is why Cholesky gets one and LDL does not. Cholesky's
    // block leaves its upper triangle as unused zeros, and `herk` is content with that: it wants no
    // stored intermediate, only C21. LDL instead materializes U := D L21^H, and this is a choice, not
    // a necessity. We could have wasted the upper triangle the same way and written a fused kernel, an
    // `oblioHerk` that forms and consumes D L21^H internally and stores nothing. We do not, because U
    // is worth keeping: in the factor kernel the same U = D L^H is needed twice (to solve for the next
    // L, and to update), so it is formed once and stored where Cholesky wasted space. Here in the
    // update that value goes to a local scratch rather than the block, but the split is the same one:
    // form U, then multiply. That is why LDL is formStaticUpper + gemmLower where Cholesky is a
    // single herk.
    //
    // The scratch is f by jjKkWidth, and it is the whole price of the D.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(jjKkWidth), Val(0));
    formStaticUpper(jjKkWidth, f, l21Val, sld, upper.data(), f, jjVal, sld, withHermitian);

    // The square part, the counterpart of Cholesky's herk: symmetric, so only its lower triangle is
    // filled. BLAS has nothing for this either (syrk does A A^T, not A B with B known to make the
    // product symmetric), so gemmLower is ours as well. It multiplies L21 against the U just formed.
    gemmLower(jjKkWidth, f, l21Val, sld, upper.data(), f, u22Val, dld);

    // The rectangle below, present only when jj reaches past this ancestor: not symmetric, so a
    // plain GEMM. Note 'N','N': the transpose is already baked into U.
    if (jjKkHeight > jjKkWidth)
        gemm('N', 'N', jjKkHeight - jjKkWidth, jjKkWidth, f,
             Val(-1), l21Val + jjKkWidth, sld,
             upper.data(), f,
             Val(1), u22Val + jjKkWidth, dld);
}

template<class Val, class Factor>
void NumFactorEngine::updateStaticMultifrontal(const Factor& nf, std::int32_t kk,
                                               UpdateMatrix<Val>& kkUpdate) const {
    const int u = static_cast<int>(nf.updateSize(kk));
    if (u == 0)
        return;   // a root, or a supernode reaching nowhere: no contribution to leave

    const int  f   = static_cast<int>(nf.frontSize(kk));
    const int  ld  = f + u;                 // kk's lu block height
    const Val* val = nf.val(kk);            // kk's factored block
    const Val* l21 = val + f;               // its update rows
    Val*       uVal = kkUpdate.val();       // the contribution block, ld == u

    if (mFactorization == Factorization::Cholesky) {
        // U -= L21 L21^H, one rank-k call. herk resolves to syrk for real and herk for complex by
        // the scalar type, so the complex case is Hermitian without the engine choosing, which is
        // where 0.9's multifrontal was silently wrong (it used syrk for both).
        herk('L', 'N', u, f, -1.0, l21, ld, 1.0, uVal, u);
        return;
    }

    // LDL: U -= L21 D L21^H. The D in the middle rules out a rank-k call, so form U := D L21^H into a
    // scratch and then multiply, exactly the square part of the left/right-looking LDL update with no
    // rectangle below (the whole contribution block is the symmetric square).
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(u), Val(0));
    formStaticUpper(u, f, l21, ld, upper.data(), f, val, ld, hermitian(mFactorization));
    gemmLower(u, f, l21, ld, upper.data(), f, uVal, u);
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
    initNumFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // For each kk, the descendants queued to update it; and how far each has got.
    std::vector<std::list<std::int32_t>> descendantUpdateQueue(snodeSize);
    std::vector<std::size_t>             nextUpdateSp(snodeSize, 0);   // sp of jj's next update row

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A's own values into kk's front, on the map just set. The static front is at its final
        // width from the start, so unlike the dynamic driver this could be a separate prepass; it
        // is folded in here to keep all four drivers one shape.
        if (!assembleFromA(A, p, gblToLcl, 0, kkFrontSize, kkNumNodeIdx, kkNodeIdx, kkVal)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Every descendant jj queued to update kk.
        while (!descendantUpdateQueue[kk].empty()) {
            const std::int32_t jj = descendantUpdateQueue[kk].front();
            descendantUpdateQueue[kk].pop_front();

            const std::size_t   jjFrontSize  = nf.frontSize(jj);
            const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
            const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);

            // How many of jj's remaining rows belong to kk. They are contiguous, because jj's
            // index set is sorted and the supernodes partition it in increasing order.
            const std::size_t jjKkUpdateSp = nextUpdateSp[jj];
            const std::size_t jjKkHeight   = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth    = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> updateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateStaticSupernode(nf, jj, jjKkUpdateSp, updateBlock);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);

            // jj has updated kk. Queue it for the next ancestor it must update.
            nextUpdateSp[jj] = jjKkUpdateSp + jjKkWidth;
            if (nextUpdateSp[jj] < jjNumNodeIdx)
                descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        }

        if (!factorStaticSupernode<Val>(nf, kk)) {
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
    initNumFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);
        Val*                jjVal        = nf.val(jj);

        setGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        if (!assembleFromA(A, p, gblToLcl, 0, jjFrontSize,
                           jjNumNodeIdx, jjNodeIdx, jjVal)) {
            clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
            return false;
        }

        clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
    }

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        const std::size_t   jjFrontSize  = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx    = nf.nodeIdx(jj);

        if (!factorStaticSupernode<Val>(nf, jj))
            return false;   // not positive definite (Cholesky only; LDL perturbs instead)

        // Walk jj's update rows. Each run of them belonging to one ancestor is one update.
        std::size_t jjKkUpdateSp = jjFrontSize;
        while (jjKkUpdateSp < jjNumNodeIdx) {
            const std::int32_t kk = nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp]);

            const std::size_t   kkFrontSize  = nf.frontSize(kk);
            const std::size_t   kkNumNodeIdx = kkFrontSize + nf.updateSize(kk);
            const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
            Val*                kkVal        = nf.val(kk);

            const std::size_t jjKkHeight = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth  = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> updateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateStaticSupernode(nf, jj, jjKkUpdateSp, updateBlock);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjKkUpdateSp += jjKkWidth;
        }
    }

    return true;
}

// =================================================================================================
// Static multifrontal. One postorder-compatible pass (supernodes in increasing order): assemble A
// and every child's contribution block into the frontal, factor, then leave this supernode's own
// contribution block on the stack for its parent. Cholesky only for now; see the header.
// =================================================================================================

template<class Val, class Factor>
bool NumFactorEngine::factorStaticMultifrontal(const SparseMatrix<Val>& A, const Permutation& p,
                                               const SymFactor& sf, Factor& nf) const {
    initNumFactor(sf, nf);

    const std::size_t snodeSize = nf.snodeSize();

    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(nf.size(), NIL);

    // One contribution block per supernode, allocated when the supernode is reached and freed once
    // its parent has assembled it. Sized once; the slots start empty.
    std::vector<UpdateMatrix<Val>> stack(snodeSize);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkUpdateSize = nf.updateSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + kkUpdateSize;
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        // Map kk's whole index set, front and update rows alike: assembleFromA fills the front
        // columns, and the extend-add routes each child entry by where its index falls in kk.
        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A's own values into kk's lu block.
        if (!assembleFromA(A, p, gblToLcl, 0, kkFrontSize, kkNumNodeIdx, kkNodeIdx, kkVal)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Allocate kk's contribution block over its update rows; allocate zeroes the values, so only
        // the indices are filled here, from kk's update-row indices.
        stack[kk].allocate(kkUpdateSize);
        {
            std::int32_t* kkUpdateIdx = stack[kk].nodeIdx();
            for (std::size_t sp = 0; sp < kkUpdateSize; ++sp)
                kkUpdateIdx[sp] = kkNodeIdx[kkFrontSize + sp];
        }

        // Extend-add each child's contribution block into kk's frontal, then free it. A child of kk
        // is factored in an earlier iteration (increasing order), so its block is present here.
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
            assembleUpdateMatrix(stack[jj], kk, nf, stack[kk], gblToLcl);
            stack[jj].discard();
        }

        // Factor kk's pivots, the same per-supernode kernel the other traversals use: Cholesky is
        // POTRF then TRSM, static LDL is the unpivoted ldl then TRSM against the upper triangle.
        if (!factorStaticSupernode<Val>(nf, kk)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;   // Cholesky non-positive-definite; LDL perturbs instead and cannot fail
        }

        // Form kk's contribution block from the factored pivots and leave it on the stack.
        updateStaticMultifrontal(nf, kk, stack[kk]);

        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
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
    std::int32_t* idx = nf.nodeIdx(jj);

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
                                             std::size_t jjKkUpdateSp, UpdateBlock<Val>& updateBlock) const {
    // `jjKkUpdateSp` is a position into jj's index set, and it is what selects the ancestor kk this
    // update reaches: kk == nf.nodeToSnode(nf.nodeIdx(jj)[jjKkUpdateSp]). The kernel never names kk
    // (it forms jj's outer product from that row down, and the caller has already sized the block to
    // the rows kk owns), but that is the position's meaning and it is worth seeing first.
    //
    // This is the same update the static twin performs: jj's update area updates kk, driven by the
    // update (descendant-to-ancestor) relationship. It is agnostic to whether jj has been contracted.
    // Left-looking calls it with jj already contracted, right-looking with jj not yet contracted, and
    // the result is identical, because the read touches only jj's update area. That area is disjoint
    // from both delay regions: jj's own delayed columns sit below its front but the walk starts past
    // them at `jjKkUpdateSp`, and kk's inbound delays arrive by a separate path (assembleDelay), never as
    // this update. The stride below counts the delayed rows, so the column layout is the same whether
    // or not the delayed-column storage past the front has been trimmed. updateSize never changes, so
    // the update area is stable ground. The delay flow (child-to-parent) runs independently of the
    // update flow (descendant-to-ancestor); they share these drivers but are not aligned.
    //
    // The update's rank is `f`, jj's *post-factor* front, the pivots actually eliminated here. A
    // delayed column is not among them, so it never updates an ancestor from jj: its row was updated
    // in place by those pivots back inside factorDynamicSupernode, and it then migrates to the parent
    // as a delayed column (assembleDelay, contractVal) to be pivoted there. So this kernel touches
    // only the update area; the delay area was finished in the factor step and leaves by the delay
    // path.
    const std::size_t jjFrontSize  = nf.frontSize(jj);
    // The one term that differs from the static twin, where it is just frontSize + updateSize: a
    // delayed column keeps its row, so the index count (and thus the stride) still counts it.
    const std::size_t jjNumNodeIdx = jjFrontSize + nf.delaySize(jj) + nf.updateSize(jj);
    const Val*        jjVal        = nf.val(jj);

    const int f = static_cast<int>(jjFrontSize);
    if (f == 0)
        return;   // every column of jj was delayed: there is no pivot here to update anyone with

    const int sld        = static_cast<int>(jjNumNodeIdx);
    const int jjKkHeight = static_cast<int>(updateBlock.mHeight);
    const int jjKkWidth  = static_cast<int>(updateBlock.mWidth);
    const int dld        = jjKkHeight;

    // As in the static twin, jjKkHeight and jjKkWidth are the (jj, kk) edge's dimensions, not either
    // supernode's: jjKkWidth is jj's rows landing in kk (change jj or kk and it moves), jjKkHeight all
    // of jj's rows from jjKkUpdateSp down. Only f is jj's alone. This kernel is the edge's work, and
    // formDynamicUpper forms only the f-by-jjKkWidth slice for this one (jj, kk) pair.

    const bool          withHermitian = hermitian(nf.factorization());
    const std::int32_t* jjNodeIdx     = nf.nodeIdx(jj);
    // jj's rows from `jjKkUpdateSp` down: the block this ancestor is about to receive.
    const Val*          l21Val        = jjVal + jjKkUpdateSp;
    Val*                u22Val        = updateBlock.mVal.data();

    // LDL, exactly as in the static twin: the update is `block -= L21 D L21^H`, no BLAS routine
    // computes it (the D in the middle rules out a rank-k call), so we form U := D L21^H into a
    // scratch and then multiply, gemmLower for the symmetric square and gemm for the rectangle below.
    // Dynamic runs only for LDL, so there is no Cholesky branch here; the whole function is the LDL
    // half of the static twin, and past forming U it is line-for-line the static code.
    //
    // The one departure is which form-upper. Static calls formStaticUpper, whose D is a plain
    // diagonal; dynamic calls formDynamicUpper, whose D is block-diagonal, 1x1 and 2x2 pivots marked
    // by mPivotType (indexed by the global node jjNodeIdx). Both fill the same f-by-jjKkWidth scratch.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(jjKkWidth), Val(0));
    formDynamicUpper(jjKkWidth, f, l21Val, sld, upper.data(), f, jjVal, sld,
                     nf.mPivotType.data(), jjNodeIdx, withHermitian);

    // From here it is the static twin exactly. The square part, the counterpart of Cholesky's herk:
    // symmetric, so gemmLower fills only its lower triangle, multiplying L21 against the U just formed.
    gemmLower(jjKkWidth, f, l21Val, sld, upper.data(), f, u22Val, dld);

    // The rectangle below, present only when jj reaches past this ancestor: not symmetric, so a plain
    // GEMM. Note 'N','N': the transpose is already baked into U.
    if (jjKkHeight > jjKkWidth)
        gemm('N', 'N', jjKkHeight - jjKkWidth, jjKkWidth, f,
             Val(-1), l21Val + jjKkWidth, sld,
             upper.data(), f,
             Val(1), u22Val + jjKkWidth, dld);
}

template<class Val>
void NumFactorEngine::updateDynamicMultifrontal(const NumFactorDynamic<Val>& nf, std::int32_t kk,
                                                UpdateMatrix<Val>& kkUpdate) const {
    const int f = static_cast<int>(nf.frontSize(kk));    // post-factor front
    const int u = static_cast<int>(nf.updateSize(kk));
    if (f == 0 || u == 0)
        return;

    // The block height counts the front, the delayed columns, and the update rows: the update rows
    // sit past the front *and* the delayed columns, which is the only place this differs from the
    // static form.
    const int           d       = static_cast<int>(nf.delaySize(kk));
    const int           ld      = f + d + u;
    const bool          herm    = hermitian(mFactorization);
    const Val*          val     = nf.val(kk);
    const std::int32_t* nodeIdx = nf.nodeIdx(kk);
    const Val*          l21     = val + f + d;

    // U -= L21 D L21^H with a block-diagonal D. formDynamicUpper walks pivotType over the front,
    // handling 1x1 and 2x2 pivots, to build D L21^H into a scratch; gemmLower then multiplies.
    std::vector<Val> upper(static_cast<std::size_t>(f) * static_cast<std::size_t>(u), Val(0));
    formDynamicUpper(u, f, l21, ld, upper.data(), f, val, ld, nf.mPivotType.data(), nodeIdx, herm);
    gemmLower(u, f, l21, ld, upper.data(), f, kkUpdate.val(), u);
}

template<class Val>
void NumFactorEngine::assembleDelay(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t kk,
                                    const std::vector<std::int32_t>& gblToLcl) const {
    const std::int32_t jjFrontSize = static_cast<std::int32_t>(nf.mFrontSize[jj]);
    const std::int32_t jjDelayed   = static_cast<std::int32_t>(nf.mDelaySize[jj]);
    const std::int32_t jjRows      = jjFrontSize + jjDelayed
                                   + static_cast<std::int32_t>(nf.mUpdateSize[jj]);
    const std::int32_t kkRows      = static_cast<std::int32_t>(nf.mFrontSize[kk] + nf.mUpdateSize[kk]);

    const std::int32_t* jjNodeIdx = nf.nodeIdx(jj);
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
//
// **The extreme case: every column of a supernode delays.** frontSize then reaches zero. The
// supernode still exists, its val still holds every row (they were reclassified front -> delayed,
// not removed), and the height identity still holds. Once its parent takes those delayed columns,
// contractVal resizes the val to frontSize * height = 0: the val goes genuinely empty, since a val
// is dimensioned by column count and there are no columns left. The nodeIdx is not shrunk to match,
// and it does not need to be. Nothing reads it afterward: assembleDelay already extracted the
// delayed globals while the columns were still present, and every later reader is gated on
// frontSize (the solve loops j < frontSize, updateDynamicSupernode returns on f == 0), so a
// zero-front supernode is uniformly skipped. The surviving rows are vestigial, kept because
// "nodeIdx never contracts" is a simpler invariant to hold than one with an emptying special case,
// not because anyone consults them. Clearing them would be equally correct and save a little
// memory; the code does not bother.
// =================================================================================================

template<class Val>
bool NumFactorEngine::factorDynamicLeftLooking(const SparseMatrix<Val>& A, const Permutation& p,
                                               const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    initNumFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    const std::vector<std::int32_t>& parent      = sf.parent();
    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // No A-assembly prepass here, and this is the essential left-vs-right-looking difference, not
    // an incidental one. Left-looking *pulls*: when kk is reached it gathers updates from its
    // descendants, all already factored, and nothing writes into kk's front until kk's own turn. So
    // A can be assembled into kk lazily, at the top of its iteration, right before the pulls begin.
    // Right-looking instead *pushes* (see factorDynamicRightLooking), which forces a prepass; the two
    // drivers mirror each other around the direction of update flow, and that mirror is the one
    // asymmetry between them that must stay. A second reason reinforces the lazy assembly: a front's
    // final shape is not known until its descendants are factored, since their delayed columns widen
    // it, so kk could not be initialized up front even if the flow allowed it.
    std::vector<std::list<std::int32_t>> descendantUpdateQueue(snodeSize);
    std::vector<std::size_t>             nextUpdateSp(snodeSize, 0);   // sp of jj's next update row

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        std::size_t kkInboundDelaySize = 0;
        for (std::int32_t ii = firstChild[kk]; ii != NIL; ii = nextSibling[ii])
            kkInboundDelaySize += nf.delaySize(ii);

        if (kkInboundDelaySize > 0) {
            // The height after expanding, which is also the height before it: the new columns come
            // from rows the val already had. Only the front/update split moves.
            const std::size_t kkPreExpandFrontSize = nf.frontSize(kk);
            const std::size_t kkPostExpandNumNodeIdx =
                kkInboundDelaySize + kkPreExpandFrontSize + nf.updateSize(kk);

            // The index set. Extend first, then shift kk's own indices right by the delayed count,
            // descending so the copy does not overwrite its own source.
            nf.expandNodeIdx(kk, kkInboundDelaySize);
            std::vector<std::int32_t>& kkPostExpandNodeIdx = nf.mNodeIdx[kk];

            for (std::size_t ssp = kkPostExpandNumNodeIdx - kkInboundDelaySize; ssp-- > 0; ) {
                const std::size_t dsp = ssp + kkInboundDelaySize;
                kkPostExpandNodeIdx[dsp] = kkPostExpandNodeIdx[ssp];
            }

            // Then the vacated slots at the left, filled from the children in sibling order. Each
            // child's delayed columns are the run just past its (already reduced) front.
            std::size_t dsp = 0;
            for (std::int32_t ii = firstChild[kk]; ii != NIL; ii = nextSibling[ii]) {
                const std::size_t   iiPostFactorFrontSize = nf.frontSize(ii);
                const std::int32_t* iiNodeIdx             = nf.nodeIdx(ii);

                for (std::size_t ssp = iiPostFactorFrontSize;
                     ssp < iiPostFactorFrontSize + nf.delaySize(ii); ++ssp, ++dsp)
                    kkPostExpandNodeIdx[dsp] = iiNodeIdx[ssp];
            }

            // And the val. The front is wider, so the old contents are discarded rather than
            // moved: nothing has been written into kk yet, which is why left-looking never needs
            // 0.9's extendEntry_.
            nf.mFrontSize[kk] += kkInboundDelaySize;
            nf.resetVal(kk);
        }

        // The full height, captured before the factorization reclassifies part of the front as
        // delayed. Every use below wants this number and not the contracted front.
        const std::size_t   kkPreFactorFrontSize = nf.frontSize(kk);
        const std::size_t   kkNumNodeIdx         = kkPreFactorFrontSize + nf.updateSize(kk);
        const std::int32_t* kkNodeIdx            = nf.nodeIdx(kk);
        Val*                kkVal                = nf.val(kk);   // stays valid across the loop below: kk was
                                                         // resized once in the migration above, and
                                                         // nothing inside the loop resizes it again

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        if (!assembleFromA(A, p, gblToLcl, kkInboundDelaySize,
                           kkPreFactorFrontSize, kkNumNodeIdx, kkNodeIdx, kkVal)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Take every child's delayed columns, then let the child reclaim their storage. Tested per
        // child, not just on the total: one child delaying does not mean its siblings did, and the
        // calls are no-ops for those that did not. A child of kk is factored in an earlier iteration
        // (ascending order), so its delays are known here, exactly as in the right-looking driver.
        if (kkInboundDelaySize > 0)
            for (std::int32_t ii = firstChild[kk]; ii != NIL; ii = nextSibling[ii]) {
                const std::size_t iiDelaySize = nf.delaySize(ii);
                if (iiDelaySize > 0) {
                    assembleDelay(nf, ii, kk, gblToLcl);
                    nf.contractVal(ii, iiDelaySize);
                }
            }

        // Every descendant jj queued to update kk.
        while (!descendantUpdateQueue[kk].empty()) {
            const std::int32_t jj = descendantUpdateQueue[kk].front();
            descendantUpdateQueue[kk].pop_front();

            const std::size_t   jjPostFactorFrontSize = nf.frontSize(jj);
            const std::size_t   jjDelaySize           = nf.delaySize(jj);
            const std::size_t   jjNumNodeIdx          = jjPostFactorFrontSize + jjDelaySize
                                                      + nf.updateSize(jj);
            const std::int32_t* jjNodeIdx             = nf.nodeIdx(jj);

            // How many of jj's remaining rows belong to kk. Contiguous, as in the static case.
            const std::size_t jjKkUpdateSp = nextUpdateSp[jj];
            const std::size_t jjKkHeight   = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth    = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> updateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateDynamicSupernode(nf, jj, jjKkUpdateSp, updateBlock);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);

            // jj has updated kk. Queue it for the next ancestor it must update.
            nextUpdateSp[jj] = jjKkUpdateSp + jjKkWidth;
            if (nextUpdateSp[jj] < jjNumNodeIdx)
                descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        }

        if (!factorDynamicSupernode(nf, kk, gblToLcl)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // kk now has updates of its own to deliver. **The advance is over the front and the delayed columns
        // together**, not the front alone: both are kk's own rows, neither updates an ancestor, and
        // the delayed ones are handed over by assembleDelay rather than as an update. Getting
        // this wrong sends the delayed rows into a temporary and corrupts the parent quietly.
        const std::size_t kkDelaySize = nf.delaySize(kk);
        nextUpdateSp[kk] = nf.frontSize(kk) + kkDelaySize;
        if (nextUpdateSp[kk] < kkNumNodeIdx)
            descendantUpdateQueue[nf.nodeToSnode(kkNodeIdx[nextUpdateSp[kk]])].push_back(kk);

        // Clear the whole height, delayed columns included. 0.9 clears only frontSize + updateSize
        // here, and since frontSize has just contracted it leaves the delayed entries stale.
        // Harmless
        // there, because the map is only ever read at indices known to be in the current
        // supernode's set, but it costs the array its stated invariant (NIL everywhere outside the
        // supernode in hand), and assembleFromA's NIL test is exactly a check that relies on it.
        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A root has no parent to take a delayed column, so a delay there is unrecoverable. 0.9
        // treats it as an error rather than a numeric failure, and so do we: it means the pivoting
        // strategy did not do its job, not that the matrix is singular.
        if (kkDelaySize > 0 && parent[kk] == NIL)
            return false;   // delayed at a root: nowhere left to put it
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
    initNumFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    const std::vector<std::int32_t>& parent      = sf.parent();
    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // A into every front, up front. This prepass is required by right-looking and is the mirror of
    // left-looking's lazy assembly (see factorDynamicLeftLooking): right-looking *pushes*, so once
    // the main loop factors a supernode it immediately writes updates into ancestors it has not yet
    // reached, and those ancestors' fronts must already hold A for the update to land on. assembleFromA
    // *assigns*, so it cannot run after any push has landed without clobbering it; hence every front
    // is filled before the loop begins. This runs while every front is still the width symbolic
    // predicted, so nothing has expanded yet and the delayed-column offset is zero everywhere.
    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        const std::size_t   jjPreFactorFrontSize = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx         = jjPreFactorFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx            = nf.nodeIdx(jj);
        Val*                jjVal                = nf.val(jj);

        setGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        if (!assembleFromA(A, p, gblToLcl, 0, jjPreFactorFrontSize,
                           jjNumNodeIdx, jjNodeIdx, jjVal)) {
            clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
            return false;
        }

        clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
    }

    for (std::int32_t jj = 0; jj < static_cast<std::int32_t>(snodeSize); ++jj) {
        std::size_t jjInboundDelaySize = 0;
        for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii])
            jjInboundDelaySize += nf.delaySize(ii);

        if (jjInboundDelaySize > 0) {
            const std::size_t jjPreExpandFrontSize = nf.frontSize(jj);
            const std::size_t jjPostExpandNumNodeIdx =
                jjInboundDelaySize + jjPreExpandFrontSize + nf.updateSize(jj);

            // The index set, exactly as in the left-looking driver: extend, shift right, prepend
            // the children's delayed globals in sibling order.
            nf.expandNodeIdx(jj, jjInboundDelaySize);
            std::vector<std::int32_t>& jjPostExpandNodeIdx = nf.mNodeIdx[jj];

            for (std::size_t ssp = jjPostExpandNumNodeIdx - jjInboundDelaySize; ssp-- > 0; ) {
                const std::size_t dsp = ssp + jjInboundDelaySize;
                jjPostExpandNodeIdx[dsp] = jjPostExpandNodeIdx[ssp];
            }

            std::size_t dsp = 0;
            for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii]) {
                const std::size_t   iiPostFactorFrontSize = nf.frontSize(ii);
                const std::int32_t* iiNodeIdx             = nf.nodeIdx(ii);

                for (std::size_t ssp = iiPostFactorFrontSize;
                     ssp < iiPostFactorFrontSize + nf.delaySize(ii); ++ssp, ++dsp)
                    jjPostExpandNodeIdx[dsp] = iiNodeIdx[ssp];
            }

            // And the val, keeping what A and the descendants already put there.
            nf.mFrontSize[jj] += jjInboundDelaySize;
            nf.expandVal(jj, jjInboundDelaySize);
        }

        // The full height, captured before the factorization reclassifies part of the front as
        // delayed. The map is computed once and serves both the delayed assembly and the factor,
        // since jj's index set does not change between them.
        const std::size_t   jjPreFactorFrontSize = nf.frontSize(jj);
        const std::size_t   jjNumNodeIdx         = jjPreFactorFrontSize + nf.updateSize(jj);
        const std::int32_t* jjNodeIdx            = nf.nodeIdx(jj);

        setGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        // Take every child's delayed columns, then let the child reclaim their storage. Tested per
        // child, not just on the total: one child delaying does not mean its siblings did, and the
        // calls are no-ops for those that did not.
        if (jjInboundDelaySize > 0)
            for (std::int32_t ii = firstChild[jj]; ii != NIL; ii = nextSibling[ii]) {
                const std::size_t iiDelaySize = nf.delaySize(ii);
                if (iiDelaySize > 0) {
                    assembleDelay(nf, ii, jj, gblToLcl);
                    nf.contractVal(ii, iiDelaySize);
                }
            }

        if (!factorDynamicSupernode(nf, jj, gblToLcl)) {
            clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);
            return false;
        }

        clearGlobalToLocal(jjNumNodeIdx, jjNodeIdx, gblToLcl);

        // From here jj is factored: part of its front was reclassified as delayed, so the front is
        // now smaller. Everything below reads this post-factor front, never nf.frontSize(jj) bare,
        // which would hide whether it is the pre- or post-factor value.
        const std::size_t jjPostFactorFrontSize = nf.frontSize(jj);
        const std::size_t jjDelaySize           = nf.delaySize(jj);

        // Push. The walk starts past the front *and* the delayed columns, for the same reason the
        // left-looking nextUpdateSp seed does: both are jj's own rows and neither updates an
        // ancestor, the delayed ones going up by assembleDelay instead.
        std::size_t jjKkUpdateSp = jjPostFactorFrontSize + jjDelaySize;
        while (jjKkUpdateSp < jjNumNodeIdx) {
            const std::int32_t kk = nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp]);

            // kk has not expanded yet, and need not have: jj's update rows are kk's own nodes, which
            // its index set already holds. When kk later expands, expandVal carries these values
            // along with the rest.
            const std::size_t   kkPreFactorFrontSize = nf.frontSize(kk);
            const std::size_t   kkNumNodeIdx         = kkPreFactorFrontSize + nf.updateSize(kk);
            const std::int32_t* kkNodeIdx            = nf.nodeIdx(kk);
            Val*                kkVal                = nf.val(kk);

            const std::size_t jjKkHeight = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth  = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> updateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, updateBlock.mRowIdx.begin());

            updateDynamicSupernode(nf, jj, jjKkUpdateSp, updateBlock);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdate(gblToLcl, updateBlock, kkNumNodeIdx, kkVal);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjKkUpdateSp += jjKkWidth;
        }

        if (jjDelaySize > 0 && parent[jj] == NIL)
            return false;   // delayed at a root: nowhere left to put it
    }

    return true;
}

// =================================================================================================
// Dynamic LDL, multifrontal. The left-looking dynamic skeleton (expand a front by its children's
// delayed columns, assemble A past them, factor, which may delay again) with the update stack in
// place of the pull queue: fold each child by both halves of the extend-add, then factor and leave
// this supernode's contribution block on the stack. This is where delayed columns meet the stack.
// =================================================================================================

template<class Val>
bool NumFactorEngine::factorDynamicMultifrontal(const SparseMatrix<Val>& A, const Permutation& p,
                                                const SymFactor& sf, NumFactorDynamic<Val>& nf) const {
    initNumFactor(sf, nf);

    const std::size_t size      = nf.size();
    const std::size_t snodeSize = nf.snodeSize();

    const std::vector<std::int32_t>& parent      = sf.parent();
    const std::vector<std::int32_t>& firstChild  = sf.firstChild();
    const std::vector<std::int32_t>& nextSibling = sf.nextSibling();

    std::vector<std::int32_t> gblToLcl(size, NIL);

    // One contribution block per supernode, allocated when the supernode is reached and freed once
    // its parent has folded it in.
    std::vector<UpdateMatrix<Val>> stack(snodeSize);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        std::size_t kkInboundDelaySize = 0;
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj])
            kkInboundDelaySize += nf.delaySize(jj);

        if (kkInboundDelaySize > 0) {
            // Expand kk's front by its children's delayed columns: extend the index set, shift kk's
            // own indices right, prepend the children's delayed globals in sibling order, then
            // discard and re-zero the val at the wider shape. Identical to the left-looking driver;
            // multifrontal also discards (nothing is in kk yet), so resetVal, not expandVal.
            const std::size_t kkPreExpandFrontSize = nf.frontSize(kk);
            const std::size_t kkPostExpandNumNodeIdx =
                kkInboundDelaySize + kkPreExpandFrontSize + nf.updateSize(kk);

            nf.expandNodeIdx(kk, kkInboundDelaySize);
            std::vector<std::int32_t>& kkPostExpandNodeIdx = nf.mNodeIdx[kk];

            for (std::size_t ssp = kkPostExpandNumNodeIdx - kkInboundDelaySize; ssp-- > 0; ) {
                const std::size_t dsp = ssp + kkInboundDelaySize;
                kkPostExpandNodeIdx[dsp] = kkPostExpandNodeIdx[ssp];
            }

            std::size_t dsp = 0;
            for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
                const std::size_t   jjPostFactorFrontSize = nf.frontSize(jj);
                const std::int32_t* jjNodeIdx             = nf.nodeIdx(jj);
                for (std::size_t ssp = jjPostFactorFrontSize;
                     ssp < jjPostFactorFrontSize + nf.delaySize(jj); ++ssp, ++dsp)
                    kkPostExpandNodeIdx[dsp] = jjNodeIdx[ssp];
            }

            nf.mFrontSize[kk] += kkInboundDelaySize;
            nf.resetVal(kk);
        }

        // The full height, captured before the factorization reclassifies part of the front as
        // delayed. As in the left-looking driver, every use below wants this and not the contracted
        // front.
        const std::size_t   kkPreFactorFrontSize = nf.frontSize(kk);
        const std::size_t   kkUpdateSize         = nf.updateSize(kk);
        const std::size_t   kkNumNodeIdx         = kkPreFactorFrontSize + kkUpdateSize;
        const std::int32_t* kkNodeIdx            = nf.nodeIdx(kk);
        Val*                kkVal                = nf.val(kk);

        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A into kk's lu block, past the delayed columns (those come from children, not A).
        if (!assembleFromA(A, p, gblToLcl, kkInboundDelaySize,
                           kkPreFactorFrontSize, kkNumNodeIdx, kkNodeIdx, kkVal)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Allocate kk's contribution block over its update rows; allocate zeroes the values, so only
        // the indices are filled here.
        stack[kk].allocate(kkUpdateSize);
        {
            std::int32_t* kkUpdateIdx = stack[kk].nodeIdx();
            for (std::size_t sp = 0; sp < kkUpdateSize; ++sp)
                kkUpdateIdx[sp] = kkNodeIdx[kkPreFactorFrontSize + sp];
        }

        // Fold in each child. Its delayed columns become kk front columns (assembleDelay) and their
        // storage is then reclaimed (contractVal); that pair is the child's delayed-column handling,
        // a no-op for a child that delayed nothing. Its contribution block extend-adds into kk's
        // frontal (assembleUpdateMatrix) and is then freed. The extend-add is independent of the
        // contraction, it reads the block on the stack rather than jj's front, so the only order that
        // matters is assembleDelay before contractVal (read the columns before reclaiming them) and
        // the extend-add before discard.
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
            const std::size_t jjDelaySize = nf.delaySize(jj);
            if (jjDelaySize > 0) {
                assembleDelay(nf, jj, kk, gblToLcl);
                nf.contractVal(jj, jjDelaySize);
            }
            assembleUpdateMatrix(stack[jj], kk, nf, stack[kk], gblToLcl);
            stack[jj].discard();
        }

        // Factor kk's pivots (dynamic threshold pivoting; may delay columns to kk's parent).
        if (!factorDynamicSupernode(nf, kk, gblToLcl)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Form kk's contribution block from the factored pivots and leave it on the stack.
        updateDynamicMultifrontal(nf, kk, stack[kk]);

        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A delay that reaches a root has nowhere to go: a correct pivoting strategy prevents it.
        if (nf.delaySize(kk) > 0 && parent[kk] == NIL)
            return false;
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
        case Traversal::Multifrontal: return factorStaticMultifrontal(A, p, sf, nf);
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
            case Traversal::Multifrontal: return factorDynamicMultifrontal(A, p, sf, nf);
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
