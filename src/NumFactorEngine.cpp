#include "oblio/NumFactorEngine.h"

#include "oblio/BlasLapack.h"

#include <algorithm>
#include <cmath>
#include <list>

namespace Oblio {

namespace {

// The 2x2 pivot val D, read out of a front at local columns k1 and k2.
//
// **This is the one place the symmetry of D is decided, and that is the point of it existing.**
// `d12 = d21` is the symmetric statement, and it is what `LDL^T` means over the reals and over the
// complex field alike. A Hermitian factorization wants the conjugate here instead, and wants `d11`
// and `d22` known real. Everything downstream, the acceptance test and the elimination, works from
// what this returns, so that change is made once rather than in four places that must not drift.
//
// The four entries are all of it. The determinant belongs to the acceptance test, which is the only
// thing that reads it, and is formed there.
template<class Val>
struct PivotBlock2x2 {
    Val d11, d22, d21, d12;
};

// The dense pivot threshold, (1 + sqrt(17)) / 8, from AGL section 2.1, where it is chosen to
// equalize the worst-case growth of one 2x2 step with that of two consecutive 1x1 steps. Unlike
// mPivotThreshold it is a constant and not a setting, because the trade a setting would express
// does not exist at a root: the front is dense, stays dense under any symmetric permutation, and no
// column can be delayed, so a larger threshold costs comparisons and never fill.
constexpr double alphaRoot = 0.64038820320220757;

// The scans behind pivot selection, one per kernel, because the two want different things from a
// column. Both walk column j of jj's block over the unfactored positions only, from nextPivot on,
// which is the row segment to the left of the diagonal and the column segment below it.

// Root: one maximum and where it sits. There are no update rows at a root, so the scan stops at the
// end of the front and every position it sees is a legal pivot partner. `r` is AGL's dense letter,
// Figure 2.1 onward; `q` would import the A11 restriction, which does not exist here.
struct RootPivotColumnScan {
    double       gamma;   // largest off-diagonal magnitude in the column, 0 if there are none
    std::int32_t r;       // where it sits, NIL if the column has no unfactored off-diagonal
};

template<class Val>
RootPivotColumnScan scanRootPivotColumn(const Val* val, std::ptrdiff_t ld, std::int32_t nextPivot,
                                        std::int32_t j, std::int32_t frontSize) {
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    // gamma starts at 0, not at a negative sentinel, because the maximum over an empty set of
    // magnitudes is what the caller's "nothing to eliminate" test means. The last column of a front
    // scans two empty ranges, and a sentinel there would slip past that test and take the pivot
    // without the zero-diagonal check, which is how it read before 2026-07-26.
    RootPivotColumnScan s{0, NIL};

    for (std::int32_t i = nextPivot; i < j; ++i) {          // j's row, to the left
        const double m = std::abs(val[at(j, i)]);
        if (s.gamma < m) { s.r = i; s.gamma = m; }
    }
    for (std::int32_t i = j + 1; i < frontSize; ++i) {      // its column, below
        const double m = std::abs(val[at(i, j)]);
        if (s.gamma < m) { s.r = i; s.gamma = m; }
    }

    return s;
}

// Non-root: two maxima and one position. `gamma` runs the full block height, update rows included,
// because the entries of L this pivot writes reach into them and a shorter maximum would not bound
// them. `frontGamma` is the same maximum stopped at the end of the front, and `q` is where it sits,
// so a 2x2 partner is always a column this front can eliminate; an update row has no column here
// and no diagonal, so it can be measured but not pivoted on. AGL section 3.1 draws exactly this
// distinction, writing the restricted one as an entry, a_q1, rather than as a gamma.
struct NonRootPivotColumnScan {
    double       gamma;        // full block height
    double       frontGamma;   // front columns only
    std::int32_t q;            // where frontGamma sits, NIL if the front holds no partner
};

// `frontGamma` is a prefix of `gamma`, so the third loop continues from it rather than starting
// over, and simply stops tracking `q`: that omission is the whole of the front-partner rule.
template<class Val>
NonRootPivotColumnScan scanNonRootPivotColumn(const Val* val, std::ptrdiff_t ld,
                                              std::int32_t nextPivot, std::int32_t j,
                                              std::int32_t frontSize, std::int32_t numNodeIdx) {
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    NonRootPivotColumnScan s{-1, -1, NIL};

    for (std::int32_t i = nextPivot; i < j; ++i) {          // j's row, to the left
        const double m = std::abs(val[at(j, i)]);
        if (s.frontGamma < m) { s.q = i; s.frontGamma = m; }
    }
    for (std::int32_t i = j + 1; i < frontSize; ++i) {      // its column, front part
        const double m = std::abs(val[at(i, j)]);
        if (s.frontGamma < m) { s.q = i; s.frontGamma = m; }
    }

    s.gamma = s.frontGamma;
    for (std::int32_t i = frontSize; i < numNodeIdx; ++i) { // its column, update rows
        const double m = std::abs(val[at(i, j)]);
        if (s.gamma < m) s.gamma = m;
    }

    return s;
}

template<class Val>
PivotBlock2x2<Val> readPivotBlock2x2(const Val* val, std::ptrdiff_t ld,
                                     std::int32_t k1, std::int32_t k2, bool withHermitian) {
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    // forceReal restores what a Hermitian diagonal is mathematically and floating point is not; see
    // its definition in Types.h. It matters here as well as at the point of use: the acceptance
    // test compares magnitudes against this block's determinant, and that determinant is real only
    // if these two are.
    PivotBlock2x2<Val> d;
    d.d11 = forceReal(val[at(k1, k1)], withHermitian);
    d.d22 = forceReal(val[at(k2, k2)], withHermitian);

    // Only the lower triangle is occupied before the front is factored, so exactly one of the two
    // off-diagonal positions is stored, and **which one it is depends on the order of k1 and k2**.
    // That distinction is invisible in the symmetric case, where the two are equal, and matters
    // here: the stored entry is d21 when k2 sits below k1 and d12 when it sits above, and the other
    // is its conjugate.
    if (k2 > k1) {
        d.d21 = val[at(k2, k1)];
        d.d12 = maybeConjugate(d.d21, withHermitian);
    } else {
        d.d12 = val[at(k1, k2)];
        d.d21 = maybeConjugate(d.d12, withHermitian);
    }

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

    // Rank starts at full and drops per zero pivot, as 0.9's rank_ = a.getSize() then rank_--.
    nf.mRank = nf.mSize;

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

template<class Val, class Factor>
void NumFactorEngine::assembleUpdateBlock(const UpdateBlock<Val>& jjKkUpdateBlock, Factor& nf,
                                          std::int32_t kk,
                                          const std::vector<std::int32_t>& gblToLcl) const {
    const std::size_t   jjKkHeight    = jjKkUpdateBlock.height();
    const std::size_t   jjKkWidth     = jjKkUpdateBlock.width();
    const std::int32_t* jjKkRowIdx    = jjKkUpdateBlock.rowIdx();
    const Val*          jjKkUpdateVal = jjKkUpdateBlock.val();

    // kk is always unfactored when this runs, so its front is at full width and the destination
    // derives cleanly, the same two values every caller used to pass.
    const std::size_t kkNumNodeIdx = nf.frontSize(kk) + nf.updateSize(kk);
    Val*              kkFrontVal   = nf.val(kk);

    // The update block's rows and columns carry global row indices; gblToLcl maps them into kk's
    // local coordinates. Only its lower triangle is meaningful (row at or below column), which is
    // exactly the part the two BLAS calls filled. Entry (si, sj) of the update sits at
    // jjKkUpdateVal[sj * jjKkHeight + si], the same column-major layout the solve reads, and
    // (di, dj) is where that entry lands in kk's front. The global (L-ordering) index in between
    // is what gblToLcl converts.
    //
    // si, sj, di, dj are column/row indices (int32_t), as in the solve; the two positions here are
    // the flat column offsets into the source and destination vals, both size_t in the cp family.

    for (std::int32_t sj = 0; sj < static_cast<std::int32_t>(jjKkWidth); ++sj) {
        const std::int32_t lj = jjKkRowIdx[sj];                // this column's global index
        const std::int32_t dj = gblToLcl[lj];                  // its column in kk's front

        const std::size_t scp = sj * jjKkHeight;     // source column
        const std::size_t dcp = dj * kkNumNodeIdx;   // dest column

        for (std::int32_t si = sj; si < static_cast<std::int32_t>(jjKkHeight); ++si) {
            const std::int32_t li = jjKkRowIdx[si];
            const std::int32_t di = gblToLcl[li];

            kkFrontVal[dcp + di] += jjKkUpdateVal[scp + si];
        }
    }
}

template<class Val, class Factor>
void NumFactorEngine::assembleUpdateMatrix(const UpdateMatrix<Val>& jjUpdateMatrix, Factor& nf,
                                           std::int32_t kk, UpdateMatrix<Val>& kkUpdateMatrix,
                                           const std::vector<std::int32_t>& gblToLcl) const {
    const std::size_t   jjUpdateSize = jjUpdateMatrix.size();
    const std::int32_t* jjRowIdx     = jjUpdateMatrix.rowIdx();
    const Val*          jjUpdateVal  = jjUpdateMatrix.val();   // jjUpdateSize square, ld the same

    const std::size_t kkFrontSize  = nf.frontSize(kk);
    const std::size_t kkUpdateSize = nf.updateSize(kk);
    const std::size_t kkNumNodeIdx = kkFrontSize + kkUpdateSize;
    Val*              kkFrontVal   = nf.val(kk);              // lu block, ld == kkNumNodeIdx
    Val*              kkUpdateVal  = kkUpdateMatrix.val();    // ld == kkUpdateSize

    // Each column of jj's update matrix carries a global row index that lies somewhere in kk's
    // index set. Where it lands decides which block receives it: a pivot column of kk (local
    // position below kkFrontSize) goes into the lu block, an update row of kk into kk's own update
    // matrix. jj's row indices are sorted and gblToLcl preserves order, so within a column the rows
    // run at or below the diagonal, and only the lower triangle is written.
    //
    // Same shape and same names as assembleUpdateBlock: sj and si index jj's columns and rows, lj
    // and li are the global (L-ordering) indices they carry, dj and di are where those land in kk,
    // and scp and dcp are the flat column positions in the source and destination vals. The sizes
    // measure, so they are size_t, and the comparisons against them carry the one cast that meets
    // the int32_t indices.
    for (std::int32_t sj = 0; sj < static_cast<std::int32_t>(jjUpdateSize); ++sj) {
        const std::int32_t lj  = jjRowIdx[sj];   // this column's global (L-ordering) index
        const std::int32_t dj  = gblToLcl[lj];   // its column in kk
        const std::size_t  scp = sj * jjUpdateSize;   // source column

        if (dj < static_cast<std::int32_t>(kkFrontSize)) {
            const std::size_t dcp = dj * kkNumNodeIdx;   // dest column

            for (std::int32_t si = sj; si < static_cast<std::int32_t>(jjUpdateSize); ++si) {
                const std::int32_t li = jjRowIdx[si];
                const std::int32_t di = gblToLcl[li];

                kkFrontVal[dcp + di] += jjUpdateVal[scp + si];
            }
        } else {
            const std::size_t dcp = (dj - kkFrontSize) * kkUpdateSize;   // dest column

            for (std::int32_t si = sj; si < static_cast<std::int32_t>(jjUpdateSize); ++si) {
                const std::int32_t li = jjRowIdx[si];
                const std::int32_t di = gblToLcl[li];

                kkUpdateVal[dcp + di - kkFrontSize] += jjUpdateVal[scp + si];
            }
        }
    }
}

template<class Val, class Factor>
bool NumFactorEngine::factorStaticSupernode(Factor& nf, std::int32_t jj) const {
    const std::size_t jjFrontSize  = nf.frontSize(jj);
    const std::size_t jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
    Val*              jjVal        = nf.val(jj);

    const int f  = static_cast<int>(jjFrontSize);
    const int u  = static_cast<int>(jjNumNodeIdx - jjFrontSize);
    const int ld = static_cast<int>(jjNumNodeIdx);

    if (mFactorization == Factorization::Cholesky) {
        // The front, which is the diagonal val: A11 = L11 L11^H.
        int info = 0;
        potrf('L', f, jjVal, ld, &info);
        if (info > 0)
            return false;   // not positive definite: the leading minor of order `info` failed

        // The update rows: L21 = A21 (L11^H)^-1, solved in place.
        //
        // Blas<Val>::conjTrans is 'T' for real and 'C' for complex, which is the whole of what the
        // scalar type decides here. 0.9 writes 'T' unconditionally, which is wrong for a complex
        // Hermitian factor, and there is nothing at its call site to reveal that.
        if (u > 0)
            trsm('R', 'L', Blas<Val>::conjTrans, 'N', u, f, Val(1), jjVal, ld, jjVal + f, ld);

        return true;
    }

    // LDL. The kernel is ours: LAPACK has no unpivoted LDL^T (its ?sytrf pivots, and pivoting is
    // what a *static* factorization refuses to do). It cannot fail, because there is no positive
    // definiteness to violate; a pivot too small to divide by is perturbed and counted.
    int numPert = 0;
    ldl(f, jjVal, ld, mPerturbation, &numPert, hermitian(mFactorization));
    nf.numPerturbations() += static_cast<std::size_t>(numPert);

    // The update rows: L21 = A21 U11^-1, where U11 = D11 L11^H sits in the front's *upper*
    // triangle. So the solve is against the upper, untransposed, which is exactly what storing U
    // buys: Cholesky would have to transpose, and does.
    if (u > 0)
        trsm('R', 'U', 'N', 'N', u, f, Val(1), jjVal, ld, jjVal + f, ld);

    return true;
}

template<class Val, class Factor>
void NumFactorEngine::updateStaticUpdateBlock(const Factor& nf, std::int32_t jj,
                                      std::size_t jjKkUpdateSp, UpdateBlock<Val>& jjKkUpdateBlock) const {
    // `jjKkUpdateSp` is a position into jj's index set, and it is what selects the ancestor kk this
    // update reaches: kk == nf.nodeToSnode(nf.nodeIdx(jj)[jjKkUpdateSp]). The kernel never names kk
    // (it forms jj's outer product from that row down, and the caller has already sized the block to
    // the rows kk owns), but that is the position's meaning and it is worth seeing first.
    const std::size_t jjFrontSize  = nf.frontSize(jj);
    const std::size_t jjNumNodeIdx = jjFrontSize + nf.updateSize(jj);
    const Val*        jjVal        = nf.val(jj);

    const int f          = static_cast<int>(jjFrontSize);
    const int sld        = static_cast<int>(jjNumNodeIdx);
    const int jjKkHeight = static_cast<int>(jjKkUpdateBlock.height());
    const int jjKkWidth  = static_cast<int>(jjKkUpdateBlock.width());
    const int dld        = jjKkHeight;

    const bool withHermitian = hermitian(mFactorization);
    const Val* l21Val        = jjVal + jjKkUpdateSp;   // the update rows that reach this ancestor, and below
    Val*       u22Val        = jjKkUpdateBlock.val();

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
    // stored intermediate, only C21. LDL has to fold the D in somewhere, and materializing
    // U12 := D L21^H into scratch is one of three places it could go, not a necessity. The
    // alternative worth naming is a fused kernel, an `oblioHerk` that forms and consumes D L21^H
    // internally and stores nothing; with the D hoisted per (row, column) it costs the same
    // arithmetic and no buffer. What it costs instead is the gemm below, which is a real library
    // call only because U12 exists in memory for it to read, and a second hand-written kernel to go
    // with gemmLower. See "The scratch upper, and why it is not stored" in docs/ARCHITECTURE.md.
    //
    // Nothing of U12 is kept. It is consumed by the two multiplies and freed on return, and it is a
    // different object from U11 = D L11^H, which does persist, in the front's own upper triangle,
    // where factorStaticSupernode's trsm solves against it for L21.
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
void NumFactorEngine::updateStaticUpdateMatrix(const Factor& nf, std::int32_t kk,
                                               UpdateMatrix<Val>& kkUpdateMatrix) const {
    const int u = static_cast<int>(nf.updateSize(kk));
    if (u == 0)
        return;   // a root, or a supernode reaching nowhere: no contribution to leave

    const int  f   = static_cast<int>(nf.frontSize(kk));
    const int  ld  = f + u;                 // kk's lu block height
    const Val* val = nf.val(kk);            // kk's factored block
    const Val* l21 = val + f;               // its update rows
    Val*       uVal = kkUpdateMatrix.val();       // the contribution block, ld == u

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

            UpdateBlock<Val> jjKkUpdateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, jjKkUpdateBlock.rowIdx());

            updateStaticUpdateBlock(nf, jj, jjKkUpdateSp, jjKkUpdateBlock);
            assembleUpdateBlock(jjKkUpdateBlock, nf, kk, gblToLcl);

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

            const std::size_t jjKkHeight = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth  = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> jjKkUpdateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, jjKkUpdateBlock.rowIdx());

            updateStaticUpdateBlock(nf, jj, jjKkUpdateSp, jjKkUpdateBlock);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdateBlock(jjKkUpdateBlock, nf, kk, gblToLcl);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjKkUpdateSp += jjKkWidth;
        }
    }

    return true;
}

// =================================================================================================
// Static multifrontal. One postorder-compatible pass (supernodes in increasing order): assemble A
// and every child's contribution block into the frontal, factor, then leave this supernode's own
// contribution block in updateMatrix for its parent. Cholesky only for now; see the header.
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
    std::vector<UpdateMatrix<Val>> updateMatrix(snodeSize);

    for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(snodeSize); ++kk) {
        const std::size_t   kkFrontSize  = nf.frontSize(kk);
        const std::size_t   kkUpdateSize = nf.updateSize(kk);
        const std::size_t   kkNumNodeIdx = kkFrontSize + kkUpdateSize;
        const std::int32_t* kkNodeIdx    = nf.nodeIdx(kk);
        Val*                kkVal        = nf.val(kk);

        // Map kk's whole index set, front and update rows alike: assembleFromA fills the front
        // columns, and the assembly routes each child entry by where its index falls in kk.
        setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

        // A's own values into kk's lu block.
        if (!assembleFromA(A, p, gblToLcl, 0, kkFrontSize, kkNumNodeIdx, kkNodeIdx, kkVal)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;
        }

        // Allocate kk's contribution block over its update rows; allocate zeroes the values, so only
        // the indices are filled here, from kk's update-row indices.
        updateMatrix[kk].allocate(kkUpdateSize);
        {
            std::int32_t* kkUpdateRowIdx = updateMatrix[kk].rowIdx();
            for (std::size_t sp = 0; sp < kkUpdateSize; ++sp)
                kkUpdateRowIdx[sp] = kkNodeIdx[kkFrontSize + sp];
        }

        // Assemble each child's contribution block into kk's frontal, then free it. A child of kk
        // is factored in an earlier iteration (increasing order), so its block is present here.
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
            assembleUpdateMatrix(updateMatrix[jj], nf, kk, updateMatrix[kk], gblToLcl);
            updateMatrix[jj].discard();
        }

        // Factor kk's pivots, the same per-supernode kernel the other traversals use: Cholesky is
        // POTRF then TRSM, static LDL is the unpivoted ldl then TRSM against the upper triangle.
        if (!factorStaticSupernode<Val>(nf, kk)) {
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            return false;   // Cholesky non-positive-definite; LDL perturbs instead and cannot fail
        }

        // Form kk's contribution block from the factored pivots and leave it for kk's parent.
        updateStaticUpdateMatrix(nf, kk, updateMatrix[kk]);

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

// The 2x2 acceptance test: AGL Figure 3.3's fourth branch, restated with the candidate as j and its
// partner as q. That restatement is the paper's own; Figure 3.6 carries the same test with generic
// indices when it embeds it in a search.
//
// Read-only by construction: judging a block and eliminating it are different jobs, so this takes
// the factor by const reference where factor2x2 takes it by non-const.
//
// 0.9 also accepted on a symmetric-maximum clause, `max == max1 && max == max2 && max != 0`, which
// is not in Figure 3.3 and was removed on 2026-07-26. It reached for the Bunch-Parlett condition of
// AGL section 2.3, which needs the off-diagonal to be a local maximum in both columns *and* both
// diagonals to be at most alpha times it. Only the candidate's diagonal is known small here, from
// the 1x1 test that failed to get us this far; the partner's was never tested, and without it the
// bound on L does not hold. Adding the missing conjunct would have made the clause redundant rather
// than correct, since under both conditions the determinant test below already accepts the block
// for any threshold at or below one half, so the clause was deleted rather than repaired. The root
// kernel keeps the same predicate, where the chase supplies the missing conjunct by control flow.
template<class Val>
bool NumFactorEngine::acceptPivot2x2(const Val* jjVal, std::int32_t jjNumNodeIdx,
                                     std::int32_t nextPivot, std::int32_t j, std::int32_t q,
                                     double jGamma, double jFrontGamma, bool withHermitian) const {
    // The block, not the factor: this is called once per 2x2 candidate, rejections included, and
    // everything it would have re-derived from nf is already in the caller's hand. jjNumNodeIdx is
    // the leading dimension as well as the scan bound, so it arrives once rather than twice.
    const std::ptrdiff_t ld = static_cast<std::ptrdiff_t>(jjNumNodeIdx);

    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    // The partner column's largest off-diagonal magnitude, over the full block height. Read here
    // rather than in the caller because this test is the only thing that wants it.
    double qGamma = -1;
    for (std::int32_t i = nextPivot; i < q; ++i) {
        const double m = std::abs(jjVal[at(q, i)]);
        if (qGamma < m) qGamma = m;
    }
    for (std::int32_t i = q + 1; i < jjNumNodeIdx; ++i) {
        const double m = std::abs(jjVal[at(i, q)]);
        if (qGamma < m) qGamma = m;
    }

    const PivotBlock2x2<Val> d = readPivotBlock2x2(jjVal, ld, j, q, withHermitian);

    // The determinant lives here rather than in the block, because it is what this test needs and
    // nothing else reads it: factor2x2 works from the four entries alone. Hermitian: d11 * d22 is
    // real and d12 * d21 is |d21|^2, so this is real, as the comparison below assumes.
    const Val det = d.d11 * d.d22 - d.d12 * d.d21;

    // The growth bound the determinant is tested against: the larger of the two ways of pairing
    // each diagonal with the other column's maximum.
    const double maxmax = std::max(std::abs(d.d22) * jGamma + jFrontGamma * qGamma,
                                   std::abs(d.d11) * qGamma + jFrontGamma * jGamma);

    return std::abs(det) > 0 && std::abs(det) >= mPivotThreshold * maxmax;
}

template<class Val>
void NumFactorEngine::factor1x1(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t nextPivot,
                                std::int32_t j1, std::int32_t lj1,
                                std::size_t jjPreFactorFrontSize, std::size_t jjNumNodeIdx,
                                std::vector<std::int32_t>& gblToLcl) const {
    Val*                 jjVal = nf.val(jj);
    const std::ptrdiff_t ld    = static_cast<std::ptrdiff_t>(jjNumNodeIdx);

    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    const bool withHermitian = hermitian(nf.factorization());

    // The position this pivot is going to. Named as in factor2x2, where there are two of them.
    const std::int32_t np0 = nextPivot;

    // Read before the swap, which is why it is not simply jjVal[at(np0, np0)]: the swap is what
    // puts the pivot there.
    //
    // forceReal restores what a Hermitian diagonal is mathematically and floating point is not; see
    // its definition in Types.h. **This is the last moment it can be applied**, the value being
    // about to be divided by below and stored for the solve to divide by again.
    const Val d = forceReal(jjVal[at(j1, j1)], withHermitian);

    if (np0 != j1) nf.swap(jj, np0, j1, gblToLcl);

    // The corrected value replaces the stored one, because the corrected one is the right one and
    // this is the last place it can be put right: from here it is read by the elimination below,
    // which divides this column by it, and by the solve, which divides again. Neither could correct
    // it, having no way to know it should be real. As a store it is redundant exactly when there
    // was no swap and the factorization is not Hermitian, and it is one store per pivot either way.
    jjVal[at(np0, np0)] = d;

    // L column: divide by the pivot
    for (std::int32_t i = np0 + 1; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
        jjVal[at(i, np0)] /= d;

    // D L^H row, in the upper part
    for (std::int32_t k = np0 + 1; k < static_cast<std::int32_t>(jjPreFactorFrontSize); ++k)
        jjVal[at(np0, k)] = d * maybeConjugate(jjVal[at(k, np0)], withHermitian);

    // rank-1 trailing update
    for (std::int32_t k = np0 + 1; k < static_cast<std::int32_t>(jjPreFactorFrontSize); ++k)
        for (std::int32_t i = k; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
            jjVal[at(i, k)] -= jjVal[at(i, np0)] * jjVal[at(np0, k)];

    nf.mPivotType[lj1] = 1;
}

template<class Val>
void NumFactorEngine::factor2x2(NumFactorDynamic<Val>& nf, std::int32_t jj, std::int32_t nextPivot,
                                std::int32_t j1, std::int32_t j2, std::int32_t lj1,
                                std::int32_t lj2, std::size_t jjPreFactorFrontSize,
                                std::size_t jjNumNodeIdx,
                                std::vector<std::int32_t>& gblToLcl) const {
    Val*                 jjVal = nf.val(jj);
    const std::ptrdiff_t ld    = static_cast<std::ptrdiff_t>(jjNumNodeIdx);

    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    const bool withHermitian = hermitian(nf.factorization());

    // Read before the swaps, as above. The caller read the same jjVal to decide on it; nothing
    // touches the front in between.
    const PivotBlock2x2<Val> d = readPivotBlock2x2(jjVal, ld, j1, j2, withHermitian);

    // The two positions this pivot is going to, nextPivot and the one after it.
    const std::int32_t np0 = nextPivot;
    const std::int32_t np1 = nextPivot + 1;

    // Bring lj1 and lj2 to the front's next two columns. Which swaps are needed depends on where
    // they already are relative to each other, and doing them in the wrong order would undo one.
    if (j1 < j2) {
        if (!(np0 == j1 && np1 == j2)) {
            if (np0 != j1) nf.swap(jj, np0, j1, gblToLcl);
            nf.swap(jj, np1, j2, gblToLcl);
        }
    } else {
        if (np0 == j2 && np1 == j1) {
            nf.swap(jj, np0, np1, gblToLcl);
        } else {
            if (np1 != j2) nf.swap(jj, np1, j2, gblToLcl);
            nf.swap(jj, np0, j1, gblToLcl);
        }
    }

    // D's own four entries, put where they belong. This is the last place the diagonal pair can be
    // put right: from here d11 and d22 are read by the elimination below, which solves this block's
    // columns against them, and by the solve, which solves against them again, and neither could
    // correct them, having no way to know they should be real. The four divide into three cases and
    // only the last is unconditional:
    //
    //   d21, at (np1, np0)    already there, carried by the swaps, conjugated on the way if
    //                       Hermitian. Nothing to write.
    //   d11 and d22         already there too, but as the raw accumulated values, which for a
    //                       Hermitian D are wrong: readPivotBlock2x2 applied forceReal and these
    //                       writes are what make that stick. See forceReal in Types.h. As stores
    //                       they are redundant exactly when there was no swap and the factorization
    //                       is not Hermitian.
    //   d12, at (np0, np1)    never redundant. Only the lower triangle is populated before a pivot
    //                       is applied, so this position has held nothing until now, and
    //                       readPivotBlock2x2 reconstructed the value from d21. This is its only
    //                       home.
    jjVal[at(np0, np0)] = d.d11;
    jjVal[at(np1, np1)] = d.d22;
    jjVal[at(np0, np1)] = d.d12;

    // L columns: **x A = b**, with A = D and one such row system per row below the block, x being
    // that row's pair of L entries and b its pair of front values. A row system becomes a column
    // system either way round, and both identities hold:
    //
    //     x A = b   <=>   A^T x^T = b^T          unknown x^T, recovered by ^T
    //     x A = b   <=>   A^H x^H = b^H          unknown x^H, recovered by ^H
    //
    // **The first is taken, so what is factored here is A^T**, not A. The second is less a rival
    // than the same equation wearing conjugations, and for a Hermitian A it collapses into the
    // first:
    //
    //     A^H x^H = b^H                       the starting point
    //     A   x^H = b^H                       A^H = A, so A is untouched and b and x carry the conj
    //     conj(A) conj(x^H) = conj(b^H)       conjugate the whole equation
    //     A^T x^T = b^T                       conj(A) = A^T, conj(x^H) = x^T, conj(b^H) = b^T
    //
    // The conjugation has moved off b and x and onto A, and on A it is not an operation at all:
    // conj(A) = A^T is the exchange of a12 and a21, two loads naming different entries. That is the
    // whole gain, and it is worth having because b and x are per row below the block while A is
    // read once. Conjugating them would also put maybeConjugate, and with it withHermitian, inside
    // the row loop.
    //
    // For a symmetric A the same form holds with nothing to move, A^T = A making the exchange a
    // no-op, so this one line covers LDL^T and LDL^H alike and needs no maybeConjugate at all.
    //
    // The ^T route is also the symmetry-blind one, which is why there is no withHermitian branch
    // here. What the symmetry does decide is whether the a12 and a21 exchange in the factorization
    // is doing any work:
    //
    //     D real or complex symmetric    A^T = A          the exchange is a bitwise no-op
    //     D Hermitian                    A^T = conj(A)    the exchange is what makes it right
    //
    // So Hermitian fronts are the only ones that can detect an error in it. For those the ^H route
    // would have factored A itself, A^H being A, paying the two conjugations instead of the
    // transpose; neither is better, and only the ^T route is uniform across all three.
    //
    // diagonalDynamic in SolveEngine solves the other half of this same system at solve time,
    // A x = b, already a column system and so factoring A itself. The two are written to look alike
    // on purpose: same pivot choice, same four names, same two-step substitution. All that
    // distinguishes them is the transpose, and it shows up in one place, the exchange of a12 and
    // a21 in the factorization below.
    //
    // **Explicit LU with partial pivoting, not Cramer's rule.** AGL page 29 rules Cramer out for
    // this family: the 2x2 acceptance test bounds the entries of L but says nothing about the
    // condition number of D, and Cramer is backward stable only under a bounded condition number.
    // Their codes use Gaussian elimination instead. This is the same scheme diagonalDynamic runs on
    // the other half of the same system at solve time, so both halves are now solved alike; before
    // 2026-07-26 the factorization used Cramer and the solve used this, which is the arrangement
    // the paper warns against. Measured on the block shapes this test accepts, small diagonals
    // against a larger off-diagonal, Cramer's worst backward error is about fifteen times this
    // one's.
    //
    // Partial rather than complete pivoting, though either would be correct here. What the proof
    // needs depends on whether L is bounded: unbounded, as in Bunch-Kaufman, it needs the 2x2 solve
    // to be *componentwise* backward stable, which scaled Cramer and partial pivoting supply and
    // complete pivoting and the SVD do not (introduction and Appendix A). Bounded, it needs only
    // *normwise*, which all of them supply (Appendix B), and that is the regime we are in on both
    // paths, non-roots by the Figure 3.3 test and roots by bounded Bunch-Kaufman. AGL's own sparse
    // codes use complete pivoting for exactly that reason. Partial is chosen only because
    // diagonalDynamic already does it, so the two halves of this system match, and because it stays
    // valid if some later path ever admits an unbounded L.
    //
    // The factorization of A^T does not depend on the row, so it is done once here and applied
    // below. Partial pivoting is on A^T's rows, which are A's columns, so the pivoted branch loads
    // A with its columns exchanged.
    Val  a11, a12, a21, a22;
    bool swapped;
    if (std::abs(d.d11) >= std::abs(d.d12)) {
        a11 = d.d11;  a12 = d.d12;  a21 = d.d21;  a22 = d.d22;  swapped = false;
    } else {                                    // pivot: A^T's second row is the larger
        a11 = d.d12;  a12 = d.d11;  a21 = d.d22;  a22 = d.d21;  swapped = true;
    }
    const Val l21 = a12 / a11;                  // A^T = L U here, against the A = L U that
    const Val u11 = a11;                        // SolveEngine's diagonalDynamic runs: a12 and
    const Val u12 = a21;                        // a21 exchange, and nothing else differs
    const Val u22 = a22 - l21 * u12;

    for (std::int32_t i = np0 + 2; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i) {
        const Val t1 = jjVal[at(i, np0)];
        const Val t2 = jjVal[at(i, np1)];
        const Val b1 = swapped ? t2 : t1;       // the swap permutes the equations, not the unknowns
        const Val b2 = swapped ? t1 : t2;

        const Val x2 = (b2 - l21 * b1) / u22;
        const Val x1 = (b1 - u12 * x2) / u11;

        jjVal[at(i, np0)] = x1;
        jjVal[at(i, np1)] = x2;
    }

    // D L^H rows, in the upper part
    for (std::int32_t k = np0 + 2; k < static_cast<std::int32_t>(jjPreFactorFrontSize); ++k) {
        const Val l1 = maybeConjugate(jjVal[at(k, np0)], withHermitian);
        const Val l2 = maybeConjugate(jjVal[at(k, np1)], withHermitian);
        jjVal[at(np0, k)] = d.d11 * l1 + d.d12 * l2;
        jjVal[at(np1, k)] = d.d21 * l1 + d.d22 * l2;
    }

    // Rank-2 trailing update, in one pass over the trailing block rather than two rank-1 sweeps.
    // The bounds already coincided, np1 + 1 and np0 + 2 being the same column since np1 = np0 + 1,
    // so this is a merge and not a restructuring. Same flop count; what halves is traffic on the
    // trailing block, which is the large object here, while the two multipliers are invariant in i
    // and are hoisted. The merge is safe because the first sweep wrote only rows at or below k, and
    // the second read column np1 and row np1, both above it, so neither fed the other.
    for (std::int32_t k = np0 + 2; k < static_cast<std::int32_t>(jjPreFactorFrontSize); ++k) {
        const Val l1k = jjVal[at(np0, k)];
        const Val l2k = jjVal[at(np1, k)];
        for (std::int32_t i = k; i < static_cast<std::int32_t>(jjNumNodeIdx); ++i)
            jjVal[at(i, k)] -= jjVal[at(i, np0)] * l1k + jjVal[at(i, np1)] * l2k;
    }

    nf.mPivotType[lj1] = 2;
    nf.mPivotType[lj2] = 3;
}


// **The strategy here is threshold pivoting, not Bunch-Kaufman, and the two are easy to confuse.**
// Both select between 1x1 and 2x2 pivots by comparing a diagonal against off-diagonal maxima, and
// both descend from Bunch and Parlett (1971) and Bunch and Kaufman (1977), but the tests differ and
// so do their guarantees. Bunch-Kaufman turns on a fixed constant, `(1 + sqrt(17)) / 8`, chosen so
// that its 2x2 case needs no test at all, and it bounds element growth while leaving the entries of
// `L` free. What runs below is Duff and Reid's (1983), used by MA27, MA57 and MUMPS: a tunable
// `threshold`, a 2x2 determinant condition that is allowed to *refuse*, and a bound on `L` rather
// than on growth. The refusal is the part a sparse solver cannot do without, because the partner
// scan is confined to fully summed columns and so can come up empty, where Bunch-Kaufman's
// unconditional fallback assumes the whole trailing submatrix is available. A refused column is
// delayed to an ancestor instead.
//
// The confusion is worth guarding against because the Duff literature writes its threshold as
// `alpha` too, in the same position of the same inequality, where it means roughly 0.01 and not
// 0.64. See Section 7 of archive/sparse_factorization.md for both algorithms side by side.

// 0.9's factorDynamicLDL_ is one function with two passes selected by jjUpdateSize. The selector is
// exactly "is this a root": the update rows are the parent edges (ElmForestEngine builds updateSize
// and parent from the same nonzeros), so no update rows means no later column reaches this
// supernode, hence no parent. A root has no ancestor waiting, so a column it cannot pivot has
// nowhere to be delayed to, and that single fact drives everything the two kernels do differently.
//
// They are two algorithms rather than two settings of one, so they are two functions, chosen by the
// caller. That costs the callers nothing: all three traversals already bind `parent` from the
// symbolic factor, and `parent[jj] == NIL` is the same fact as `mUpdateSize[jj] == 0` read off the
// forest instead of off the storage.
//
// **Both are selection only.** The eliminations are factor1x1 and factor2x2, shared by
// both, because once a pivot is accepted the arithmetic is identical; 0.9's two copies of it were
// an artifact of writing the passes out separately. The column scans are not shared: the two want
// different things from a column, so there is one apiece.

// **Bounded Bunch-Kaufman, AGL Figure 2.4.** This is the one place the port is not a transcription.
// 0.9 ran a weakened copy of the non-root pass here, forcing its last remaining candidate as a 1x1
// whatever it looked like and accepting a 2x2 on `max1 == max2` without reading the block's
// determinant. Neither bounds the entries of L; a zero diagonal beside a large one drove max|L| to
// 1e6 and cost five digits of the solution. Replaced 2026-07-26. See archive/pivoting.md.
//
// A root front is dense, `A11` is the whole matrix, and section 3's central difficulty cannot
// arise: there is no A22 to exclude a partner from, and nothing to postpone into. What is needed
// instead is an algorithm that cannot refuse, and the chase is one. Take the next unfactored
// column, follow its largest off-diagonal to the column holding it, and repeat until either a
// diagonal passes its own 1x1 test or the maximum becomes mutual, at which point the 2x2 is taken.
// It terminates because gamma never decreases along the chase, that entry lying in both columns,
// and strictly increases except in the case that stops it, so no column is visited twice. AGL's
// Appendix C bounds the expected number of column searches at about e.
//
// On letters: the partner here is r, not q. AGL introduce q in section 3.1 for the partner
// *restricted to* A11, a distinction that exists only because the unrestricted one may sit in A21
// and be unavailable; the dense figures have no such restriction and write r from Figure 2.1 on. A
// root has no A21, so q would import a distinction this algorithm does not make, which is why the
// root has its own scan returning r rather than borrowing the non-root's q. The candidate stays j
// where Figure 2.4 writes i: eliminating a column is the mindset, and i is a row index by house
// convention.
//
// Two guarantees follow from the shape rather than from tests. The 2x2 is reached only after both
// diagonals have failed, so |det| > gamma^2 (1 - alpha^2) and the block cannot be singular; and
// kappa(D) < (1 + alpha) / (1 - alpha), which is the condition under which AGL permit explicit
// inversion, so Cramer's rule in factor2x2 is sound for blocks accepted here. There is no
// queue, no forced acceptance and no failing branch, which is why this returns void and why
// mDelaySize needs no epilogue: it keeps the zero setSymFactor gave it.

template<class Val>
void NumFactorEngine::factorDynamicRootSupernode(NumFactorDynamic<Val>& nf, std::int32_t jj,
                                                 std::vector<std::int32_t>& gblToLcl) const {
    // No pre/post distinction here: a root delays nothing, so the front never contracts.
    const std::size_t    jjFrontSize  = nf.mFrontSize[jj];
    const std::size_t    jjNumNodeIdx = jjFrontSize + nf.mUpdateSize[jj];
    std::int32_t*        jjNodeIdx    = nf.nodeIdx(jj);
    Val*                 jjVal        = nf.val(jj);
    const std::ptrdiff_t ld           = static_cast<std::ptrdiff_t>(jjNumNodeIdx);

    // Column-major position of (row r, column c), in ptrdiff_t to avoid overflow.
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    for (std::int32_t nextPivot = 0;
         nextPivot < static_cast<std::int32_t>(jjFrontSize); ) {
        std::int32_t              j      = nextPivot;
        const RootPivotColumnScan jScan  =
            scanRootPivotColumn(jjVal, ld, nextPivot, j,
                                static_cast<std::int32_t>(jjFrontSize));
        double                    jGamma = jScan.gamma;
        std::int32_t              r      = jScan.r;

        if (jGamma == 0) {                          // isolated column: 1x1, nothing to eliminate
            if (std::abs(jjVal[at(j, j)]) == 0) --nf.rank();   // 0.9: zero pivot drops rank
            nf.mPivotType[jjNodeIdx[j]] = 1;
            ++nextPivot;
            continue;
        }

        if (std::abs(jjVal[at(j, j)]) >= alphaRoot * jGamma) {   // 1x1 on this column's diagonal
            factor1x1(nf, jj, nextPivot, j, jjNodeIdx[j], jjFrontSize, jjNumNodeIdx,
                          gblToLcl);
            ++nextPivot;
            continue;
        }

        // The chase. rGamma >= jGamma always, because r holds j's largest off-diagonal and that
        // same entry lies in column r, so it is a lower bound on column r's maximum. Equality is
        // the stopping case; otherwise the maximum strictly increases and no column is visited
        // twice, which is why this always reaches a pivot and a root never delays.
        for (;;) {
            const RootPivotColumnScan rScan =
                scanRootPivotColumn(jjVal, ld, nextPivot, r,
                                    static_cast<std::int32_t>(jjFrontSize));
            const double rGamma = rScan.gamma;

            if (std::abs(jjVal[at(r, r)]) >= alphaRoot * rGamma) {   // 1x1 on the column reached
                factor1x1(nf, jj, nextPivot, r, jjNodeIdx[r], jjFrontSize,
                              jjNumNodeIdx, gblToLcl);
                ++nextPivot;
                break;
            }

            if (jGamma == rGamma) {                 // the maximum is mutual: accept the 2x2
                factor2x2(nf, jj, nextPivot, j, r, jjNodeIdx[j], jjNodeIdx[r],
                              jjFrontSize, jjNumNodeIdx, gblToLcl);
                nextPivot += 2;
                break;
            }

            j = r;  jGamma = rGamma;  r = rScan.r;
        }
    }

    // A root cannot delay: every path out of the chase accepts a pivot, so mDelaySize[jj] keeps the
    // zero setSymFactor gave it.
}

// **Threshold pivoting, AGL Figures 3.4 and 3.3.** The port of 0.9's jjUpdateSize != 0 pass,
// reshaped to the figure but choosing the same pivots, which the pinned counts confirmed.
//
// An ancestor exists, so nothing is forced: a column that cannot be pivoted safely here is delayed
// to the parent, where the rows it needs will have been assembled. The cost is fill, since a
// delayed column widens the parent's front, and that trade is why the threshold is tunable here
// (mPivotThreshold) where the root kernel's is a fixed constant.
//
// Figure 3.4 is the search. Candidates sit in a queue in front order; pop one, and either accept a
// pivot or push it to the rear and try the next, with no elimination in between, so a later
// candidate sees exactly the values an earlier one saw. A sweep in which every candidate has been
// tried since the last acceptance ends the loop and delays whatever is left. Figure 3.3 is the
// acceptance test, in acceptPivot2x2.
//
// The one thing sparsity forces on the search: a 2x2 partner must come from the front. An update
// row has no column here and no diagonal, so it can be measured but not eliminated, which is why
// scanNonRootPivotColumn stops tracking `q` once it passes jjFrontSize while `gamma` goes
// to the full height. The two quantities are different on purpose, and 7.4 of
// archive/sparse_factorization.md is where that distinction is argued.
template<class Val>
void NumFactorEngine::factorDynamicNonRootSupernode(NumFactorDynamic<Val>& nf, std::int32_t jj,
                                                    std::vector<std::int32_t>& gblToLcl) const {
    const std::size_t    jjPreFactorFrontSize = nf.mFrontSize[jj];
    const std::size_t    jjNumNodeIdx         = jjPreFactorFrontSize + nf.mUpdateSize[jj];
    std::int32_t*        jjNodeIdx            = nf.nodeIdx(jj);
    Val*                 jjVal                = nf.val(jj);
    const std::ptrdiff_t ld                   = static_cast<std::ptrdiff_t>(jjNumNodeIdx);
    const double         threshold            = mPivotThreshold;
    const bool           withHermitian        = hermitian(nf.factorization());

    // Column-major position of (row r, column c), in ptrdiff_t to avoid overflow.
    const auto at = [ld](std::int32_t r, std::int32_t c) {
        return static_cast<std::ptrdiff_t>(c) * ld + static_cast<std::ptrdiff_t>(r);
    };

    // The candidate pivot columns, by global index, in front order.
    std::list<std::int32_t> pivotList;
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(jjPreFactorFrontSize); ++j)
        pivotList.push_back(jjNodeIdx[j]);

    std::int32_t nextPivot = 0;


    // A non-root has update rows, so an ancestor exists to take a column this front cannot pivot.
    // Ported from 0.9 factorDynamicLDL_, the jjUpdateSize != 0 branch, and reshaped to AGL Figure
    // 3.4: pop a candidate, choose its partner in the same breath, then one chain of four outcomes,
    // 1x1 on an isolated column, 1x1 on the threshold test, 2x2 on the acceptance test, or delay.
    // The test itself is behind acceptPivot2x2, which is the hook the figure leaves open.
    while (!pivotList.empty()) {
        bool         pivotFound = false;
        std::int32_t trials     = static_cast<std::int32_t>(pivotList.size());

        while (trials > 0) {
            const std::int32_t lj = pivotList.front(); pivotList.pop_front();
            const std::int32_t j  = gblToLcl[lj];

            const NonRootPivotColumnScan jScan =
                scanNonRootPivotColumn(jjVal, ld, nextPivot, j,
                                       static_cast<std::int32_t>(jjPreFactorFrontSize),
                                       static_cast<std::int32_t>(jjNumNodeIdx));
            const double       jGamma      = jScan.gamma;
            const double       jFrontGamma = jScan.frontGamma;
            const std::int32_t q           = jScan.q;

            const Val jDiagonal = jjVal[at(j, j)];

            if (jGamma == 0) {                      // isolated column: 1x1, nothing to eliminate
                pivotFound = true;
                if (nextPivot != j) nf.swap(jj, nextPivot, j, gblToLcl);
                if (std::abs(jDiagonal) == 0) --nf.rank();   // 0.9: zero pivot drops rank
                nf.mPivotType[lj] = 1;
                ++nextPivot;
                break;
            }

            if (std::abs(jDiagonal) > 0 && std::abs(jDiagonal) >= threshold * jGamma) {   // 1x1
                pivotFound = true;
                factor1x1(nf, jj, nextPivot, j, lj, jjPreFactorFrontSize, jjNumNodeIdx,
                              gblToLcl);
                ++nextPivot;
                break;
            }

            // q == NIL means the front holds no partner, which happens only when lj is the last
            // unfactored column. That is the same condition as pivotList being empty after the pop,
            // asked of the scan rather than of the queue.
            if (q != NIL && acceptPivot2x2(jjVal, static_cast<std::int32_t>(jjNumNodeIdx),
                                           nextPivot, j, q, jGamma, jFrontGamma,
                                           withHermitian)) {                              // 2x2
                pivotFound = true;
                const std::int32_t lq = jjNodeIdx[q];
                pivotList.remove(lq);
                factor2x2(nf, jj, nextPivot, j, q, lj, lq, jjPreFactorFrontSize, jjNumNodeIdx,
                              gblToLcl);
                nextPivot += 2;
                break;
            }

            pivotList.push_back(lj);                    // nothing acceptable: delay lj
            --trials;
        }

        if (!pivotFound)
            break;
    }

    // Whatever is left could not be pivoted here; delay it to an ancestor. frontSize contracts by
    // that count, and the val height (frontSize + delaySize + updateSize) is preserved.
    const std::size_t delaySize = pivotList.size();
    nf.mDelaySize[jj]  = delaySize;
    nf.mFrontSize[jj] -= delaySize;
}

template<class Val>
void NumFactorEngine::updateDynamicUpdateBlock(const NumFactorDynamic<Val>& nf, std::int32_t jj,
                                             std::size_t jjKkUpdateSp, UpdateBlock<Val>& jjKkUpdateBlock) const {
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
    // in place by those pivots back inside factorDynamicNonRootSupernode, and it then migrates
    // to the parent
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
    const int jjKkHeight = static_cast<int>(jjKkUpdateBlock.height());
    const int jjKkWidth  = static_cast<int>(jjKkUpdateBlock.width());
    const int dld        = jjKkHeight;

    // As in the static twin, jjKkHeight and jjKkWidth are the (jj, kk) edge's dimensions, not either
    // supernode's: jjKkWidth is jj's rows landing in kk (change jj or kk and it moves), jjKkHeight all
    // of jj's rows from jjKkUpdateSp down. Only f is jj's alone. This kernel is the edge's work, and
    // formDynamicUpper forms only the f-by-jjKkWidth slice for this one (jj, kk) pair.

    const bool          withHermitian = hermitian(nf.factorization());
    const std::int32_t* jjNodeIdx     = nf.nodeIdx(jj);
    // jj's rows from `jjKkUpdateSp` down: the block this ancestor is about to receive.
    const Val*          l21Val        = jjVal + jjKkUpdateSp;
    Val*                u22Val        = jjKkUpdateBlock.val();

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
void NumFactorEngine::updateDynamicUpdateMatrix(const NumFactorDynamic<Val>& nf, std::int32_t kk,
                                                UpdateMatrix<Val>& kkUpdateMatrix) const {
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
    gemmLower(u, f, l21, ld, upper.data(), f, kkUpdateMatrix.val(), u);
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
//   Update.     For each jj owing kk: if kk is jj's parent, assemble jj's delayed columns and only
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
// frontSize (the solve loops j < frontSize, updateDynamicUpdateBlock returns on f == 0), so a
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
        Val*                kkVal                = nf.val(kk);   // for assembleFromA just below

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

            UpdateBlock<Val> jjKkUpdateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, jjKkUpdateBlock.rowIdx());

            updateDynamicUpdateBlock(nf, jj, jjKkUpdateSp, jjKkUpdateBlock);
            assembleUpdateBlock(jjKkUpdateBlock, nf, kk, gblToLcl);

            // jj has updated kk. Queue it for the next ancestor it must update.
            nextUpdateSp[jj] = jjKkUpdateSp + jjKkWidth;
            if (nextUpdateSp[jj] < jjNumNodeIdx)
                descendantUpdateQueue[nf.nodeToSnode(jjNodeIdx[nextUpdateSp[jj]])].push_back(jj);
        }

        // A root has no parent to delay a column to, and the two kernels differ throughout on
        // that one fact, so the choice is made here rather than inside.
        if (parent[kk] == NIL) factorDynamicRootSupernode(nf, kk, gblToLcl);
        else                   factorDynamicNonRootSupernode(nf, kk, gblToLcl);

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
// The second difference is smaller and follows from the direction. Left-looking assembles a child's
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

        // A root has no parent to delay a column to, and the two kernels differ throughout on
        // that one fact, so the choice is made here rather than inside.
        if (parent[jj] == NIL) factorDynamicRootSupernode(nf, jj, gblToLcl);
        else                   factorDynamicNonRootSupernode(nf, jj, gblToLcl);

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

            const std::size_t jjKkHeight = jjNumNodeIdx - jjKkUpdateSp;
            std::size_t       jjKkWidth  = 0;
            while (jjKkUpdateSp + jjKkWidth < jjNumNodeIdx
                   && nf.nodeToSnode(jjNodeIdx[jjKkUpdateSp + jjKkWidth]) == kk)
                ++jjKkWidth;

            UpdateBlock<Val> jjKkUpdateBlock(jjKkHeight, jjKkWidth);
            std::copy(jjNodeIdx + jjKkUpdateSp, jjNodeIdx + jjNumNodeIdx, jjKkUpdateBlock.rowIdx());

            updateDynamicUpdateBlock(nf, jj, jjKkUpdateSp, jjKkUpdateBlock);

            setGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
            assembleUpdateBlock(jjKkUpdateBlock, nf, kk, gblToLcl);
            clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);

            jjKkUpdateSp += jjKkWidth;
        }
    }

    return true;
}

// =================================================================================================
// Dynamic LDL, multifrontal. The left-looking dynamic skeleton (expand a front by its children's
// delayed columns, assemble A past them, factor, which may delay again) with the per-supernode
// update matrices in place of the pull queue: assemble each child by both halves of the assembly, then
// factor and leave this supernode's contribution block for its parent. This is where delayed columns
// meet the update matrices.
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
    // its parent has assembled it.
    std::vector<UpdateMatrix<Val>> updateMatrix(snodeSize);

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
        updateMatrix[kk].allocate(kkUpdateSize);
        {
            std::int32_t* kkUpdateRowIdx = updateMatrix[kk].rowIdx();
            for (std::size_t sp = 0; sp < kkUpdateSize; ++sp)
                kkUpdateRowIdx[sp] = kkNodeIdx[kkPreFactorFrontSize + sp];
        }

        // Assemble each child. Its delayed columns become kk front columns (assembleDelay) and their
        // storage is then reclaimed (contractVal); that pair is the child's delayed-column handling,
        // a no-op for a child that delayed nothing. Its contribution block is assembled into kk's
        // frontal (assembleUpdateMatrix) and is then freed. The assembly is independent of the
        // contraction, it reads jj's update matrix rather than jj's front, so the only order that
        // matters is assembleDelay before contractVal (read the columns before reclaiming them) and
        // the assembly before discard.
        for (std::int32_t jj = firstChild[kk]; jj != NIL; jj = nextSibling[jj]) {
            const std::size_t jjDelaySize = nf.delaySize(jj);
            if (jjDelaySize > 0) {
                assembleDelay(nf, jj, kk, gblToLcl);
                nf.contractVal(jj, jjDelaySize);
            }
            assembleUpdateMatrix(updateMatrix[jj], nf, kk, updateMatrix[kk], gblToLcl);
            updateMatrix[jj].discard();
        }

        // Factor kk's pivots (dynamic threshold pivoting; may delay columns to kk's parent).
        // A root has no parent to delay a column to, and the two kernels differ throughout on
        // that one fact, so the choice is made here rather than inside.
        if (parent[kk] == NIL) factorDynamicRootSupernode(nf, kk, gblToLcl);
        else                   factorDynamicNonRootSupernode(nf, kk, gblToLcl);

        // Form kk's contribution block from the factored pivots and leave it for kk's parent.
        updateDynamicUpdateMatrix(nf, kk, updateMatrix[kk]);

        clearGlobalToLocal(kkNumNodeIdx, kkNodeIdx, gblToLcl);
    }

    return true;
}

template<class Val>
bool NumFactorEngine::compute(const SparseMatrix<Val>& A, const Permutation& p, const SymFactor& sf,
                              NumFactorStatic<Val>& nf) const {
    if (A.size() != p.size() || A.size() != sf.size())
        return false;

    // **Dynamic pivoting cannot go into this storage, and never will.** A delayed column expands
    // its parent's front, and this factor's value buffer is one flat array sized once from the
    // symbolic factorization; expanding a front in the middle of it would mean moving everything
    // after it. So this is a design refusal rather than a missing feature, and the combination is
    // asserted to be refused in test_pipeline. Callers wanting a dynamic factorization pass
    // NumFactorDynamic, which the overload below takes.
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

    // Dynamic LDL is the reason this storage exists. Everything runs except complex `LDL^H`. The
    // static factorizations run unchanged, below.
    //
    // **The kernels needed nothing to become complex.** 0.9's complex `factorDynamicLDL_` differs
    // from its real one in six lines, all the same edit: the pivot *magnitudes* (`max1`, `max2`,
    // `maxmax`) are declared real rather than scalar. This port declared them `double` from the
    // start, so it was already the complex form; `updateDynamicLDL_` is byte-identical between
    // 0.9's two engines to begin with. Everything else is `Val` arithmetic and `std::abs`, which
    // means modulus for complex and is the right comparison either way.
    //
    // **Complex `LDL^H` is an extension rather than a port, and it runs.** 0.9's complex LDL is
    // symmetric only, so there was nothing to transcribe and nothing to check against: the oracles
    // are the residual in `test_pipeline`, on a genuinely Hermitian band through all three
    // traversals, and reconstruction of `L D L^H` on dense fronts in `test_numfactor`. What the
    // Hermitian form needs beyond the symmetric one is concentrated in two places,
    // `readPivotBlock2x2`, which conjugates the off-diagonal and forces both diagonals real, and
    // `forceReal` on the 1x1 pivot; everything downstream works from what those return.
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
        case Traversal::Multifrontal: return factorStaticMultifrontal(A, p, sf, nf);
    }
    return false;
}

template bool NumFactorEngine::compute(const SparseMatrix<double>&, const Permutation&,
                                       const SymFactor&, NumFactorDynamic<double>&) const;
template bool NumFactorEngine::compute(const SparseMatrix<std::complex<double>>&,
                                       const Permutation&, const SymFactor&,
                                       NumFactorDynamic<std::complex<double>>&) const;

} // namespace Oblio
