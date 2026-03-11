#include "oblio/Matrix.h"
#include <complex>
#include <stdexcept>

namespace Oblio {

template<class Val>
Matrix<Val> Matrix<Val>::fromCOO(Size n,
                                  const std::vector<Size>& rows,
                                  const std::vector<Size>& cols,
                                  const std::vector<Val>&  vals) {
    if (rows.size() != cols.size() || rows.size() != vals.size())
        throw std::invalid_argument("COO arrays must have equal length");
    Matrix m;
    m.mSize   = n;
    m.mIsSymm = true;
    // Phase 1: count entries per column (lower triangle: row >= col).
    m.mColPtr.assign(n + 1, 0);
    for (Size k = 0; k < rows.size(); ++k) {
        Size lo = (rows[k] >= cols[k]) ? cols[k] : rows[k];
        m.mColPtr[lo + 1]++;
    }
    for (Size j = 0; j < n; ++j) m.mColPtr[j + 1] += m.mColPtr[j];
    Size nnz = m.mColPtr[n];
    m.mRowIdx.resize(nnz);
    m.mVal.resize(nnz);
    std::vector<Size> cur(m.mColPtr.begin(), m.mColPtr.end());
    for (Size k = 0; k < rows.size(); ++k) {
        Size r = rows[k], c = cols[k];
        Size lo = (r >= c) ? c : r;
        Size hi = (r >= c) ? r : c;
        Size pos = cur[lo]++;
        m.mRowIdx[pos] = hi;
        m.mVal[pos]    = vals[k];
    }
    // Phase 2: sort each column by row index and merge duplicates.
    for (Size j = 0; j < n; ++j) {
        Size beg = m.mColPtr[j], end = m.mColPtr[j + 1];
        // Insertion sort (columns are small).
        for (Size a = beg + 1; a < end; ++a) {
            Size ri = m.mRowIdx[a]; Val vi = m.mVal[a];
            Size b = a;
            while (b > beg && m.mRowIdx[b - 1] > ri) {
                m.mRowIdx[b] = m.mRowIdx[b - 1];
                m.mVal[b]    = m.mVal[b - 1];
                --b;
            }
            m.mRowIdx[b] = ri; m.mVal[b] = vi;
        }
        // Merge duplicates in-place (last value wins for symmetric entries).
        Size w = beg;
        for (Size r = beg; r < end; ) {
            Size ri = m.mRowIdx[r]; Val vi = m.mVal[r]; ++r;
            while (r < end && m.mRowIdx[r] == ri) { vi = m.mVal[r]; ++r; }
            m.mRowIdx[w] = ri; m.mVal[w] = vi; ++w;
        }
        // Shift remaining columns if entries were removed.
        Size removed = end - w;
        if (removed > 0) {
            Size tail = m.mRowIdx.size() - end;
            for (Size t = 0; t < tail; ++t) {
                m.mRowIdx[w + t] = m.mRowIdx[end + t];
                m.mVal[w + t]    = m.mVal[end + t];
            }
            m.mRowIdx.resize(m.mRowIdx.size() - removed);
            m.mVal.resize(m.mVal.size() - removed);
            for (Size jj = j + 1; jj <= n; ++jj) m.mColPtr[jj] -= removed;
        }
    }
    return m;
}

template<class Val>
Size Matrix<Val>::numOffDiag() const {
    Size cnt = 0;
    for (Size j = 0; j < mSize; ++j)
        for (Size p = mColPtr[j]; p < mColPtr[j + 1]; ++p)
            if (mRowIdx[p] != j) ++cnt;
    return cnt;
}

// ── Explicit instantiations ──────────────────────────────────────────────────
template struct Matrix<double>;
template struct Matrix<std::complex<double>>;

} // namespace Oblio
