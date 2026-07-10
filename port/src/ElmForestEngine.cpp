#include "oblio/ElmForestEngine.h"

namespace Oblio {

template<class Val>
bool ElmForestEngine::compute(const SparseMatrix<Val>& A, const Permutation& p,
                              ElmForest& f) const {
    const std::size_t n = A.numCols();
    if (p.size() != n)
        return false;

    // A is full-symmetric: each column already holds its complete neighbour list,
    // so the etree reads it directly (diagonal entries map to lc1 == lc2 and are
    // skipped by the lc1 < lc2 test).
    computeParent(n, A.colPtr(), A.rowIdx(), p.oldToNew(), p.newToOld(), f.mParent);
    return true;
}

void ElmForestEngine::computeParent(std::size_t n,
        const std::vector<std::size_t>&  colPtr,
        const std::vector<std::int32_t>& rowIdx,
        const std::vector<std::int32_t>& oldToNew,
        const std::vector<std::int32_t>& newToOld,
        std::vector<std::int32_t>& parent) const {
    parent.assign(n, NIL);
    std::vector<std::int32_t> ancestor(n, NIL);

    // For each factor column lc2 (increasing), look at its neighbours mapping to
    // earlier columns lc1 < lc2; path-compress to attach lc1's subtree under lc2.
    for (std::size_t lc2 = 0; lc2 < n; ++lc2) {
        const std::size_t ac2 = static_cast<std::size_t>(newToOld[lc2]);
        for (std::size_t sp = colPtr[ac2]; sp < colPtr[ac2 + 1]; ++sp) {
            const std::size_t lc1 = static_cast<std::size_t>(oldToNew[static_cast<std::size_t>(rowIdx[sp])]);
            if (lc1 >= lc2)
                continue;   // later column or the diagonal itself
            std::size_t lc3 = lc1;
            while (ancestor[lc3] != NIL && static_cast<std::size_t>(ancestor[lc3]) != lc2) {
                const std::size_t lc4 = static_cast<std::size_t>(ancestor[lc3]);
                ancestor[lc3] = static_cast<std::int32_t>(lc2);
                lc3 = lc4;
            }
            if (ancestor[lc3] == NIL) {
                ancestor[lc3] = static_cast<std::int32_t>(lc2);
                parent[lc3] = static_cast<std::int32_t>(lc2);
            }
        }
    }
}

template bool ElmForestEngine::compute(const SparseMatrix<double>&, const Permutation&, ElmForest&) const;
template bool ElmForestEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, ElmForest&) const;

} // namespace Oblio
