#pragma once

// SymFactorEngine.h - computes the symbolic factorization of a SparseMatrix under a
// given Permutation, from an already-computed ElmForest. Copies the forest
// attributes the factor needs, then computes each supernode's index set as a
// union: the sparsity patterns of the supernode's own front columns, plus the
// update indices of its children. Returns true on success.
//
// Corollary 2 of the theory is a recurrence, in the mathematical sense: Idx(kk) is defined in
// terms of Idx(jj) for its children jj. But the implementation is a plain loop, not a
// recursion, and that is not a choice. A parent's index exceeds its children's, so increasing
// supernode order is already a valid bottom-up order of the forest: children are complete by
// the time their parent is reached, with no stack and no tree traversal.
//
// Reads the matrix, permutation and forest via their public accessors and fills
// the symbolic factorization via friend access.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/ElmForest.h"
#include "oblio/SymFactor.h"

#include <vector>
#include <complex>
#include <cstddef>
#include <cstdint>

namespace Oblio {

class SymFactorEngine {
public:
    SymFactorEngine() = default;

    // Compute the symbolic factorization of A under p, from the forest f.
    //
    // The factor's pattern depends on A's pattern alone, never on its values, so the
    // implementation is the non-templated overload, taking colPtr and rowIdx, compiled
    // once. The templated overload is an adapter over it. The pattern overload is public:
    // a caller holding a graph with no numbers attached can compute the symbolic
    // factorization without inventing a scalar type to satisfy the signature.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Permutation& p,
                 const ElmForest& f, SymFactor& s) const;

    bool compute(const std::vector<std::size_t>&  colPtr,
                 const std::vector<std::int32_t>& rowIdx,
                 const Permutation& p, const ElmForest& f, SymFactor& s) const;

private:
    // The front indices of each supernode, laid out contiguously. The map runs column to
    // supernode, so a supernode's columns are generally scattered through it (and after
    // amalgamation they need not even be consecutive), while the union loop wants them
    // together. Counting sort by supernode: walking columns in increasing order leaves each
    // supernode's front indices sorted.
    //
    // A pass over the map is needed in *both* regimes, and cannot be avoided in either: the map
    // runs column to supernode and we want the reverse, and a supernode's columns need not be
    // contiguous (the forest is topological, not postordered), so its lowest column cannot be
    // found by inspection. That pass is the price of admission.
    //
    // What amalgamation changes is how much the pass must produce. With fundamental supernodes
    // the columns share a pattern exactly, so Struct(k'') is contained in Struct(k') for any two
    // front columns k' < k'' of a supernode, and reading the first (lowest) front column's
    // A-pattern is enough: everything the later ones would contribute is already there. (Nor do
    // the front indices need seeding. The first column's Struct already holds it and every later
    // column of the supernode, since those are rows below it, so the block diagonal emerges from
    // the union rather than being put into it.) Amalgamation merges columns whose patterns are
    // only *nearly* identical, paying explicit zeros for the difference, so the containment
    // fails and every front column must genuinely be read.
    //
    // So the alternative (getFirstFrontalIndex, below) differs in *size*, not in passes: it
    // produces the lowest column of each supernode, supSize entries, where this produces all
    // front columns, size entries. Both references take this general path unconditionally, and
    // 0.9's comment says exactly why one must when amalgamation is on. See Section 4.6 of the
    // sparse-factorization notes.
    void getFrontalIndices(const SymFactor& s,
                              std::vector<std::size_t>&  frontSupPtr,
                              std::vector<std::int32_t>& frontRowIdx) const;

    // The specialization: the lowest front column of each supernode, and nothing else. One
    // vector of supSize entries, no offsets, since each supernode has exactly one and the
    // position is the supernode.
    //
    // Correct exactly when the forest reports exactPatterns(): a supernode's columns then share
    // one pattern *in L*, so L(k'') is contained in L(k') for any two front columns k' < k''.
    // That covers nodal forests, fundamental supernodes, and amalgamation at threshold zero
    // (which merges only where the merge is free, hence only where the patterns already agree).
    // It is wrong the moment amalgamation stores a zero: the later front columns then carry rows
    // the first does not, and skipping them would silently lose indices.
    //
    // Be careful about *what* is sufficient, since the obvious gloss is false. It is **not**
    // that one column of A suffices: A's columns of a supernode are not nested, and a later one
    // can hold indices the first does not. What is sufficient is the *L* column. The union
    // computes Idx(K) = L(k') from two inputs, A's column k' and the children, and anything a
    // later front column of A could contribute is already in L(k'), so it arrives through the
    // children instead. Nothing is lost; the reason is containment in L, not in A. See Section
    // 4.6 of the notes, which works this through on a three-column example.
    //
    // What it saves is a great deal of A. getFrontalIndices makes the union read every column of
    // A, all `size` of them; this one makes it read `supSize`, one per supernode. On a forest
    // that compresses well that is most of the traversal of A.
    void getFirstFrontalIndex(const SymFactor& s,
                                   std::vector<std::int32_t>& firstFrontRowIdx) const;

    void sortIndices(SymFactor& s) const;
};

extern template bool SymFactorEngine::compute(const SparseMatrix<double>&, const Permutation&, const ElmForest&, SymFactor&) const;
extern template bool SymFactorEngine::compute(const SparseMatrix<std::complex<double>>&, const Permutation&, const ElmForest&, SymFactor&) const;

} // namespace Oblio
