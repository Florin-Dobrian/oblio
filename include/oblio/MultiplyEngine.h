#pragma once

// MultiplyEngine.h - y = A x, for a sparse A.
//
// Small, and it exists for one reason: **the residual.** Everything before this point is checked
// against an oracle of some kind (a dense Cholesky, a reconstruction), which tells us whether one
// phase computed what it should. None of it tells us whether the *pipeline* does. With A x and a
// solve we can finally ask the only question that matters end to end:
//
//     || A x - b ||  /  || b ||
//
// through ordering, forest, symbolic factorization, numeric factorization and solve, in one
// number. Nothing else exercises all six at once.
//
// A is stored full-symmetric (both triangles), so this is an ordinary sparse matvec: for each
// column, scatter its entries into y. No symmetry trick, no half-storage bookkeeping. That is one
// of the things full storage buys, and the flat-vs-VV experiment measured this exact kernel.

#include "oblio/SparseMatrix.h"
#include "oblio/Vector.h"

#include <complex>

namespace Oblio {

class MultiplyEngine {
public:
    MultiplyEngine() = default;

    // y := A x. Returns false if the sizes disagree.
    template<class Val>
    bool compute(const SparseMatrix<Val>& A, const Vector<Val>& x, Vector<Val>& y) const;

    // y := A x - b, the residual vector. The one call a test actually wants.
    template<class Val>
    bool residual(const SparseMatrix<Val>& A, const Vector<Val>& x, const Vector<Val>& b,
                  Vector<Val>& r) const;
};

} // namespace Oblio
