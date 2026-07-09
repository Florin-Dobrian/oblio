// Smoke test for the parallel port tree: SparseMatrix (CSC) + Permutation.
#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include <iostream>

using Oblio::SparseMatrix;
using Oblio::Permutation;

int main() {
    int pass = 0, fail = 0;
    auto ck = [&](bool ok, const char* n){ std::cout << "  " << n << (ok?"  PASS":"  FAIL") << "\n"; ok?++pass:++fail; };

    // --- SparseMatrix: 4x4 symmetric, lower triangle (diagonal included) ---
    //   col0: rows 0,1,2   col1: rows 1,3   col2: row 2   col3: row 3
    SparseMatrix<double> A(4, {0,3,5,6,7}, {0,1,2, 1,3, 2, 3}, {4,1,1, 4,1, 4, 4});
    ck(A.numCols()==4, "matrix numCols == 4    ");
    ck(A.nnz()==7,     "matrix nnz == 7        ");

    // --- Permutation: identity + inverse round-trip ---
    Permutation p(A.numCols());
    ck(p.size()==4 && p.validate(), "identity perm valid    ");
    bool rt = true;
    for (std::size_t i = 0; i < p.size(); ++i) rt = rt && p.newToOld()[p.oldToNew()[i]] == i;
    ck(rt, "perm inverse round-trip");

    std::cout << "port smoke: " << pass << "/" << (pass+fail) << " passed\n";
    return fail==0 ? 0 : 1;
}
