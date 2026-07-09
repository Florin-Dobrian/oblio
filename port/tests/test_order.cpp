// test_order.cpp — tests for OrderEngine (AMD and MMD) on SparseMatrix,
// producing a Permutation. Focus: every ordering must be a valid permutation
// (bijection) of the right size, on a range of structures.
//
// A "valid ordering" here means: order() succeeds, p.size() == n, and
// p.validate() (genuine bijection, consistent inverse maps). We don't assert a
// *specific* permutation — AMD and MMD may legitimately differ — only that each
// produces a correct permutation of the matrix.

#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"

#include <iostream>
#include <vector>

using namespace Oblio;

namespace {

int gPass = 0;
int gFail = 0;

void ck(bool ok, const std::string& name) {
    std::cout << "  " << (ok ? "PASS  " : "FAIL  ") << name << "\n";
    ok ? ++gPass : ++gFail;
}

// Run one method on one matrix; check the result is a valid permutation.
template<class Val>
void checkOrder(const SparseMatrix<Val>& A, OrderMethod method,
                const std::string& label) {
    OrderEngine eng(method);
    Permutation p;
    const bool ok = eng.order(A, p);
    ck(ok && p.size() == A.numCols() && p.validate(), label);
}

// Build a symmetric tridiagonal matrix (lower triangle, diagonal included) of
// size n: each column j has the diagonal j and, if j+1 < n, the sub-diagonal j+1.
SparseMatrix<double> tridiag(std::size_t n) {
    std::vector<std::size_t> colPtr(n + 1, 0);
    std::vector<std::size_t> rowIdx;
    std::vector<double>      val;
    for (std::size_t j = 0; j < n; ++j) {
        rowIdx.push_back(j); val.push_back(2.0);          // diagonal
        if (j + 1 < n) { rowIdx.push_back(j + 1); val.push_back(-1.0); } // sub-diag
        colPtr[j + 1] = rowIdx.size();
    }
    return SparseMatrix<double>(n, colPtr, rowIdx, val);
}

} // namespace

int main() {
    std::cout << "=== OrderEngine tests (AMD / MMD) ===\n";

    // 1) Small hand-built "arrow" matrix (6x6). Arrow matrices are the classic
    //    case where ordering matters: natural order fills in, a good order doesn't.
    //    col0: 0,1,2,5  col1: 1  col2: 2,3  col3: 3,4  col4: 4,5  col5: 5
    {
        std::vector<std::size_t> colPtr = {0,4,5,7,9,11,12};
        std::vector<std::size_t> rowIdx = {0,1,2,5, 1, 2,3, 3,4, 4,5, 5};
        std::vector<double>      val(rowIdx.size(), 1.0);
        SparseMatrix<double> A(6, colPtr, rowIdx, val);
        checkOrder(A, OrderMethod::AMD, "arrow 6x6           : AMD valid");
        checkOrder(A, OrderMethod::MMD, "arrow 6x6           : MMD valid");
    }

    // 2) Tridiagonal matrices of a few sizes.
    for (std::size_t n : {1u, 2u, 10u, 100u}) {
        SparseMatrix<double> A = tridiag(n);
        checkOrder(A, OrderMethod::AMD, "tridiag n=" + std::to_string(n) +
                   std::string(n < 10 ? "          " : (n < 100 ? "         " : "        ")) + ": AMD valid");
        checkOrder(A, OrderMethod::MMD, "tridiag n=" + std::to_string(n) +
                   std::string(n < 10 ? "          " : (n < 100 ? "         " : "        ")) + ": MMD valid");
    }

    // 3) Diagonal-only matrix (no off-diagonals): ordering is trivial but must
    //    still yield a valid permutation.
    {
        std::size_t n = 5;
        std::vector<std::size_t> colPtr(n + 1);
        std::vector<std::size_t> rowIdx(n);
        std::vector<double>      val(n, 1.0);
        for (std::size_t j = 0; j < n; ++j) { colPtr[j] = j; rowIdx[j] = j; }
        colPtr[n] = n;
        SparseMatrix<double> A(n, colPtr, rowIdx, val);
        checkOrder(A, OrderMethod::AMD, "diagonal 5x5        : AMD valid");
        checkOrder(A, OrderMethod::MMD, "diagonal 5x5        : MMD valid");
    }

    // 4) Complex matrix — ordering uses structure only, so it must work identically.
    {
        std::vector<std::size_t> colPtr = {0,4,5,7,9,11,12};
        std::vector<std::size_t> rowIdx = {0,1,2,5, 1, 2,3, 3,4, 4,5, 5};
        std::vector<std::complex<double>> val(rowIdx.size(), {1.0, 0.0});
        SparseMatrix<std::complex<double>> C(6, colPtr, rowIdx, val);
        checkOrder(C, OrderMethod::AMD, "arrow 6x6 complex   : AMD valid");
        checkOrder(C, OrderMethod::MMD, "arrow 6x6 complex   : MMD valid");
    }

    // 5) AMD and MMD on the same matrix each give a valid (if possibly different)
    //    permutation — sanity that both drive the same matrix without interfering.
    {
        SparseMatrix<double> A = tridiag(20);
        OrderEngine amd(OrderMethod::AMD), mmd(OrderMethod::MMD);
        Permutation pa, pm;
        bool ok = amd.order(A, pa) && mmd.order(A, pm);
        ck(ok && pa.validate() && pm.validate() && pa.size() == 20 && pm.size() == 20,
           "tridiag n=20        : AMD & MMD both valid");
    }

    std::cout << "\nOrderEngine tests: " << gPass << "/" << (gPass + gFail)
              << " passed\n";
    return gFail == 0 ? 0 : 1;
}
