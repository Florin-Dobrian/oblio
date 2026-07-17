// examples/pipeline.cpp
// The solve pipeline by hand, without the OblioEngine facade (which does not exist yet). This is
// what that facade will one day encapsulate: order the matrix, build the elimination forest,
// factor it symbolically, factor it numerically, and solve. Here each phase is a separate object,
// wired in the order the data flows, so the seams are visible.
//
// It runs the same 4x4 tridiagonal matrix through every factorization and every traversal, under
// both ordering algorithms, and prints the residual ||Ax - b|| / ||b|| for each. Where a
// configuration is not implemented yet (dynamic LDL, multifrontal), the numeric factorization
// returns false and the row says so: the sweep doubles as a map of what is wired today.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/pipeline.cpp src/*.cpp \
//       -framework Accelerate -o pipeline
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <cstdio>
#include <vector>

using namespace Oblio;

namespace {

const char* name(OrderMethod m) {
    switch (m) {
        case OrderMethod::Natural: return "Natural";
        case OrderMethod::MMD:     return "MMD";
        case OrderMethod::AMD:     return "AMD";
    }
    return "?";
}
const char* name(Factorization f) {
    switch (f) {
        case Factorization::Cholesky:    return "Cholesky";
        case Factorization::StaticLDLT:  return "StaticLDLT";
        case Factorization::StaticLDLH:  return "StaticLDLH";
        case Factorization::DynamicLDLT: return "DynamicLDLT";
        case Factorization::DynamicLDLH: return "DynamicLDLH";
    }
    return "?";
}
const char* name(Traversal t) {
    switch (t) {
        case Traversal::LeftLooking:  return "LeftLooking";
        case Traversal::RightLooking: return "RightLooking";
        case Traversal::Multifrontal: return "Multifrontal";
    }
    return "?";
}

} // namespace

int main() {
    // The 4x4 tridiagonal, diagonal 4, off-diagonal -1, stored full (both triangles) in CSC:
    //
    //   [ 4 -1  0  0 ]
    //   [-1  4 -1  0 ]
    //   [ 0 -1  4 -1 ]
    //   [ 0  0 -1  4 ]
    //
    const std::size_t n = 4;
    const std::vector<std::size_t>  colPtr = {0, 2, 5, 8, 10};
    const std::vector<std::int32_t> rowIdx = {0, 1,  0, 1, 2,  1, 2, 3,  2, 3};
    const std::vector<double>       val    = {4, -1, -1, 4, -1, -1, 4, -1, -1, 4};
    const SparseMatrix<double> A(n, colPtr, rowIdx, val);

    // Right-hand side: all ones.
    Vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = 1.0;
    const double normB = b.norm();

    MultiplyEngine mul;
    SolveEngine    sol;

    printf("4x4 tridiagonal, b = ones. Residual ||Ax - b|| / ||b|| per configuration:\n\n");
    printf("  %-8s  %-13s  %-13s  %s\n", "order", "factorization", "traversal", "residual");
    printf("  %-8s  %-13s  %-13s  %s\n", "-----", "-------------", "---------", "--------");

    // The forest and symbolic factor depend only on the pattern and the ordering, not on the
    // factorization or traversal, so they are computed once per ordering. Only the numeric factor
    // varies inside the inner loops.
    for (OrderMethod method : {OrderMethod::MMD, OrderMethod::AMD}) {
        OrderEngine ord(method);
        Permutation p;
        if (!ord.compute(A, p)) { printf("  %-8s  ordering failed\n", name(method)); continue; }

        ElmForest f;
        ElmForestEngine fe;
        if (!fe.compute(A, p, f)) { printf("  %-8s  forest failed\n", name(method)); continue; }

        SymFactor s;
        SymFactorEngine se;
        if (!se.compute(A, p, f, s)) { printf("  %-8s  symbolic failed\n", name(method)); continue; }

        for (Factorization fact : {Factorization::Cholesky, Factorization::StaticLDLT,
                                   Factorization::StaticLDLH, Factorization::DynamicLDLT,
                                   Factorization::DynamicLDLH}) {
            for (Traversal trav : {Traversal::LeftLooking, Traversal::RightLooking,
                                   Traversal::Multifrontal}) {
                NumFactorStatic<double> nf;
                NumFactorEngine ne(fact, trav);
                if (!ne.compute(A, p, s, nf)) {
                    printf("  %-8s  %-13s  %-13s  %s\n",
                           name(method), name(fact), name(trav), "not implemented");
                    continue;
                }

                Vector<double> x(n);
                if (!sol.compute(p, nf, b, x)) {
                    printf("  %-8s  %-13s  %-13s  %s\n",
                           name(method), name(fact), name(trav), "solve failed");
                    continue;
                }

                Vector<double> r(n);
                if (!mul.residual(A, x, b, r)) {
                    printf("  %-8s  %-13s  %-13s  %s\n",
                           name(method), name(fact), name(trav), "residual failed");
                    continue;
                }

                const double rel = normB > 0 ? r.norm() / normB : r.norm();
                printf("  %-8s  %-13s  %-13s  %.3e\n",
                       name(method), name(fact), name(trav), rel);
            }
        }
    }

    // One actual solution, so the example ends on a number a reader can check by hand.
    {
        OrderEngine ord(OrderMethod::AMD);
        Permutation p;      ord.compute(A, p);
        ElmForest f;        ElmForestEngine fe; fe.compute(A, p, f);
        SymFactor s;        SymFactorEngine se; se.compute(A, p, f, s);
        NumFactorStatic<double> nf;
        NumFactorEngine ne(Factorization::Cholesky, Traversal::LeftLooking);
        ne.compute(A, p, s, nf);
        Vector<double> x(n);
        sol.compute(p, nf, b, x);

        printf("\nSolution (AMD, Cholesky, LeftLooking):\n");
        for (std::size_t i = 0; i < n; ++i)
            printf("  x[%zu] = %.10f\n", i, x[i]);
    }

    return 0;
}
