// examples/example_pipeline_real.cpp
// The solve pipeline by hand, without the DirectSolver facade. This is what that facade
// encapsulates (see examples/example_basic.cpp): order the matrix, build the elimination forest,
// factor it symbolically, factor it numerically, and solve. Here each phase is a separate object,
// wired in the order the data flows, so the seams are visible.
//
// It runs the same 4x4 tridiagonal matrix through every factorization and every traversal, under
// every ordering method, and prints the residual ||Ax - b|| / ||b|| for each. Every phase can
// refuse rather than answer, so each returns a bool and the row names the ENGINE that said no,
// "elimination forest failed" for efEng and so on. The three engines named for the object they
// fill say that object; the three named for what they do say the verb, "order", "solve",
// "multiply". Today every cell resolves to a residual.
//
// **The storage is chosen from the factorization, not fixed.** Dynamic pivoting delays a column up
// to an ancestor, which grows that ancestor's front, and only NumFactorDynamic can grow; asking for
// dynamic LDL in the flat storage is refused by design rather than by omission. The predicate that
// says so is dynamicPivoting() in Types.h, and the sweep consults it, so this example never asks
// for that combination, and a "numeric factor failed" row here would mean something else had gone
// wrong.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_pipeline_real.cpp src/*.cpp -framework Accelerate -o example_pipeline_real_cpp
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

// Spelled as an alias rather than used directly, so that this file and its complex twin differ in
// one line rather than throughout. Every scalar below is Val; the two quantities that stay `double`
// are norms, which are real whatever Val is.
using Val = double;

const char* name(Ordering m) {
    switch (m) {
        case Ordering::Natural: return "Natural";
        case Ordering::MMD:     return "MMD";
        case Ordering::MMD1:    return "MMD1";
        case Ordering::MMD2:    return "MMD2";
        case Ordering::AMD:     return "AMD";
        case Ordering::AMD1:    return "AMD1";
        case Ordering::AMD2:    return "AMD2";
        case Ordering::AMD1B:   return "AMD1B";
        case Ordering::AMD2B:   return "AMD2B";
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
    const std::vector<Val>          val    = {4, -1, -1, 4, -1, -1, 4, -1, -1, 4};
    const SparseMatrix<Val> A(n, colPtr, rowIdx, val);

    // Right-hand side: all ones.
    Vector<Val> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = 1.0;
    const double normB = b.norm();

    MultiplyEngine mulEng;
    SolveEngine    solEng;

    printf("4x4 tridiagonal, b = ones. Residual ||Ax - b|| / ||b|| per configuration:\n\n");
    printf("  %-8s  %-13s  %-13s  %s\n", "order", "factorization", "traversal", "residual");
    printf("  %-8s  %-13s  %-13s  %s\n", "-----", "-------------", "---------", "--------");

    // The forest and symbolic factor depend only on the pattern and the ordering, not on the
    // factorization or traversal, so they are computed once per ordering. Only the numeric factor
    // varies inside the inner loops.
    for (Ordering ordering : {Ordering::Natural, Ordering::MMD, Ordering::MMD1,
                              Ordering::MMD2, Ordering::AMD, Ordering::AMD1,
                              Ordering::AMD2, Ordering::AMD1B, Ordering::AMD2B}) {
        OrderEngine ordEng(ordering);
        Permutation P;
        if (!ordEng.compute(A, P))
            { printf("  %-8s  order failed\n", name(ordering)); continue; }

        ElmForest ef;
        ElmForestEngine efEng;
        if (!efEng.compute(A, P, ef))
            { printf("  %-8s  elimination forest failed\n", name(ordering)); continue; }

        SymFactor sf;
        SymFactorEngine sfEng;
        if (!sfEng.compute(A, P, ef, sf))
            { printf("  %-8s  symbolic factor failed\n", name(ordering)); continue; }

        for (Factorization factorization : {Factorization::Cholesky, Factorization::StaticLDLT,
                                            Factorization::StaticLDLH, Factorization::DynamicLDLT,
                                            Factorization::DynamicLDLH}) {
            for (Traversal traversal : {Traversal::LeftLooking, Traversal::RightLooking,
                                        Traversal::Multifrontal}) {
                // Factor, solve and take the residual, in whichever storage this factorization
                // needs. The two branches differ only in the type of nf, so the body is a lambda
                // templated on it rather than written twice.
                const auto attempt = [&](auto& nf) -> const char* {
                    NumFactorEngine nfEng(factorization, traversal);
                    if (!nfEng.compute(A, P, sf, nf))  return "numeric factor failed";

                    Vector<Val> x(n);
                    if (!solEng.compute(P, nf, b, x))  return "solve failed";

                    Vector<Val> r(n);
                    if (!mulEng.residual(A, x, b, r))  return "multiply failed";

                    const double relRes = normB > 0 ? r.norm() / normB : r.norm();
                    printf("  %-8s  %-13s  %-13s  %.3e\n",
                           name(ordering), name(factorization), name(traversal), relRes);
                    return nullptr;
                };

                const char* why = nullptr;
                if (dynamicPivoting(factorization)) {
                    NumFactorDynamic<Val> nf;
                    why = attempt(nf);
                } else {
                    NumFactorStatic<Val> nf;
                    why = attempt(nf);
                }

                if (why)
                    printf("  %-8s  %-13s  %-13s  %s\n",
                           name(ordering), name(factorization), name(traversal), why);
            }
        }
    }

    // One actual solution, so the example ends on a number a reader can check by hand.
    {
        OrderEngine ordEng(Ordering::MMD2);
        Permutation P;  ordEng.compute(A, P);
        ElmForest ef;   ElmForestEngine efEng;  efEng.compute(A, P, ef);
        SymFactor sf;   SymFactorEngine sfEng;  sfEng.compute(A, P, ef, sf);
        NumFactorStatic<Val> nf;
        NumFactorEngine nfEng(Factorization::Cholesky, Traversal::LeftLooking);
        nfEng.compute(A, P, sf, nf);
        Vector<Val> x(n);
        solEng.compute(P, nf, b, x);

        printf("\nSolution (MMD2, Cholesky, LeftLooking):\n");
        for (std::size_t i = 0; i < n; ++i)
            printf("  x[%zu] = %.10f\n", i, x[i]);
    }

    return 0;
}
