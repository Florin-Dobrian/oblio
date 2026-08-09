// examples/example_pipeline_complex.cpp
// The complex counterpart of examples/example_pipeline_real.cpp: the same pipeline wired by hand,
// the same sweep over every ordering, factorization and traversal, with Val = complex<double>.
// Read the real one first; only what complex changes is explained here.
//
// **One matrix will not do, and that is the whole difference from the real sweep.** Over the reals
// symmetric and Hermitian are the same condition, so one matrix serves all five factorizations.
// Over the complex field they are different conditions and the five split into two groups:
//
//   Cholesky, StaticLDLH, DynamicLDLH    need a HERMITIAN matrix,        A = A^H, diagonal real
//   StaticLDLT, DynamicLDLT              need a COMPLEX SYMMETRIC one,   A = A^T, diagonal complex
//
// So this example builds both and hands each factorization the one it is entitled to. The
// predicate that decides is hermitian() in Types.h, the library's own, the same one the numeric
// kernels and the solve consult to decide whether they conjugate. The real example uses
// dynamicPivoting() the same way to pick the storage; here two predicates are in play, one
// choosing the container and one choosing the input.
//
// **Handing a factorization the other group's matrix is not caught, and would not look wrong.**
// Nothing validates the input's symmetry today (docs/TODO.md, "Validate the input matrix"), and
// zpotrf reads only the lower triangle and assumes the upper is its conjugate. So complex Cholesky
// on a complex symmetric matrix runs, succeeds, and returns a plausible wrong answer: it has
// factored the Hermitian matrix agreeing with the lower triangle, which is a different matrix. The
// pairing below is therefore a precondition the caller owes the library, not a check the library
// performs, which is exactly why it is written out rather than left implicit.
//
// Complex Cholesky on a complex symmetric matrix is not merely unchecked, it does not exist:
// positive definiteness requires x* A x to be real, which requires Hermitian. See the comment
// above Factorization in Types.h.
//
// Compile (macOS, from repo root):
//   g++ -std=c++17 -O3 -DOBLIO_BLAS_UNDERSCORE -Iinclude examples/example_pipeline_complex.cpp src/*.cpp -framework Accelerate -o example_pipeline_complex_cpp
// Linux: replace `-framework Accelerate` with `-llapack -lblas`.

#include "oblio/ElmForestEngine.h"
#include "oblio/MultiplyEngine.h"
#include "oblio/NumFactorEngine.h"
#include "oblio/OrderEngine.h"
#include "oblio/SolveEngine.h"
#include "oblio/SymFactorEngine.h"

#include <complex>
#include <cstdio>
#include <vector>

using namespace Oblio;

namespace {

// The one line that differs from example_pipeline_real.cpp, which declares `using Val = double;`
// here and is otherwise written in Val throughout. The two quantities that stay `double` below are
// norms, which are real whatever Val is.
using Val = std::complex<double>;

const char* name(Ordering m) {
    switch (m) {
        case Ordering::Natural: return "Natural";
        case Ordering::MMD:     return "MMD";
        case Ordering::MMD1:    return "MMD1";
        case Ordering::MMD2:    return "MMD2";
        case Ordering::MMD3:    return "MMD3";
        case Ordering::AMD:     return "AMD";
        case Ordering::AMD1:    return "AMD1";
        case Ordering::AMD2:    return "AMD2";
        case Ordering::AMD3:    return "AMD3";
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
    // Both matrices are 4x4 tridiagonal and share one sparsity pattern, stored full (both
    // triangles) in CSC. Sharing the pattern is what lets the ordering, the elimination forest and
    // the symbolic factorization be computed once per ordering and serve both: those three phases
    // read the pattern and never a value.
    const std::size_t n = 4;
    const std::vector<std::size_t>  colPtr = {0, 2, 5, 8, 10};
    const std::vector<std::int32_t> rowIdx = {0, 1,  0, 1, 2,  1, 2, 3,  2, 3};

    // Hermitian positive definite: A = A^H, so the diagonal is real and the entry above the
    // diagonal is the conjugate of the one below. Diagonally dominant (4 against 2*sqrt(2)), which
    // is what makes it positive definite and so a valid input for Cholesky.
    //
    //   [   4     -1+i     0       0   ]
    //   [ -1-i      4    -1+i      0   ]
    //   [   0     -1-i     4     -1+i  ]
    //   [   0       0    -1-i      4   ]
    //
    const Val d(4.0, 0.0), lower(-1.0, -1.0), upper(-1.0, 1.0);
    const std::vector<Val> hermitianVal = {d, lower,  upper, d, lower,  upper, d, lower,  upper, d};
    const SparseMatrix<Val> AH(n, colPtr, rowIdx, hermitianVal);

    // Complex symmetric: A = A^T, the same entry on both sides of the diagonal and the diagonal
    // itself complex. NOT Hermitian, and not positive definite in any sense, which is why Cholesky
    // is never handed this one.
    //
    //   [ 4+i    -1+i     0       0   ]
    //   [-1+i    4+i    -1+i      0   ]
    //   [   0    -1+i    4+i    -1+i  ]
    //   [   0      0     -1+i    4+i  ]
    //
    const Val ds(4.0, 1.0), off(-1.0, 1.0);
    const std::vector<Val> symmetricVal = {ds, off,  off, ds, off,  off, ds, off,  off, ds};
    const SparseMatrix<Val> AS(n, colPtr, rowIdx, symmetricVal);

    // Right-hand side: 1 + i throughout, so that a lost conjugation shows up rather than cancelling
    // against a real vector.
    Vector<Val> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = Val(1.0, 1.0);
    const double normB = b.norm();

    MultiplyEngine mulEng;
    SolveEngine    solEng;

    printf("4x4 tridiagonal, complex, b = 1+i. Residual ||Ax - b|| / ||b|| per configuration.\n");
    printf("Input: H = Hermitian positive definite, S = complex symmetric.\n\n");
    printf("  %-8s  %-13s  %-13s  %-5s  %s\n", "order", "factorization", "traversal", "input",
           "residual");
    printf("  %-8s  %-13s  %-13s  %-5s  %s\n", "-----", "-------------", "---------", "-----",
           "--------");

    for (Ordering ordering : {Ordering::Natural, Ordering::MMD, Ordering::MMD1,
                              Ordering::MMD2, Ordering::MMD3,
                              Ordering::AMD, Ordering::AMD1,
                              Ordering::AMD2, Ordering::AMD3, Ordering::AMD1B, Ordering::AMD2B}) {
        // Structural only, so one permutation, forest and symbolic factor serve both matrices.
        OrderEngine ordEng(ordering);
        Permutation P;
        if (!ordEng.compute(AH, P))
            { printf("  %-8s  order failed\n", name(ordering)); continue; }

        ElmForest ef;
        ElmForestEngine efEng;
        if (!efEng.compute(AH, P, ef))
            { printf("  %-8s  elimination forest failed\n", name(ordering)); continue; }

        SymFactor sf;
        SymFactorEngine sfEng;
        if (!sfEng.compute(AH, P, ef, sf))
            { printf("  %-8s  symbolic factor failed\n", name(ordering)); continue; }

        for (Factorization factorization : {Factorization::Cholesky, Factorization::StaticLDLT,
                                            Factorization::StaticLDLH, Factorization::DynamicLDLT,
                                            Factorization::DynamicLDLH}) {
            // The matrix this factorization is entitled to. One predicate, asked once.
            const bool                 withHermitian = hermitian(factorization);
            const SparseMatrix<Val>&   A             = withHermitian ? AH : AS;
            const char*                which         = withHermitian ? "H" : "S";

            for (Traversal traversal : {Traversal::LeftLooking, Traversal::RightLooking,
                                        Traversal::Multifrontal}) {
                // As in the real example, the two storages differ only in the type of nf, so the
                // body is a lambda templated on it rather than written twice.
                const auto attempt = [&](auto& nf) -> const char* {
                    NumFactorEngine nfEng(factorization, traversal);
                    if (!nfEng.compute(A, P, sf, nf))  return "numeric factor failed";

                    Vector<Val> x(n);
                    if (!solEng.compute(P, nf, b, x))  return "solve failed";

                    Vector<Val> r(n);
                    if (!mulEng.residual(A, x, b, r))  return "multiply failed";

                    const double relRes = normB > 0 ? r.norm() / normB : r.norm();
                    printf("  %-8s  %-13s  %-13s  %-5s  %.3e\n",
                           name(ordering), name(factorization), name(traversal), which, relRes);
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
                    printf("  %-8s  %-13s  %-13s  %-5s  %s\n",
                           name(ordering), name(factorization), name(traversal), which, why);
            }
        }
    }

    // One actual solution, so the example ends on numbers a reader can check. Cholesky needs the
    // Hermitian matrix, so that is the one solved here.
    {
        OrderEngine ordEng(Ordering::MMD2);
        Permutation P;  ordEng.compute(AH, P);
        ElmForest ef;   ElmForestEngine efEng;  efEng.compute(AH, P, ef);
        SymFactor sf;   SymFactorEngine sfEng;  sfEng.compute(AH, P, ef, sf);
        NumFactorStatic<Val> nf;
        NumFactorEngine nfEng(Factorization::Cholesky, Traversal::LeftLooking);
        nfEng.compute(AH, P, sf, nf);
        Vector<Val> x(n);
        solEng.compute(P, nf, b, x);

        printf("\nSolution (Hermitian input, MMD2, Cholesky, LeftLooking):\n");
        for (std::size_t i = 0; i < n; ++i)
            printf("  x[%zu] = %+.10f %+.10fi\n", i, x[i].real(), x[i].imag());
    }

    return 0;
}
