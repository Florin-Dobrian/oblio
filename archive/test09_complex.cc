// test09_complex.cc — test oblio 0.9 complex factorization
// Uses the same 6×6 complex symmetric tridiagonal as the modern test.
#include <iostream>
#include <cmath>
using namespace std;

#include "Oblio.h"

#undef __CLASS__
#define __CLASS__ "test09c"

static const char* algName(FactorizationAlgorithm a) {
    switch(a) {
        case LeftLooking:  return "LL ";
        case RightLooking: return "RL ";
        case Multifrontal: return "MF ";
    } return "?";
}
static const char* typeName(FactorizationType t) {
    switch(t) {
        case StaticCC:  return "CC  ";
        case StaticLDL: return "SLDL";
        case DynamicLDL:return "DLDL";
    } return "?";
}

// Build n×n complex symmetric tridiagonal (full columns):
// diag = diagVal (real), off-diag = offVal (complex)
void buildComplexTridiag(MatrixComplex& a, int n, Real diagVal,
                         Complex offVal) {
    int colSz[128];
    for (int j = 0; j < n; j++) {
        colSz[j] = 1;
        if (j > 0)   colSz[j]++;
        if (j < n-1)  colSz[j]++;
    }
    a.setColumnSizes(colSz);
    Complex diag(diagVal, 0.0);
    for (int j = 0; j < n; j++) {
        int ri[3]; Complex ev[3]; int k = 0;
        if (j > 0) { ri[k] = j-1; ev[k] = offVal; k++; }
        ri[k] = j; ev[k] = diag; k++;
        if (j < n-1) { ri[k] = j+1; ev[k] = offVal; k++; }
        a.setColumn(j, ri, ev);
    }
    a.validate();
}

#undef __FUNC__
#define __FUNC__ "main"
int main()
{
    BEGIN_FUNCTION();
    OblioSetNewHandler();

    FactorizationAlgorithm algs[] = {LeftLooking, RightLooking, Multifrontal};
    FactorizationType types[] = {StaticCC, StaticLDL, DynamicLDL};

    int pass = 0, fail = 0;
    const int n = 6;
    Complex offVal(-1.0, 0.5);

    cout << "=== 6x6 complex symmetric tridiag (d=8, off=-1+0.5i) ===" << endl;
    for (int ai = 0; ai < 3; ai++) {
        for (int ti = 0; ti < 3; ti++) {
            MatrixComplex a(n);
            buildComplexTridiag(a, n, 8.0, offVal);

            // Set x = {1, 2, ..., n} (complex with zero imag part)
            SingleVectorComplex xtrue(n);
            xtrue.initialize(); // sets to (i+1, 0)

            // Compute b = A * xtrue
            SingleVectorComplex b(n);
            MultiplyEngineComplex mul;
            mul.run(a, xtrue, b);

            SingleVectorComplex x;
            Real error = -1;

            OblioEngineComplex eng;
            eng.setOrderingAlgorithm(MultipleMinimumDegree);
            eng.setFactorizationType(types[ti]);
            eng.setFactorizationAlgorithm(algs[ai]);
            eng.setFactorsDataType(DynamicDataType);
            eng.setStackDataType(DynamicDataType);
            eng.setSupernodal(true);
            eng.setPerturbation(1e-15);
            eng.setPivotingThreshold(1e-1);
            eng.setTrace(false);

            eng.run(a, b, x, &error);

            bool ok = (error >= 0 && error < 1e-10);
            cout << "  " << algName(algs[ai]) << typeName(types[ti])
                 << "  err=" << error
                 << (ok ? "  PASS" : "  FAIL") << endl;
            if (ok) pass++; else fail++;
        }
    }

    // Also test with a larger complex matrix: 4x4 2D Laplacian (n=16)
    // with complex off-diagonal
    cout << "\n=== 4x4 2D complex Laplacian (n=16, off=-1+0.3i) ===" << endl;
    {
        int m = 4, nn = 16;
        Complex coff(-1.0, 0.3);
        Complex cdiag(4.0, 0.0);

        for (int ai = 0; ai < 3; ai++) {
            for (int ti = 0; ti < 3; ti++) {
                MatrixComplex a(nn);

                int colSz[16];
                for (int j = 0; j < nn; j++) {
                    int row = j / m, col = j % m;
                    int sz = 1;
                    if (row > 0)   sz++;
                    if (row < m-1) sz++;
                    if (col > 0)   sz++;
                    if (col < m-1) sz++;
                    colSz[j] = sz;
                }
                a.setColumnSizes(colSz);
                for (int j = 0; j < nn; j++) {
                    int row = j / m, col = j % m;
                    int ri[5]; Complex ev[5]; int k = 0;
                    if (col > 0)   { ri[k] = j-1; ev[k] = coff; k++; }
                    if (row > 0)   { ri[k] = j-m; ev[k] = coff; k++; }
                    ri[k] = j; ev[k] = cdiag; k++;
                    if (row < m-1) { ri[k] = j+m; ev[k] = coff; k++; }
                    if (col < m-1) { ri[k] = j+1; ev[k] = coff; k++; }
                    // Sort by row index
                    for (int a_ = 0; a_ < k-1; a_++)
                        for (int b_ = a_+1; b_ < k; b_++)
                            if (ri[a_] > ri[b_]) {
                                int tr = ri[a_]; ri[a_] = ri[b_]; ri[b_] = tr;
                                Complex te = ev[a_]; ev[a_] = ev[b_]; ev[b_] = te;
                            }
                    a.setColumn(j, ri, ev);
                }
                a.validate();

                SingleVectorComplex xtrue(nn);
                xtrue.initialize();
                SingleVectorComplex b(nn);
                MultiplyEngineComplex mul;
                mul.run(a, xtrue, b);

                SingleVectorComplex x;
                Real error = -1;

                OblioEngineComplex eng;
                eng.setOrderingAlgorithm(MultipleMinimumDegree);
                eng.setFactorizationType(types[ti]);
                eng.setFactorizationAlgorithm(algs[ai]);
                eng.setFactorsDataType(DynamicDataType);
                eng.setStackDataType(DynamicDataType);
                eng.setSupernodal(true);
                eng.setPerturbation(1e-15);
                eng.setPivotingThreshold(1e-1);
                eng.setTrace(false);

                eng.run(a, b, x, &error);

                bool ok = (error >= 0 && error < 1e-10);
                cout << "  " << algName(algs[ai]) << typeName(types[ti])
                     << "  err=" << error
                     << (ok ? "  PASS" : "  FAIL") << endl;
                if (ok) pass++; else fail++;
            }
        }
    }

    cout << "\nTotal: " << pass << "/" << (pass+fail) << " passed" << endl;

    OblioResetNewHandler();
    END_FUNCTION_RETURN(0);
}
