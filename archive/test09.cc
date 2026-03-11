// test09.cc — test oblio 0.9 with all 9 factorization combinations
// on a 5×5 and 6×6 tridiagonal, plus 5×5 2D Laplacian (n=25).
#include <iostream>
#include <cmath>
using namespace std;

#include "Oblio.h"

#undef __CLASS__
#define __CLASS__ "test09"

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

void buildTridiag(MatrixReal& a, int n, Real diagVal, Real offVal) {
    int colSz[128];
    for (int j = 0; j < n; j++) {
        colSz[j] = 1;
        if (j > 0)   colSz[j]++;
        if (j < n-1)  colSz[j]++;
    }
    a.setColumnSizes(colSz);
    for (int j = 0; j < n; j++) {
        int ri[3]; Real ev[3]; int k = 0;
        if (j > 0) { ri[k] = j-1; ev[k] = offVal; k++; }
        ri[k] = j; ev[k] = diagVal; k++;
        if (j < n-1) { ri[k] = j+1; ev[k] = offVal; k++; }
        a.setColumn(j, ri, ev);
    }
    a.validate();
}

void build2DLaplacian(MatrixReal& a, int m) {
    int n = m * m;
    int colSz[1024];
    for (int j = 0; j < n; j++) {
        int row = j / m, col = j % m;
        int sz = 1;
        if (row > 0)   sz++;
        if (row < m-1) sz++;
        if (col > 0)   sz++;
        if (col < m-1) sz++;
        colSz[j] = sz;
    }
    a.setColumnSizes(colSz);
    for (int j = 0; j < n; j++) {
        int row = j / m, col = j % m;
        int ri[5]; Real ev[5]; int k = 0;
        if (col > 0)   { ri[k] = j - 1; ev[k] = -1.0; k++; }
        if (row > 0)   { ri[k] = j - m; ev[k] = -1.0; k++; }
        ri[k] = j; ev[k] = 4.0; k++;
        if (row < m-1) { ri[k] = j + m; ev[k] = -1.0; k++; }
        if (col < m-1) { ri[k] = j + 1; ev[k] = -1.0; k++; }
        // Sort by row index (bubble sort on small array)
        for (int a_ = 0; a_ < k-1; a_++)
            for (int b_ = a_+1; b_ < k; b_++)
                if (ri[a_] > ri[b_]) {
                    int tr = ri[a_]; ri[a_] = ri[b_]; ri[b_] = tr;
                    Real te = ev[a_]; ev[a_] = ev[b_]; ev[b_] = te;
                }
        a.setColumn(j, ri, ev);
    }
    a.validate();
}

int runTest(const char* label, int n,
            void (*builder)(MatrixReal&, int, Real, Real),
            int barg, Real d, Real o,
            void (*builder2d)(MatrixReal&, int) = 0, int m2d = 0)
{
    FactorizationAlgorithm algs[] = {LeftLooking, RightLooking, Multifrontal};
    FactorizationType types[] = {StaticCC, StaticLDL, DynamicLDL};
    int pass = 0, fail = 0;

    cout << "=== " << label << " ===" << endl;
    for (int ai = 0; ai < 3; ai++) {
        for (int ti = 0; ti < 3; ti++) {
            MatrixReal a(n);
            if (builder2d)
                builder2d(a, m2d);
            else
                builder(a, n, d, o);

            SingleVectorReal xtrue(n);
            xtrue.initialize();
            SingleVectorReal b(n);
            MultiplyEngineReal mul;
            mul.run(a, xtrue, b);

            SingleVectorReal x;
            Real error = -1;

            OblioEngineReal eng;
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
    return pass;
}

#undef __FUNC__
#define __FUNC__ "main"
int main()
{
    BEGIN_FUNCTION();
    OblioSetNewHandler();

    int total_pass = 0;
    int total_tests = 0;

    // T1: 5×5 tridiagonal
    total_pass += runTest("5x5 tridiag (d=4 o=-1)", 5, buildTridiag, 5, 4.0, -1.0);
    total_tests += 9;

    // T2: 6×6 tridiagonal (same size as our complex test)
    total_pass += runTest("6x6 tridiag (d=8 o=-1)", 6, buildTridiag, 6, 8.0, -1.0);
    total_tests += 9;

    // T3: 5×5 2D Laplacian (n=25)
    {
        FactorizationAlgorithm algs[] = {LeftLooking, RightLooking, Multifrontal};
        FactorizationType types[] = {StaticCC, StaticLDL, DynamicLDL};
        int pass = 0;

        cout << "\n=== 5x5 2D Laplacian (n=25) ===" << endl;
        for (int ai = 0; ai < 3; ai++) {
            for (int ti = 0; ti < 3; ti++) {
                MatrixReal a(25);
                build2DLaplacian(a, 5);

                SingleVectorReal xtrue(25);
                xtrue.initialize();
                SingleVectorReal b(25);
                MultiplyEngineReal mul;
                mul.run(a, xtrue, b);

                SingleVectorReal x;
                Real error = -1;

                OblioEngineReal eng;
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
                if (ok) pass++;
            }
        }
        total_pass += pass;
        total_tests += 9;
    }

    cout << "\nTotal: " << total_pass << "/" << total_tests << " passed" << endl;

    OblioResetNewHandler();
    END_FUNCTION_RETURN(0);
}
