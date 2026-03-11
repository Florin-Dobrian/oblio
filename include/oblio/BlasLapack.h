#pragma once
#include "oblio/Types.h"
#include <complex>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace Oblio {
inline int toBI(Size n) {
    if (n > (Size)std::numeric_limits<int>::max())
        throw std::overflow_error("dimension exceeds BLAS int");
    return (int)n;
}
} // namespace Oblio

// ---- Fortran symbol convention -------------------------------------------
#ifdef OBLIO_BLAS_UNDERSCORE
extern "C" {
void dgemm_ (char*,char*,int*,int*,int*,double*,const double*,int*,const double*,int*,double*,double*,int*);
void dsyrk_ (char*,char*,int*,int*,double*,const double*,int*,double*,double*,int*);
void dtrsm_ (char*,char*,char*,char*,int*,int*,double*,const double*,int*,double*,int*);
void dpotrf_(char*,int*,double*,int*,int*);
void zgemm_ (char*,char*,int*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,const std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
void zsyrk_ (char*,char*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
void zherk_ (char*,char*,int*,int*,double*,const std::complex<double>*,int*,double*,std::complex<double>*,int*);
void ztrsm_ (char*,char*,char*,char*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,std::complex<double>*,int*);
void zpotrf_(char*,int*,std::complex<double>*,int*,int*);
}
#else
extern "C" {
void dgemm (char*,char*,int*,int*,int*,double*,const double*,int*,const double*,int*,double*,double*,int*);
void dsyrk (char*,char*,int*,int*,double*,const double*,int*,double*,double*,int*);
void dtrsm (char*,char*,char*,char*,int*,int*,double*,const double*,int*,double*,int*);
void dpotrf(char*,int*,double*,int*,int*);
void zgemm (char*,char*,int*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,const std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
void zsyrk (char*,char*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
void zherk (char*,char*,int*,int*,double*,const std::complex<double>*,int*,double*,std::complex<double>*,int*);
void ztrsm (char*,char*,char*,char*,int*,int*,std::complex<double>*,const std::complex<double>*,int*,std::complex<double>*,int*);
void zpotrf(char*,int*,std::complex<double>*,int*,int*);
}
#endif

namespace Oblio {

// ---- BlasTraits<Val> -------------------------------------------------------
template<class Val> struct BT;

template<> struct BT<double> {
    static void gemm(char ta,char tb,int m,int n,int k,double a,const double*A,int lda,const double*B,int ldb,double b,double*C,int ldc){
#ifdef OBLIO_BLAS_UNDERSCORE
        dgemm_(&ta,&tb,&m,&n,&k,&a,A,&lda,B,&ldb,&b,C,&ldc);
#else
        dgemm (&ta,&tb,&m,&n,&k,&a,A,&lda,B,&ldb,&b,C,&ldc);
#endif
    }
    static void syrk(char up,char tr,int n,int k,double a,const double*A,int lda,double b,double*C,int ldc){
#ifdef OBLIO_BLAS_UNDERSCORE
        dsyrk_(&up,&tr,&n,&k,&a,A,&lda,&b,C,&ldc);
#else
        dsyrk (&up,&tr,&n,&k,&a,A,&lda,&b,C,&ldc);
#endif
    }
    // rankUpdate: A*A^T for real (syrk). Uniform interface used by Cholesky factor.
    static void rankUpdate(char up,int n,int k,double a,const double*A,int lda,double b,double*C,int ldc){
        syrk(up,'N',n,k,a,A,lda,b,C,ldc);
    }
    static void trsm(char s,char up,char ta,char d,int m,int n,double a,const double*A,int lda,double*B,int ldb){
#ifdef OBLIO_BLAS_UNDERSCORE
        dtrsm_(&s,&up,&ta,&d,&m,&n,&a,A,&lda,B,&ldb);
#else
        dtrsm (&s,&up,&ta,&d,&m,&n,&a,A,&lda,B,&ldb);
#endif
    }
    static int potrf(char up,int n,double*A,int lda){
        int info=0;
#ifdef OBLIO_BLAS_UNDERSCORE
        dpotrf_(&up,&n,A,&lda,&info);
#else
        dpotrf (&up,&n,A,&lda,&info);
#endif
        return info;
    }
};

template<> struct BT<std::complex<double>> {
    using C=std::complex<double>;
    static void gemm(char ta,char tb,int m,int n,int k,C a,const C*A,int lda,const C*B,int ldb,C b,C*Cv,int ldc){
#ifdef OBLIO_BLAS_UNDERSCORE
        zgemm_(&ta,&tb,&m,&n,&k,&a,A,&lda,B,&ldb,&b,Cv,&ldc);
#else
        zgemm (&ta,&tb,&m,&n,&k,&a,A,&lda,B,&ldb,&b,Cv,&ldc);
#endif
    }
    static void syrk(char up,char tr,int n,int k,C a,const C*A,int lda,C b,C*Cv,int ldc){
#ifdef OBLIO_BLAS_UNDERSCORE
        zsyrk_(&up,&tr,&n,&k,&a,A,&lda,&b,Cv,&ldc);
#else
        zsyrk (&up,&tr,&n,&k,&a,A,&lda,&b,Cv,&ldc);
#endif
    }
    // rankUpdate: A*A^H for complex (herk). Uniform interface used by Cholesky factor.
    // Note: zherk takes real alpha/beta, not complex.
    static void rankUpdate(char up,int n,int k,C a,const C*A,int lda,C b,C*Cv,int ldc){
        double ra=std::real(a), rb=std::real(b); char nt='N';
#ifdef OBLIO_BLAS_UNDERSCORE
        zherk_(&up,&nt,&n,&k,&ra,A,&lda,&rb,Cv,&ldc);
#else
        zherk (&up,&nt,&n,&k,&ra,A,&lda,&rb,Cv,&ldc);
#endif
    }
    static void trsm(char s,char up,char ta,char d,int m,int n,C a,const C*A,int lda,C*B,int ldb){
#ifdef OBLIO_BLAS_UNDERSCORE
        ztrsm_(&s,&up,&ta,&d,&m,&n,&a,A,&lda,B,&ldb);
#else
        ztrsm (&s,&up,&ta,&d,&m,&n,&a,A,&lda,B,&ldb);
#endif
    }
    static int potrf(char up,int n,C*A,int lda){
        int info=0;
#ifdef OBLIO_BLAS_UNDERSCORE
        zpotrf_(&up,&n,A,&lda,&info);
#else
        zpotrf (&up,&n,A,&lda,&info);
#endif
        return info;
    }
};

// ---- Custom oblio kernels (ported from 0.9 BlasLapack.cc) -----------------

// oblioComputeU: U(n x k) = D(diag k x k) * L^T(k x n)
// L is n x k col-major (ldl), D is k x k col-major (ldd, only diagonal used),
// U is k x n col-major (ldu).
template<class Val>
inline void oblioComputeU(int n,int k,
    const Val*l,int ldl, Val*u,int ldu, const Val*d,int ldd)
{
    for(int i=0,lp=0,up=0,dp=0;i<k;i++,lp+=ldl,up++,dp+=ldd+1)
        for(int j=0,lq=lp,uq=up;j<n;j++,lq++,uq+=ldu)
            u[uq]=d[dp]*l[lq];
}

// oblioGemm: A(n x n, lower tri) -= L(n x k) * U(k x n)
template<class Val>
inline void oblioGemm(int n,int k,
    const Val*l,int ldl, const Val*u,int ldu, Val*a,int lda)
{
    if(n==1){Val a1=-1,b1=1; BT<Val>::gemm('N','N',1,1,k,a1,l,ldl,u,ldu,b1,a,lda);return;}
    int n1=n/2,n2=n-n1;
    const Val*l11=l,*l21=l+n1,*u11=u,*u12=u+ldu*n1;
    Val*a11=a,*a21=a+n1,*a22=a21+lda*n1;
    oblioGemm(n1,k,l11,ldl,u11,ldu,a11,lda);
    Val am1=-1,b1=1;
    BT<Val>::gemm('N','N',n2,n1,k,am1,l21,ldl,u11,ldu,b1,a21,lda);
    oblioGemm(n2,k,l21,ldl,u12,ldu,a22,lda);
}

// oblioPotrfCC: recursive blocked Cholesky (no LAPACK required).
template<class Val>
inline int oblioPotrfCC(int n,Val*a,int lda)
{
    if(n==1){if(std::real(*a)<=0)return 1;return 0;}
    int n1=n/2,n2=n-n1;
    Val*a11=a,*a12=a+lda*n1,*a21=a+n1,*a22=a21+lda*n1;
    int info=oblioPotrfCC(n1,a11,lda); if(info>0)return info;
    Val one=1; BT<Val>::trsm('R','U','N','N',n2,n1,one,a11,lda,a21,lda);
    oblioComputeU(n2,n1,a21,lda,a12,lda,a11,lda);
    oblioGemm(n2,n1,a21,lda,a12,lda,a22,lda);
    return oblioPotrfCC(n2,a22,lda);
}

// oblioPotrfLDL: same but with diagonal perturbation (StaticLDL).
template<class Val>
inline int oblioPotrfLDL(int n,Val*a,int lda,RVal eps,int*npert)
{
    if(n==1){if(std::abs(*a)<std::abs(eps)){*a=Val{eps};(*npert)++;}return 0;}
    int n1=n/2,n2=n-n1;
    Val*a11=a,*a12=a+lda*n1,*a21=a+n1,*a22=a21+lda*n1;
    int info=oblioPotrfLDL(n1,a11,lda,eps,npert); if(info>0)return info;
    Val one=1; BT<Val>::trsm('R','U','N','N',n2,n1,one,a11,lda,a21,lda);
    oblioComputeU(n2,n1,a21,lda,a12,lda,a11,lda);
    oblioGemm(n2,n1,a21,lda,a12,lda,a22,lda);
    return oblioPotrfLDL(n2,a22,lda,eps,npert);
}

template<class Val>
inline RVal absv(Val v){return std::abs(v);}

// conjTrans<Val>() — returns the BLAS transpose character for the
// conjugate-transpose of a Val matrix:
//   real:    'T'  (transpose == conjugate-transpose)
//   complex: 'C'  (conjugate-transpose)
// Needed because zpotrf produces L s.t. A = L*L^H, so factor and
// backward solve must use L^H (not L^T) for complex types.
template<class Val> inline char conjTrans() { return 'T'; }
template<> inline char conjTrans<std::complex<double>>() { return 'C'; }

// conjv<Val>(x) — conjugate that stays in type Val.
// For real Val, conj is identity. For complex, uses std::conj.
template<class Val> inline Val conjv(Val x) { return x; }
template<> inline std::complex<double> conjv(std::complex<double> x) { return std::conj(x); }

} // namespace Oblio
