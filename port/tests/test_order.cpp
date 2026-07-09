#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include <iostream>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
template<class Val> static void checkOrder(const SparseMatrix<Val>& A, OrderMethod m, const std::string& lbl){
    OrderEngine e(m); Permutation p; bool ok=e.order(A,p);
    ck(ok && p.size()==A.numCols() && p.validate(), lbl); }
static SparseMatrix<double> tridiagFull(std::size_t n){
    std::vector<std::size_t> cp(n+1,0), ri; std::vector<double> v;
    for(std::size_t j=0;j<n;++j){ if(j>0){ri.push_back(j-1);v.push_back(-1);} ri.push_back(j);v.push_back(2);
        if(j+1<n){ri.push_back(j+1);v.push_back(-1);} cp[j+1]=ri.size(); }
    return SparseMatrix<double>(n,cp,ri,v); }
int main(){
    std::cout<<"=== OrderEngine tests (AMD / MMD), full-symmetric A ===\n";
    // full arrow
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::size_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      checkOrder(A,OrderMethod::AMD,"arrow 6x6      : AMD valid");
      checkOrder(A,OrderMethod::MMD,"arrow 6x6      : MMD valid"); }
    for(std::size_t n : {1u,2u,10u,100u}){ auto A=tridiagFull(n);
      checkOrder(A,OrderMethod::AMD,"tridiag n="+std::to_string(n)+" : AMD valid");
      checkOrder(A,OrderMethod::MMD,"tridiag n="+std::to_string(n)+" : MMD valid"); }
    // diagonal only
    { std::size_t n=5; std::vector<std::size_t> cp(n+1),ri(n); std::vector<double> v(n,1.0);
      for(std::size_t j=0;j<n;++j){cp[j]=j;ri[j]=j;} cp[n]=n;
      SparseMatrix<double> A(n,cp,ri,v);
      checkOrder(A,OrderMethod::AMD,"diagonal 5x5   : AMD valid");
      checkOrder(A,OrderMethod::MMD,"diagonal 5x5   : MMD valid"); }
    // complex full arrow
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::size_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<std::complex<double>> v(ri.size(),{1,0}); SparseMatrix<std::complex<double>> C(6,cp,ri,v);
      checkOrder(C,OrderMethod::AMD,"arrow complex  : AMD valid");
      checkOrder(C,OrderMethod::MMD,"arrow complex  : MMD valid"); }
    std::cout<<"\nOrderEngine tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
