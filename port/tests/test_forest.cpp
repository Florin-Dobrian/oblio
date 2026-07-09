#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include <iostream>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
static const std::size_t NP = ElmForest::cNoParent;

// FULL-symmetric tridiagonal of size n: col j = {j-1, j, j+1} (clipped).
static SparseMatrix<double> tridiagFull(std::size_t n){
    std::vector<std::size_t> cp(n+1,0), ri; std::vector<double> v;
    for(std::size_t j=0;j<n;++j){
        if(j>0){ ri.push_back(j-1); v.push_back(-1.0); }
        ri.push_back(j); v.push_back(2.0);
        if(j+1<n){ ri.push_back(j+1); v.push_back(-1.0); }
        cp[j+1]=ri.size();
    }
    return SparseMatrix<double>(n,cp,ri,v);
}
static bool eq(const std::vector<std::size_t>& a, std::initializer_list<std::size_t> b){
    std::vector<std::size_t> bb(b); if(a.size()!=bb.size())return false;
    for(std::size_t i=0;i<a.size();++i) if(a[i]!=bb[i])return false; return true; }
static bool validEtree(const std::vector<std::size_t>& par, std::size_t& roots){
    roots=0; for(std::size_t i=0;i<par.size();++i){
        if(par[i]==NP){++roots;continue;} if(par[i]<=i) return false; if(par[i]>=par.size()) return false; }
    return true; }

int main(){
    ElmForestEngine eng;
    { auto A=tridiagFull(4); Permutation p(4); ElmForest f;
      ck(eng.compute(A,p,f) && eq(f.parent(),{1,2,3,NP}), "tridiag n=4 natural : path"); }
    { auto A=tridiagFull(6); Permutation p(6); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.parent(),{1,2,3,4,5,NP}), "tridiag n=6 natural : path"); }
    // FULL arrow: col0={0..5}; colj={0,j}
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::size_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0);
      SparseMatrix<double> A(6,cp,ri,v); Permutation p(6); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.parent(),{1,2,3,4,5,NP}), "arrow 6x6 natural   : path"); }
    // FULL two disconnected tridiag blocks 0-1-2 and 3-4-5
    { std::vector<std::size_t> cp={0,2,5,7,9,12,14};
      std::vector<std::size_t> ri={0,1, 0,1,2, 1,2, 3,4, 3,4,5, 4,5};
      std::vector<double> v(ri.size(),1.0);
      SparseMatrix<double> A(6,cp,ri,v); Permutation p(6); ElmForest f; eng.compute(A,p,f);
      std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      ck(inv && roots==2 && eq(f.parent(),{1,2,NP,4,5,NP}), "two blocks natural  : 2 trees"); }
    { auto A=tridiagFull(30); OrderEngine ord(OrderMethod::AMD); Permutation p; ord.order(A,p);
      ElmForest f; bool ok=eng.compute(A,p,f); std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      ck(ok && inv && f.size()==30 && roots>=1, "tridiag n=30 + AMD  : valid etree"); }
    std::cout<<"\nElmForest tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
