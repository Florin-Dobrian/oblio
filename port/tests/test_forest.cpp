#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/Types.h"
#include "test_util.h"
#include <cstdint>
#include <iostream>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
static const std::int32_t NP = NIL;
template<class Val> static void reqSym(const SparseMatrix<Val>& A, const std::string& lbl){
    ck(OblioTest::isStructurallySymmetric(A), lbl); }
static SparseMatrix<double> tridiagFull(std::size_t n){
    std::vector<std::size_t> cp(n+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
    for(std::size_t j=0;j<n;++j){ if(j>0){ri.push_back(static_cast<std::int32_t>(j-1));v.push_back(-1);} ri.push_back(static_cast<std::int32_t>(j));v.push_back(2);
        if(j+1<n){ri.push_back(static_cast<std::int32_t>(j+1));v.push_back(-1);} cp[j+1]=ri.size(); }
    return SparseMatrix<double>(n,cp,ri,v); }
static bool eq(const std::vector<std::int32_t>& a, std::initializer_list<std::int32_t> b){
    std::vector<std::int32_t> bb(b); if(a.size()!=bb.size())return false;
    for(std::size_t i=0;i<a.size();++i){ if(a[i]!=bb[i]) return false; }
    return true; }
static bool validEtree(const std::vector<std::int32_t>& par, std::size_t& roots){
    roots=0; for(std::size_t i=0;i<par.size();++i){
        if(par[i]==NP){++roots;continue;} std::size_t pv=static_cast<std::size_t>(par[i]); if(pv<=i) return false; if(pv>=par.size()) return false; }
    return true; }
int main(){
    ElmForestEngine eng;
    { auto A=tridiagFull(4); reqSym(A,"tridiag n=4         : symmetric");
      Permutation p(4); ElmForest f;
      ck(eng.compute(A,p,f) && eq(f.parent(),{1,2,3,NP}), "tridiag n=4 natural : path"); }
    { auto A=tridiagFull(6); reqSym(A,"tridiag n=6         : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.parent(),{1,2,3,4,5,NP}), "tridiag n=6 natural : path"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"arrow 6x6           : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.parent(),{1,2,3,4,5,NP}), "arrow 6x6 natural   : path"); }
    { std::vector<std::size_t> cp={0,2,5,7,9,12,14};
      std::vector<std::int32_t> ri={0,1, 0,1,2, 1,2, 3,4, 3,4,5, 4,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"two blocks          : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      ck(inv && roots==2 && eq(f.parent(),{1,2,NP,4,5,NP}), "two blocks natural  : 2 trees"); }
    { auto A=tridiagFull(30); reqSym(A,"tridiag n=30        : symmetric");
      OrderEngine ord(OrderMethod::AMD); Permutation p; ord.order(A,p);
      ElmForest f; bool ok=eng.compute(A,p,f); std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      ck(ok && inv && f.size()==30 && roots>=1, "tridiag n=30 + AMD  : valid etree"); }
    std::cout<<"\nElmForest tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
