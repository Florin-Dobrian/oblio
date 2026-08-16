#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/Mmd3B.h"
#include "test_util.h"
#include <cstdint>
#include <iostream>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
template<class Val> static void checkOrder(const SparseMatrix<Val>& A, Ordering m, const std::string& lbl){
    OrderEngine e(m); Permutation P; bool ok=e.compute(A,P);
    ck(ok && P.size()==A.size() && P.validate(), lbl); }
// THE B LAYERS ARE REACHED AS FREE FUNCTIONS, not through the Ordering enum, because a B is not
// an ordering a caller should choose: it is the SAME ordering computed on a different schedule,
// so it exists to be measured against its original. See OrderEngine.h. Every check they had
// through the enum is kept here unchanged in what it asserts; only the call changes.
//
// And these are the strongest oracles in the ordering suite. The vendored pairs and the
// digit-suffixed pairs are different orderings and can only be compared on fill, where a B pair
// must agree ENTRY FOR ENTRY, and a difference is a defect rather than a tie-break. That matters
// more than usual for these two: they have no prototype in experiments/ordering and are not in
// its PORTED list, so this is the only thing checking them at all.
using OrderFn = std::vector<std::int32_t>(*)(const std::vector<std::size_t>&,
                                             const std::vector<std::int32_t>&);
// `orderMmd3B` takes a third argument, `delta`, with a default, so its type is not OrderFn and it
// cannot be named directly where one is wanted. Forwarded rather than widening OrderFn: the helpers
// exist to compare a free function against an enum member, and delta is Mmd3B's business.
static std::vector<std::int32_t> mmd3bDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmd3B(colPtr, rowIdx);
}

template<class Val> static void checkOrderFn(const SparseMatrix<Val>& A, OrderFn f,
                                             const std::string& lbl){
    const std::vector<std::int32_t> order = f(A.colPtr(), A.rowIdx());
    Permutation P(A.size());
    const bool set = (order.size()==A.size()) && P.setNewToOld(order);
    ck(set && P.validate(), lbl); }
template<class Val> static void checkSameOrderFn(const SparseMatrix<Val>& A, Ordering m,
                                                 OrderFn f, const std::string& lbl){
    OrderEngine a(m); Permutation pa;
    const bool ok = a.compute(A,pa);
    const std::vector<std::int32_t> order = f(A.colPtr(), A.rowIdx());
    bool same = ok && order.size()==A.size();
    for (std::size_t k = 0; same && k < order.size(); ++k)
        same = (pa.newToOld()[k] == order[k]);
    ck(same, lbl); }
template<class Val> static void reqSym(const SparseMatrix<Val>& A, const std::string& lbl){
    ck(OblioTest::isStructurallySymmetric(A), lbl); }
static SparseMatrix<double> tridiagFull(std::size_t size){
    std::vector<std::size_t> cp(size+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
    for(std::size_t j=0;j<size;++j){ if(j>0){ri.push_back(static_cast<std::int32_t>(j-1));v.push_back(-1);} ri.push_back(static_cast<std::int32_t>(j));v.push_back(2);
        if(j+1<size){ri.push_back(static_cast<std::int32_t>(j+1));v.push_back(-1);} cp[j+1]=ri.size(); }
    return SparseMatrix<double>(size,cp,ri,v); }
int main(){
        // The vendored MMD and AMD are checked only when private/ supplies them; see
    // docs/TESTING_SPECIFICATION.md. Fourteen assertions here are theirs, one pair on each of the
    // arrow, the diagonal and the complex arrow and one pair per size in the tridiagonal loop, so
    // the total is 77 with that directory and 63 without. Nothing else changes.
    std::cout<<"=== OrderEngine tests (AMD / MMD lineages, full-symmetric A) ===\n";
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"arrow 6x6      : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"arrow 6x6      : AMD valid");
      checkOrder(A,Ordering::MMD,"arrow 6x6      : MMD valid");
#endif
      checkOrder(A,Ordering::MMD1,"arrow 6x6      : MMD1 valid");
      checkOrder(A,Ordering::MMD2,"arrow 6x6      : MMD2 valid");
      checkOrder(A,Ordering::MMD3,"arrow 6x6      : MMD3 valid");
      checkOrder(A,Ordering::AMD1,"arrow 6x6      : AMD1 valid");
      checkOrder(A,Ordering::AMD2,"arrow 6x6      : AMD2 valid");
      checkOrder(A,Ordering::AMD3,"arrow 6x6      : AMD3 valid");
      // THE ONE REMAINING B LAYER. Mmd3B is Mmd3 on the vendored clique storage scheme and must
      // reproduce its permutation entry for entry, which is the whole of what makes it a
      // measurement rather than a second ordering. The two AMD B layers that used to be checked
      // here were retired on 2026-08-16 when their schedule moved into their originals.
      checkOrderFn(A,mmd3bDefault,"arrow 6x6      : MMD3B valid");
      checkSameOrderFn(A,Ordering::MMD3,mmd3bDefault,"arrow 6x6      : MMD3B == MMD3"); }
    for(std::size_t size : {1u,2u,10u,100u}){ auto A=tridiagFull(size);
      reqSym(A,"tridiag n="+std::to_string(size)+" : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"tridiag n="+std::to_string(size)+" : AMD valid");
      checkOrder(A,Ordering::MMD,"tridiag n="+std::to_string(size)+" : MMD valid");
#endif
      checkOrder(A,Ordering::MMD1,"tridiag n="+std::to_string(size)+" : MMD1 valid");
      checkOrder(A,Ordering::MMD2,"tridiag n="+std::to_string(size)+" : MMD2 valid");
      checkOrder(A,Ordering::MMD3,"tridiag n="+std::to_string(size)+" : MMD3 valid");
      checkOrder(A,Ordering::AMD1,"tridiag n="+std::to_string(size)+" : AMD1 valid");
      checkOrder(A,Ordering::AMD2,"tridiag n="+std::to_string(size)+" : AMD2 valid");
      checkOrder(A,Ordering::AMD3,"tridiag n="+std::to_string(size)+" : AMD3 valid");
      checkSameOrderFn(A,Ordering::MMD3,mmd3bDefault,
                       "tridiag n="+std::to_string(size)+" : MMD3B == MMD3"); }
    { std::size_t size=5; std::vector<std::size_t> cp(size+1); std::vector<std::int32_t> ri(size); std::vector<double> v(size,1.0);
      for(std::size_t j=0;j<size;++j){cp[j]=j; ri[j]=static_cast<std::int32_t>(j);} cp[size]=size;
      SparseMatrix<double> A(size,cp,ri,v);
      reqSym(A,"diagonal 5x5   : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"diagonal 5x5   : AMD valid");
      checkOrder(A,Ordering::MMD,"diagonal 5x5   : MMD valid");
#endif
      checkOrder(A,Ordering::MMD1,"diagonal 5x5   : MMD1 valid");
      checkOrder(A,Ordering::MMD2,"diagonal 5x5   : MMD2 valid");
      checkOrder(A,Ordering::MMD3,"diagonal 5x5   : MMD3 valid");
      checkOrder(A,Ordering::AMD1,"diagonal 5x5   : AMD1 valid");
      checkOrder(A,Ordering::AMD2,"diagonal 5x5   : AMD2 valid");
      checkOrder(A,Ordering::AMD3,"diagonal 5x5   : AMD3 valid"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<std::complex<double>> v(ri.size(),{1,0}); SparseMatrix<std::complex<double>> C(6,cp,ri,v);
      reqSym(C,"arrow complex  : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(C,Ordering::AMD,"arrow complex  : AMD valid");
      checkOrder(C,Ordering::MMD,"arrow complex  : MMD valid");
#endif
      checkOrder(C,Ordering::MMD1,"arrow complex  : MMD1 valid");
      checkOrder(C,Ordering::MMD2,"arrow complex  : MMD2 valid");
      checkOrder(C,Ordering::MMD3,"arrow complex  : MMD3 valid");
      checkOrder(C,Ordering::AMD1,"arrow complex  : AMD1 valid");
      checkOrder(C,Ordering::AMD2,"arrow complex  : AMD2 valid");
      checkOrder(C,Ordering::AMD3,"arrow complex  : AMD3 valid"); }
    std::cout<<"\nOrderEngine tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
