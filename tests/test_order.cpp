#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/Mmd3B.h"
#include "oblio/Mmd3C.h"
#include "oblio/Amd3B.h"
#include "oblio/QuotientGraph.h"
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
// THE NON-ENUM LAYERS ARE REACHED AS FREE FUNCTIONS, not through the Ordering enum, because none
// of them is an ordering a caller should choose: each is the SAME ordering computed differently,
// so each exists to be measured against its original. See OrderEngine.h.
//
// And these are the strongest oracles in the ordering suite. The vendored pairs and the
// digit-suffixed pairs are different orderings and can only be compared on fill, where one of these
// must agree with its original ENTRY FOR ENTRY, and a difference is a defect rather than a
// tie-break. That matters more than usual here: none of the three has a prototype in
// experiments/ordering or appears in its PORTED list, so this is the only thing checking them at
// all, and `make digest` in benchmarks/ordering is the only other place any of them is exercised.
//
// ALL THREE ARE CHECKED ON EVERY MATRIX, deliberately and uniformly. Mmd3B ran on five of the seven
// until 2026-08-17, missing the 5x5 diagonal and the complex arrow, and Amd3B had no assertion
// anywhere. The two absent matrices are the ones the experiment's own graphs cannot produce: n
// isolated vertices, where every degree is zero and nothing ever merges, and the complex scalar
// type, which exercises the structural overloads on a second instantiation.
using OrderFn = std::vector<std::int32_t>(*)(const std::vector<std::size_t>&,
                                             const std::vector<std::int32_t>&);
// `orderMmd3B` and `orderMmd3C` take a third argument, `delta`, with a default, so their type is
// not OrderFn and neither can be named directly where one is wanted. Forwarded rather than widening
// OrderFn: the helpers exist to compare a free function against an enum member, and delta is the
// mmd layers' business. `orderAmd3B` needs no forwarder, the amd branch having no delta.
static std::vector<std::int32_t> mmd3bDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmd3B(colPtr, rowIdx);
}

static std::vector<std::int32_t> mmd3cDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmd3C(colPtr, rowIdx);
}

template<class Val> static void checkOrderFn(const SparseMatrix<Val>& A, OrderFn f,
                                             const std::string& lbl){
    const std::vector<std::int32_t> order = f(A.colPtr(), A.rowIdx());
    Permutation P(A.size());
    const bool set = (order.size()==A.size()) && P.setNewToOld(order);
    ck(set && P.validate(), lbl); }
// SAME PERMUTATION, AND SAME WORK BEHIND IT. The order comparison is the older half. The second
// half compares PEAK LIVE CLIQUE MEMBERS, which is a property of the algorithm rather than of the
// layout: two drivers running the same method form the same cliques and merge the same vertices at
// the same moments whatever their storage, so the figure must agree exactly. Two orderings that
// agree on output while doing different work is the failure this catches and the permutation
// comparison cannot, and it is the shape of defect this tree has repeatedly found by hand.
//
// ZERO MEANS NOT TRACKED, and then only the order is compared. `Mmd3B` does not track it: chained
// storage ends a clique at a terminator and keeps no size, so subtracting a length on death would
// need a per-vertex array in the one file whose purpose is genmmd's array economy.
template<class Val> static void checkSameOrderFn(const SparseMatrix<Val>& A, Ordering m,
                                                 OrderFn f, const std::string& lbl){
    OrderEngine a(m); Permutation pa;
    gPeakCliqueMembers = 0;
    const bool ok = a.compute(A,pa);
    const std::size_t peakEngine = gPeakCliqueMembers;
    gPeakCliqueMembers = 0;
    const std::vector<std::int32_t> order = f(A.colPtr(), A.rowIdx());
    const std::size_t peakFn = gPeakCliqueMembers;
    bool same = ok && order.size()==A.size();
    for (std::size_t k = 0; same && k < order.size(); ++k)
        same = (pa.newToOld()[k] == order[k]);
    if (same && peakEngine != 0 && peakFn != 0) same = (peakEngine == peakFn);
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
    // the total is 59 with that directory and 45 without. Nothing else changes.
    //
    // IT WAS 87 AND 73 UNTIL 2026-08-21, when MMD1, MMD2, AMD1 and AMD2 were retired to retired/
    // and their sixteen validity assertions went with them. Twelve more went with the enumerators
    // in other suites; see retired/README.md.
    //
    // COUNTS CORRECTED 2026-08-17, having been wrong here and in three other files, every figure
    // internally consistent and none of them the number the suite ran. A count in a comment is
    // worth nothing unless it is read off a run, which is how both of the figures above were
    // obtained.
    std::cout<<"=== OrderEngine tests (AMD / MMD lineages, full-symmetric A) ===\n";
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"arrow 6x6      : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"arrow 6x6      : AMD valid");
      checkOrder(A,Ordering::MMD,"arrow 6x6      : MMD valid");
#endif
      checkOrder(A,Ordering::MMD3,"arrow 6x6      : MMD3 valid");
      checkOrder(A,Ordering::AMD3,"arrow 6x6      : AMD3 valid");
      // THE THREE NON-ENUM LAYERS, each of which must reproduce its original entry for entry,
      // which is the whole of what makes it a measurement rather than a second ordering. Mmd3B is
      // Mmd3 on genmmd's clique storage and Amd3B is Amd3 on AMD_2's, both permanent; Mmd3C is
      // Mmd3 on the production layout and is transitional, carrying the amd folds onto the mmd
      // side. The two AMD B layers that used to be checked here were retired on 2026-08-16 when
      // their schedule moved into their originals.
      //
      // Validity is asserted once per layer, here, and sameness on every matrix below. Sameness
      // against an original already checked valid implies validity, so the arrow's three are what
      // exercise the free-function path through setNewToOld rather than the comparison.
      checkOrderFn(A,mmd3bDefault,"arrow 6x6      : MMD3B valid");
      checkOrderFn(A,mmd3cDefault,"arrow 6x6      : MMD3C valid");
      checkOrderFn(A,orderAmd3B,  "arrow 6x6      : AMD3B valid");
      checkSameOrderFn(A,Ordering::MMD3,mmd3bDefault,"arrow 6x6      : MMD3B == MMD3");
      checkSameOrderFn(A,Ordering::MMD3,mmd3cDefault,"arrow 6x6      : MMD3C == MMD3");
      checkSameOrderFn(A,Ordering::AMD3,orderAmd3B,  "arrow 6x6      : AMD3B == AMD3"); }
    for(std::size_t size : {1u,2u,10u,100u}){ auto A=tridiagFull(size);
      reqSym(A,"tridiag n="+std::to_string(size)+" : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"tridiag n="+std::to_string(size)+" : AMD valid");
      checkOrder(A,Ordering::MMD,"tridiag n="+std::to_string(size)+" : MMD valid");
#endif
      checkOrder(A,Ordering::MMD3,"tridiag n="+std::to_string(size)+" : MMD3 valid");
      checkOrder(A,Ordering::AMD3,"tridiag n="+std::to_string(size)+" : AMD3 valid");
      checkSameOrderFn(A,Ordering::MMD3,mmd3bDefault,
                       "tridiag n="+std::to_string(size)+" : MMD3B == MMD3");
      checkSameOrderFn(A,Ordering::MMD3,mmd3cDefault,
                       "tridiag n="+std::to_string(size)+" : MMD3C == MMD3");
      checkSameOrderFn(A,Ordering::AMD3,orderAmd3B,
                       "tridiag n="+std::to_string(size)+" : AMD3B == AMD3"); }
    { std::size_t size=5; std::vector<std::size_t> cp(size+1); std::vector<std::int32_t> ri(size); std::vector<double> v(size,1.0);
      for(std::size_t j=0;j<size;++j){cp[j]=j; ri[j]=static_cast<std::int32_t>(j);} cp[size]=size;
      SparseMatrix<double> A(size,cp,ri,v);
      reqSym(A,"diagonal 5x5   : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AMD,"diagonal 5x5   : AMD valid");
      checkOrder(A,Ordering::MMD,"diagonal 5x5   : MMD valid");
#endif
      checkOrder(A,Ordering::MMD3,"diagonal 5x5   : MMD3 valid");
      checkOrder(A,Ordering::AMD3,"diagonal 5x5   : AMD3 valid");
      // n ISOLATED VERTICES, where every degree is zero and nothing ever merges. None of the three
      // layers was checked on this shape until 2026-08-17, and it is the one the experiment's own
      // graphs cannot produce: all seven of them are connected and none is trivial.
      checkSameOrderFn(A,Ordering::MMD3,mmd3bDefault,"diagonal 5x5   : MMD3B == MMD3");
      checkSameOrderFn(A,Ordering::MMD3,mmd3cDefault,"diagonal 5x5   : MMD3C == MMD3");
      checkSameOrderFn(A,Ordering::AMD3,orderAmd3B,  "diagonal 5x5   : AMD3B == AMD3"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<std::complex<double>> v(ri.size(),{1,0}); SparseMatrix<std::complex<double>> C(6,cp,ri,v);
      reqSym(C,"arrow complex  : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(C,Ordering::AMD,"arrow complex  : AMD valid");
      checkOrder(C,Ordering::MMD,"arrow complex  : MMD valid");
#endif
      checkOrder(C,Ordering::MMD3,"arrow complex  : MMD3 valid");
      checkOrder(C,Ordering::AMD3,"arrow complex  : AMD3 valid");
      // THE SECOND SCALAR TYPE. An ordering reads only the pattern, so the permutation must be the
      // real arrow's; what this exercises is the structural overloads through a second
      // instantiation of the templated helpers.
      checkSameOrderFn(C,Ordering::MMD3,mmd3bDefault,"arrow complex  : MMD3B == MMD3");
      checkSameOrderFn(C,Ordering::MMD3,mmd3cDefault,"arrow complex  : MMD3C == MMD3");
      checkSameOrderFn(C,Ordering::AMD3,orderAmd3B,  "arrow complex  : AMD3B == AMD3"); }
    std::cout<<"\nOrderEngine tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
