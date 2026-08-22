#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/MmdChained.h"
#include "oblio/MmdCompacted.h"
#include "oblio/AmdCompacted.h"
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
// ALL THREE ARE CHECKED ON EVERY MATRIX, deliberately and uniformly. MmdChained ran on five of the
// seven until 2026-08-17, missing the 5x5 diagonal and the complex arrow, and AmdCompacted had
// no assertion anywhere. The two absent matrices are the ones the experiment's own graphs cannot produce: n
// isolated vertices, where every degree is zero and nothing ever merges, and the complex scalar
// type, which exercises the structural overloads on a second instantiation.
using OrderFn = std::vector<std::int32_t>(*)(const std::vector<std::size_t>&,
                                             const std::vector<std::int32_t>&);
// `orderMmdChained` and `orderMmdCompacted` take a third argument, `delta`, with a default, so
// their type is not OrderFn and neither can be named directly where one is wanted. Forwarded rather than widening
// OrderFn: the helpers exist to compare a free function against an enum member, and delta is the
// mmd layers' business. `orderAmdCompacted` needs no forwarder, the amd branch having no delta.
static std::vector<std::int32_t> mmdChainedDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmdChained(colPtr, rowIdx);
}

static std::vector<std::int32_t> mmdCompactedDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmdCompacted(colPtr, rowIdx);
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
// ZERO MEANS NOT TRACKED, and then only the order is compared. `MmdChained` does not track it:
// chained storage ends a clique at a terminator and keeps no size, so subtracting a length on death would
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
      checkOrder(A,Ordering::AmdVendored,"arrow 6x6      : AMD valid");
      checkOrder(A,Ordering::MmdVendored,"arrow 6x6      : MMD valid");
#endif
      checkOrder(A,Ordering::MmdFlat,"arrow 6x6      : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"arrow 6x6      : AmdFlat valid");
      // THE THREE NON-ENUM LAYERS, each of which must reproduce its original entry for entry,
      // which is the whole of what makes it a measurement rather than a second ordering. MmdChained
      // is MmdFlat on genmmd's clique storage and AmdCompacted is AmdFlat on AMD_2's, both
      // permanent; MmdCompacted is MmdFlat on the production layout and is transitional,
      // carrying the amd folds onto
      // the mmd side. The two AMD B layers that used to be checked here were retired on 2026-08-16 when
      // their schedule moved into their originals.
      //
      // Validity is asserted once per layer, here, and sameness on every matrix below. Sameness
      // against an original already checked valid implies validity, so the arrow's three are what
      // exercise the free-function path through setNewToOld rather than the comparison.
      checkOrderFn(A,mmdChainedDefault,"arrow 6x6      : MmdChained valid");
      checkOrderFn(A,mmdCompactedDefault,"arrow 6x6      : MmdCompacted valid");
      checkOrderFn(A,orderAmdCompacted,  "arrow 6x6      : AmdCompacted valid");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdChainedDefault,
                       "arrow 6x6      : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdCompactedDefault,
                       "arrow 6x6      : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdCompacted,
                       "arrow 6x6      : AmdCompacted == AmdFlat"); }
    for(std::size_t size : {1u,2u,10u,100u}){ auto A=tridiagFull(size);
      reqSym(A,"tridiag n="+std::to_string(size)+" : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AmdVendored,"tridiag n="+std::to_string(size)+" : AMD valid");
      checkOrder(A,Ordering::MmdVendored,"tridiag n="+std::to_string(size)+" : MMD valid");
#endif
      checkOrder(A,Ordering::MmdFlat,"tridiag n="+std::to_string(size)+" : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"tridiag n="+std::to_string(size)+" : AmdFlat valid");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdChainedDefault,
                       "tridiag n="+std::to_string(size)+" : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdCompactedDefault,
                       "tridiag n="+std::to_string(size)+" : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdCompacted,
                       "tridiag n="+std::to_string(size)+" : AmdCompacted == AmdFlat"); }
    { std::size_t size=5; std::vector<std::size_t> cp(size+1); std::vector<std::int32_t> ri(size); std::vector<double> v(size,1.0);
      for(std::size_t j=0;j<size;++j){cp[j]=j; ri[j]=static_cast<std::int32_t>(j);} cp[size]=size;
      SparseMatrix<double> A(size,cp,ri,v);
      reqSym(A,"diagonal 5x5   : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AmdVendored,"diagonal 5x5   : AMD valid");
      checkOrder(A,Ordering::MmdVendored,"diagonal 5x5   : MMD valid");
#endif
      checkOrder(A,Ordering::MmdFlat,"diagonal 5x5   : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"diagonal 5x5   : AmdFlat valid");
      // n ISOLATED VERTICES, where every degree is zero and nothing ever merges. None of the three
      // layers was checked on this shape until 2026-08-17, and it is the one the experiment's own
      // graphs cannot produce: all seven of them are connected and none is trivial.
      checkSameOrderFn(A,Ordering::MmdFlat,mmdChainedDefault,
                       "diagonal 5x5   : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdCompactedDefault,
                       "diagonal 5x5   : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdCompacted,
                       "diagonal 5x5   : AmdCompacted == AmdFlat"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<std::complex<double>> v(ri.size(),{1,0}); SparseMatrix<std::complex<double>> C(6,cp,ri,v);
      reqSym(C,"arrow complex  : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(C,Ordering::AmdVendored,"arrow complex  : AMD valid");
      checkOrder(C,Ordering::MmdVendored,"arrow complex  : MMD valid");
#endif
      checkOrder(C,Ordering::MmdFlat,"arrow complex  : MmdFlat valid");
      checkOrder(C,Ordering::AmdFlat,"arrow complex  : AmdFlat valid");
      // THE SECOND SCALAR TYPE. An ordering reads only the pattern, so the permutation must be the
      // real arrow's; what this exercises is the structural overloads through a second
      // instantiation of the templated helpers.
      checkSameOrderFn(C,Ordering::MmdFlat,mmdChainedDefault,
                       "arrow complex  : MmdChained == MmdFlat");
      checkSameOrderFn(C,Ordering::MmdFlat,mmdCompactedDefault,
                       "arrow complex  : MmdCompacted == MmdFlat");
      checkSameOrderFn(C,Ordering::AmdFlat,orderAmdCompacted,
                       "arrow complex  : AmdCompacted == AmdFlat"); }
    std::cout<<"\nOrderEngine tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
