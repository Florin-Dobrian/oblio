#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/AmdFlat.h"
#include "oblio/MmdChained.h"
#include "oblio/MmdFlat.h"
#include "oblio/MmdCompacted.h"
#include "oblio/AmdCompacted.h"
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
using OrderFn = ElmOrder(*)(const std::vector<std::size_t>&,
                               const std::vector<std::int32_t>&);
// `orderMmdChained` and `orderMmdCompacted` take a third argument, `delta`, with a default, so
// their type is not OrderFn and neither can be named directly where one is wanted. Forwarded rather than widening
// OrderFn: the helpers exist to compare a free function against an enum member, and delta is the
// mmd layers' business. `orderAmdCompacted` needs no forwarder, the amd branch having no delta.
static ElmOrder mmdChainedDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmdChained(colPtr, rowIdx);
}

static ElmOrder mmdCompactedDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmdCompacted(colPtr, rowIdx);
}

static ElmOrder mmdFlatDefault(const std::vector<std::size_t>&  colPtr,
                                              const std::vector<std::int32_t>& rowIdx) {
    return orderMmdFlat(colPtr, rowIdx);
}

template<class Val> static void checkOrderFn(const SparseMatrix<Val>& A, OrderFn f,
                                             const std::string& lbl){
    const std::vector<std::int32_t> order = f(A.colPtr(), A.rowIdx()).order();
    Permutation P(A.size());
    const bool set = (order.size()==A.size()) && P.setNewToOld(order);
    ck(set && P.validate(), lbl); }
// SAME PERMUTATION, AND SAME WORK BEHIND IT. Three comparisons in one assertion, and they answer
// three different questions.
//
// THE ENUM REACHES THE DRIVER IT NAMES. `m` is the branch's flat driver and `ref` is that same
// driver as a free function, so comparing the engine's permutation against `f`'s covers the
// dispatch as well as the two drivers agreeing.
//
// AND THE TWO DRIVERS DID THE SAME WORK, which the permutation cannot say. Peak live clique
// members and members born are properties of the ALGORITHM rather than of the layout: two drivers
// running the same method form the same cliques and merge the same vertices at the same moments
// whatever their storage, so both figures must agree exactly. Two orderings that agree on output
// while doing different work is the failure this catches, and it is the shape of defect this tree
// has repeatedly found by hand.
//
// THE WORK COMPARISON IS BETWEEN THE TWO FREE FUNCTIONS, not against the engine, because the engine
// calls one of them: comparing its figures with that driver's would compare a call with itself.
//
// ZERO MEANS NOT TRACKED, and then that figure alone is skipped. `MmdChained` tracks neither:
// chained storage ends a clique at a terminator and keeps no size, so subtracting a length on death would
// need a per-vertex array in the one file whose purpose is the mmd oracle's array economy.
template<class Val> static void checkSameOrderFn(const SparseMatrix<Val>& A, Ordering m,
                                                 OrderFn ref, OrderFn f, const std::string& lbl){
    OrderEngine a(m); Permutation pa;
    const bool ok = a.compute(A,pa);
    const ElmOrder r = ref(A.colPtr(), A.rowIdx());
    const ElmOrder s = f(A.colPtr(), A.rowIdx());
    bool same = ok && s.order().size()==A.size();
    for (std::size_t k = 0; same && k < s.order().size(); ++k)
        same = (pa.newToOld()[k] == s.order()[k]);
    if (same && r.numPeakCliqueMembers() != 0 && s.numPeakCliqueMembers() != 0)
        same = (r.numPeakCliqueMembers() == s.numPeakCliqueMembers());
    if (same && r.numBornCliqueMembers() != 0 && s.numBornCliqueMembers() != 0)
        same = (r.numBornCliqueMembers() == s.numBornCliqueMembers());
    ck(same, lbl); }
// TWO ENUMERATORS, ONE PERMUTATION. The mmd branch's oracle: our drivers must return exactly what
// `MmdCorrected` returns, which is `MmdVendored` with its degree scale repaired. `MmdVendored`
// keeps the original and is reference only, so nothing here compares against it.
template<class Val> static void checkSameOrder(const SparseMatrix<Val>& A, Ordering m, Ordering r,
                                               const std::string& lbl){
    OrderEngine em(m); Permutation pm;
    OrderEngine er(r); Permutation pr;
    bool same = em.compute(A,pm) && er.compute(A,pr) && pm.size()==pr.size();
    for (std::size_t k = 0; same && k < pm.size(); ++k)
        same = (pm.newToOld()[k] == pr.newToOld()[k]);
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
    // notes/TESTING_SPECIFICATION.md. Fourteen assertions here are theirs, one pair on each of the
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
      checkOrder(A,Ordering::MmdCorrected,"arrow 6x6      : MmdCorrected valid");
      checkSameOrder(A,Ordering::MmdFlat,Ordering::MmdCorrected,
                     "arrow 6x6      : MmdFlat == MmdCorrected");
#endif
      checkOrder(A,Ordering::MmdFlat,"arrow 6x6      : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"arrow 6x6      : AmdFlat valid");
      // THE THREE NON-ENUM LAYERS, each of which must reproduce its original entry for entry,
      // which is the whole of what makes it a measurement rather than a second ordering. MmdChained
      // is MmdFlat on the mmd oracle's clique storage and AmdCompacted is AmdFlat on the amd
      // oracle's, both permanent; MmdCompacted is MmdFlat on the production layout and is
      // transitional, carrying the amd folds onto the mmd side. The two AMD B layers that used to
      // be checked here were retired on 2026-08-16 when
      // their schedule moved into their originals.
      //
      // Validity is asserted once per layer, here, and sameness on every matrix below. Sameness
      // against an original already checked valid implies validity, so the arrow's three are what
      // exercise the free-function path through setNewToOld rather than the comparison.
      checkOrderFn(A,mmdChainedDefault,"arrow 6x6      : MmdChained valid");
      checkOrderFn(A,mmdCompactedDefault,"arrow 6x6      : MmdCompacted valid");
      checkOrderFn(A,orderAmdCompacted,  "arrow 6x6      : AmdCompacted valid");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdChainedDefault,
                       "arrow 6x6      : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdCompactedDefault,
                       "arrow 6x6      : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdFlat,orderAmdCompacted,
                       "arrow 6x6      : AmdCompacted == AmdFlat"); }
    for(std::size_t size : {1u,2u,10u,100u}){ auto A=tridiagFull(size);
      reqSym(A,"tridiag n="+std::to_string(size)+" : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AmdVendored,"tridiag n="+std::to_string(size)+" : AMD valid");
      checkOrder(A,Ordering::MmdVendored,"tridiag n="+std::to_string(size)+" : MMD valid");
      checkOrder(A,Ordering::MmdCorrected,
                 "tridiag n="+std::to_string(size)+" : MmdCorrected valid");
      checkSameOrder(A,Ordering::MmdFlat,Ordering::MmdCorrected,
                     "tridiag n="+std::to_string(size)+" : MmdFlat == MmdCorrected");
#endif
      checkOrder(A,Ordering::MmdFlat,"tridiag n="+std::to_string(size)+" : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"tridiag n="+std::to_string(size)+" : AmdFlat valid");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdChainedDefault,
                       "tridiag n="+std::to_string(size)+" : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdCompactedDefault,
                       "tridiag n="+std::to_string(size)+" : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdFlat,orderAmdCompacted,
                       "tridiag n="+std::to_string(size)+" : AmdCompacted == AmdFlat"); }
    { std::size_t size=5; std::vector<std::size_t> cp(size+1); std::vector<std::int32_t> ri(size); std::vector<double> v(size,1.0);
      for(std::size_t j=0;j<size;++j){cp[j]=j; ri[j]=static_cast<std::int32_t>(j);} cp[size]=size;
      SparseMatrix<double> A(size,cp,ri,v);
      reqSym(A,"diagonal 5x5   : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(A,Ordering::AmdVendored,"diagonal 5x5   : AMD valid");
      checkOrder(A,Ordering::MmdVendored,"diagonal 5x5   : MMD valid");
      checkOrder(A,Ordering::MmdCorrected,"diagonal 5x5   : MmdCorrected valid");
      checkSameOrder(A,Ordering::MmdFlat,Ordering::MmdCorrected,
                     "diagonal 5x5   : MmdFlat == MmdCorrected");
#endif
      checkOrder(A,Ordering::MmdFlat,"diagonal 5x5   : MmdFlat valid");
      checkOrder(A,Ordering::AmdFlat,"diagonal 5x5   : AmdFlat valid");
      // n ISOLATED VERTICES, where every degree is zero and nothing ever merges. None of the three
      // layers was checked on this shape until 2026-08-17, and it is the one the experiment's own
      // graphs cannot produce: all seven of them are connected and none is trivial.
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdChainedDefault,
                       "diagonal 5x5   : MmdChained == MmdFlat");
      checkSameOrderFn(A,Ordering::MmdFlat,mmdFlatDefault,mmdCompactedDefault,
                       "diagonal 5x5   : MmdCompacted == MmdFlat");
      checkSameOrderFn(A,Ordering::AmdFlat,orderAmdFlat,orderAmdCompacted,
                       "diagonal 5x5   : AmdCompacted == AmdFlat"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<std::complex<double>> v(ri.size(),{1,0}); SparseMatrix<std::complex<double>> C(6,cp,ri,v);
      reqSym(C,"arrow complex  : symmetric");
#ifdef OBLIO_VENDORED_ORDERINGS
      checkOrder(C,Ordering::AmdVendored,"arrow complex  : AMD valid");
      checkOrder(C,Ordering::MmdVendored,"arrow complex  : MMD valid");
      checkOrder(C,Ordering::MmdCorrected,"arrow complex  : MmdCorrected valid");
      checkSameOrder(C,Ordering::MmdFlat,Ordering::MmdCorrected,
                     "arrow complex  : MmdFlat == MmdCorrected");
#endif
      checkOrder(C,Ordering::MmdFlat,"arrow complex  : MmdFlat valid");
      checkOrder(C,Ordering::AmdFlat,"arrow complex  : AmdFlat valid");
      // THE SECOND SCALAR TYPE. An ordering reads only the pattern, so the permutation must be the
      // real arrow's; what this exercises is the structural overloads through a second
      // instantiation of the templated helpers.
      checkSameOrderFn(C,Ordering::MmdFlat,mmdFlatDefault,mmdChainedDefault,
                       "arrow complex  : MmdChained == MmdFlat");
      checkSameOrderFn(C,Ordering::MmdFlat,mmdFlatDefault,mmdCompactedDefault,
                       "arrow complex  : MmdCompacted == MmdFlat");
      checkSameOrderFn(C,Ordering::AmdFlat,orderAmdFlat,orderAmdCompacted,
                       "arrow complex  : AmdCompacted == AmdFlat"); }
    std::cout<<"\nOrderEngine tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
