#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/SymFact.h"
#include "oblio/SymFactEngine.h"
#include "oblio/Types.h"
#include "test_util.h"
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
template<class Val> static void reqSym(const SparseMatrix<Val>& A, const std::string& lbl){
    ck(OblioTest::isStructurallySymmetric(A), lbl); }

// Build a full-symmetric matrix from an edge list (diagonal added automatically).
static SparseMatrix<double> fromEdges(std::size_t size,
        const std::vector<std::pair<std::size_t,std::size_t>>& edges){
    std::vector<std::vector<bool>> M(size, std::vector<bool>(size,false));
    for(std::size_t i=0;i<size;++i) M[i][i]=true;
    for(const auto& e : edges){ M[e.first][e.second]=true; M[e.second][e.first]=true; }
    std::vector<std::size_t> cp(size+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
    for(std::size_t j=0;j<size;++j){
        for(std::size_t i=0;i<size;++i) if(M[i][j]){ ri.push_back(static_cast<std::int32_t>(i));
            v.push_back(i==j ? 10.0 : -1.0); }
        cp[j+1]=ri.size(); }
    return SparseMatrix<double>(size,cp,ri,v);
}

// Compare a computed SymFact against the dense oracle. A fundamental supernode's columns
// share one sparsity pattern, so the supernode's index set is exactly the pattern of its
// lowest column: its front indices are the supernode's own columns, its update indices
// the rest. This holds in both regimes, a trivial supernode being the one-column case.
static bool matchesOracle(const SparseMatrix<double>& A, const Permutation& p,
                          const SymFact& s, std::size_t& fillOut){
    const auto pattern = OblioTest::denseFactorPattern(A, p);
    const std::size_t size = A.size();
    if(s.size()!=size) return false;

    // Gather each supernode's columns from the map (increasing, so front() is lowest).
    std::vector<std::vector<std::size_t>> cols(s.supSize());
    for(std::size_t lc=0; lc<size; ++lc)
        cols[static_cast<std::size_t>(s.idxToSupIdx()[lc])].push_back(lc);

    for(std::size_t k=0; k<s.supSize(); ++k){
        if(cols[k].empty()) return false;
        const std::size_t lowest = cols[k].front();
        std::vector<std::int32_t> got(s.idx().begin()+static_cast<std::ptrdiff_t>(s.idxPtr()[k]),
                                      s.idx().begin()+static_cast<std::ptrdiff_t>(s.idxPtr()[k+1]));
        if(got != pattern[lowest]) return false;                 // index set, sorted
        if(s.frontSize()[k] != cols[k].size()) return false;     // front == own columns
        if(s.frontSize()[k]+s.updateSize()[k] != pattern[lowest].size()) return false;
        // The front indices are the supernode's own columns, in order.
        for(std::size_t t=0; t<cols[k].size(); ++t)
            if(got[t] != static_cast<std::int32_t>(cols[k][t])) return false;
    }
    // Fill: entries of L absent from the permuted A.
    fillOut = 0;
    const auto& cp=A.colPtr(); const auto& ri=A.rowIdx(); const auto& o2n=p.oldToNew();
    std::vector<std::vector<bool>> inA(size, std::vector<bool>(size,false));
    for(std::size_t ac=0;ac<size;++ac)
        for(std::size_t t=cp[ac];t<cp[ac+1];++t)
            inA[static_cast<std::size_t>(o2n[static_cast<std::size_t>(ri[t])])]
               [static_cast<std::size_t>(o2n[ac])]=true;
    for(std::size_t j=0;j<size;++j)
        for(std::int32_t i : pattern[j])
            if(!inA[static_cast<std::size_t>(i)][j]) ++fillOut;
    return true;
}

// Reconstruct the pattern of each factor column from a SymFact, in either regime. Column
// lc sits at some position t among its supernode's front indices; since a supernode's
// index set is sorted and its front indices come first, the pattern of column lc is the
// index set from position t onward. This is how numeric factorization will read a column,
// and it must agree with the dense oracle whether or not supernodes were compressed.
static bool columnPatternsMatch(const SparseMatrix<double>& A, const Permutation& p,
                                const SymFact& s){
    const auto pattern = OblioTest::denseFactorPattern(A,p);
    std::vector<std::size_t> posInFront(s.size(), 0);
    std::vector<std::size_t> cursor(s.supSize(), 0);
    for(std::size_t lc=0; lc<s.size(); ++lc){
        const std::size_t k = static_cast<std::size_t>(s.idxToSupIdx()[lc]);
        posInFront[lc] = cursor[k]++;
    }
    for(std::size_t lc=0; lc<s.size(); ++lc){
        const std::size_t k = static_cast<std::size_t>(s.idxToSupIdx()[lc]);
        const std::size_t from = s.idxPtr()[k] + posInFront[lc];
        std::vector<std::int32_t> got(s.idx().begin()+static_cast<std::ptrdiff_t>(from),
                                      s.idx().begin()+static_cast<std::ptrdiff_t>(s.idxPtr()[k+1]));
        if(got != pattern[lc]) return false;
    }
    return true;
}

static bool runCase(const SparseMatrix<double>& A, const Permutation& p, std::size_t& fillOut,
                    Supernodes mode = Supernodes::Fundamental){
    ElmForestEngine feng(mode); SymFactEngine seng;
    ElmForest f; SymFact s;
    if(!feng.compute(A,p,f)) return false;
    if(!seng.compute(A,p,f,s)) return false;
    // The offsets must bracket the index array exactly.
    if(s.idxPtr().size()!=s.supSize()+1) return false;
    if(s.idxPtr()[s.supSize()]!=s.numIdx()) return false;
    if(s.idx().size()!=s.numIdx()) return false;
    if(!columnPatternsMatch(A,p,s)) return false;
    return matchesOracle(A,p,s,fillOut);
}

int main(){
    std::size_t fill=0;

    // No fill: a path. Every column touches only the next one.
    { auto A=fromEdges(5,{{0,1},{1,2},{2,3},{3,4}}); reqSym(A,"path n=5            : symmetric");
      Permutation p(5);
      ck(runCase(A,p,fill) && fill==0, "path n=5 natural    : no fill"); }

    // Smallest genuine fill: 0 adjacent to 1 and 2, which are not adjacent.
    // Eliminating 0 fills (2,1).
    { auto A=fromEdges(3,{{0,1},{0,2}}); reqSym(A,"vee n=3             : symmetric");
      Permutation p(3);
      ck(runCase(A,p,fill) && fill==1, "vee n=3 natural     : 1 fill"); }

    // A cycle fills: the path 0-1-2-3-4 closed by the edge 0-4.
    { auto A=fromEdges(5,{{0,1},{1,2},{2,3},{3,4},{0,4}}); reqSym(A,"cycle n=5           : symmetric");
      Permutation p(5);
      ck(runCase(A,p,fill) && fill==2, "cycle n=5 natural   : 2 fill"); }

    // 3x3 grid, the classic fill generator.
    { std::vector<std::pair<std::size_t,std::size_t>> e;
      for(std::size_t r=0;r<3;++r) for(std::size_t c=0;c<3;++c){ std::size_t u=r*3+c;
          if(c+1<3){ e.push_back({u,u+1}); }
          if(r+1<3){ e.push_back({u,u+3}); } }
      auto A=fromEdges(9,e); reqSym(A,"grid 3x3            : symmetric");
      Permutation p(9);
      ck(runCase(A,p,fill) && fill==8, "grid 3x3 natural    : 8 fill"); }

    // Two disconnected blocks: two trees, so the union must not leak across them.
    { auto A=fromEdges(6,{{0,1},{1,2},{3,4},{4,5}}); reqSym(A,"two blocks          : symmetric");
      Permutation p(6);
      ck(runCase(A,p,fill) && fill==0, "two blocks natural  : 2 trees, no fill"); }

    // A chosen permutation, checked against fill counts computable by hand. The arrow
    // has a hub adjacent to all five leaves. Eliminate the hub first and its five
    // leaves become mutually adjacent: 5*4/2 = 10 fill entries. Eliminate it last and
    // each leaf is a leaf of the tree, touching only the hub: no fill at all. This
    // pins the permutation maps to arithmetic we can check on paper, unlike the AMD
    // cases below, which check only against the oracle.
    { auto A=fromEdges(6,{{0,1},{0,2},{0,3},{0,4},{0,5}});   // hub is column 0
      reqSym(A,"arrow n=6           : symmetric");
      Permutation pNat(6); std::size_t fillNat=0;
      ck(runCase(A,pNat,fillNat) && fillNat==10, "arrow n=6 hub first : 10 fill");
      // Move the hub last: old 0 -> new 5, old i -> new i-1.
      Permutation pHubLast;
      bool set = pHubLast.setOldToNew({5,0,1,2,3,4});
      std::size_t fillLast=0;
      ck(set && pHubLast.validate(), "arrow n=6 setOldToNew: accepted, inverse built");
      ck(runCase(A,pHubLast,fillLast) && fillLast==0, "arrow n=6 hub last  : 0 fill");
      // The same ordering given the other way round must produce the same permutation.
      Permutation pOther;
      bool set2 = pOther.setNewToOld({1,2,3,4,5,0});
      ck(set2 && pOther.oldToNew()==pHubLast.oldToNew()
              && pOther.newToOld()==pHubLast.newToOld(),
         "arrow n=6 setNewToOld: same permutation");
      // A map that is not a bijection must be refused, leaving the object untouched.
      Permutation pBad(3);
      ck(!pBad.setOldToNew({0,1,1}) && !pBad.setOldToNew({0,1,3})
         && !pBad.setOldToNew({0,1,-1}) && pBad.validate()
         && pBad.oldToNew()==std::vector<std::int32_t>({0,1,2}),
         "setOldToNew rejects  : duplicate, out of range, negative"); }

    // Supernodal cases. Until compression these could not exist: every supernode had one
    // front column, so the union's multi-front-column path never ran. Now it does.
    { // A dense 4x4: all four columns share a pattern, so one supernode spans them all.
      auto A=fromEdges(4,{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});
      reqSym(A,"dense n=4           : symmetric");
      ElmForestEngine feng; SymFactEngine seng;
      Permutation p(4); ElmForest f; SymFact s;
      bool ok = feng.compute(A,p,f) && seng.compute(A,p,f,s);
      std::vector<std::int32_t> idx(s.idx().begin(), s.idx().end());
      ck(ok && s.supSize()==1 && s.frontSize()[0]==4 && s.updateSize()[0]==0
         && s.numIdx()==4 && idx==std::vector<std::int32_t>({0,1,2,3}),
         "dense n=4 natural   : one supernode, 4 front columns"); }

    { // A single column hanging off a dense block: supernodes {0} and {1,2,3,4}.
      auto A=fromEdges(5,{{0,1},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}});
      reqSym(A,"tail n=5            : symmetric");
      ElmForestEngine feng; SymFactEngine seng;
      Permutation p(5); ElmForest f; SymFact s;
      bool ok = feng.compute(A,p,f) && seng.compute(A,p,f,s);
      ck(ok && s.supSize()==2 && s.frontSize()[0]==1 && s.updateSize()[0]==1
         && s.frontSize()[1]==4 && s.updateSize()[1]==0,
         "tail n=5 natural    : supernodes {0} and {1,2,3,4}"); }

    { auto A=fromEdges(5,{{0,1},{1,2},{2,3},{3,4}}); Permutation p(5);
      std::size_t f2=0;
      ck(runCase(A,p,f2), "path n=5 natural    : supernodal index sets match oracle"); }

    // The nodal regime: compression turned off. Every supernode is one column, so the
    // forest and factor degenerate to the classic per-column form. The factor itself must
    // be identical to the compressed one; only the grouping into supernodes differs.
    { auto A=fromEdges(4,{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});   // dense: 1 supernode when compressed
      ElmForestEngine feng(Supernodes::Nodal); SymFactEngine seng;
      Permutation p(4); ElmForest f; SymFact s;
      bool ok = feng.compute(A,p,f) && seng.compute(A,p,f,s);
      bool allOne = true;
      for(std::size_t k=0;k<s.supSize();++k) if(s.frontSize()[k]!=1) allOne=false;
      ck(ok && s.supSize()==4 && allOne && s.idxToSupIdx()==std::vector<std::int32_t>({0,1,2,3})
         && columnPatternsMatch(A,p,s),
         "dense n=4 nodal     : no merging, identity map, same factor"); }

    // Same matrices, both regimes: the supernode counts differ, the factor does not.
    { std::mt19937 rng(20260711);
      int bad=0; std::size_t nodalSnodes=0, fundSnodes=0;
      for(int trial=0; trial<200; ++trial){
          std::size_t size = 3 + rng()%10;
          std::vector<std::pair<std::size_t,std::size_t>> e;
          for(std::size_t i=1;i<size;++i) e.push_back({i-1,i});
          for(std::size_t i=0;i<size;++i)
              for(std::size_t j=i+2;j<size;++j)
                  if(rng()%100 < 25) e.push_back({i,j});
          auto A=fromEdges(size,e);
          Permutation p(size);
          std::size_t f1=0, f2=0;
          if(!runCase(A,p,f1,Supernodes::Nodal)) ++bad;
          if(!runCase(A,p,f2,Supernodes::Fundamental)) ++bad;
          if(f1!=f2) ++bad;                          // same fill either way

          ElmForestEngine nodal(Supernodes::Nodal), fund(Supernodes::Fundamental);
          ElmForest fn, ff; nodal.compute(A,p,fn); fund.compute(A,p,ff);
          if(fn.supSize()!=size) ++bad;              // nodal: one supernode per column
          if(ff.supSize()>fn.supSize()) ++bad;       // compression never grows the count
          nodalSnodes += fn.supSize(); fundSnodes += ff.supSize();
      }
      ck(bad==0, "random x200 both    : nodal and fundamental give the same factor");
      ck(fundSnodes<nodalSnodes,
         "random x200 both    : compression merges ("+std::to_string(nodalSnodes)+" -> "
         +std::to_string(fundSnodes)+" supernodes)"); }

    // Random connected patterns against the oracle, with AMD as well as natural order.
    { std::mt19937 rng(20260711);
      OrderEngine ord(OrderMethod::AMD);
      int bad=0; std::size_t totalFill=0;
      for(int trial=0; trial<200; ++trial){
          std::size_t size = 3 + rng()%10;
          std::vector<std::pair<std::size_t,std::size_t>> e;
          for(std::size_t i=1;i<size;++i) e.push_back({i-1,i});        // spanning path
          for(std::size_t i=0;i<size;++i)
              for(std::size_t j=i+2;j<size;++j)
                  if(rng()%100 < 25) e.push_back({i,j});
          auto A=fromEdges(size,e);
          std::size_t f1=0, f2=0;
          Permutation pNat(size);
          if(!runCase(A,pNat,f1)) ++bad;
          Permutation pAmd; ord.compute(A,pAmd);
          if(!runCase(A,pAmd,f2)) ++bad;
          totalFill += f1 + f2;
      }
      ck(bad==0, "random x200         : natural and AMD match oracle");
      ck(totalFill>0, "random x200         : fill exercised ("+std::to_string(totalFill)+" entries)"); }

    std::cout<<"\nSymFact tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
