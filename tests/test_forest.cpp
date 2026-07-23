#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "oblio/OrderEngine.h"
#include "oblio/ElmForest.h"
#include "oblio/ElmForestEngine.h"
#include "oblio/Types.h"
#include "test_util.h"
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
static const std::int32_t NP = NIL;
template<class Val> static void reqSym(const SparseMatrix<Val>& A, const std::string& lbl){
    ck(OblioTest::isStructurallySymmetric(A), lbl); }
static SparseMatrix<double> tridiagFull(std::size_t size){
    std::vector<std::size_t> cp(size+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
    for(std::size_t j=0;j<size;++j){ if(j>0){ri.push_back(static_cast<std::int32_t>(j-1));v.push_back(-1);} ri.push_back(static_cast<std::int32_t>(j));v.push_back(2);
        if(j+1<size){ri.push_back(static_cast<std::int32_t>(j+1));v.push_back(-1);} cp[j+1]=ri.size(); }
    return SparseMatrix<double>(size,cp,ri,v); }
static bool eq(const std::vector<std::int32_t>& a, std::initializer_list<std::int32_t> b){
    std::vector<std::int32_t> bb(b); if(a.size()!=bb.size())return false;
    for(std::size_t i=0;i<a.size();++i){ if(a[i]!=bb[i]) return false; }
    return true; }
static bool eq2(const std::vector<std::size_t>& a, std::initializer_list<std::size_t> b){
    std::vector<std::size_t> bb(b); if(a.size()!=bb.size())return false;
    for(std::size_t i=0;i<a.size();++i){ if(a[i]!=bb[i]) return false; }
    return true; }
static bool validEtree(const std::vector<std::int32_t>& par, std::size_t& roots){
    roots=0; for(std::size_t i=0;i<par.size();++i){
        if(par[i]==NP){++roots;continue;} std::size_t pv=static_cast<std::size_t>(par[i]); if(pv<=i) return false; if(pv>=par.size()) return false; }
    return true; }

// The doubly-linked child/sibling structure must agree with parent in every direction:
// each node's children are exactly the nodes naming it as parent, the forward and
// backward sibling chains are inverses, and the root list is itself a sibling chain.
static bool validLinks(const ElmForest& f){
    const std::size_t n=f.snodeSize();
    const auto& par=f.parent(); const auto& fc=f.firstChild(); const auto& lc=f.lastChild();
    const auto& ns=f.nextSibling(); const auto& ps=f.previousSibling();
    std::vector<std::size_t> seen(n,0);
    std::size_t roots=0;
    // walk the root chain
    for(std::int32_t r=f.firstRoot(); r!=NP; r=ns[static_cast<std::size_t>(r)]){
        if(par[static_cast<std::size_t>(r)]!=NP) return false;
        ++roots; if(roots>n) return false; }
    if(roots!=f.numTrees()) return false;
    if(roots>0 && par[static_cast<std::size_t>(f.lastRoot())]!=NP) return false;
    // walk every child list
    for(std::size_t s=0;s<n;++s){
        std::int32_t prev=NP; std::int32_t c=fc[s]; std::size_t count=0;
        for(; c!=NP; c=ns[static_cast<std::size_t>(c)]){
            const std::size_t uc=static_cast<std::size_t>(c);
            if(par[uc]!=static_cast<std::int32_t>(s)) return false;   // child names s as parent
            if(ps[uc]!=prev) return false;                            // backward chain inverts forward
            ++seen[uc]; prev=c; ++count; if(count>n) return false; }
        if(count==0){ if(fc[s]!=NP||lc[s]!=NP) return false; }
        else if(lc[s]!=prev) return false;                            // lastChild ends the chain
    }
    // every non-root appears in exactly one child list; roots in none
    for(std::size_t s=0;s<n;++s){
        const std::size_t want = (par[s]==NP) ? 0u : 1u;
        if(seen[s]!=want) return false; }
    return true; }

// Height, recomputed independently by climbing each node to its root.
static bool validHeight(const ElmForest& f){
    const std::size_t n=f.snodeSize(); const auto& par=f.parent();
    std::size_t h=0;
    for(std::size_t s=0;s<n;++s){ std::size_t d=1;
        for(std::int32_t a=par[s]; a!=NP; a=par[static_cast<std::size_t>(a)]) ++d;
        if(d>h) h=d; }
    return f.height()==h; }

// The supernodes, the map and the sizes, all against the dense oracle. A fundamental
// supernode's columns share one pattern, so the supernode's index set is exactly the
// pattern of its lowest column: hence frontSize is the column count and updateSize is
// what remains of that pattern below them.
static bool validSupernodes(const SparseMatrix<double>& A, const Permutation& p,
                            const ElmForest& f, Supernodes mode){
    const auto pattern = OblioTest::denseFactorPattern(A,p);

    if(mode==Supernodes::Fundamental){
        std::size_t wantSnodeSize=0;
        const auto wantMap = OblioTest::fundamentalSupernodes(A,p,wantSnodeSize);
        if(f.snodeSize()!=wantSnodeSize) return false;
        if(f.nodeToSnode()!=wantMap) return false;      // same merges, same labels
    } else {
        // Nodal: one supernode per column, identity map.
        if(f.snodeSize()!=f.size()) return false;
        for(std::size_t lc=0; lc<f.size(); ++lc)
            if(f.nodeToSnode(lc)!=static_cast<std::int32_t>(lc)) return false;
    }

    // Collect each supernode's columns from the map.
    std::vector<std::vector<std::size_t>> cols(f.snodeSize());
    for(std::size_t lc=0; lc<f.size(); ++lc)
        cols[static_cast<std::size_t>(f.nodeToSnode(lc))].push_back(lc);

    for(std::size_t s=0; s<f.snodeSize(); ++s){
        if(cols[s].empty()) return false;
        const std::size_t lowest = cols[s].front();     // collected in increasing order
        if(f.frontSize(s)!=cols[s].size()) return false;
        if(f.frontSize(s)+f.updateSize(s)!=pattern[lowest].size()) return false;
        // Every column of the supernode must appear in the lowest column's pattern.
        for(std::size_t c : cols[s])
            if(std::find(pattern[lowest].begin(), pattern[lowest].end(),
                         static_cast<std::int32_t>(c))==pattern[lowest].end()) return false;
    }
    // Supernode labels must stay topological: a parent's label exceeds its children's.
    for(std::size_t s=0; s<f.snodeSize(); ++s)
        if(f.parent()[s]!=NP && f.parent()[s]<=static_cast<std::int32_t>(s)) return false;
    return true; }

// Everything the forest asserts about itself, checked against independent recomputation.
static bool validForest(const SparseMatrix<double>& A, const Permutation& p, const ElmForest& f,
                        Supernodes mode){
    std::size_t roots=0;
    return validEtree(f.parent(),roots) && roots==f.numTrees()
        && validLinks(f) && validHeight(f) && validSupernodes(A,p,f,mode); }
int main(){
    ElmForestEngine eng;
    { auto A=tridiagFull(4); reqSym(A,"tridiag n=4         : symmetric");
      Permutation p(4); ElmForest f;
      // Columns 2 and 3 share a pattern ({2,3} and {3}), so they merge: 3 supernodes.
      ck(eng.compute(A,p,f) && eq(f.parent(),{1,2,NP}) && f.snodeSize()==3
         && eq(f.nodeToSnode(),{0,1,2,2}), "tridiag n=4 natural : path, 2+3 merged"); }
    { auto A=tridiagFull(6); reqSym(A,"tridiag n=6         : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      // Only the last two columns merge, as in the n=4 case: 5 supernodes.
      ck(eq(f.parent(),{1,2,3,4,NP}) && f.snodeSize()==5
         && eq(f.nodeToSnode(),{0,1,2,3,4,4}), "tridiag n=6 natural : path, 4+5 merged"); }
    { std::vector<std::size_t> cp={0,6,8,10,12,14,16};
      std::vector<std::int32_t> ri={0,1,2,3,4,5, 0,1, 0,2, 0,3, 0,4, 0,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"arrow 6x6           : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      // The hub is eliminated first, so every later column is dense: one supernode.
      ck(eq(f.parent(),{NP}) && f.snodeSize()==1 && f.frontSize(0)==6
         && f.updateSize(0)==0 && f.height()==1, "arrow 6x6 natural   : one supernode"); }
    { std::vector<std::size_t> cp={0,2,5,7,9,12,14};
      std::vector<std::int32_t> ri={0,1, 0,1,2, 1,2, 3,4, 3,4,5, 4,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      reqSym(A,"two blocks          : symmetric");
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      // Each block compresses {0},{1,2} and {3},{4,5}: 4 supernodes, still 2 trees.
      ck(inv && roots==2 && eq(f.parent(),{1,NP,3,NP}) && f.snodeSize()==4
         && eq(f.nodeToSnode(),{0,1,1,2,3,3}), "two blocks natural  : 2 trees, compressed"); }
    { auto A=tridiagFull(30); reqSym(A,"tridiag n=30        : symmetric");
      OrderEngine ord(OrderMethod::AMD); Permutation p; ord.compute(A,p);
      ElmForest f; bool ok=eng.compute(A,p,f); std::size_t roots=0; bool inv=validEtree(f.parent(),roots);
      ck(ok && inv && f.size()==30 && roots>=1, "tridiag n=30 + AMD  : valid etree"); }

    // The links, height and sizes, none of which the etree cases above touch.
    { auto A=tridiagFull(4); Permutation p(4); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.firstChild(),{NP,0,1}) && eq(f.lastChild(),{NP,0,1})
         && eq(f.nextSibling(),{NP,NP,NP}) && eq(f.previousSibling(),{NP,NP,NP})
         && f.numTrees()==1 && f.firstRoot()==2 && f.lastRoot()==2 && f.height()==3
         && eq2(f.frontSize(),{1,1,2}) && eq2(f.updateSize(),{1,1,0}),
         "tridiag n=4 natural : chain links, height 3"); }

    // A star: three leaves under one root, so the sibling chain is non-trivial.
    { std::vector<std::size_t> cp={0,2,4,6,10};
      std::vector<std::int32_t> ri={0,3, 1,3, 2,3, 0,1,2,3};
      std::vector<double> v={2,-1, 2,-1, 2,-1, -1,-1,-1,4};
      SparseMatrix<double> A(4,cp,ri,v); reqSym(A,"star n=4            : symmetric");
      Permutation p(4); ElmForest f; eng.compute(A,p,f);
      ck(eq(f.parent(),{3,3,3,NP}) && eq(f.firstChild(),{NP,NP,NP,0}) && eq(f.lastChild(),{NP,NP,NP,2})
         && eq(f.nextSibling(),{1,2,NP,NP}) && eq(f.previousSibling(),{NP,0,1,NP})
         && f.numTrees()==1 && f.height()==2,
         "star n=4 natural    : sibling chain, height 2"); }

    // Two disconnected blocks: two roots, so the root list is a real sibling chain.
    { std::vector<std::size_t> cp={0,2,5,7,9,12,14};
      std::vector<std::int32_t> ri={0,1, 0,1,2, 1,2, 3,4, 3,4,5, 4,5};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(6,cp,ri,v);
      Permutation p(6); ElmForest f; eng.compute(A,p,f);
      ck(f.numTrees()==2 && f.firstRoot()==1 && f.lastRoot()==3
         && f.nextSibling()[1]==3 && f.previousSibling()[3]==1 && f.height()==2,
         "two blocks natural  : root chain, height 2"); }

    // The full invariant set against independent recomputation, over random patterns,
    // natural and AMD ordered. Covers the links, the height, and the front/update sizes.
    { std::mt19937 rng(20260711);
      OrderEngine ord(OrderMethod::AMD);
      int bad=0;
      for(int trial=0; trial<200; ++trial){
          std::size_t size = 3 + rng()%10;
          std::vector<std::vector<bool>> M(size, std::vector<bool>(size,false));
          for(std::size_t i=0;i<size;++i) M[i][i]=true;
          for(std::size_t i=1;i<size;++i){ M[i][i-1]=M[i-1][i]=true; }
          for(std::size_t i=0;i<size;++i)
              for(std::size_t j=i+2;j<size;++j)
                  if(rng()%100<25){ M[i][j]=M[j][i]=true; }
          std::vector<std::size_t> cp(size+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
          for(std::size_t j=0;j<size;++j){
              for(std::size_t i=0;i<size;++i) if(M[i][j]){ ri.push_back(static_cast<std::int32_t>(i));
                  v.push_back(i==j?10.0:-1.0); }
              cp[j+1]=ri.size(); }
          SparseMatrix<double> A(size,cp,ri,v);
          ElmForestEngine nodal(Supernodes::Nodal);
          Permutation pNat(size); ElmForest fNat, fNodal;
          if(!eng.compute(A,pNat,fNat) || !validForest(A,pNat,fNat,Supernodes::Fundamental)) ++bad;
          if(!nodal.compute(A,pNat,fNodal) || !validForest(A,pNat,fNodal,Supernodes::Nodal)) ++bad;
          Permutation pAmd; ord.compute(A,pAmd); ElmForest fAmd;
          if(!eng.compute(A,pAmd,fAmd) || !validForest(A,pAmd,fAmd,Supernodes::Fundamental)) ++bad;
      }
      ck(bad==0, "random x200         : links, height, sizes valid, both regimes"); }


    // Amalgamation. Fundamental supernodes are free but not maximal: the only-child condition
    // refuses a parent with two children even when one shares its pattern exactly and could
    // be absorbed for nothing. Amalgamation at threshold zero drops that condition.
    { // The smallest case: a star, 0-2 and 1-2. Column 2 has two children, so fundamental
      // merges nothing. Both children are free to merge, but only one can be: absorbing the
      // first widens the front, which prices out the second.
      std::vector<std::size_t> cp={0,2,4,7};
      std::vector<std::int32_t> ri={0,2, 1,2, 0,1,2};
      std::vector<double> v(ri.size(),1.0); SparseMatrix<double> A(3,cp,ri,v);
      Permutation p(3);
      ElmForest f0, f1; ElmForestEngine e1; e1.setThreshold(0);
      eng.compute(A,p,f0); e1.compute(A,p,f1);
      ck(f0.snodeSize()==3 && f1.snodeSize()==2 && eq2(f1.frontSize(),{1,2}),
         "star n=3 amalg(0)   : 3 supernodes -> 2, no fill"); }

    { // The grid of the notes. Fundamental gives 7 supernodes (the separator {7,8,9}).
      // At threshold zero, column 3 joins it for free: update(3) == indexSet(separator).
      std::vector<std::pair<std::size_t,std::size_t>> e =
          {{0,1},{0,6},{1,2},{1,7},{2,8},{3,4},{3,6},{4,5},{4,7},{5,8},{6,7},{7,8}};
      std::vector<std::vector<bool>> M(9, std::vector<bool>(9,false));
      for(std::size_t i=0;i<9;++i) M[i][i]=true;
      for(const auto& x : e){ M[x.first][x.second]=M[x.second][x.first]=true; }
      std::vector<std::size_t> cp(10,0); std::vector<std::int32_t> ri; std::vector<double> v;
      for(std::size_t j=0;j<9;++j){
          for(std::size_t i=0;i<9;++i) if(M[i][j]){ ri.push_back(static_cast<std::int32_t>(i)); v.push_back(i==j?10.0:-1.0); }
          cp[j+1]=ri.size(); }
      SparseMatrix<double> A(9,cp,ri,v); Permutation p(9);
      ElmForest f0, f1; ElmForestEngine e1; e1.setThreshold(0);
      eng.compute(A,p,f0); e1.compute(A,p,f1);
      ck(f0.snodeSize()==7 && f1.snodeSize()==6 && f1.frontSize(5)==4,
         "grid amalg(0)       : 7 supernodes -> 6, 4-column front"); }

    // The invariants that survive tie-breaking. Amalgamation is greedy and not canonical, so
    // the partition it produces depends on how ties are broken; these hold regardless.
    { std::mt19937 rng(20260712);
      OrderEngine ord(OrderMethod::AMD);
      int bad=0; std::size_t fund=0, at8=0;
      for(int trial=0; trial<200; ++trial){
          std::size_t size = 3 + rng()%10;
          std::vector<std::vector<bool>> M(size, std::vector<bool>(size,false));
          for(std::size_t i=0;i<size;++i) M[i][i]=true;
          for(std::size_t i=1;i<size;++i){ M[i][i-1]=M[i-1][i]=true; }
          for(std::size_t i=0;i<size;++i)
              for(std::size_t j=i+2;j<size;++j)
                  if(rng()%100<25){ M[i][j]=M[j][i]=true; }
          std::vector<std::size_t> cp(size+1,0); std::vector<std::int32_t> ri; std::vector<double> v;
          for(std::size_t j=0;j<size;++j){
              for(std::size_t i=0;i<size;++i) if(M[i][j]){ ri.push_back(static_cast<std::int32_t>(i)); v.push_back(i==j?10.0:-1.0); }
              cp[j+1]=ri.size(); }
          SparseMatrix<double> A(size,cp,ri,v);
          Permutation pNat(size), pAmd; ord.compute(A,pAmd);

          const Permutation* perms[2] = {&pNat, &pAmd};
          for(int which=0; which<2; ++which){
              const Permutation& p = *perms[which];
              ElmForest f0, f1, f2;
              ElmForestEngine a0;
              ElmForestEngine a1; a1.setThreshold(0);
              ElmForestEngine a2; a2.setThreshold(8);
              if(!a0.compute(A,p,f0) || !a1.compute(A,p,f1) || !a2.compute(A,p,f2)) { ++bad; continue; }

              // The links, height and labels must stay valid after amalgamation.
              if(!validLinks(f1) || !validHeight(f1)) ++bad;
              if(!validLinks(f2) || !validHeight(f2)) ++bad;
              std::size_t roots=0;
              if(!validEtree(f1.parent(),roots) || roots!=f1.numTrees()) ++bad;

              // Amalgamation never increases the supernode count, and a larger budget never
              // increases it either.
              if(f1.snodeSize() > f0.snodeSize()) ++bad;
              if(f2.snodeSize() > f1.snodeSize()) ++bad;

              // Every column belongs to exactly one supernode, and front sizes sum to size.
              std::size_t totalFront=0;
              for(std::size_t k=0;k<f2.snodeSize();++k) totalFront += f2.frontSize(k);
              if(totalFront != size) ++bad;

              // Supernode labels stay topological.
              for(std::size_t k=0;k<f2.snodeSize();++k)
                  if(f2.parent()[k]!=NP && f2.parent()[k] <= static_cast<std::int32_t>(k)) ++bad;

              if(which==0){ fund += f0.snodeSize(); at8 += f2.snodeSize(); }
          }
      }
      ck(bad==0, "random x200 amalg   : links, height, labels valid; count never grows");
      ck(at8 < fund, "random x200 amalg   : threshold 8 merges ("+std::to_string(fund)
                     +" -> "+std::to_string(at8)+" supernodes)"); }

    // exactPatterns(): do a supernode's columns share one sparsity pattern? It decides whether
    // symbolic factorization may read one front column per supernode or must read them all, so
    // it has to be exactly right, not merely conservative.
    { std::vector<std::vector<bool>> M(6, std::vector<bool>(6,false));
      for(std::size_t i=0;i<6;++i) M[i][i]=true;
      auto edge=[&](std::size_t a,std::size_t b){ M[a][b]=M[b][a]=true; };
      edge(0,1); edge(1,2); edge(2,3); edge(3,4); edge(4,5); edge(0,3); edge(1,4);
      std::vector<std::size_t> cp(7,0); std::vector<std::int32_t> ri; std::vector<double> v;
      for(std::size_t j=0;j<6;++j){
          for(std::size_t i=0;i<6;++i) if(M[i][j]){ ri.push_back(static_cast<std::int32_t>(i)); v.push_back(i==j?10.0:-1.0); }
          cp[j+1]=ri.size(); }
      SparseMatrix<double> A(6,cp,ri,v);
      Permutation p(6);
      ElmForest fNodal, fFund, fAmal0, fAmal8;
      ElmForestEngine eNodal(Supernodes::Nodal);
      ElmForestEngine eFund;
      ElmForestEngine eAmal0; eAmal0.setThreshold(0);
      ElmForestEngine eAmal8; eAmal8.setThreshold(64);
      eNodal.compute(A,p,fNodal); eFund.compute(A,p,fFund);
      eAmal0.compute(A,p,fAmal0); eAmal8.compute(A,p,fAmal8);
      ck(fNodal.exactPatterns(), "exactPatterns       : nodal, trivially true");
      ck(fFund.exactPatterns(),  "exactPatterns       : fundamental, true by definition");
      ck(fAmal0.exactPatterns(), "exactPatterns       : amalgamation at 0 buys nothing, still true");
      ck(!fAmal8.exactPatterns() || fAmal8.snodeSize()==fFund.snodeSize(),
         "exactPatterns       : false once a zero is stored"); }

    { // ---------------------------------------------------------------------------------------
      // sortForOptimalMultifrontal: children reordered by decreasing (maximumStorage - block).
      //
      // The forest is built twice from the same matrix, once with the option off and once on, and
      // the two are compared. Off, nothing may move. On, the invariants must hold and every child
      // list must come out in non-increasing key order, which is the property the sort exists to
      // establish. A grid is used because its forest genuinely branches; a chain has no sibling
      // choice to make and would assert nothing.
      const std::size_t g = 8, n = g * g;
      std::vector<std::size_t>  cp(n + 1, 0);
      std::vector<std::int32_t> ri;
      std::vector<double>       v;
      for (std::size_t j = 0; j < n; ++j) {
          const std::size_t r = j / g, c = j % g;
          std::vector<std::int32_t> col;
          col.push_back(static_cast<std::int32_t>(j));
          if (r) col.push_back(static_cast<std::int32_t>(j - g));
          if (c) col.push_back(static_cast<std::int32_t>(j - 1));
          if (r + 1 < g) col.push_back(static_cast<std::int32_t>(j + g));
          if (c + 1 < g) col.push_back(static_cast<std::int32_t>(j + 1));
          std::sort(col.begin(), col.end());
          for (std::int32_t i : col) {
              ri.push_back(i);
              v.push_back(i == static_cast<std::int32_t>(j) ? 8.0 : -1.0);
          }
          cp[j + 1] = ri.size();
      }
      SparseMatrix<double> A(n, cp, ri, v);

      Permutation p;
      OrderEngine oe; oe.setMethod(OrderMethod::AMD);
      oe.compute(A, p);

      ElmForest fOff, fOn;
      ElmForestEngine eOff;
      ElmForestEngine eOn;  eOn.setOptimizeMultifrontal(true);
      eOff.compute(A, p, fOff);
      eOn.compute(A, p, fOn);

      ck(!eOff.optimizeMultifrontal() && eOn.optimizeMultifrontal(),
         "optimizeMultifrontal: off by default, settable");

      // The pair rewrites links and relabels, so the supernodes themselves survive but their labels
      // move. What must hold is that the forest is the same forest: same count, same height, same
      // number of trees, and the same multiset of front and update sizes.
      std::vector<std::size_t> fsOff(fOff.frontSize()), fsOn(fOn.frontSize());
      std::vector<std::size_t> usOff(fOff.updateSize()), usOn(fOn.updateSize());
      std::sort(fsOff.begin(), fsOff.end()); std::sort(fsOn.begin(), fsOn.end());
      std::sort(usOff.begin(), usOff.end()); std::sort(usOn.begin(), usOn.end());
      ck(fOn.snodeSize() == fOff.snodeSize() && fOn.height() == fOff.height()
             && fOn.numTrees() == fOff.numTrees() && fsOn == fsOff && usOn == usOff,
         "optimizeMultifrontal: same forest, relabeled");

      ck(validLinks(fOn), "optimizeMultifrontal: links still consistent both ways");

      // Every child list in non-increasing key order, the key being 0.9's
      // maximumStorage(jj) - updateSize(jj)^2, recomputed here from the sorted forest.
      std::vector<std::size_t> ms(fOn.snodeSize(), 0);
      bool sorted = true;
      for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(fOn.snodeSize()); ++kk) {
          const std::size_t kkBlock = fOn.updateSize(kk) * fOn.updateSize(kk);
          std::size_t running = 0, best = 0, prevKey = 0;
          bool first = true;
          for (std::int32_t jj = fOn.firstChild()[kk]; jj != NIL; jj = fOn.nextSibling()[jj]) {
              const std::size_t key = ms[jj] - fOn.updateSize(jj) * fOn.updateSize(jj);
              if (!first && key > prevKey) sorted = false;
              prevKey = key; first = false;
              best = std::max(best, running + ms[jj]);
              running += fOn.updateSize(jj) * fOn.updateSize(jj);
          }
          ms[kk] = std::max(best, running + kkBlock);
      }
      ck(sorted, "optimizeMultifrontal: children in non-increasing key order");

      // The labels are now a postorder: every subtree holds a contiguous run ending at its root.
      // This is the property the multifrontal drivers need, since they walk labels, not links.
      // Asserted for the sorted forest and denied for nothing: the unsorted one may or may not be
      // a postorder already, which is exactly why the relabeling cannot be skipped.
      std::vector<std::size_t> sub(fOn.snodeSize(), 1);
      bool postorder = true;
      for (std::int32_t kk = 0; kk < static_cast<std::int32_t>(fOn.snodeSize()); ++kk) {
          std::size_t total = 1;
          for (std::int32_t jj = fOn.firstChild()[kk]; jj != NIL; jj = fOn.nextSibling()[jj]) {
              if (jj > kk) postorder = false;          // a child must be labeled below its parent
              total += sub[jj];
          }
          sub[kk] = total;
          // the subtree of kk is exactly [kk - sub(kk) + 1, kk]
          if (static_cast<std::size_t>(kk) + 1 < sub[kk]) postorder = false;
      }
      for (std::int32_t kk = 0; postorder && kk < static_cast<std::int32_t>(fOn.snodeSize()); ++kk) {
          const std::int32_t lo = kk - static_cast<std::int32_t>(sub[kk]) + 1;
          for (std::int32_t jj = lo; jj <= kk; ++jj) {
              // every label in the range must have kk as an ancestor
              std::int32_t aa = jj;
              while (aa != NIL && aa != kk) aa = fOn.parent()[aa];
              if (aa != kk) { postorder = false; break; }
          }
      }
      ck(postorder, "optimizeMultifrontal: labels are a postorder, subtrees contiguous");
    }

    std::cout<<"\nElmForest tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
