#include "oblio/Permutation.h"
#include "oblio/Types.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace Oblio;
static int pass=0, fail=0;
static void ck(bool ok,const std::string& n){ std::cout<<"  "<<(ok?"PASS  ":"FAIL  ")<<n<<"\n"; ok?++pass:++fail; }
using Map = std::vector<std::int32_t>;

// A permutation is well formed if the two maps are consistent inverses (validate),
// and, applied to indices, oldToNew and newToOld undo each other.
static bool wellFormed(const Permutation& p){
    if(!p.validate()) return false;
    for(std::size_t i=0;i<p.size();++i)
        if(p.newToOld()[static_cast<std::size_t>(p.oldToNew()[i])]!=static_cast<std::int32_t>(i))
            return false;
    return true; }

// The composition, computed the obvious way, independently of compose().
static Map composedBy(const Map& first, const Map& second){
    Map out(first.size());
    for(std::size_t i=0;i<first.size();++i)
        out[i]=second[static_cast<std::size_t>(first[i])];
    return out; }

int main(){
    // Setters: accept a bijection, build the inverse, refuse anything else.
    { Permutation p;
      ck(p.setOldToNew({1,2,0}) && p.oldToNew()==Map({1,2,0}) && p.newToOld()==Map({2,0,1})
         && wellFormed(p), "setOldToNew        : adopts map, builds inverse"); }
    { Permutation p;
      ck(p.setNewToOld({2,0,1}) && p.oldToNew()==Map({1,2,0}) && p.newToOld()==Map({2,0,1})
         && wellFormed(p), "setNewToOld        : same permutation, other direction"); }
    { Permutation p(3);
      ck(!p.setOldToNew({0,1,1}) && !p.setOldToNew({0,1,3}) && !p.setOldToNew({0,1,-1})
         && !p.setNewToOld({2,2,0}) && p.oldToNew()==Map({0,1,2}) && wellFormed(p),
         "setters reject     : duplicate, out of range, negative; state untouched"); }

    // Compose: this runs first, p second. Worked by hand: [1,2,0] then [2,0,1] is the
    // identity, because [2,0,1] is the inverse of [1,2,0].
    { Permutation a; a.setOldToNew({1,2,0});
      Permutation b; b.setOldToNew({2,0,1});
      ck(a.compose(b) && a.oldToNew()==Map({0,1,2}) && wellFormed(a),
         "compose            : with its own inverse gives the identity"); }

    // Order matters, and it is "this first, then p", not the reverse.
    { Permutation a; a.setOldToNew({1,0,2});     // swap 0 and 1
      Permutation b; b.setOldToNew({0,2,1});     // swap 1 and 2
      Permutation ab=a, ba=b;
      ck(ab.compose(b) && ab.oldToNew()==composedBy({1,0,2},{0,2,1}), "compose            : this first, p second");
      ck(ba.compose(a) && ba.oldToNew()==composedBy({0,2,1},{1,0,2}), "compose            : reversed gives the other product");
      ck(ab.oldToNew()!=ba.oldToNew(), "compose            : not commutative (as expected)"); }

    // The identity is neutral on both sides.
    { Permutation a; a.setOldToNew({3,0,2,1});
      Permutation id(4);
      Permutation ai=a, ia=id;
      ck(ai.compose(id) && ai.oldToNew()==a.oldToNew() && wellFormed(ai),
         "compose            : identity on the right leaves it unchanged");
      ck(ia.compose(a) && ia.oldToNew()==a.oldToNew() && wellFormed(ia),
         "compose            : identity on the left leaves it unchanged"); }

    // Size mismatch is refused and changes nothing.
    { Permutation a; a.setOldToNew({1,0});
      Permutation b(3);
      ck(!a.compose(b) && a.oldToNew()==Map({1,0}) && wellFormed(a),
         "compose            : size mismatch refused, state untouched"); }

    // Random: compose must agree with the direct computation, stay well formed, and
    // composing with the inverse must return the identity.
    { std::mt19937 rng(20260711);
      int bad=0;
      for(int trial=0; trial<500; ++trial){
          const std::size_t size = 1 + rng()%12;
          Map ma(size), mb(size);
          std::iota(ma.begin(), ma.end(), 0);
          std::iota(mb.begin(), mb.end(), 0);
          std::shuffle(ma.begin(), ma.end(), rng);
          std::shuffle(mb.begin(), mb.end(), rng);

          Permutation a; if(!a.setOldToNew(ma)){ ++bad; continue; }
          Permutation b; if(!b.setOldToNew(mb)){ ++bad; continue; }
          Permutation ab=a;
          if(!ab.compose(b) || ab.oldToNew()!=composedBy(ma,mb) || !wellFormed(ab)) ++bad;

          // The inverse of a is a's newToOld, read as an oldToNew map.
          Permutation inv; if(!inv.setOldToNew(a.newToOld())){ ++bad; continue; }
          Permutation ident=a;
          Map want(size); std::iota(want.begin(), want.end(), 0);
          if(!ident.compose(inv) || ident.oldToNew()!=want || !wellFormed(ident)) ++bad;
      }
      ck(bad==0, "random x500        : matches direct composition; inverse gives identity"); }

    std::cout<<"\nPermutation tests: "<<pass<<"/"<<(pass+fail)<<" passed\n";
    return fail==0?0:1;
}
