// amd3.cpp -- approximate minimum degree, iteration 3: aligned to the vendored AMD.
//
// **amd3 adds no mechanism. It is amd2 with the vendored routine's list order, and
// it exists to return AMD_2's pivot sequence.** mmd3 is its counterpart on the
// other branch and the digit means the same thing on both: 3 is the layer aligned
// to the vendored code. What amd2 has, this file has, and nothing else has been
// added.
//
// WHY A WHOLE LAYER FOR A TIE-BREAK. Minimum degree is a tie-break algorithm: at
// almost every iteration several vertices share the least degree and the winner is
// whichever the data structure hands over first. So two codes can agree on every
// rule and part company on the first tie, and from there they are ordering
// different graphs. While that is true, every comparison against the vendored
// routine measures two things at once, a difference of MECHANISM and a difference
// of ARBITRARY CHOICE, and a fill gap cannot be attributed to either. Alignment
// turns the comparison into an equality test.
//
// WHAT IS COMPARED, and it is not the permutation. AMD_2 ends with a postorder of
// the assembly tree, which Oblio does not want and this file does not do: a
// postorder relabels the output after elimination has finished and cannot touch
// which pivot was chosen when. So the acceptance test is the PIVOT SEQUENCE, which
// is the algorithm, and the permutations legitimately differ by that relabeling.
// mmd learned this the expensive way: with fill exact at every size and the
// permutation still diverging at pivot 700 of 1024, the pivot sequences were
// identical, 788 of 788.
//
// THE ORACLE is `amd_order` on a scratch copy of private/Amd.cpp with
// `Control[AMD_DENSE]` raised above sqrt(n), which drives the dense threshold into
// its `MIN(n, dense)` clamp so `deg > dense` is unreachable. Dense-row removal is
// the one thing amd2 lacks that AMD_2 does and that the Control array can switch
// off; at the default it never fires here anyway, the floor being MAX(16, dense)
// against a grid's maximum degree of four. `Control[AMD_AGGRESSIVE]` stays at its
// default, aggressive absorption being a mechanism amd2 already has.
//
// THE LEDGER
//
// Append only. A row is never edited once closed, so the sequence stays a record
// of what was wrong rather than a summary written afterwards. The nature column is
// the one to read first: CONVENTION means neither side is wrong and only the
// winner among equals moves, DEFECT means wrong on its own terms and must be fixed
// wherever the same code sits, COSMETIC means it cannot change the answer at all.
//
//     #  what diverged              where in ours      AMD_2                 nature
//     -  --------------------------  -----------------  --------------------  ----------
//     1  hash bucket walk            the hash pass      the head push at      convention
//                                                       its line 1940
//     2  reachable set layout        amd3Neighbors     construct new         convention
//        (cliques before explicit)                      element, knt1 loop
//     3  mass elimination ran        the eliminator,    scan 2, after the     convention
//        before absorption           now the driver     aggressive absorb
//     4  the vertex's own weight     the bound loop,    the fourth pass,      DEFECT
//        subtracted before the       now a fourth       `deg = Degree[i]
//        hash merge that grows it    pass               + degme - nvi`
//     5  the new clique appended     the eliminator     `Iw [p1] = me` and     convention
//        to I[u] instead of                             the two moves above
//        prepended, with a rotation                     it
//     6  the exact test walked the    the hash pass      `for (p = Pe[j] + 1`,  COST
//        new clique, which entry 5                      "skip the first
//        had just made a guaranteed                     element in the list
//        match at position zero                         (me)"
//     7  the stored clique degree    beginElimination,  `Degree [me] = degme`,  DEFECT
//        not rewritten after mass    PRODUCTION AmdFlat    written TWICE, at its
//        elimination trimmed the     alone              lines 1676 and 1940
//        clique
//     8  the hash key's incidence    the hash pass,     `hval += e` and        COST
//        half annihilated by         EVERY amd layer    `hval += j`, one sum
//        taking the modulus to be    and Amd2, Amd2B    with no stride, then
//        the stride                                     `hval % n`
//
// **Entry 8 is the largest of the eight and the only one that was never a
// divergence.** Entries 1 to 5 are lines of `AMD_2` we failed to reproduce; this is
// a defect of our own that it does not have. The key added the incidence half as
// `(c + 1) * (n + 1)`, so that a vertex and a clique of the same index could not
// cancel, and then reduced modulo the same number, which annihilates that half
// exactly: the bucket was a function of the ADJACENCY ALONE, and `A[u]` empties as
// the elimination proceeds, so the key carried less and less and the buckets grew.
// The invariant the two lines have to hold TOGETHER is that the modulus must not
// divide the stride, and `AMD_2` holds it by having no stride, letting a vertex and
// a clique collide on purpose because the hash is a filter and never the decision.
//
// It changed NO OUTPUT, which is why nothing here could see it: twins collide under
// any function of the pattern, so the merges found were the vendored routine's
// throughout. Against it on the same graphs, for the same merges, we tested 19.0
// pairs per pivot at 140 a side where it tests 0.333, and 155.3 at 26 cubed where it
// tests 0.484. Fixed on alpamayo, `AMD2` at 26 a side falls from 14.88 ms to 5.45
// and `AmdFlat` from 12.30 to 5.83, with the vendored routine and `AMD1` unmoved.
// `hash pairs tested` is the standing witness and should stay near one per merge.
// `docs/DESIGN_DECISIONS.md` (2026-08-09) has why five separate oracles were blind.
//
//
// **Entry 4's nature said `convention` here until 2026-08-09, and it is a DEFECT.**
// The README, `AMD3.md` and `docs/DESIGN_DECISIONS.md` have all said DEFECT since
// the day it was closed, and it is one by this ledger's own definition: it filed
// every supervariable one bucket too high per vertex a hash merge absorbed, which
// is wrong on its own terms with no appeal to AMD_2, and it was fixed in `amd2`,
// `Amd2` and `Amd2B` where it had been costing 3 to 9 percent of fill. The column
// is corrected rather than left, and the correction dated rather than made
// silently: append-only protects the record from being rewritten as a tidy summary
// afterwards, not from being wrong about itself. Worth knowing that the copy called
// authoritative had drifted from its mirrors in exactly the column the ledger tells
// a reader to look at first, and that nothing compares them.
//
// **Entry 7 is PRODUCTION'S ALONE, which is why its middle column says so.** These
// prototypes cannot have it: they obtain `|C[c] - C[p]|` by walking the members of
// C[c] and counting the live ones outside C[p], so they recompute it from the truth
// at every step and there is nothing to go stale. Production maintains a clique
// degree and obtains the same quantity by subtraction, which is amd2's pass 3 and
// is the encoding the prototypes deliberately do not carry.
//
// So the twin check could not have caught it, and not merely for want of a bigger
// case: a prototype written to read as the algorithm cannot model a hazard that
// lives in an optimization only production has, which makes it blind to exactly the
// class of defect that optimization introduces. That is the divergence
// `REPORT.md` parked as its fifth lead, and this is the first time it has cost
// anything.
//
// What the entry is. AMD_2 writes `Degree [me] = degme` twice, before scan 1 and
// again after supervariable detection, and the second write is the durable one:
// by then scan 2 has run `degme -= nvi` for every vertex mass elimination took, so
// what a later step reads as |C[me]| is the post-merge size. Production wrote it
// once, with the pre-merge value, so any pivot that mass-eliminated left a clique
// degree permanently too large by the merged weight and every later bound taken
// through that clique inherited it. It is half a mechanism again, as entry 6 was:
// ledger entry 3 moved mass elimination out of the eliminator and did not carry the
// second write that the move is the whole reason for. `Amd1` and `Amd2` are
// unaffected, mass-eliminating inside the eliminator, so their single write already
// sees a trimmed clique.
//
// It moves an ordering only when an inflated bound moves the head of the minimum
// bucket, which no 2D grid does at any size to 140 a side. A 3D grid at 16 finds
// it. That is what widened the acceptance test to four shapes.
//
//
// **Entry 6 is a fourth NATURE, and the ledger needs the word.** It is neither a
// convention nor a defect nor cosmetic: it changes no ordering, no fill and no
// permutation, and it changes the COST. Entry 5 put the new clique at the front of
// every I[u], which is right; it did not also skip that entry in the exact test,
// which Amd.cpp does in the same breath. So a guaranteed match sat at the head of a
// short-circuiting walk and every failing pair paid one extra iteration of the
// hottest line in the run. Measured at 2.08 incidence iterations per pair against
// amd2's 1.21, and 1.08 after, with Instruments putting that one line at 6.22 s of
// a 14.90 s run.
//
// The lesson is about porting rather than about AMD: **half a mechanism can be
// correct and still be wrong.** Entry 5 alone gives the right answer and pays for
// it, and nothing in the output could ever have shown that. Only the profile could.
//
// **Entry 2 is the deep one and entry 1 cannot be judged without it.** Entry 1
// alone took the examples from 2 of 7 to 2 of 7, closing graph1 and opening
// graph7; with entry 2 beneath it they went to 4 of 7. That is the same shape mmd
// had, where three reversed walks stalled at 4 of 7 until the element expansion
// was fixed under them, and the reason is the same: entry 2 fixes the content
// order of C[pivot], which is the order entry 1's buckets are built from.
//
// **Entry 3 closed two graphs at once**, which is what its shape predicted: graph4
// and graph6 each matched the oracle for their whole prefix and then took one
// extra pivot, and one cause was likelier than two. It was the same cause.
//
//     after entry 1        2 of 7      graph1 graph5
//     after entry 2        4 of 7      graph1 graph2 graph5 graph7
//     after entry 3        6 of 7      all but graph3
//     after entry 4        7 of 7      every example
//     after entry 5        7 of 7      and every grid tested
//
// nnz(L) is exact against the oracle throughout, at every entry and every size, so
// all five entries moved the tie-break and nothing else.
//
// ## Where it stands: ALIGNED
//
// **amd3 returns AMD_2's permutation exactly, up to the postorder.** Not merely
// its pivot sequence: the full expanded order agrees, member order within each
// supervariable included, on the seven examples and on every square grid tested
// from 3 a side to 40, `n = 1600`. The mechanism counters agree too, hash merges
// and aggressive absorptions alike, and `bound below exact` stays 0 throughout,
// which is the invariant a wrong bound would break.
//
// The strongest single check is against numbers this experiment did not produce:
// nnz(L) comes out 206332 at 100 a side and 474995 at 140, which is what
// `benchmarks/ordering/README.md` records for the vendored AMD, digit for digit.
//
// **How that is checked without the postorder, since AMD_2 always postorders.**
// Its `PIVOT` selection and its supervariable finalization both sit inside the
// main loop, and AMD_postorder runs after it, so the raw elimination order can be
// reconstructed upstream of the relabeling: track membership alongside, a hash
// merge moving j's members to i and a mass elimination moving i's to me, and emit
// each pivot's supervariable as its iteration closes. Concatenating those is what
// the vendored routine would return if it stopped at the end of its main loop, and
// it is what this layer returns. The scratch probe that does it is not a
// repository artifact; the README's method section describes rebuilding it.
//
// **The postorder itself is not a gap on our side.** It NEARLY cannot change the
// fill, and the qualifier was added on 2026-08-09. A postorder of the ELIMINATION
// tree cannot: every node is numbered after all its descendants either way, so the
// tree and its fill are unchanged. AMD postorders its ASSEMBLY tree, which its own
// header says need not be that tree, mass elimination under an approximate degree
// merging vertices that were never adjacent, so its relabeling is not guaranteed
// fill-neutral. The two orders give identical nnz(L) on every square grid and on
// cubic grids from 7 a side up, and differ by one to three entries at 4^3, 5^3 and
// 6^3. Nor is the postorder needed for CORRECTNESS by any traversal: what those require is
// a TOPOLOGICAL order, children before parents, and that holds whatever permutation
// arrives, since ElmForestEngine builds parent links from the permuted matrix and a
// parent's column is numbered above its child's. Left-looking pulls from
// descendants already factored, right-looking pushes to ancestors not yet reached,
// multifrontal consumes children at the parent, and a raw elimination order serves
// all three exactly as well.
//
// What a postorder buys is a smaller multifrontal STACK PEAK, and nothing else.
// The raw order is consistent with the assembly tree but its subtrees interleave,
// because minimum degree jumps around the graph, so a contribution block stays live
// across more of the numbering than it needs to. Left- and right-looking hold one
// update block at a time and never keep a standing Schur complement, so they have
// no peak for contiguity to shrink; and left-looking does not even gain locality
// from it, its relay hopping to whichever supernode owns a descendant's next
// remaining row rather than climbing the tree.
//
// AMD has to bake it into the permutation because the permutation is all it
// returns, having built the assembly tree and thrown it away. Oblio does not,
// because ElmForest holds the forest and reorders it directly, which
// `setOptimizeMultifrontal` does by Liu's rule on the supernodal tree with real
// front and update sizes.
//
// And that option is not something a multifrontal run can miss on the ordinary
// path: `DirectSolver` constructs the forest engine with
// `mTraversal == Traversal::Multifrontal`, so choosing multifrontal turns it on.
// `labelDepthFirst` then relabels the SUPERNODES into Liu's postorder, and the
// drivers loop over those labels, so the peak is set by that relabeling and not by
// the column permutation the ordering produced. AMD's postorder is redone in full
// and this layer's absence of one costs nothing. The engine default of off reaches
// only a caller wiring the engines by hand who intends multifrontal, which is that
// caller's obligation and is unaffected by which ordering ran.
//
// One second-order effect is left, and it is not special to this layer. A
// relabeling can change which columns are consecutive, and amalgamation is greedy
// with a tie-break on list position, so two permutations differing only by a
// postorder could in principle amalgamate differently. That would be a different
// supernode partition rather than a different peak, it applies to any pair of
// orderings related that way, and it is unmeasured.
//
// **What alignment bought, which is the reason for the whole layer.** Any future
// divergence from the vendored routine is now a named pivot in a small grid rather
// than a fill number somebody has to interpret. And the gap against amd2 can now be
// attributed: whatever it is, it is not a difference of arbitrary choice.
//
// **What it turned up about the layers below, recorded and NOT acted on.** Two of
// these entries look like defects in amd2 and in production Amd2 and Amd2B rather
// than conventions, entry 4 above all, which is mmd's entry 5 in a different array.
// Nothing outside this file has been touched: the mmd work fixed Mmd2 only after
// mmd3 was fully aligned, and the same order applies. `docs/TODO.md` carries the
// two items and two corrections this work owes the README.
//
//
// Run: ./amd3_cpp            every example
//      ./amd3_cpp 3          just the third
//      ./amd3_cpp grid 20
//
// What follows is amd2's own description, which stands unchanged: this file adds
// no mechanism to it.
//
// amd1 has the idea, the degree bound. This file adds the two mechanisms the
// vendored amd_2 carries beyond it, and nothing else. Both ride along with the
// bound rather than being about the degree at all, and both are cheap only because
// the bound's work has already been done. Section 5.13 of
// archive/sparse_factorization.md.
//
//   1. AGGRESSIVE ABSORPTION, which kills a clique the moment the bound work shows
//      it lies inside the new one. |C[c] - C[pivot]| has just been computed for
//      every touched clique; if it is zero, C[c] lies entirely inside C[pivot], so
//      that clique is dead. Ordinary absorption kills only the cliques the PIVOT
//      touched; this kills cliques that ANY reached vertex touched.
//
//   2. HASH SUPERVARIABLE DETECTION, which finds vertices indistinguishable from
//      EACH OTHER rather than from the pivot. Mass elimination merges into the
//      pivot, and two vertices can match each other while neither matches the
//      pivot. The hash is a filter and never the decision: a collision costs a
//      comparison, and no merge is ever missed, since equal sets give equal keys.
//
// **This file is production's Amd2, layer for layer**, and `make test` checks it by
// PERMUTATION rather than by fill. What amd3 adds beyond here is not ordering ideas
// at all: dense row detection, the pattern of A + A', the postorder and the
// Control/Info interface, none of which production has or wants.
//
// **One thing comes from amd3 and not from amd1**, and it is not optional. A hash
// merge leaves the merged vertex in place with weight zero rather than removing it
// from every list, so the walks must skip eliminated vertices. amd1 has no such
// vertices and its core takes no `eliminated` argument; this file's does. That is
// the prototype's version of what production calls live merges.
//
// The other fork from md5. Section 5.13 of archive/sparse_factorization.md.
//
// The other fork from md5. Section 5.13 of archive/sparse_factorization.md.
//
// md5 has the quotient graph, supervariables, maintained degrees and buckets, and
// returns exactly md1's ordering. What is left costing anything is the refresh
// itself, which for each reached vertex u unites the members of every clique in
// I[u] and counts the result. That union is the expensive object.
//
// MMD made the refresh RARE. AMD makes each one CHEAP, and the two are the same
// answer reached from opposite ends: do the expensive thing less.
//
// THE BOUND. Rather than uniting the cliques, sum their separate contributions:
//
//   degree(u) <= min( n - k - weight(u),                nothing exceeds what remains
//                     degree_old[u] + |C[p] - {u}|,     it can only grow by the new clique
//                     |A[u] - C[p]| + |C[p] - {u}|
//                                   + sum |C[c] - C[p]| )   over c in I[u] - {p}
//
// where p is the pivot, so C[p] is the new clique, k is the count of original vertices
// eliminated so far, and weight(u) is the size of u's supervariable. The third line
// OVERCOUNTS, because two cliques may overlap outside C[p] and the overlap is counted
// twice. So it is an upper bound, not the degree.
//
// WHY THAT IS FAST, which is the entire point and is easy to miss. The quantity
// |C[c] - C[p]| depends only on the clique c, not on the vertex u, so it is
// computed ONCE PER CLIQUE and then read by every vertex whose incidence list
// holds c. The exact degree costs, per vertex, a walk over the members of all its
// cliques. The bound costs, per vertex, one addition per clique. Both are counted
// below, and the gap widens with the size of the cliques, which is to say with the
// amount of fill, which is to say exactly where it matters.
//
// WHAT IS GIVEN UP, and it is a different kind of loss from mmd1's. Every layer up
// to here picks a true minimum-degree vertex and differs only in how it finds one
// or how ties fall. This one can pick the WRONG vertex outright, because an
// overcounted bound can hide the true minimum. It is the first layer whose
// heuristic changes rather than its implementation, and the first whose pivot is
// not guaranteed to be minimal at all.
//
// The trace prints the exact degree beside the bound, so the gap is visible at
// every iteration, and the closing lines count how often the bound was loose.
//
// THIS FILE IS THE IDEA ALONE. Aggressive absorption, hash supervariable
// detection, the two-pass update, dense row handling and the rest of amd_1 and
// amd_2 are amd3's business, exactly as mmd1 held only the batching and mmd2 took
// the rest of genmmd.
//
// Build:  g++ -std=c++17 -O3 amd3.cpp -o amd3_cpp  (or: make)
// Run:    ./amd3_cpp
//         ./amd3_cpp 3      just the third example

#include <cstdlib>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Plain vectors, UNSORTED, and a vector indexed by clique id, not std::map. A set
// costs O(log d) per membership test and per insertion; keeping a vector sorted
// costs a merge per union. Neither is needed: membership comes from a MARK array
// stamped with a tag, so "is v in the new clique" is one comparison, and every
// pass is linear in what it touches. That is what the vendored codes and Oblio's
// own SymFactorEngine do. See the README section on complexity.
//
// Types follow Oblio's rule: an INDEX names a vertex or a clique and is a
// std::int32_t, with NIL for "none"; a POSITION locates something inside a vector
// and is a std::size_t.
constexpr std::int32_t NIL = -1;

// The mark array is a set and the tag names it, so a tag must never repeat: a
// repeat makes a stale stamp read as a match, which is wrong silently. The tag
// only ever climbs, so the ceiling is where it has to be swept back. Half the
// positive range of std::int32_t, which is a pragmatic choice and not a derived
// one: nothing here stores anything but a tag, so the true ceiling is the type's
// own maximum, and the room left over is against a later layer wanting some of it.
constexpr std::int32_t TAG_CEILING = (1 << 30) - 1;

// Above this n, nothing is printed from inside the run: no initial state, no
// per-iteration trace. That output is for reading a small example by eye, and at any
// size worth calling large it is O(n) lines of O(n) each, so it is unreadable and slow
// to produce. What still prints at every size is the end of the run, the counters and
// the order, since each is O(1) lines and that is what the twin comparison comes down
// to. To watch a larger run, raise this.
constexpr std::size_t SHOW_THRESHOLD = 32;

using Graph = std::vector<std::vector<std::int32_t>>;

// C[c] holds the members of clique c, and cliqueLive[c] says whether c exists.
// A clique id is the pivot that created it, so the id space is the vertex space.
struct Cliques {
    std::vector<std::vector<std::int32_t>> members;
    std::vector<bool> live;
    std::uint32_t count = 0;

    explicit Cliques(std::size_t n) : members(n), live(n, false) {}
    const std::vector<std::int32_t>& at(std::int32_t c) const { return members[c]; }
    std::vector<std::int32_t>& operator[](std::int32_t c) { return members[c]; }
    void create(std::int32_t c, std::vector<std::int32_t> m) {
        if (!live[c]) ++count;
        live[c] = true;
        members[c] = std::move(m);
    }
    void erase(std::int32_t c) {
        if (live[c]) --count;
        live[c] = false;
        members[c].clear();
    }
    std::size_t size() const { return count; }
};

// I[u] cliques that contain u
// C[c] vertices that c contains

std::vector<std::int32_t> amd3Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u);

// The degree buckets, as the vendored codes hold them: one doubly linked list per
// degree, threaded through arrays of size n. Push, pop and splice are all O(1),
// which an ordered container cannot give. MMD spells these fwd/bwd and AMD
// Next/Last. The Python twin mirrors the same sequence with a list whose position
// 0 is the head, so both pick the same pivot.
struct Buckets {
    std::vector<std::int32_t> head;   // head[d], the first live vertex of degree d
    std::vector<std::int32_t> next;   // next[u], toward the tail
    std::vector<std::int32_t> prev;   // prev[u], toward the head
    std::vector<bool> filed;          // whether u is in a bucket at all

    explicit Buckets(std::size_t n)
        : head(n, NIL), next(n, NIL), prev(n, NIL), filed(n, false) {}

    void file(std::size_t d, std::int32_t u) {          // buckets[d].add(u), at the head
        next[u] = head[d];
        prev[u] = NIL;
        if (head[d] != NIL) prev[head[d]] = u;
        head[d] = u;
        filed[u] = true;
    }
    void unfile(std::size_t d, std::int32_t u) {        // buckets[d].discard(u)
        if (!filed[u]) return;                          // idempotent, as set.discard was
        if (prev[u] != NIL) next[prev[u]] = next[u];
        else head[d] = next[u];
        if (next[u] != NIL) prev[next[u]] = prev[u];
        next[u] = NIL;
        prev[u] = NIL;
        filed[u] = false;
    }
    bool empty(std::size_t d) const { return head[d] == NIL; }
};

// Print a quotient graph: adjacency, incidence, cliques, in the order the
// structure holds them.
void amd3Show(const Graph& A, const Graph& I, const Cliques& C,
              const std::vector<std::uint32_t>& degrees,
              const std::vector<std::uint32_t>& exact, const std::string& title = "",
              const std::vector<bool>* eliminated = nullptr) {
    const std::size_t n = A.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    std::vector<std::int32_t> aliveVertices;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        if (eliminated == nullptr || !(*eliminated)[u]) aliveVertices.push_back(u);
    std::size_t numAliveEdges = 0;
    for (std::int32_t u : aliveVertices) numAliveEdges += A[u].size();
    numAliveEdges /= 2;
    std::size_t numAliveIncidences = 0;
    for (std::int32_t u : aliveVertices) numAliveIncidences += I[u].size();
    std::size_t numAliveCliques = C.size();
    if (!title.empty()) std::cout << title << "\n";
    std::ostringstream aliveVerticesText;
    if (eliminated == nullptr) aliveVerticesText << n;
    else aliveVerticesText << aliveVertices.size() << " of " << n;
    std::cout << "num alive vertices = " << aliveVerticesText.str()
              << ", num alive edges = " << numAliveEdges
              << ", num alive cliques = " << numAliveCliques
              << ", storage = " << 2 * numAliveEdges << " + " << 2 * numAliveIncidences
              << " = " << 2 * (numAliveEdges + numAliveIncidences) << "\n";
    for (std::int32_t u : aliveVertices) {
        std::ostringstream adjacencyText;
        bool first = true;
        for (std::int32_t v : A[u]) {
            adjacencyText << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::ostringstream incidenceText;
        first = true;
        for (std::int32_t c : I[u]) {
            incidenceText << (first ? "" : " ") << "c" << c;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": {" << adjacencyText.str()
                  << "} {" << incidenceText.str() << "} bound " << degrees[u]
                  << " exact " << exact[u] << "\n";
    }
    for (std::int32_t c = 0; c < static_cast<std::int32_t>(n); ++c) {
        if (!C.live[c]) continue;
        std::ostringstream cliqueMembersText;
        bool first = true;
        for (std::int32_t u : C.at(c)) {
            cliqueMembersText << (first ? "" : " ") << std::setw(width) << u;
            first = false;
        }
        std::cout << "  c" << c << ": {" << cliqueMembersText.str() << "}\n";
    }
    std::cout << "\n";
}

// Print the state arrays: degrees, buckets, min degree, members, eliminated,
// and the order so far.
void amd3ShowState(const std::vector<std::uint32_t>& degrees, const Buckets& buckets,
                  std::uint32_t minDegree,
                  const std::vector<std::vector<std::int32_t>>& superMembers,
                  const std::vector<bool>& eliminated,
                  const std::vector<std::int32_t>& pivots, const std::string& title = "") {
    const std::size_t n = superMembers.size();
    int width = static_cast<int>(std::to_string(n > 0 ? n - 1 : 0).size());
    if (!title.empty()) std::cout << title << "\n";
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        std::string status;
        if (!eliminated[u]) status = "live";
        else if (!superMembers[u].empty()) status = "done";
        else status = "merged";
        std::ostringstream superMemberList;
        bool first = true;
        for (std::int32_t v : superMembers[u]) {
            superMemberList << (first ? "" : " ") << std::setw(width) << v;
            first = false;
        }
        std::cout << "  " << std::setw(width) << u << ": members ["
                  << superMemberList.str() << "] " << status << "\n";
    }
    std::ostringstream degreesText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        degreesText << (u == 0 ? "" : " ") << std::setw(width) << degrees[u];
    std::ostringstream superMembersText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) {
        superMembersText << (u == 0 ? "" : " ") << "[";
        bool firstMember = true;
        for (std::int32_t v : superMembers[u]) {
            superMembersText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        superMembersText << "]";
    }
    std::ostringstream eliminatedText;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        eliminatedText << (u == 0 ? "" : " ") << std::setw(width) << (eliminated[u] ? 1 : 0);
    std::ostringstream bucketsText;
    bool firstBucket = true;
    for (std::size_t d = 0; d < buckets.head.size(); ++d) {
        if (buckets.empty(d)) continue;
        bucketsText << (firstBucket ? "" : "  ") << d << ": [";
        bool firstMember = true;
        for (std::int32_t v = buckets.head[d]; v != NIL; v = buckets.next[v]) {
            bucketsText << (firstMember ? "" : " ") << v;
            firstMember = false;
        }
        bucketsText << "]";
        firstBucket = false;
    }
    std::cout << "  degrees: [" << degreesText.str() << "]\n";
    std::cout << "  buckets: " << (firstBucket ? "all empty" : bucketsText.str()) << "\n";
    std::cout << "  min degree: " << minDegree << "\n";
    std::cout << "  members: [" << superMembersText.str() << "]\n";
    std::cout << "  eliminated: [" << eliminatedText.str() << "]\n";
    std::cout << "  pivots: [";
    for (std::size_t k = 0; k < pivots.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << pivots[k];
    std::cout << "]\n";
    std::cout << "  order: [";
    bool firstOrder = true;
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) {
            std::cout << (firstOrder ? "" : ", ") << u;
            firstOrder = false;
        }
    std::cout << "]\n\n";
}

// Entries actually stored. Each edge costs two, one per endpoint in A. Each
// incidence costs two as well, the clique id in I and the member in C. Watch
// the total fall monotonically; the naive graph's only rises.
std::size_t amd3Storage(const Graph& A, const Graph& I, const Cliques& C) {
    std::size_t total = 0;
    for (const std::vector<std::int32_t>& adjacency : A) total += adjacency.size();
    for (const std::vector<std::int32_t>& incidence : I) total += incidence.size();
    for (std::size_t c = 0; c < C.members.size(); ++c)
        if (C.live[c]) total += C.members[c].size();
    return total;
}

// The neighbors of live vertex u: its explicit adjacency A[u] together with the
// members of every clique that contains u, minus u itself, which the cliques
// always carry. This is George and Liu's reachable set, and it is what the
// elimination graph would hold explicitly.
std::vector<std::int32_t> amd3Neighbors(const Graph& A, const Graph& I, const Cliques& C,
                                       const std::vector<bool>& eliminated,
                                       std::vector<std::int32_t>& mark, std::int32_t& tag,
                                       std::int32_t u) {
    // In set terms this is one line, and it is worth keeping in view because the
    // code below is that line with the set taken away:
    //
    //     reach(u) = ( A[u] | C[c] for every c in I[u] ) - {u}
    //
    // The mark array IS the set. mark[v] == tag is the membership test, one
    // comparison; mark[v] = tag is the insertion, one store. So the union costs one
    // pass per source rather than a hash per member, and nothing is allocated.
    //
    // One pass per source, with the mark array doing the deduplication, so the
    // cost is linear in what is touched. Nothing is sorted: the order is the order
    // the sources were walked in.
    //
    // LEDGER ENTRY 2, convention, and it is the DEEP one. The cliques are walked
    // BEFORE the explicit adjacency, which is the opposite of every other layer
    // here and is what AMD_2 does: its construct-new-element loop runs
    // `for (knt1 = 1; knt1 <= elenme + 1; knt1++)`, taking the ELEMENTS of me for
    // knt1 <= elenme and the supervariables only on the last pass, so Lme comes
    // out cliques-then-explicit where ours came out explicit-then-cliques. genmmd
    // is the other way round, variables first and elements last, which is why md1
    // through mmd3 are laid out as they are and why this is an amd convention
    // rather than a correction to the ladder.
    //
    // This fixes the CONTENT order of C[pivot], which is the order every later
    // walk reads: the hash buckets of entry 1, the pair comparisons, the degree
    // update. So entry 1 cannot be judged without it, exactly as mmd's three
    // reversed walks stalled at 4 of 7 until the element expansion was fixed
    // beneath them.
    ++tag;
    std::vector<std::int32_t> neighbors;
    mark[u] = tag;                          // never its own neighbor
    for (std::int32_t c : I[u])             // elements first, as AMD_2 lays out Lme
        for (std::int32_t v : C.at(c))
            if (mark[v] != tag && !eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    for (std::int32_t v : A[u])             // then the supervariables, its last pass
        if (mark[v] != tag && !eliminated[v]) { mark[v] = tag; neighbors.push_back(v); }
    return neighbors;
}

// The degree md5 would have computed: the union of A[u] with the members of every
// clique in I[u], counted in original vertices. Kept only so the trace can show
// the bound beside the truth and count how often the bound is loose.
//
// Set view: sum of |superMembers[v]| over v in reach(u). It is the union the bound
// exists to avoid, so this function is instrumentation and nothing more.
std::size_t amd3ExactDegree(const Graph& A, const Graph& I, const Cliques& C,
                            const std::vector<bool>& eliminated,
                            const std::vector<std::vector<std::int32_t>>& superMembers,
                            std::vector<std::int32_t>& mark, std::int32_t& tag,
                            std::int32_t u) {
    std::uint32_t degree = 0;
    for (std::int32_t v : amd3Neighbors(A, I, C, eliminated, mark, tag, u))
        degree += superMembers[v].size();
    return degree;
}

// Turn the pivot into a clique, then merge in every member it makes
// indistinguishable. Identical to md5Eliminate: this layer changes how a degree is
// estimated afterwards, not what an elimination does.
//
// Returns (neighbors, absorbedCliques, prunedEdges), as in md2. The last two are
// reported for display; only neighbors is used by the caller.//
// Set view of the whole function, in the order the code does it:
//
//     C[pivot] = reach(pivot)                    absorb into C[pivot]
//     C        = C - I[pivot]                    reclaim I[pivot]
//     for u in C[pivot]:
//         A[u] = A[u] - C[pivot] - {pivot}       prune
//         I[u] = ( I[u] - I[pivot] ) | {pivot}   absorb into C[pivot], reclaim I[pivot]
//
// The new clique is C[pivot] and gets no name of its own, so the first line reads
// as what an elimination IS: the pivot stops being a vertex with a reachable set
// and becomes a clique holding that same set. The last line is the first two
// written on the I side, since u is in C[c] exactly when c is in I[u].
//
// Three set differences, and not one of them builds a set. Each is a single stamp
// of the subtrahend followed by one compaction pass over the minuend, which turns
// |A[u]| * |C[pivot]| comparisons into |A[u]| + |C[pivot]|.
//
// LEDGER ENTRY 3, convention. Mass elimination used to be the last thing this
// function did, and it is now the driver's, run after aggressive absorption. The
// TEST is unchanged and is textually AMD_2's `Elen[i] == 1 && p3 == pn`, one
// element left and no supervariables; what moved is WHERE it runs. AMD_2 makes it
// in scan 2, over an element list from which aggressive absorption has already
// dropped every clique lying inside the new one, and says as much in its own
// comment: with aggressive absorption, `deg == 0` is identical to that structural
// test. Running it first, as this function did, asks the question of an I[u] that
// still holds cliques the absorption is about to remove, so the cheap test
// declines vertices AMD merges. graph6 is the case: at its third pivot the
// vendored routine mass-eliminates all three members of the new clique and
// finishes in three pivots, where we merged two, left vertex 2 behind and took a
// fourth.
//
// So this function now stops at the prune, and C[pivot] is reach(pivot) exactly,
// which is the md2 identity restored. The driver trims it after the merges.
std::tuple<std::vector<std::int32_t>, std::vector<std::int32_t>,
           std::vector<std::pair<std::int32_t, std::int32_t>>>
amd3Eliminate(Graph& A, Graph& I, Cliques& C, std::vector<bool>& eliminated,
             std::vector<std::int32_t>& mark, std::int32_t& tag, std::int32_t pivot) {
    const std::vector<std::int32_t> neighbors =
        amd3Neighbors(A, I, C, eliminated, mark, tag, pivot);
    const std::vector<std::int32_t> absorbedCliques = I[pivot];
    for (std::int32_t c : absorbedCliques)
        C.erase(c);
    C.create(pivot, neighbors);     // becomes the column pattern of the pivot

    // Stamp the new clique once, and the absorbed cliques once. Membership is then
    // a comparison, and both loops below are compactions in place. cliqueTag is the
    // set C[pivot] and absorbedTag is the set I[pivot], each built in one pass and
    // then queried for free.
    ++tag;
    const std::int32_t cliqueTag = tag;
    for (std::int32_t v : neighbors) mark[v] = cliqueTag;
    ++tag;
    const std::int32_t absorbedTag = tag;
    for (std::int32_t c : absorbedCliques) mark[c] = absorbedTag;

    std::vector<std::pair<std::int32_t, std::int32_t>> prunedEdges;
    std::vector<std::int32_t> kept;
    for (std::int32_t u : neighbors) {
        kept.clear();
        for (std::int32_t v : A[u]) {
            if (v == pivot) continue;            // the pivot is no longer a variable
            if (mark[v] == cliqueTag) {          // both ends inside the new clique
                if (u < v) prunedEdges.push_back({u, v});
                continue;                        // implicit now: drop the explicit copy
            }
            // A hash merge leaves the merged vertex in place with weight zero rather
            // than removing it from every list, so a dead vertex can still sit in
            // A[u]. Dropping it here matters beyond tidiness: the mass elimination
            // test below asks whether A[u] is EMPTY, and a stale entry makes it
            // answer no. Without this line the pivot count is higher and the
            // permutation the same, since a vertex missed by mass elimination is
            // simply eliminated on its own later.
            if (eliminated[v]) continue;
            kept.push_back(v);
        }
        A[u].swap(kept);                         // what survives is A[u] - C[pivot] - {pivot}

        kept.clear();                            // I[u] loses the absorbed cliques
        for (std::int32_t c : I[u])
            if (mark[c] != absorbedTag) kept.push_back(c);

        // LEDGER ENTRY 5, convention. The new clique goes to the FRONT of the
        // incidence list, not the back, and the displaced entries ROTATE rather
        // than shift. AMD_2 makes room for me in one move at the end of its scan 2:
        //
        //     Iw [pn] = Iw [p3] ;   /* move first supervariable to end of list */
        //     Iw [p3] = Iw [p1] ;   /* move first element to end of element part */
        //     Iw [p1] = me ;        /* add new element, me, to front of list */
        //
        // Its two lists live back to back in one run, elements then supervariables,
        // so inserting at the front of the elements would mean shifting everything
        // right. Instead it lifts the two entries sitting at the boundaries to the
        // two ends: the elements come out [me, e2, ..., ek, e1] and the
        // supervariables [j2, ..., jm, j1]. That is a data-structure trick with no
        // meaning of its own, and it decides the order every later walk reads.
        //
        // Both rotations are load-bearing here, because entry 2 made the reachable
        // set walk the cliques and then the adjacency, so both feed C[pivot]'s
        // content order, and that order decides the hash survivor of entry 1. The
        // empty cases fall out as they do there: with no surviving element the
        // first two moves are self-assignments and the list is just [pivot].
        if (!A[u].empty())
            std::rotate(A[u].begin(), A[u].begin() + 1, A[u].end());
        I[u].clear();
        I[u].push_back(pivot);
        if (!kept.empty()) {
            I[u].insert(I[u].end(), kept.begin() + 1, kept.end());
            I[u].push_back(kept.front());
        }
    }

    // Mass elimination is NOT here. It is the driver's, run after aggressive
    // absorption, which is where AMD_2's scan 2 makes the same test. See the
    // comment block above, ledger entry 3.

    A[pivot].clear();
    I[pivot].clear();
    eliminated[pivot] = true;
    return {neighbors, absorbedCliques, prunedEdges};
}

// Move u from the bucket for its old degree to the one for newDegree. Filing
// pushes at the head, which is the O(1) end of the list.
//
// Set view: buckets[d] is the set of live vertices whose current degree is d, and
// filing, unfiling and refiling are add, discard and move between two of them. A
// linked list gives all three in O(1) and gives the head in O(1) too, which is
// everything the picker asks of it. What it does not give is a minimum, which is
// why minDegree walks. A sorted container would hand over the minimum directly and
// charge a log on every file, and files outnumber picks.
void amd3Refile(Buckets& buckets, std::vector<std::uint32_t>& degrees,
               std::int32_t u, std::size_t newDegree) {
    buckets.unfile(degrees[u], u);
    degrees[u] = newDegree;
    buckets.file(newDegree, u);
}

// Same as md5, with the exact refresh replaced by the approximate bound.
// Everything else, the quotient graph, mass elimination, the buckets, is md5's.
std::vector<std::int32_t> amd3MinimumDegree(const Graph& G) {
    const std::size_t n = G.size();
    std::size_t nnzTrilA = 0;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) nnzTrilA += G[u].size();
    nnzTrilA = nnzTrilA / 2 + n;
    Graph A = G;                                  // explicit vertex neighbors
    Graph I(n);                                   // cliques that contain each vertex
    Cliques C(n);      // clique id -> member list
    // Twice n: vertices are stamped below n and cliques at c + n, so the exact
    // comparison in the hash pass cannot confuse a vertex with a clique.
    std::vector<std::int32_t> mark(2 * n, NIL);   // scratch for membership
    // The stride separating the two halves of `mark`, and the one place a COUNT becomes an offset
    // in the INDEX space: the hash comparison stamps a clique at c + cliqueStamp so that vertices
    // and cliques share one stamp. Named once rather than cast at each site that makes the
    // crossing. It exists only because there is no type for a count, which is one dimensional and
    // bounded by a SIDE where a position is bounded by an AREA; see REPORT.md.
    const std::int32_t cliqueStamp = static_cast<std::int32_t>(n);
    std::int32_t tag = 0;
    // Calls to the eliminate procedure, one per pivot. Not the count of vertices
    // removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    // an iteration batches several eliminations before one degree update pass. The three
    // counts coincide only where both of those are absent.
    std::uint32_t numEliminations = 0;
    // Summed over the eliminations, |C[p]| being the new clique AFTER the trim, so
    // in supernodal terms the update rather than the front. It is the raw reach of
    // the eliminations, undeduplicated: where a layer deduplicates, the degree
    // update count comes out below this, and the gap is what the batching saved.
    // In md2 it is nnz(L) - n, there being no mass elimination to shrink a clique.
    std::size_t numCliqueEntries = 0;
    // Passes of the outer loop, each one a batch of eliminations followed by one
    // degree update pass. Here the batch is always a single elimination, so this
    // equals numEliminations; from mmd1 up the two come apart.
    std::uint32_t numIterations = 0;
    std::vector<std::vector<std::int32_t>> superMembers(n);   // for the expansion
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        superMembers[u].push_back(u);
    std::vector<bool> eliminated(n, false);
    std::vector<std::int32_t> pivots;             // the order over supervariables
    std::uint32_t numEliminatedVertices = 0;                // a counter, not a scan of eliminated
    // Live ORIGINAL vertices, which is not n - numEliminatedVertices. numEliminatedVertices counts
    // what has left the SELECTION, and a hash merge folds v into a LIVE u, so v
    // stops being selectable while the vertices it stands for are still live inside
    // u. The first cap of the bound needs the second reading, so it gets its own
    // counter: only an elimination reduces it.
    //
    // amd1 has no such counter and does not need one, since it has no hash merges
    // and every increment of numEliminatedVertices really is an original leaving. This file
    // inherited that line along with the driver, which made the cap too tight and
    // drove the bound BELOW the true degree, 22 times on a 10 by 10 grid. The
    // vendored amd_2 is the oracle here: its nel is advanced only at the pivot and
    // at mass elimination, never at the hash merge, which moves the weight with
    // Nv[i] += Nv[j] and leaves nel alone.
    std::uint32_t numLive = n;
    std::size_t nnzL = 0;

    // The cache, and the count of degree computations, which is what this layer
    // exists to reduce. Built once, then touched only where it can be wrong.
    // The cache, as in md5, except that from the first elimination it holds a
    // BOUND rather than a degree. exact[] is carried alongside for the trace only.
    std::vector<std::uint32_t> degrees(n);          // a degree counts, so it measures
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u) degrees[u] = A[u].size();
    std::vector<std::uint32_t> exact = degrees;
    // Only the updates are counted. The total, including the initial pass over all
    // n vertices, is that plus n, so the report derives it. That first pass finds
    // |A[u]| with no clique yet formed, which is the bound formula on an empty
    // clique set and so is exact; the bound becomes a bound from the first
    // elimination on.
    std::size_t numBoundUpdates = 0;
    // Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    // is here as the witness that the guard is inert rather than as a statistic.
    std::size_t numTagSweeps = 0;
    std::size_t numMemberVisits = 0;              // what an exact refresh would cost
    std::size_t numCliqueReads = 0;               // what the bound costs instead
    // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    std::size_t numAbsorbed = 0;                 // cliques killed aggressively
    std::size_t numHashMerges = 0;               // pairs found by the hash
    // THE STANDING WITNESS FOR LEDGER ENTRY 8, and it is here for the reason `tag sweeps` and
    // `bound below exact` are: to make a claim checkable rather than to measure anything. The key
    // spread badly for months without a single output moving, because twins collide under any
    // function of the pattern, so the merges found were right the whole time and only the pairs
    // tested to find them were absurd. This is that number. It should read about half a pair per
    // pivot; before the fix it read 19 on a 140x140 grid and 155 at 26 cubed, against the vendored
    // routine's 0.33 and 0.48 for the same merges.
    std::size_t numHashPairs = 0;                // pairs the exact test was run on
    std::vector<std::vector<std::int32_t>> hashBucket(n);   // Amd.cpp's Head[hval]
    std::size_t numBoundChecks = 0;
    std::size_t numLooseBounds = 0;
    std::size_t numBoundsBelowExact = 0;   // an invariant, not a measurement

    // The buckets, and minDegree, a LOWER BOUND on the current minimum degree.
    // The search starts at minDegree rather than at 0, so it never looks at
    // buckets known to be empty. The bound may lag, and the walk corrects it; what
    // it must never do is overshoot, since a vertex below it would never be seen.
    //
    // n buckets is exactly right. A live vertex counts only live neighbors, so its
    // degree is at most n - 1, and the walk stops at the first non-empty bucket,
    // which exists while anything is live.
    Buckets buckets(n);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(n); ++u)
        buckets.file(degrees[u], u);
    std::uint32_t minDegree = n > 0 ? *std::min_element(degrees.begin(), degrees.end()) : 0;
    std::size_t numBucketProbes = 0;

    // |C[c] - C[pivot]| per clique, indexed by clique id, hoisted out of the loop it is used in.
    // Allocating and zeroing it per pivot reads better and is O(n) per iteration, hence O(n * n) over
    // the run in bookkeeping alone, independent of the graph, which would swamp the very cost the
    // bound exists to save. Only the entries an iteration writes are touched, and they are exactly the
    // ones it reads, so the iteration clears what it wrote rather than the array being rebuilt. The
    // Python twin has no such line, its outside being a dict over the cliques the iteration touched,
    // which is already the right shape.
    std::vector<std::uint32_t> outside(n, 0);

    // NOT PRODUCTION: display only. The trace is what makes these files teachable and
    // is the whole reason they exist; nothing downstream reads it.
    if (n <= SHOW_THRESHOLD) {
        amd3Show(A, I, C, degrees, exact,
                 "start: every edge explicit, no clique yet, degrees exact", &eliminated);
        amd3ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
    }
    int iteration = 0;
    while (numEliminatedVertices < n) {
        ++numIterations;
        while (buckets.empty(minDegree)) {      // walk up to the first live bucket
            ++minDegree;
            ++numBucketProbes;
        }
        ++numBucketProbes;
        std::int32_t pivot = buckets.head[minDegree];   // whatever was filed last
        // Sweep the tag back before it can wrap. THREE sites in this layer, unlike
        // the two everywhere else, and each placed where nothing in mark is live.
        // Note mark is 2n long here: the hash pass stamps cliques at c + n, so this
        // is the one layer whose array is not n. Not inside amd3Eliminate, which
        // holds three stamps live in turn: cliqueTag and absorbedTag across the
        // prune loop, then the merged set across the C[pivot] compaction. Never
        // observed to fire.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        auto [neighbors, absorbedCliques, prunedEdges] =
            amd3Eliminate(A, I, C, eliminated, mark, tag, pivot);
        ++numEliminations;
        std::uint32_t degree = neighbors.size();
        pivots.push_back(pivot);

        buckets.unfile(degrees[pivot], pivot);  // the pivot has left the graph
        degrees[pivot] = 0;

        // ---- the BOUND, in place of md5's exact refresh --------------------
        // C[pivot] is the new clique. Everything it reached needs a new degree, and
        // the bound replaces the union with a sum of separate contributions.
        //
        // Set view of the three quantities, none of which is built as a set:
        //
        //     pivotClique = C[pivot]
        //     degme       = |C[pivot]|             weighted, original vertices
        //     outside[c]  = |C[c] - C[pivot]|      ONE value per CLIQUE
        //
        // The last line is the whole idea. |C[c] - C[pivot]| depends on c and not
        // on the vertex asking, so it is computed once and read many times, where
        // the exact degree recomputes a union per vertex. mark[v] == inClique is
        // the membership test for C[pivot]; a second tag makes touchedCliques a set
        // too, so a clique is listed once however many vertices reach it.
        // The second site, guarding one contiguous region: inClique, seenClique, the
        // outside[c] loop, aggressive absorption's deadTag, and the bound loop's
        // per-vertex amd3ExactDegree calls. inClique is stamped here and still read
        // inside the outside[c] loop, so no sweep may land between them. The whole
        // region advances the tag by about n.
        if (tag > TAG_CEILING) {
            std::fill(mark.begin(), mark.end(), NIL);
            tag = 0;
            ++numTagSweeps;
        }
        std::vector<std::int32_t> pivotClique = C[pivot];
        ++tag;
        const std::int32_t inClique = tag;      // membership of C[pivot], one test
        for (std::int32_t v : pivotClique) mark[v] = inClique;
        // degme is NOT taken here. It is |C[pivot]| after mass elimination, and mass
        // elimination now runs below, after the absorption; see entry 3. The scan
        // over the cliques that follows is deliberately over the UNTRIMMED clique,
        // which is what AMD_2's scan 1 walks: it runs over the whole of Lme before
        // any of it has been mass eliminated. Nothing is lost by that, since a
        // vertex the merge will take belongs to no clique but the new one, so it
        // cannot appear in any touched clique's member list either way.

        // |C[c] - C[pivot]| ONCE PER CLIQUE. This is the whole reason the bound is
        // cheap: the quantity depends on c alone, so every vertex whose incidence
        // list holds c reads it rather than recomputing it.
        std::vector<std::int32_t> touchedCliques;
        ++tag;
        const std::int32_t seenClique = tag;
        for (std::int32_t u : pivotClique)
            for (std::int32_t c : I[u])
                if (c != pivot && mark[c] != seenClique) {
                    mark[c] = seenClique;
                    touchedCliques.push_back(c);
                }
        for (std::int32_t c : touchedCliques) {
            std::size_t total = 0;
            for (std::int32_t v : C[c])
                if (mark[v] != inClique && !eliminated[v]) total += superMembers[v].size();
            outside[c] = total;
            numMemberVisits += C[c].size();     // what an exact degree pays PER VERTEX
        }

        // AGGRESSIVE ABSORPTION, the first of amd3's two extras. Set view:
        // dead = { c : C[c] <= C[pivot] }, the
        // containment decided by the count already computed for the bound, since
        // |C[c] - C[pivot]| == 0 IS C[c] <= C[pivot]. Then I[u] = I[u] - dead for
        // every u in C[pivot], one stamp and one compaction pass, the same shape as
        // the absorption in the eliminator. Ordinary absorption kills only what the
        // pivot touched; this kills what any reached vertex touched, and it is free
        // because the quantity was computed for the bound anyway.
        std::vector<std::int32_t> deadCliques;
        for (std::int32_t c : touchedCliques)
            if (outside[c] == 0) deadCliques.push_back(c);
        if (!deadCliques.empty()) {
            ++tag;
            const std::int32_t deadTag = tag;
            for (std::int32_t c : deadCliques) { C.erase(c); mark[c] = deadTag; }
            std::vector<std::int32_t> keptCliques;
            for (std::int32_t u : pivotClique) {
                keptCliques.clear();
                for (std::int32_t c : I[u])
                    if (mark[c] != deadTag) keptCliques.push_back(c);
                I[u].swap(keptCliques);         // I[u] - dead
            }
            numAbsorbed += deadCliques.size();
        }

        // MASS ELIMINATION, and it runs HERE rather than in the eliminator, which is
        // ledger entry 3. u is INDISTINGUISHABLE from the pivot when the two have
        // the same closed neighborhood, so that everything u can still reach lies
        // inside the new clique. The test is a cheap sufficient condition for that
        // and is textually AMD_2's `Elen[i] == 1 && p3 == pn`: nothing explicit left
        // and no clique but the new one.
        //
        // It has to come AFTER the absorption above, because absorption is what
        // makes the cheap test agree with the true one. A clique whose members all
        // lie inside C[pivot] contributes nothing to what u can reach, yet its
        // presence in I[u] makes `I[u] == [pivot]` false, so the test declines a
        // genuine merge. AMD_2 says this itself: with aggressive absorption,
        // `deg == 0` is identical to the structural test. Asking first, as this
        // layer's parent does, is asking of an I[u] that still holds cliques about
        // to be removed.
        //
        //     merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
        //     C[pivot] = C[pivot] - merged
        std::vector<std::int32_t> mergedVertices;
        for (std::int32_t u : pivotClique) {
            if (A[u].empty() && I[u].size() == 1 && I[u][0] == pivot) {
                I[u].clear();
                eliminated[u] = true;
                mergedVertices.push_back(u);
            }
        }
        if (!mergedVertices.empty()) {      // C[pivot] - merged, one compaction pass
            ++tag;
            const std::int32_t mergedTag = tag;
            for (std::int32_t u : mergedVertices) mark[u] = mergedTag;
            std::vector<std::int32_t> keptClique;
            for (std::int32_t v : C[pivot])
                if (mark[v] != mergedTag) keptClique.push_back(v);
            C[pivot].swap(keptClique);
            pivotClique = C[pivot];
        }

        numCliqueEntries += C[pivot].size();
        numEliminatedVertices += 1 + mergedVertices.size();
        for (std::int32_t u : mergedVertices) {   // the pivot now stands for them too
            superMembers[pivot].insert(superMembers[pivot].end(),
                                       superMembers[u].begin(), superMembers[u].end());
            superMembers[u].clear();
        }
        numLive -= superMembers[pivot].size();  // every original the pivot stands for
        for (std::int32_t u : mergedVertices) { // and so have the merged vertices
            buckets.unfile(degrees[u], u);
            degrees[u] = 0;
        }

        // |C[pivot]| weighted, taken AFTER the merges, which is the value AMD_2
        // reads. Its `degme` is decremented inside scan 2 as each vertex is mass
        // eliminated, but it is not consumed there: the term enters a survivor's
        // degree only in the later pass that restores the degree lists,
        // `deg = Degree[i] + degme - nvi`, by which point degme is final. So every
        // survivor sees the same number, which is what this line gives.
        std::uint32_t degme = 0;
        for (std::int32_t v : pivotClique) degme += superMembers[v].size();

        const std::uint32_t numLeft = numLive;
        const std::vector<std::int32_t>& refreshedVertices = pivotClique;
        // LEDGER ENTRY 4, convention. This loop no longer finishes a bound. It
        // computes the part that does not involve u's own weight and stores it, and
        // the fourth pass below adds the new clique's size and subtracts the weight.
        // That split is AMD_2's: scan 2 forms `deg = sum dext(e) + sum nvj` and keeps
        // `Degree[i] = MIN(Degree[i], deg)`, and only the later pass that restores
        // the degree lists finishes it with `deg = Degree[i] + degme - nvi` and
        // `deg = MIN(deg, nleft - nvi)`.
        //
        // Why the split is load-bearing rather than a rearrangement: supervariable
        // detection runs BETWEEN the two, and a hash merge does `Nv[i] += Nv[j]`, so
        // `nvi` in the fourth pass is the weight AFTER the merge. Subtracting it here
        // would use the weight u had before absorbing v, which files a supervariable
        // one bucket too high per vertex merged into it, and it is never picked as
        // early as its size has earned. graph3 is the case: 8 absorbs 5, the vendored
        // routine files it at 2 and we filed 3, and from that pivot the two runs
        // order different graphs.
        //
        // This is mmd's ledger entry 5 in a different array, which the README says
        // cannot happen on the amd branch because AMD files at an external degree
        // that does not move when a weight changes. The external degree does not; the
        // `- nvi` term does, and it is the term that decides the bucket.
        std::vector<std::uint32_t> partial(n, 0);
        for (std::int32_t u : refreshedVertices) {
            // bound = |A[u]| + |C[pivot] - {u}| + sum |C[c] - C[pivot]| over the
            // cliques in I[u] - {pivot}, against the exact
            // |( A[u] | C[c] for c in I[u] ) - {u}|. The bound replaces the union
            // by a sum, so an overlap outside C[pivot] is counted once per clique
            // that holds it, which is exactly where it overcounts. The middle term
            // and the caps are the fourth pass's; what is formed here is the rest.
            std::size_t explicitPart = 0;
            for (std::int32_t v : A[u]) explicitPart += superMembers[v].size();
            std::size_t deg = explicitPart;
            for (std::int32_t c : I[u]) {
                if (c == pivot) continue;
                deg += outside[c];
                ++numCliqueReads;               // what the bound pays instead
            }
            // AMD_2's `Degree[i] = MIN (Degree[i], deg)`. The stored degree is a full
            // one from a previous iteration and this is a partial, so the two are
            // comparable only once the fourth pass adds `degme - nvi` to whichever
            // won. That is why the minimum is taken here and the common term added
            // there, rather than the other way round.
            // `deg` accumulates `outside[c]` over I[u], each term up to n and O(n) of them, so
            // it stays WIDE and the minimum is what brings it back into range.
            partial[u] = static_cast<std::uint32_t>(std::min<std::size_t>(deg, degrees[u]));
        }
        // HASH SUPERVARIABLE DETECTION, the second extra. Vertices indistinguishable from EACH
        // OTHER, which the pivot test cannot see. Hash first so the exact
        // comparison runs only within a bucket; the hash is a filter, never the
        // decision.
        // The buckets are an array indexed by the hash value, allocated once and
        // cleared only where it was used, which is Amd.cpp's Head[hval]. A map keyed
        // by the hash would cost a log per insertion and a node per group, for a
        // quantity that is already an index into 0 .. n.
        std::vector<std::uint32_t> usedKeys;
        for (std::int32_t u : pivotClique) {
            if (eliminated[u]) continue;
            // The hash stands for the PAIR of sets (A[u], I[u]), so equal sets
            // always collide and unequal ones usually do not. It is a filter and
            // never the decision: a collision costs a comparison, not a wrong merge.
            //
            // A SUM, because addition has no order and neither do the sets. Sorting
            // to build a key would be a log factor for nothing; Amd.cpp sums the
            // indices and reduces modulo n for the same reason.
            std::size_t key = 0;
            for (std::int32_t v : A[u])
                if (!eliminated[v]) key += static_cast<std::size_t>(v) + 1;
            // ONE SUM, WITH NO STRIDE, and that is ledger entry 8. This added the
            // incidence half as `(c + 1) * (n + 1)` until 2026-08-09, so that a vertex
            // and a clique of the same index could not cancel. True of the KEY and false
            // of the BUCKET: the modulus below is the same number as the stride, so the
            // incidence term is annihilated exactly and the hash came out a function of
            // the ADJACENCY ALONE. Amd.cpp lets a vertex and a clique collide on purpose,
            // the hash being a filter and never the decision. The invariant the two lines
            // hold TOGETHER is that the modulus must not divide the stride, and having no
            // stride is the cheapest way to hold it.
            //
            // THE MODULUS IS n, AS Amd.cpp'S IS AND AS PRODUCTION'S IS. It was n + 1 until
            // 2026-08-23, left behind when the stride went on 2026-08-09 and never
            // followed. The hash is a filter, so the two moduli give the same merges and
            // the same permutation and differ only in how many exact comparisons they
            // cost.
            for (std::int32_t c : I[u])
                key += static_cast<std::size_t>(c) + 1;
            const std::size_t k = key % n;
            if (hashBucket[k].empty()) usedKeys.push_back(k);
            hashBucket[k].push_back(u);
        }
        std::vector<std::pair<std::int32_t, std::int32_t>> hashPairs;
        for (std::size_t k : usedKeys) {
            std::vector<std::int32_t>& group = hashBucket[k];
            if (group.size() < 2) continue;
            // LEDGER ENTRY 1, convention. The bucket is walked BACKWARD, because
            // Amd.cpp builds it by pushing at the head while scan 2 walks Lme
            // FORWARD, `Next[i] = Head[hval]; Head[hval] = i` at its line 1940,
            // so its chain comes out reversed against C[p] and the pair loop then
            // takes the head. The survivor of an indistinguishable pair is
            // therefore the member of C[p] seen LAST, where appending and walking
            // forward keeps the one seen FIRST. Same set of merges either way;
            // only which of the two absorbs the other moves, and minimum degree
            // is settled by exactly that.
            //
            // The walk is reversed rather than the fill, which is the same order
            // at O(1) per insertion instead of O(bucket). `usedKeys` keeps its
            // forward order deliberately: Amd.cpp reaches a bucket by rescanning
            // Lme forward and taking the first member it meets whose bucket is
            // still full, so buckets are PROCESSED in C[p] order while their
            // CONTENTS are reversed. The two orders are independent and only the
            // second decides a merge.
            for (std::size_t x = group.size(); x-- > 0; ) {
                std::int32_t u = group[x];
                if (eliminated[u]) continue;
                for (std::size_t y = x; y-- > 0; ) {
                    std::int32_t v = group[y];
                    if (eliminated[v]) continue;
                    ++numHashPairs;              // the witness; see its declaration
                    // The exact test, which the hash only filters for:
                    //     A[u] - {v} == A[v] - {u}  and  I[u] == I[v]
                    // Decided by stamping one side and counting matches on the
                    // other, as every other membership test in this file is, so it
                    // costs one pass and no sort.
                    // The third site, and the reason this layer has one more than the
                    // others. `other` advances once per PAIR TESTED rather than once
                    // per pass, and the pair count is quadratic in the bucket sizes
                    // with no clean bound, so a check before the pass would leave the
                    // gap between checks unbounded. Safe at the top of a pair because
                    // the previous pair's stamps are spent and hashBucket, usedKeys
                    // and eliminated are separate structures.
                    if (tag > TAG_CEILING) {
                        std::fill(mark.begin(), mark.end(), NIL);
                        tag = 0;
                        ++numTagSweeps;
                    }
                    ++tag;
                    const std::int32_t other = tag;
                    std::uint32_t sizeV = 0;
                    for (std::int32_t w : A[v])
                        if (w != u && !eliminated[w]) { mark[w] = other; ++sizeV; }
                    // LEDGER ENTRY 6, and it is the other half of entry 5 rather than a new
                    // idea. Both walks skip I[..][0], the NEW CLIQUE. Entry 5 put it at the front
                    // of every I[u], and u and v are both members of C[pivot], so both lists begin
                    // with the same entry and it can never discriminate. Amd.cpp skips it outright,
                    // `for (p = Pe[j] + 1; ...)` with the comment "skip the first element in the
                    // list (me)", and that skip is part of the same mechanism as putting me first:
                    // porting one without the other leaves a guaranteed match at the head of a
                    // short-circuiting walk, so every FAILING pair pays one extra iteration before
                    // it can fail.
                    //
                    // Measured in production on a 140x140 grid: 2.08 incidence iterations per pair
                    // against amd2's 1.21, on the line Instruments put at 6.22 s of a 14.90 s run,
                    // and 1.08 after. Exactly one extra iteration per pair, which is what a
                    // guaranteed match at position zero predicts.
                    //
                    // Skipped on BOTH sides, since the two walks feed sizeV and sizeU and those are
                    // compared: dropping it from one alone would make every pair fail on the count.
                    for (std::size_t ci = 1; ci < I[v].size(); ++ci) {   // stamped past the vertices
                        mark[I[v][ci] + cliqueStamp] = other;
                        ++sizeV;
                    }
                    std::uint32_t sizeU = 0;
                    bool same = true;
                    for (std::int32_t w : A[u]) {
                        if (w == v || eliminated[w]) continue;
                        ++sizeU;
                        if (mark[w] != other) { same = false; break; }
                    }
                    if (same)
                        for (std::size_t ci = 1; ci < I[u].size(); ++ci) {   // skip the new clique
                            ++sizeU;
                            if (mark[I[u][ci] + cliqueStamp] != other) {
                                same = false;
                                break;
                            }
                        }
                    if (!same || sizeU != sizeV) continue;
                    // v is folded into u and left exactly where it lies, with a
                    // weight of zero, which is Amd.cpp's Nv[v] = 0 and the same move
                    // the dense prepass makes. Removing it from every clique and
                    // every adjacency would cost a pass over the whole structure per
                    // merge, against O(1) for this.
                    //
                    // Nothing is lost by leaving it. The merge required
                    // A[u] - {v} == A[v] - {u}, so every list holding v holds u as
                    // well: v is redundant wherever it appears, never the only way
                    // to reach anything. The walks skip it.
                    superMembers[u].insert(superMembers[u].end(),
                                           superMembers[v].begin(), superMembers[v].end());
                    superMembers[v].clear();
                    buckets.unfile(degrees[v], v);
                    A[v].clear();
                    I[v].clear();
                    eliminated[v] = true;
                    ++numEliminatedVertices;
                    hashPairs.push_back({u, v});
                    ++numHashMerges;
                }
            }
        }

        for (std::size_t k : usedKeys) hashBucket[k].clear();  // only what was used

        // THE FOURTH PASS, which finishes the bounds and files them. AMD_2 spells it
        // `deg = Degree[i] + degme - nvi` and `deg = MIN (deg, nleft - nvi)`, under
        // its RESTORE DEGREE LISTS heading, and it runs here for the reason in entry
        // 4: `nvi` is read after supervariable detection, so a vertex that absorbed
        // another subtracts the combined weight.
        //
        // `degme` is |C[pivot]| after the merges and is the same for everyone, and
        // `numLeft` is the live original count, both settled before this loop, so the
        // pass is one addition and two comparisons per survivor.
        for (std::int32_t u : refreshedVertices) {
            if (eliminated[u]) continue;        // absorbed by the hash a moment ago
            const std::uint32_t weightU = superMembers[u].size();  // POST-merge, the point
            // One operand widened, so the sum is formed in std::size_t: `partial[u]` and `degme`
            // each reach n, so the sum reaches 2n. Widening cannot be done after the addition
            // the way narrowing is done after a subtraction.
            std::size_t bound = static_cast<std::size_t>(partial[u]) + degme - weightU;
            bound = std::min<std::size_t>(bound, numLeft - weightU);
            // NOT PRODUCTION: instrumentation. This computes the very union the bound exists to
            // avoid, and its only purpose is to show the truth beside the estimate. It is taken
            // here rather than in the pass above so that it sees the same graph the bound does,
            // supervariable detection having run in between.
            exact[u] = amd3ExactDegree(A, I, C, eliminated, superMembers, mark, tag, u);
            ++numBoundChecks;
            if (bound > exact[u]) ++numLooseBounds;
            // The other direction, which is not a quality signal but an INVARIANT. A bound
            // may exceed the degree by any amount and still be a bound; falling below it is
            // the one thing it must never do, since the picker would then be told a vertex is
            // cheaper than it is. Counted because it was not: the layer measured looseness
            // only, and a cap taken from the wrong counter drove this negative 22 times on a
            // 10 by 10 grid while every example stayed green. Anything but zero here is a
            // defect. See the amd3 subsection of README.md.
            if (bound < exact[u]) ++numBoundsBelowExact;
            amd3Refile(buckets, degrees, u, bound);
        }

        numBoundUpdates += refreshedVertices.size();
        for (std::int32_t u : refreshedVertices) minDegree = std::min(minDegree, degrees[u]);

        for (std::int32_t c : touchedCliques) outside[c] = 0;   // clear what this iteration wrote

        // A supervariable of size w is w consecutive columns of L. Its external
        // degree is what remains of the clique after the merges, since a merged
        // vertex joins the supervariable instead of neighboring it, and every
        // member left there is a live vertex standing for itself alone. The first
        // column then holds ext + w - 1 entries below its diagonal, the next
        // ext + w - 2, down to ext, and each column contributes its own diagonal.
        std::uint32_t superSize = superMembers[pivot].size();
        // Weighted, because hash detection folds a vertex into a LIVE one, so a
        // member of the new clique can stand for several original vertices, and a
        // merged one stands for none and is still lying there. amd1 can use the
        // plain count and this file cannot, which is the same inheritance the
        // numLive counter above records.
        std::uint32_t externalDegree = 0;
        for (std::int32_t v : C[pivot])
            if (!eliminated[v]) externalDegree += superMembers[v].size();
        // ONE FACTOR WIDENED: a product of two one-dimensional quantities is TWO dimensional, so
        // it is formed in std::size_t. Widening cannot be done after the multiply the way
        // narrowing is done after a subtraction; one operand is enough, the other promoting to
        // meet it. Both factors are bounded by n, so the product reaches n^2.
        nnzL += static_cast<std::size_t>(superSize) * externalDegree
              + static_cast<std::size_t>(superSize) * (superSize - 1) / 2 + superSize;

        std::ostringstream absorbedCliquesText;
        if (absorbedCliques.empty()) {
            absorbedCliquesText << "none";
        } else {
            bool first = true;
            for (std::int32_t c : absorbedCliques) {
                absorbedCliquesText << (first ? "" : ", ") << "c" << c;
                first = false;
            }
        }
        std::ostringstream prunedEdgesText;
        if (prunedEdges.empty()) {
            prunedEdgesText << "none";
        } else {
            bool first = true;
            for (auto [u, v] : prunedEdges) {
                prunedEdgesText << (first ? "" : ", ") << u << "-" << v;
                first = false;
            }
        }
        std::ostringstream mergedVerticesText;
        if (mergedVertices.empty()) {
            mergedVerticesText << "none";
        } else {
            bool first = true;
            for (std::int32_t u : mergedVertices) {
                mergedVerticesText << (first ? "" : ", ") << u;
                first = false;
            }
        }
        std::ostringstream refreshedVerticesText;
        if (refreshedVertices.empty()) {
            refreshedVerticesText << "none";
        } else {
            bool first = true;
            for (std::int32_t u : refreshedVertices) {
                refreshedVerticesText << (first ? "" : ", ") << u;
                first = false;
            }
        }
        std::ostringstream title;
        title << "iteration " << iteration << ": eliminate " << pivot << " (degree " << degree
              << ", size " << superSize << ", external degree " << externalDegree
              << "), absorbed cliques: " << absorbedCliquesText.str()
              << ", pruned edges: " << prunedEdgesText.str()
              << ", merged vertices: " << mergedVerticesText.str()
              << ", refreshed vertices: " << refreshedVerticesText.str();
        // NOT PRODUCTION: display only. The trace is what makes these files teachable and
        // is the whole reason they exist; nothing downstream reads it.
        if (n <= SHOW_THRESHOLD) {
            amd3Show(A, I, C, degrees, exact, title.str(), &eliminated);
            amd3ShowState(degrees, buckets, minDegree, superMembers, eliminated, pivots);
        }
        ++iteration;
    }

    std::vector<std::int32_t> order;
    for (std::int32_t pivot : pivots)
        for (std::int32_t u : superMembers[pivot]) order.push_back(u);
    std::cout << "n = " << n << ", nnz(L) = " << nnzL
              << " against nnz(tril A) = " << nnzTrilA
              << ", fill = " << (nnzL - nnzTrilA) << "\n";
    std::cout << "iterations: " << numIterations << "\n";
    std::cout << "eliminations: " << numEliminations << "\n";
    std::cout << "sum of |C[p]|: " << numCliqueEntries << "\n";
    std::cout << "bound computations: " << (numBoundUpdates + n)
              << ", bound updates: " << numBoundUpdates
              << ", bucket probes: " << numBucketProbes << "\n";
    std::cout << "clique-member visits an exact degree would need: "
              << numMemberVisits << "\n";
    std::cout << "clique reads the bound needed:                    "
              << numCliqueReads << "\n";
    // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    std::cout << "aggressively absorbed: " << numAbsorbed
              << ", hash merges: " << numHashMerges
              << ", hash pairs tested: " << numHashPairs << "\n";
    // NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    std::cout << "bound below exact " << numBoundsBelowExact
              << " times, which must be zero\n";
    std::cout << "bound was loose " << numLooseBounds << " times out of "
              << numBoundChecks << "\n";
    std::cout << "tag sweeps: " << numTagSweeps << "\n";
    std::cout << "order: [";
    for (std::size_t k = 0; k < order.size(); ++k)
        std::cout << (k == 0 ? "" : ", ") << order[k];
    std::cout << "]\n";
    return order;
}


// A square grid graph, four-neighbor, for running the counters at a size the seven examples
// cannot reach. It is here rather than among them because it is not an example: nothing about it
// illustrates a mechanism, and above SHOW_THRESHOLD its trace is not printed at all.
//
// It must match the Python twin's grid_graph exactly, vertex for vertex, or `make test` would be
// diffing two different problems.
static Graph gridGraph(int side) {
    const int n = side * side;
    Graph graph(n);
    for (int r = 0; r < side; ++r)
        for (int c = 0; c < side; ++c) {
            const int u = r * side + c;
            if (r > 0)        graph[u].push_back(u - side);
            if (c > 0)        graph[u].push_back(u - 1);
            if (c + 1 < side) graph[u].push_back(u + 1);
            if (r + 1 < side) graph[u].push_back(u + side);
        }
    return graph;
}

void run(const std::string& name, const Graph& G) {
    std::cout << "=== " << name << " ===\n";
    amd3MinimumDegree(G);
    std::cout << "\n";
}

int main(int argc, char** argv) {
    // Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff
    // at a size the seven examples cannot reach. Nothing is filtered: the run is silent above
    // SHOW_THRESHOLD and prints its closing lines as always.
    //
    //   ./amd3_cpp grid 22
    if (argc > 2 && std::string(argv[1]) == "grid") {
        const int side = std::atoi(argv[2]);
        std::cout << "=== grid " << side << "x" << side << " (n = " << side * side << ") ===\n";
        amd3MinimumDegree(gridGraph(side));
        return 0;
    }

    // The same three graphs as md1 and md2.
    //
    //   graph1, a 4-cycle: eliminating any vertex forces its two neighbors
    //   together, so it is the smallest graph that fills (one fill edge).
    //
    //      0---1          edges: 0-1 1-2 2-3 3-0
    //      |   |
    //      3---2
    //
    //   graph2, uneven degrees so the picker actually chooses; it fills twice.
    //
    //        0            edges: 0-1 0-2 1-3 2-4
    //       / \                  3-4 3-5 4-5
    //      1   2
    //      |   |
    //      3---4
    //       \ /
    //        5
    //
    //   graph3, twelve vertices: a path 0-1-...-11 with eight extra edges. Big
    //   enough that cliques grow past two members, which is where the quotient
    //   graph starts to pay, and its elimination order is not the identity.
    //
    //      edges: 0-1 0-3 0-8 1-2 1-6 1-8 2-3 2-5 3-4 4-5
    //             5-6 5-9 6-7 6-10 7-8 8-9 9-10 10-11
    Graph graph1 = {
        {1, 3}, {0, 2}, {1, 3}, {0, 2},
    };
    Graph graph2 = {
        {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
    };
    Graph graph3 = {
        {1, 3, 8},        // 0
        {0, 2, 6, 8},     // 1
        {1, 3, 5},        // 2
        {0, 2, 4},        // 3
        {3, 5},           // 4
        {2, 4, 6, 9},     // 5
        {1, 5, 7, 10},    // 6
        {6, 8},           // 7
        {0, 1, 7, 9},     // 8
        {5, 8, 10},       // 9
        {6, 9, 11},       // 10
        {10},             // 11
    };

    // graph4, eight vertices and fourteen edges. Denser than the others, and here
    // for one specific reason: it is the smallest graph we could find on which
    // AMD's degree BOUND is ever loose. The bound overcounts only when a vertex
    // belongs to two cliques that overlap outside the new one, which needs enough
    // eliminations to have made several cliques and enough fill for them to
    // intersect. Every connected graph on five or six vertices is exact (checked
    // exhaustively), and so are graph1 to graph3, so without this one the amd
    // trace would never show the approximation approximating. The other layers use
    // it as an ordinary denser test.
    //
    //   edges: 0-2 0-3 0-4 0-7 1-3 1-4 1-6 1-7 2-3 2-5 3-6 3-7 4-5 5-6
    Graph graph4 = {
        {2, 3, 4, 7},     // 0
        {3, 4, 6, 7},     // 1
        {0, 3, 5},        // 2
        {0, 1, 2, 6, 7},  // 3
        {0, 1, 5},        // 4
        {2, 4, 6},        // 5
        {1, 3, 5},        // 6
        {0, 1, 3},        // 7
    };

    // graph5, five vertices and four edges, two paths joined at 4: 2-1-4-0-3.
    // Small and fill free, and here for one reason: it is the smallest graph on
    // which amd3's merge test declines a genuine supervariable. At the iteration whose
    // pivot is 0 and whose clique is {4}, vertex 4 has nothing explicit left but
    // belongs to c1 as well as to the new clique, so I[4] == {pivot} fails even
    // though c1's only member is 4 itself and everything 4 reaches lies inside
    // the new clique. The exact test amd3Neighbors(A, I, C, u) contained in
    // C[pivot] would merge it. See the README section on mass elimination.
    //
    //   edges: 0-3 0-4 1-2 1-4
    Graph graph5 = {
        {3, 4},           // 0
        {2, 4},           // 1
        {1},              // 2
        {0},              // 3
        {0, 1},           // 4
    };

    // graph6, six vertices and eight edges. Here because one small graph carries
    // three things at once. Its supervariable {0, 4} is a supernode but NOT a
    // fundamental one: the elimination forest is 2 -> 1 -> 4 and 3 -> 0 -> 4, so
    // 4 already has 1 as a child when 0 merges into it. The merge happens at iteration
    // 2 of 5, so the run continues afterwards and the selection degree, 3 over
    // {2, 3, 4}, differs from the external degree, 2 over {2, 3}, with the
    // difference being the size of what merged. And superMembers ends with a hole
    // in the middle, slot 4 empty between two used ones, while no pivot equals
    // its own iteration number. See the README sections on mass elimination and on
    // external degree.
    //
    //   edges: 0-2 0-3 0-4 1-3 2-3 2-4 2-5 3-4
    Graph graph6 = {
        {2, 3, 4},        // 0
        {3},              // 1
        {0, 3, 4, 5},     // 2
        {0, 1, 2, 4},     // 3
        {0, 2, 3},        // 4
        {2},              // 5
    };

    // graph7, five vertices and six edges. The pairwise case: at the iteration whose
    // pivot is 0 and whose clique is {2, 4}, vertices 2 and 4 are
    // indistinguishable FROM EACH OTHER, both reaching the same closed
    // neighborhood, yet neither is absorbable into the pivot, since each still
    // reaches 3 from outside the clique. No test framed against the pivot finds
    // them, and the exact test does not help either: both orders are 1 0 (2 3 4).
    // Catching such pairs needs a comparison between candidates, which is what
    // amd's hashing does. See the README section on detecting supervariables
    // against each other.
    //
    //   edges: 0-1 0-2 0-4 1-4 2-3 2-4 3-4
    Graph graph7 = {
        {1, 2, 4},        // 0
        {0, 4},           // 1
        {0, 3, 4},        // 2
        {2, 4},           // 3
        {0, 1, 2, 3},     // 4
    };

    std::vector<std::pair<std::string, Graph>> examples = {
        {"graph1", graph1}, {"graph2", graph2},
        {"graph3", graph3}, {"graph4", graph4},
        {"graph5", graph5}, {"graph6", graph6},
        {"graph7", graph7},
    };

    // All of them by default. To run just one, pass its number: ./amd3_cpp 3
    int selected = (argc > 1) ? std::atoi(argv[1]) : 0;
    for (int number = 1; number <= static_cast<int>(examples.size()); ++number) {
        if (selected != 0 && number != selected) continue;
        run(examples[number - 1].first, examples[number - 1].second);
    }
    return 0;
}
