#pragma once

// Buckets.h - the degree buckets Oblio's own minimum-degree orderings pick from, and the scan
// descriptor the amd branch passes alongside them. Shared by every driver and by all three
// quotient graph classes, none of which owns either type.
//
// They are here rather than in QuotientGraphFlat.h because they belong to no one graph: the
// compacted and chained classes used to include the flat class's header to reach them, which made
// two siblings depend on a third for something none of them defines.

#include "oblio/Types.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace Oblio {

// The degree buckets: one doubly linked list per degree, threaded through arrays of size n, so
// that filing, unfiling and taking the head are all O(1). MMD spells these fwd/bwd and AMD
// Next/Last. An ordered container cannot give O(1) removal from the middle, which is what a
// degree change needs and what happens far more often than a pick.
//
// The bodies are inline because they are single-statement pointer splices on the hot path, the
// same exception the tree makes for trivial accessors.
class Buckets {
public:
    // ALL THREE ARE SIZE n. The heads are indexed by DEGREE and the links by VERTEX, two ranges of
    // the same extent: every branch files at the true degree, which reaches n - 1. A caller that
    // reads a bucket must bound its own index; the mmd prepass does.

    // FOUR FACTS IN `mPrev[u]`, told apart by sign, which is what lets three arrays be three
    // rather than four:
    //
    //     mPrev[u] in [0, n)       filed; the vertex before it is mPrev[u]
    //     mPrev[u] <= -2           filed AT THE HEAD of bucket -mPrev[u] - 2
    //     mPrev[u] == UNFILED      not in any list                            (-1)
    //     mPrev[u] == OUTMATCHED   withheld from the buckets; see outmatch()  (INT32_MAX)
    //
    // THE HEAD CASE IS WHAT DELETES A DEGREE ARRAY. `unfile` reads the bucket back out of the
    // slot and takes no degree argument, so no caller has to keep one to unfile with. The amd
    // drivers keep theirs for a different reason: their bound READS the value.
    //
    // A PREDECESSOR IS STORED RAW AND THE HEAD FORM IS SHIFTED BY TWO. The shift is two and not
    // one because DEGREE 0 IS REACHABLE, an isolated vertex filing at 0, and -1 is taken by
    // UNFILED.
    //
    // BOTH SENTINELS ARE EQUALITY-ONLY. Nothing here orders `mPrev` against either, so they sit at
    // the two values a predecessor and a head cannot take, leaving the whole range below -1 to the
    // head form.
    //
    // IT FITS EXACTLY AT THE LARGEST ADMISSIBLE n, with nothing to spare at either end: a head
    // encodes a degree reaching n - 1, so the head form runs to -(n + 1) and lands on the type's
    // minimum at n == MAX_IDX, while a predecessor reaches n - 1 and leaves INT32_MAX free above
    // it. Filing at anything above the true degree does not fit.
    //
    // A DEGREE IS ONE DIMENSIONAL and so `std::uint32_t` in every signature here. The links stay
    // `std::int32_t`, carrying NIL and UNFILED; only the key does not.
    static constexpr std::int32_t UNFILED    = -1;
    static constexpr std::int32_t OUTMATCHED = std::numeric_limits<std::int32_t>::max();

    explicit Buckets(std::size_t size)
        : mHead(size, NIL), mNext(size, NIL), mPrev(size, UNFILED) {}

    // buckets[degree].add(u), at the head. The head is the only O(1) end of a singly reachable
    // list, so the winner among equal degrees is whatever was filed last rather than the lowest
    // index. That tie-break is what the oracles use, and it is why an ordering differs from an
    // exact scan's in its ties.
    void file(std::uint32_t degree, std::int32_t u) {
        mNext[u] = mHead[degree];
        mPrev[u] = -static_cast<std::int32_t>(degree) - 2;   // head of `degree`; see the encoding
        if (mHead[degree] != NIL) mPrev[mHead[degree]] = u;
        mHead[degree] = u;
    }

    // buckets[degree].discard(u). Idempotent, which matters during a batch: a vertex evicted
    // early can be merged away by a later pivot in the same round, and unfiling it twice must
    // not splice a list it is no longer in.
    void unfile(std::int32_t u) {
        const std::int32_t prev = mPrev[u];
        if (prev == UNFILED || prev == OUTMATCHED) return;   // not in a list; see the encoding
        if (prev >= 0) mNext[prev]      = mNext[u];
        else           mHead[-prev - 2] = mNext[u];          // u headed bucket -prev - 2
        if (mNext[u] != NIL) mPrev[mNext[u]] = prev;
        mNext[u] = NIL;
        mPrev[u] = UNFILED;
    }

    // Withhold u from the buckets without filing it anywhere. It stays live and reachable; it
    // simply cannot be the minimum before the vertex that outmatched it, so it is not a candidate
    // until an elimination reaches it and `restore` puts it back. Unfiling and restoring are one
    // operation in two calls, which is why they sit together at every call site here.
    void outmatch(std::int32_t u)         { unfile(u); mPrev[u] = OUTMATCHED; }
    void restore(std::int32_t u)          { if (mPrev[u] == OUTMATCHED) mPrev[u] = UNFILED; }
    bool outmatched(std::int32_t u) const { return mPrev[u] == OUTMATCHED; }

    // The next vertex in the same bucket, and whether u is filed at all. Only mmd reads either.
    // A VERTEX OUT OF EVERY LIST HAS BOTH LINKS FREE, and the amd branch parks the hash key in one
    // and the hash chain in the other for the middle of an elimination step, which is why it needs
    // no key array and no hash-head links of its own.
    //
    // TWO FORMS IN ONE SLOT, and a pair of names for each so a call site says which it means. The
    // prune parks the UNREDUCED sum, `setHashKey`/`hashKey`; the driver reduces it modulo the
    // bucket count and writes the result back over it, `setHashBucket`/`hashBucket`, which is what
    // lets detection find a vertex's bucket from the vertex. Both pairs are the same slot, and the
    // SIGNATURES are what tell them apart: a key is any bit pattern and cannot index anything, a
    // bucket is in [0, n) and indexes the hash heads directly.
    //
    // THE KEY PAIR CARRIES THE CONVERSION, so no caller writes a cast. A key is accumulated
    // unsigned, wrapping being defined there and `% n` being well defined only there, and the slot
    // is `std::int32_t` because it also has to hold NIL and the two sentinels. The round trip is
    // the two's complement reinterpretation; the uint32-to-int32 half of it is
    // implementation-defined before C++20.
    //
    // LEGAL ONLY BETWEEN unfile() AND file(). Either form is an arbitrary int32 and can look like
    // any of the encodings mPrev carries, so `unfile()`, `filed()` and `outmatched()` MUST NOT be
    // called on a vertex holding one; they would splice a list on garbage. No mmd driver may call
    // these: mmd leaves its candidates filed and asks `filed()` about them.
    void          setHashKey(std::int32_t u, std::uint32_t uHashKey)
                                       { mPrev[u] = static_cast<std::int32_t>(uHashKey); }
    std::uint32_t hashKey(std::int32_t u) const
                                       { return static_cast<std::uint32_t>(mPrev[u]); }
    void          setHashBucket(std::int32_t u, std::int32_t uHashBucket)
                                       { mPrev[u] = uHashBucket; }
    std::int32_t  hashBucket(std::int32_t u) const
                                       { return mPrev[u]; }
    void         setChain(std::int32_t u, std::int32_t v)                { mNext[u] = v; }
    std::int32_t chain(std::int32_t u) const                             { return mNext[u]; }

    std::int32_t next(std::int32_t u) const      { return mNext[u]; }
    bool         filed(std::int32_t u) const     { return mPrev[u] != UNFILED
                                                          && mPrev[u] != OUTMATCHED; }

    std::int32_t head(std::uint32_t degree) const  { return mHead[degree]; }
    bool         empty(std::uint32_t degree) const { return mHead[degree] == NIL; }

private:
    std::vector<std::int32_t> mHead;   // mHead[d], the first live vertex of degree d
    std::vector<std::int32_t> mNext;   // mNext[u], toward the tail
    std::vector<std::int32_t> mPrev;   // mPrev[u], toward the head
    // A byte per vertex, not std::vector<bool>. That specialization packs one bit per entry, and
    // what it costs here is CONSTRUCTION: a vector of n false is built through the bit-reference
    // machinery word by word where a vector of n zero bytes is a memset, and this is built once
    // per ordering.
};


// The same, for a driver carrying a TAGGED array instead of a value array and a separate
// seen-this-step mark. One array holds three facts: `markAmd[c] == 0` is absorbed, `0 < markAmd[c]
// < tagAmd` is alive but stale, and `markAmd[c] >= tagAmd` is seen this step with `markAmd[c] -
// tagAmd` the value. So there is no mark to carry and no clearing pass.
//
// `key` carries the ADJACENCY HALF of the hash key alone, and the reason is a phase boundary
// rather than a preference: aggressive absorption runs between this prune and the driver's bound
// pass and COMPACTS I[u] in place, so the list the other half must sum over does not exist yet
// here.
//
// REDUCED AS IT ACCUMULATES, modulo the driver's bucket count, which is why it is an int32. The
// key is only ever used modulo that number, so reducing early cannot change which bucket a vertex
// lands in, and the running value stays narrow enough to ride in an array the driver already has.
struct TaggedScan {
    // THE BUCKETS, OR NOT. A driver that parks its hash key in the links passes them: the prune
    // then takes every member of C[pivot] out of the degree lists and parks the adjacency half of
    // the hash key in the predecessor link it has just freed. A driver that refiles inside its own
    // bound pass, which is Amd2 and Amd2B, cannot have either: its links are still degree links
    // when the hash runs. Those pass NULL and build their key in a pass of their own.
    //
    // Null therefore means "leave the degree lists alone and store no key". It does not change
    // what the scan computes, only where the by-products go.
    Buckets*                          buckets;
    std::vector<std::int32_t>&        markAmd;       // per clique, the tagged workspace
    // Serves a LIVE vertex's degree and a DEAD one's clique weight from
    // one array, the two being disjoint because a clique id is the id of the pivot that formed it.
    // The scan reads only the clique half.
    const std::vector<std::uint32_t>& degree;
    std::vector<std::int32_t>&        touchedCliques;// the cliques this step reached, once each
    std::int32_t                      tagAmd;        // the tag for this elimination
    std::int32_t                      modulus;       // the driver's bucket count, n + 1
};

} // namespace Oblio
