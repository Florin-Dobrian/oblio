#pragma once

// ElmOrder.h - the elimination order of a sparse matrix, and what computing it cost in clique
// storage. Produced by all five of Oblio's own ordering drivers.
//
// HALF A PERMUTATION, and deliberately. `Permutation` holds both directions and computes the
// inverse when it is built, which is an O(n) pass most callers of an ordering never want: the
// digest hashes the order, the timing and profile harnesses discard it, and the two checkers
// compare it against a vendored routine. `OrderEngine` is the one caller that wants a permutation,
// and it completes this into one.
//
// THE STATISTICS RIDE WITH THE RESULT because they are a property of what was computed rather than
// of whatever computed it, which is the arrangement the rest of the pipeline already uses:
// `NumFactorStatic` carries `numPerturbations` and its delayed columns, both written by an engine
// that keeps no such state of its own.
//
// The quotient graph a driver builds is a local, destroyed when the driver returns, so a caller
// that wants these figures has to be handed them. They travel here rather than in the signature.
//
// BUILT BY ITS CONSTRUCTOR AND NOT WRITTEN INTO AFTERWARDS, which is why there is no friend and no
// setter: a driver creates the whole thing at once, where an engine fills a factor the caller
// declared. The accessors carry the quotient graph's own names, so a driver's last statement reads
// straight across.

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace Oblio {

template<class QuotientGraph> class AmdEngine;
template<class QuotientGraph> class MmdEngine;

class ElmOrder {
public:
    ElmOrder() = default;

    // TEMPORARY, for the one driver not yet converted to an engine: `MmdChained` still builds the
    // whole object and returns it. It goes when that one follows.
    ElmOrder(std::vector<std::int32_t> order, std::size_t numPeakCliqueMembers,
             std::size_t numBornCliqueMembers, std::size_t numCompactions)
        : mOrder(std::move(order)),
          mNumPeakCliqueMembers(numPeakCliqueMembers),
          mNumBornCliqueMembers(numBornCliqueMembers),
          mNumCompactions(numCompactions) {}

    // The elimination order over the ORIGINAL vertices: `order()[k]` is the vertex eliminated k-th.
    const std::vector<std::int32_t>& order() const { return mOrder; }

    std::size_t size() const { return mOrder.size(); }
    bool        empty() const { return mOrder.empty(); }

    // PEAK LIVE CLIQUE MEMBERS. A member is a vertex in a live clique at this instant, against an
    // ENTRY, which is a slot in the store. It is a property of the ALGORITHM and not of the layout,
    // so a branch's drivers must agree on it however they store their cliques, which is what makes
    // it worth comparing: two drivers agreeing on a permutation can still be caught doing different
    // work, and have been.
    std::size_t numPeakCliqueMembers() const { return mNumPeakCliqueMembers; }

    // MEMBERS BORN, the same quantity summed over the whole run rather than taken at an instant.
    // Also a property of the algorithm. Zero from a driver that does not track it.
    std::size_t numBornCliqueMembers() const { return mNumBornCliqueMembers; }

    // HOW OFTEN THE CLIQUE STORE RAN OUT AND WAS COMPACTED. Zero for a store that never compacts,
    // which is an answer rather than the absence of one: an arena only appends and a chained store
    // needs no room.
    std::size_t numCompactions() const { return mNumCompactions; }

    // HOW OFTEN THE TAG ARRAY HAD TO BE RESET. A tag scheme is cheap only while the tag can keep
    // advancing; when it approaches its type's ceiling the array has to be swept back to the
    // alive-and-unseen state and the tag restarted, which is O(n) each time. Zero on everything the
    // suite runs today, and that is what makes it worth reporting: the day it is not, the scheme has
    // started paying for itself and nothing else would say so.
    std::size_t numTagResets() const { return mNumTagResets; }

private:
    template<class QuotientGraph> friend class AmdEngine;   // fills the order via the engine
    template<class QuotientGraph> friend class MmdEngine;

    std::vector<std::int32_t> mOrder;
    std::size_t               mNumPeakCliqueMembers = 0;
    std::size_t               mNumBornCliqueMembers = 0;
    std::size_t               mNumCompactions       = 0;
    std::size_t               mNumTagResets         = 0;
};

} // namespace Oblio
