#include "oblio/Mmd1.h"

#include <algorithm>
#include <utility>

namespace Oblio {

std::vector<std::int32_t> orderMmd1(const std::vector<std::size_t>&  colPtr,
                                    const std::vector<std::int32_t>& rowIdx,
                                    std::int32_t delta) {
    if (colPtr.empty()) return std::vector<std::int32_t>();
    const std::size_t size = colPtr.size() - 1;
    if (size == 0) return std::vector<std::int32_t>();

    QuotientGraph qg(colPtr, rowIdx);
    std::vector<std::int32_t> pivots;              // the order over supervariables
    pivots.reserve(size);
    std::uint32_t numEliminated = 0;               // a counter, not a scan of the flags

    // The degree cache, exact throughout: this branch keeps the honest set union and pays for it
    // by refreshing rarely. No weight array, because mass elimination merges only into the
    // pivot, which is eliminated in the same call, so no live vertex ever stands for more than
    // one original vertex.
    std::vector<std::uint32_t> degrees(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        degrees[u] = qg.adjacencySize(u);

    // minDegree is a LOWER BOUND on the current minimum, not the minimum itself. It is allowed
    // to lag and the walk corrects it; what it must never do is overshoot, since a vertex filed
    // below it would never be seen. So the walk raises it and a refresh lowers it, and nothing
    // else touches it.
    Buckets buckets(size);
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(size); ++u)
        buckets.file(degrees[u], u);
    std::uint32_t minDegree = *std::min_element(degrees.begin(), degrees.end());

    // The round a vertex was last evicted in, which accumulates the refresh set without a set
    // and without a sort: a vertex reached by two pivots in the same round is listed once.
    std::vector<std::int32_t> touchedRound(size, NIL);
    std::vector<std::int32_t> touched;              // hoisted, cleared per round
    std::int32_t numRounds = 0;

    while (numEliminated < size) {
        while (buckets.empty(minDegree)) ++minDegree;   // walk up to the first live bucket

        // ---- one batch, no degree refreshed inside it ---------------------------------
        // Clamped, since a live vertex's degree cannot exceed size - 1 and a wider window would
        // walk the bucket array off its end.
        std::uint32_t batchLimit = minDegree;
        if (delta > 0)
            batchLimit = std::min(minDegree + static_cast<std::uint32_t>(delta),
                                  static_cast<std::uint32_t>(size) - 1);

        touched.clear();
        while (true) {
            if (buckets.empty(minDegree)) {             // this degree is drained
                if (minDegree >= batchLimit) break;
                ++minDegree;
                continue;
            }
            // The head, whatever was filed last. On entry the walk above has left this bucket
            // non-empty, so the first pass always takes a pivot whatever the limit says; only
            // then does delta decide whether the round continues.
            const std::int32_t pivot = buckets.head(minDegree);
            buckets.unfile(degrees[pivot], pivot);

            const std::vector<std::int32_t>& merged = qg.eliminate(pivot);
            pivots.push_back(pivot);
            numEliminated += 1 + static_cast<std::uint32_t>(merged.size());
            for (std::int32_t u : merged) {
                buckets.unfile(degrees[u], u);          // unfile before zeroing: the bucket
                degrees[u] = 0;                         //   index is read from the degree
            }
            degrees[pivot] = 0;

            // Evict everything the pivot reached, with a stale degree. This is what keeps the
            // batch independent, and it is also what makes the deferred refresh safe: a vertex
            // with a stale degree is not a candidate, because it is not in a bucket to be found
            // in.
            const std::int32_t* pivotClique = qg.clique(pivot);
            const std::uint32_t pivotCliqueSize = qg.cliqueSize(pivot);
            for (std::uint32_t k = 0; k < pivotCliqueSize; ++k) {
                const std::int32_t u = pivotClique[k];
                buckets.unfile(degrees[u], u);
                if (touchedRound[u] != numRounds) {
                    touchedRound[u] = numRounds;
                    touched.push_back(u);
                }
            }

            if (delta < 0) break;                       // one pivot per round
        }

        // ---- one refresh, for everything the batch reached ------------------------------
        // The vertices are already out of their buckets, so this writes the degree and files,
        // where a one-pivot-per-step layer would refile. The eliminated filter matters because a
        // vertex evicted early in the round can be merged away by a later pivot in the same one.
        for (std::int32_t u : touched) {
            if (qg.eliminated(u)) continue;
            degrees[u] = qg.reachableSize(u);       // the size, without materializing the set
            buckets.file(degrees[u], u);
            minDegree = std::min(minDegree, degrees[u]);
        }
        ++numRounds;
    }

    return qg.order(pivots);
}


} // namespace Oblio
