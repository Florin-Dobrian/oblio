# %% [markdown]
# # Multiple minimum degree
#
# md6 finished the cheap wins. It has the quotient graph, supervariables,
# maintained degrees and degree buckets, and it returns exactly the ordering md1
# returns, only far faster. Everything left costs something.
#
# This is the first layer that changes the ANSWER. Section 5.11 of
# archive/sparse_factorization.md.
#
# The idea, from Liu (1985), is the M in MMD. Refreshing degrees is the expensive
# step, so do it less often: eliminate a whole INDEPENDENT SET of least-degree
# vertices before refreshing anything. Non-adjacent pivots cannot disturb each
# other's degrees, so every pivot in a batch is still a true minimum-degree
# vertex when it is taken.
#
# We never search for the independent set. It falls out of the bookkeeping:
# eliminating a pivot EVICTS every vertex it reached from the degree buckets, so
# whatever is still sitting in the bucket was not reached, hence is non-adjacent
# to everything already taken this round. Draining the bucket drains an
# independent set.
#
# WHAT THIS GIVES UP, and it is not what one would guess. The pivots are exact,
# but the vertices the batch evicted are invisible for the rest of the round, so
# the choice is made among the untouched remainder rather than among all
# candidates. The batch does not pick a worse vertex, it picks a different vertex
# OF THE SAME DEGREE. Minimum degree is famously sensitive to tie-breaks, so the
# fill moves by a fraction of a percent, in either direction.
#
# Batching across connected components is free (5.4); batching within one is the
# wager. This code does not distinguish them, exactly as the vendored MMD does
# not.
#
# HOW THIS DIFFERS FROM THE VENDORED Mmd.cpp. Two algorithmic gaps, to fill later;
# neither changes the shape of what is here.
#
#   - the PREPASS that numbers degree 0 and 1 vertices before the main loop,
#     leaving their neighbors' degrees stale (genmmd, the loop over head[1])
#   - mmdupd's q2h path. mmdelm stashes each reached vertex's pruned adjacency
#     count as fwd[rn] = nq+1, and mmdupd routes the nq==1 cases into a separate
#     list where it merges indistinguishable PAIRS. Our merge test only catches
#     vertices indistinguishable from the pivot, so MMD's supervariables are at
#     least as coarse as ours and sometimes coarser.
#
# And two bookkeeping conventions, which matter only if one tries to match bucket
# indices against the vendored code, not to the ordering:
#
#   - MMD never uses bucket 0. mmdint maps degree 0 to 1, and mmdupd floors with
#     `if(dg<1)dg=1`, so its least bucket is 1 where ours is 0.
#   - MMD files at `dg - qsize[en] + 1`, an offset convention; we file at the
#     plain external degree. A monotone shift, so the pivot choice is the same,
#     except that flooring merges degree 0 and 1 into one bucket and can therefore
#     break a tie differently.
#
# The tag/marker machinery with its overflow reset, and the ncsub statistic, are
# implementation detail with no counterpart here (we use Python sets).

# %%
def reachable(A, C, cliques, i):
    """The true neighbors of live variable i."""
    neighbors = set(A[i])
    for c in C[i]:
        neighbors |= cliques[c]
    neighbors.discard(i)
    return neighbors


def weighted_degree(A, C, cliques, weight, i):
    return sum(weight[j] for j in reachable(A, C, cliques, i))


def storage(A, cliques):
    return sum(len(a) for a in A) + sum(len(m) for m in cliques.values())


def eliminate(A, C, cliques, weight, eliminated, pivot):
    """Unchanged from md6: absorb, prune, merge. The batching lives in the
    caller, not here, which is also true of the vendored code (mmdelm knows
    nothing about batches)."""
    neighbors = reachable(A, C, cliques, pivot)
    absorbed = set(C[pivot])
    for c in absorbed:
        del cliques[c]
    cliques[pivot] = set(neighbors)

    pruned = []
    for i in neighbors:
        redundant = A[i] & neighbors
        for j in redundant:
            if i < j:
                pruned.append((i, j))
        A[i] -= redundant
        A[i].discard(pivot)
        C[i] -= absorbed
        C[i].add(pivot)

    merged = []
    for i in sorted(neighbors):
        if not A[i] and C[i] == {pivot}:
            weight[pivot] += weight[i]
            weight[i] = 0
            C[i].clear()
            eliminated[i] = True
            merged.append(i)
    for i in merged:
        cliques[pivot].discard(i)
        for c in cliques.values():
            c.discard(i)

    A[pivot].clear()
    C[pivot].clear()
    eliminated[pivot] = True
    return neighbors, absorbed, pruned, merged


def mmd(graph, delta=0, verbose=True):
    """Multiple minimum degree.

    delta widens the batch to vertices within delta of the minimum degree, which
    buys still fewer refreshes for a real concession: those vertices are not
    minimal, so taking them is a worse choice and not merely a different one.
    delta = 0 keeps the batch to true minima and is what we pass.
    """
    n = len(graph)
    nnz_tril_A = sum(len(neigh) for neigh in graph) // 2 + n
    A = [set(neigh) for neigh in graph]
    C = [set() for _ in range(n)]
    cliques = {}
    weight = [1] * n
    members = [[v] for v in range(n)]
    eliminated = [False] * n
    pivots = []
    nnz_L = 0

    degrees = [len(A[v]) for v in range(n)]
    buckets = [set() for _ in range(n + 1)]
    for v in range(n):
        buckets[degrees[v]].add(v)
    mdeg = min(degrees) if n else 0
    rounds = 0
    refreshes = n

    def show(label, indent="         "):
        if not verbose:
            return
        a = "{" + ", ".join(f"{v}: {sorted(A[v])}" for v in range(n)) + "}"
        c = "{" + ", ".join(
            f"{v}: {[f'c{x}' for x in sorted(C[v])]}" for v in range(n)) + "}"
        cl = "{" + ", ".join(
            f"c{x}: {sorted(m)}" for x, m in sorted(cliques.items())) + "}"
        w = "{" + ", ".join(f"{v}: {weight[v]}" for v in range(n)) + "}"
        dg = "{" + ", ".join(
            f"{v}: {degrees[v] if not eliminated[v] else '-'}" for v in range(n)) + "}"
        bk = "{" + ", ".join(f"{d}: {sorted(b)}"
                             for d, b in enumerate(buckets) if b) + "}"
        print(f"{indent}A       = {a}")
        print(f"{indent}C       = {c}")
        print(f"{indent}cliques = {cl}")
        print(f"{indent}weights = {w}   storage {storage(A, cliques)}")
        print(f"{indent}degrees = {dg}   {label}")
        print(f"{indent}buckets = {bk}   mdeg = {mdeg}")

    if verbose:
        print("start: no cliques yet, degrees from A, vertices filed by degree")
    show("(all fresh)")

    while any(not eliminated[v] for v in range(n)):
        while mdeg <= n and not buckets[mdeg]:
            mdeg += 1
        if mdeg > n:
            break

        # ---- one BATCH ----------------------------------------------------
        # Take pivots from buckets [mdeg, mdeg+delta] without refreshing any
        # degree. Eviction inside eliminate() is what keeps them independent.
        mdlmt = mdeg + delta
        batch = []
        touched = set()
        if verbose:
            print(f"round {rounds}: mdeg = {mdeg}, batch limit = {mdlmt}")
        while True:
            while mdeg <= mdlmt and not buckets[mdeg]:
                mdeg += 1
            if mdeg > mdlmt:
                break
            pivot = min(buckets[mdeg])
            d = degrees[pivot]
            w = weight[pivot]
            buckets[d].discard(pivot)

            neighbors, absorbed, pruned, merged = eliminate(
                A, C, cliques, weight, eliminated, pivot)
            for i in merged:
                members[pivot] += members[i]
                buckets[degrees[i]].discard(i)

            # EVICT. Every reached vertex leaves its bucket with a stale degree,
            # which is what makes the rest of this bucket an independent set.
            evicted = []
            for i in neighbors:
                if not eliminated[i]:
                    buckets[degrees[i]].discard(i)
                    touched.add(i)
                    evicted.append(i)

            batch.append(pivot)
            wp = weight[pivot]
            ext = sum(weight[j] for j in cliques[pivot])
            nnz_L += wp * ext + wp * (wp - 1) // 2 + wp

            if verbose:
                mg = ", ".join(map(str, merged)) if merged else "none"
                ev = ", ".join(map(str, sorted(evicted))) if evicted else "none"
                print(f"  eliminate {pivot} (degree {d}, weight {w} -> {wp})"
                      f"   merged = {mg}   evicted = {ev}")

        # ---- one REFRESH, for everything the batch reached -----------------
        refreshed = sorted(i for i in touched if not eliminated[i])
        for i in refreshed:
            new_degree = weighted_degree(A, C, cliques, weight, i)
            buckets[degrees[i]].discard(i)
            degrees[i] = new_degree
            buckets[new_degree].add(i)
            if new_degree < mdeg:
                mdeg = new_degree
        refreshes += len(refreshed)
        pivots += batch
        rounds += 1

        if verbose:
            print(f"  batch of {len(batch)}: {batch}")
            print(f"  refreshed = {', '.join(map(str, refreshed)) if refreshed else 'none'}")
            show(f"(refreshed {len(refreshed)})", indent="  ")

    order = [v for p in pivots for v in members[p]]
    if verbose:
        print(f"rounds = {rounds}, pivots = {len(pivots)}, "
              f"average batch = {len(pivots)/max(rounds,1):.1f}")
        print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
              f"fill = {nnz_L - nnz_tril_A}")
        print(f"degree computations: {refreshes}")
    return order


# %%
# The same three graphs as md1 through md6.
#
#   graph1, a 4-cycle:        graph2, uneven degrees:
#
#      0---1                       0            edges: 0-1 0-2 1-3 2-4
#      |   |                      / \                  3-4 3-5 4-5
#      3---2                     1   2
#                                |   |
#   edges: 0-1 1-2 2-3 3-0       3---4
#                                 \ /
#                                  5
#
#   graph3, twelve vertices: a path 0-1-...-11 with eight extra edges.
#
#      edges: 0-1 0-3 0-8 1-2 1-6 1-8 2-3 2-5 3-4 4-5
#             5-6 5-9 6-7 6-10 7-8 8-9 9-10 10-11
graph1 = [
    {1, 3}, {0, 2}, {1, 3}, {0, 2},
]
graph2 = [
    {1, 2}, {0, 3}, {0, 4}, {1, 4, 5}, {2, 3, 5}, {3, 4},
]
graph3 = [
    {1, 3, 8},        # 0
    {0, 2, 6, 8},     # 1
    {1, 3, 5},        # 2
    {0, 2, 4},        # 3
    {3, 5},           # 4
    {2, 4, 6, 9},     # 5
    {1, 5, 7, 10},    # 6
    {6, 8},           # 7
    {0, 1, 7, 9},     # 8
    {5, 8, 10},       # 9
    {6, 9, 11},       # 10
    {10},             # 11
]

# graph4, eight vertices and fourteen edges. Denser than the others, and here
# for one specific reason: it is the smallest graph we could find on which AMD's
# degree BOUND is ever loose. The bound overcounts only when a vertex belongs to
# two elements that overlap outside the new one, which needs enough eliminations
# to have made several elements and enough fill for them to intersect. Every
# connected graph on five or six vertices is exact (checked exhaustively), and so
# are graph1 to graph3, so without this one the amd trace would never show the
# approximation approximating. The other layers use it as an ordinary denser test.
#
#   edges: 0-2 0-3 0-4 0-7 1-3 1-4 1-6 1-7 2-3 2-5 3-6 3-7 4-5 5-6
graph4 = [
    {2, 3, 4, 7},     # 0
    {3, 4, 6, 7},     # 1
    {0, 3, 5},        # 2
    {0, 1, 2, 6, 7},  # 3
    {0, 1, 5},        # 4
    {2, 4, 6},        # 5
    {1, 3, 5},        # 6
    {0, 1, 3},        # 7
]

for name, g in [("graph1", graph1), ("graph2", graph2),
                ("graph3", graph3), ("graph4", graph4)]:
    print(f"=== {name} ===")
    order = mmd(g)
    print("order:", order)
    print()
