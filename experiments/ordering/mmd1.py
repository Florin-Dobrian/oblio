# %% [markdown]
# # Multiple minimum degree
#
# md5 finished the cheap wins. It has the quotient graph, supervariables,
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
# WHAT IS HERE, AND WHAT MMD2 ADDS. This file is the idea alone. Everything else
# genmmd does is deliberately left to mmd1's successor, which completes it:
#
#   - the PREPASS that numbers degree 0 and 1 vertices before the main loop,
#     leaving their neighbors' degrees stale (genmmd, the loop over head[1])
#   - mmdupd's q2h path. mmdelm stashes each reached vertex's pruned adjacency
#     count as fwd[rn] = nq+1, and mmdupd routes the nq==1 cases into a separate
#     list where it merges indistinguishable PAIRS. The merge test here catches
#     only vertices indistinguishable from the pivot, so MMD's supervariables are
#     at least as coarse as ours and sometimes coarser.
#   - OUTMATCHED marking, bwd[nd] = -maxint, which takes a vertex out of the
#     degree lists without merging it.
#   - the filing convention: MMD files at `dg - qsize[en] + 1` floored at 1, so
#     its least bucket is 1 where ours is 0, and it never uses bucket 0. Plus the
#     ncsub subscript statistic.
#
# NO WEIGHT ARRAY, for the same reason md3 through md5 have none: mass elimination
# merges only into the PIVOT, which is eliminated in the same call, so no live
# vertex ever stands for more than one original vertex, and a supervariable's size
# is len(super_members[pivot]) whenever it is wanted. mmd2 needs one, because its
# q2h merge folds a vertex into a LIVE one.
#
# TIE-BREAKS. Our buckets are index-ordered, min(buckets[min_degree]), which is
# md5's convention and the reason md1 through md5 agree. MMD's degree lists are
# linked chains prepended at head[dg], so its bucket is a stack and the winner is
# whatever was pushed last, which after construction is the highest-numbered
# vertex of that degree. There is no quality claim behind it: prepending is the
# cheap end of a linked list. We keep our convention and the orderings differ in
# ties; see the README.
#
# The tag/marker machinery with its maxint overflow reset is not modeled at all.
# It exists because the marks live in reusable integer arrays; we use sets.

#
# COMPLEXITY, AND ONE PLACE THE PYTHON PAYS MORE THAN THE C++. The goal is the
# same asymptotic cost as the vendored routines, without their coding style. Two
# things were wrong and are fixed: the driver loop counts eliminations rather than
# scanning `eliminated` (O(n) per step before, O(1) now), and the mass elimination
# block strips a merged vertex from C[pivot] alone rather than from every clique,
# which is sound because I[u] was {pivot}. On a 20 by 20 grid those two cost 14800
# and 4247 elementary steps before, against 34 and 47 after, with the real
# neighbor work at 26408.
#
# What remains is min(buckets[min_degree]), which is O(bucket size) because a
# Python set is unordered. The C++ twin does not pay it: std::set is ordered, so
# *buckets[minDegree].begin() is O(1) and matches the vendored head[dg] in cost
# while keeping our index-ordered tie-break. Closing the gap in Python would need
# a heap per bucket with lazy deletion, plus a membership set to skip stale
# entries, which is more machinery than a prototype should carry. It is the one
# documented place where the Python is asymptotically worse than the C++.

# %%
import sys

# I[u] cliques that contain u
# C[c] vertices that c contains

def mmd1_show(A, I, C, degrees, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, degrees. The degree
    shown is the stored one, never recomputed, which is the point of this layer."""
    n = len(A)
    width = len(str(max(n - 1, 0)))
    alive_vertices = [u for u in range(n) if eliminated is None or not eliminated[u]]
    num_alive_edges = sum(len(A[u]) for u in alive_vertices) // 2
    num_alive_incidences = sum(len(I[u]) for u in alive_vertices)
    num_alive_cliques = len(C)
    if title:
        print(title)
    alive_vertices_text = f"{n}" if eliminated is None else f"{len(alive_vertices)} of {n}"
    print(f"num alive vertices = {alive_vertices_text}, "
          f"num alive edges = {num_alive_edges}, "
          f"num alive cliques = {num_alive_cliques}, "
          f"storage = {2 * num_alive_edges} + {2 * num_alive_incidences} = {2 * (num_alive_edges + num_alive_incidences)}")
    for u in alive_vertices:
        adjacency_text = " ".join(f"{v:>{width}}" for v in sorted(A[u]))
        incidence_text = " ".join(f"c{c}" for c in sorted(I[u]))
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} "
              f"degree {degrees[u]}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in sorted(C[c]))
        print(f"  c{c}: {{{clique_members_text}}}")
    print()

def mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots,
                   title=None):
    """Print the state arrays: degrees, buckets, min degree, members, eliminated,
    and the order so far."""
    n = len(super_members)
    width = len(str(max(n - 1, 0)))
    if title:
        print(title)
    for u in range(n):
        if not eliminated[u]:
            status = "live"
        elif super_members[u]:
            status = "done"
        else:
            status = "merged"
        super_members_text = " ".join(f"{v:>{width}}" for v in super_members[u])
        print(f"  {u:>{width}}: members [{super_members_text}] {status}")
    degrees_text = " ".join(f"{degrees[u]:>{width}}" for u in range(n))
    super_members_text = " ".join("[" + " ".join(str(v) for v in super_members[u]) + "]"
                                  for u in range(n))
    eliminated_text = " ".join(f"{int(e):>{width}}" for e in eliminated)
    buckets_text = "  ".join(f"{d}: {{{' '.join(str(v) for v in sorted(buckets[d]))}}}"
                             for d in range(len(buckets)) if buckets[d])
    print(f"  degrees: [{degrees_text}]")
    print(f"  buckets: {buckets_text if buckets_text else 'all empty'}")
    print(f"  min degree: {min_degree}")
    print(f"  members: [{super_members_text}]")
    print(f"  eliminated: [{eliminated_text}]")
    print(f"  pivots: {pivots}")
    print(f"  order: {[u for pivot in pivots for u in super_members[pivot]]}")
    print()

def mmd1_storage(A, I, C):
    """Entries actually stored, as in md5. Batching changes when degrees are
    refreshed, not what the quotient graph holds."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def mmd1_neighbors(A, I, C, u):
    """The neighbors of live vertex u, exactly as in md5: its explicit adjacency
    A[u] together with the members of every clique that contains u, minus u."""
    neighbors = set(A[u])
    for c in I[u]:
        neighbors |= C[c]
    neighbors.discard(u)
    return neighbors

def mmd1_eliminate(A, I, C, eliminated, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md5_eliminate: this layer changes how often
    degrees are refreshed, not what an elimination does.

    Returns (neighbors, absorbed_cliques, pruned_edges, merged_vertices).
    """
    neighbors = mmd1_neighbors(A, I, C, pivot)
    absorbed_cliques = set(I[pivot])
    for c in absorbed_cliques:
        del C[c]
    C[pivot] = set(neighbors)     # becomes L_pivot, the column pattern

    pruned_edges = []
    for u in neighbors:
        redundant = A[u] & neighbors    # both ends inside the new clique
        for v in redundant:
            if u < v:
                pruned_edges.append((u, v))
        A[u] -= redundant               # implicit now: delete the explicit copy
        A[u].discard(pivot)             # the pivot is no longer a variable
        I[u] -= absorbed_cliques        # its absorbed cliques are gone
        I[u].add(pivot)                 # u joins the new clique, whose id is the pivot

    # Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    # the same closed neighborhood, mmd1_neighbors(u) | {u} == mmd1_neighbors(pivot)
    # | {pivot}, as it stood before the step. Equivalently, now that the clique is
    # formed, when everything u can still reach lies inside it. The test below is
    # a cheap sufficient condition for that: nothing explicit left and no clique
    # but the new one means u sees exactly what the pivot sees, so eliminating it
    # next would cost no fill. Fold it into the pivot now and strip it from the
    # cliques, since it is no longer a vertex.
    merged_vertices = []
    for u in sorted(neighbors):
        if not A[u] and I[u] == {pivot}:
            I[u].clear()
            eliminated[u] = True
            merged_vertices.append(u)
    for u in merged_vertices:
        C[pivot].discard(u)     # I[u] was {pivot}, so no other clique holds u

    A[pivot].clear()
    I[pivot].clear()
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges, merged_vertices

def mmd1_refile(buckets, degrees, u, new_degree):
    """Move u from the bucket for its old degree to the one for new_degree.
    Removal from the middle of a bucket must be O(1), which is why a bucket is a
    set here; the vendored codes use doubly linked lists for the same reason."""
    buckets[degrees[u]].discard(u)
    degrees[u] = new_degree
    buckets[new_degree].add(u)

def mmd1_minimum_degree(G, delta=0):
    """Multiple elimination: a batch of independent pivots per degree refresh.

    delta widens the batch to vertices within delta of the minimum degree, which
    buys still fewer refreshes for a real concession, since those vertices are not
    minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
    pivot per round, which is md5's behavior reached through this code path.
    """
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n
    A = [set(adjacency) for adjacency in G]    # explicit vertex neighbors
    I = [set() for _ in range(n)]              # cliques that contain each vertex
    C = {}                                     # clique id -> member set
    super_members = [[u] for u in range(n)]    # the vertices each pivot stands for
    eliminated = [False] * n
    pivots = []                                # the order over supervariables
    num_eliminated = 0                         # a counter, not a scan of eliminated
    nnz_L = 0

    degrees = [len(A[u]) for u in range(n)]
    num_degree_computations = n

    buckets = [set() for _ in range(n)]        # buckets[d] holds the live degree-d
    for u in range(n):
        buckets[degrees[u]].add(u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0
    num_rounds = 0                             # batches, the metric this layer adds

    mmd1_show(A, I, C, degrees, "start: every edge explicit, no clique yet",
              eliminated=eliminated)
    mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
    while num_eliminated < n:
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1

        # ---- one BATCH, no degree refreshed inside it ----------------------
        # Take pivots from buckets [min_degree, min_degree + delta]. Eviction is
        # what keeps them independent: eliminating a pivot pulls every vertex it
        # reached out of the buckets, so whatever is still filed was not reached,
        # hence is not adjacent to anything taken this round.
        batch_limit = min_degree + delta
        batch = []
        touched = set()
        while True:
            if not buckets[min_degree]:        # this degree is drained
                if min_degree >= batch_limit:
                    break
                min_degree += 1
                num_bucket_probes += 1
                continue
            pivot = min(buckets[min_degree])
            degree = degrees[pivot]
            buckets[degree].discard(pivot)

            neighbors, absorbed_cliques, pruned_edges, merged_vertices = mmd1_eliminate(
                A, I, C, eliminated, pivot)
            batch.append(pivot)
            pivots.append(pivot)
            num_eliminated += 1 + len(merged_vertices)
            for u in merged_vertices:          # the pivot now stands for them too
                super_members[pivot] += super_members[u]
                super_members[u] = []
                buckets[degrees[u]].discard(u)
                degrees[u] = 0
            degrees[pivot] = 0

            for u in C[pivot]:                 # EVICT, with a stale degree
                buckets[degrees[u]].discard(u)
                touched.add(u)

            super_size = len(super_members[pivot])
            external_degree = len(C[pivot])
            nnz_L += (super_size * external_degree
                      + super_size * (super_size - 1) // 2
                      + super_size)

            absorbed_cliques_text = ", ".join(f"c{c}" for c in sorted(absorbed_cliques)) if absorbed_cliques else "none"
            pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
            merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
            evicted_text = ", ".join(str(u) for u in sorted(C[pivot])) if C[pivot] else "none"
            print(f"round {num_rounds}: eliminate {pivot} (degree {degree}, size {super_size}, "
                  f"external degree {external_degree}), "
                  f"absorbed cliques: {absorbed_cliques_text}, pruned edges: {pruned_edges_text}, "
                  f"merged vertices: {merged_vertices_text}, evicted: {evicted_text}")
            if delta < 0:                      # one pivot per round, as md5 does
                break

        # ---- one REFRESH, for everything the batch reached -----------------
        refreshed_vertices = sorted(u for u in touched if not eliminated[u])
        for u in refreshed_vertices:
            degrees[u] = len(mmd1_neighbors(A, I, C, u))
            buckets[degrees[u]].add(u)
        num_degree_computations += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
        num_rounds += 1

        batch_text = ", ".join(str(u) for u in batch)
        refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
        mmd1_show(A, I, C, degrees,
                  (f"round {num_rounds - 1} done: batch of {len(batch)}: {batch_text}, "
                   f"refreshed: {refreshed_vertices_text}"),
                  eliminated=eliminated)
        mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)

    order = [u for pivot in pivots for u in super_members[pivot]]
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"degree computations: {num_degree_computations}, "
          f"bucket probes: {num_bucket_probes}, rounds: {num_rounds}")
    print(f"order: {order}")
    return order


# %%
# The same three graphs as md1 and md2.
#
#   graph1, a 4-cycle: eliminating any vertex forces its two neighbors
#   together, so it is the smallest graph that fills (one fill edge).
#
#      0---1          edges: 0-1 1-2 2-3 3-0
#      |   |
#      3---2
#
#   graph2, uneven degrees so the picker actually chooses; it fills twice.
#
#        0            edges: 0-1 0-2 1-3 2-4
#       / \                  3-4 3-5 4-5
#      1   2
#      |   |
#      3---4
#       \ /
#        5
#
#   graph3, twelve vertices: a path 0-1-...-11 with eight extra edges. Big
#   enough that cliques grow past two members, which is where the quotient
#   graph starts to pay, and its elimination order is not the identity.
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

# graph5, five vertices and four edges, two paths joined at 4: 2-1-4-0-3. Small
# and fill free, and here for one reason: it is the smallest graph on which md3's
# merge test declines a genuine supervariable. At the step whose pivot is 0 and
# whose clique is {4}, vertex 4 has nothing explicit left but belongs to c1 as
# well as to the new clique, so I[4] == {pivot} fails even though c1's only
# member is 4 itself and everything 4 reaches lies inside the new clique. The
# exact test md3_neighbors(A, I, C, u) <= C[pivot] would merge it. See the README
# section on mass elimination.
#
#   edges: 0-3 0-4 1-2 1-4
graph5 = [
    {3, 4},           # 0
    {2, 4},           # 1
    {1},              # 2
    {0},              # 3
    {0, 1},           # 4
]

# graph6, six vertices and eight edges. Here because one small graph carries
# three things at once. Its supervariable {0, 4} is a supernode but NOT a
# fundamental one: the elimination forest is 2 -> 1 -> 4 and 3 -> 0 -> 4, so 4
# already has 1 as a child when 0 merges into it. The merge happens at step 2 of
# 5, so the run continues afterwards and the selection degree, 3 over {2, 3, 4},
# differs from the external degree, 2 over {2, 3}, with the difference being the
# size of what merged. And super_members ends with a hole in the middle, slot 4
# empty between two used ones, while no pivot equals its own step number. See the
# README sections on mass elimination and on external degree.
#
#   edges: 0-2 0-3 0-4 1-3 2-3 2-4 2-5 3-4
graph6 = [
    {2, 3, 4},        # 0
    {3},              # 1
    {0, 3, 4, 5},     # 2
    {0, 1, 2, 4},     # 3
    {0, 2, 3},        # 4
    {2},              # 5
]

# graph7, five vertices and six edges. The pairwise case: at the step whose pivot
# is 0 and whose clique is {2, 4}, vertices 2 and 4 are indistinguishable FROM
# EACH OTHER, both reaching the same closed neighborhood, yet neither is
# absorbable into the pivot, since each still reaches 3 from outside the clique.
# No test framed against the pivot finds them, and the exact test does not help
# either: both orders are 1 0 (2 3 4). Catching such pairs needs a comparison
# between candidates, which is what amd's hashing does. See the README section on
# detecting supervariables against each other.
#
#   edges: 0-1 0-2 0-4 1-4 2-3 2-4 3-4
graph7 = [
    {1, 2, 4},        # 0
    {0, 4},           # 1
    {0, 3, 4},        # 2
    {2, 4},           # 3
    {0, 1, 2, 3},     # 4
]

examples = [("graph1", graph1), ("graph2", graph2),
            ("graph3", graph3), ("graph4", graph4),
            ("graph5", graph5), ("graph6", graph6),
            ("graph7", graph7)]

# All of them by default. To run just one, pass its number: python3 mmd1.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    mmd1_minimum_degree(g)
    print()
