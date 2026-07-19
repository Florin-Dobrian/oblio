# %% [markdown]
# # Minimum degree, step 6: degree buckets
#
# md5 stopped recomputing degrees that could not have changed. What it left in
# place is the scan: the picker still walks every live vertex to find the
# smallest cached degree, O(n) per step, now over integers rather than set
# unions. Cheap, but still the only remaining O(n) per pivot.
#
# The fix is to file each supervariable in a bucket indexed by its degree, so
# the minimum can be found by walking UP from the last known minimum rather than
# looking at everything. Section 5.9 of archive/sparse_factorization.md describes
# this as common ground: both MMD and AMD do it, neither invented it.
#
# Two things make the walk cheap:
#
#   - mdeg, a LOWER BOUND on the current minimum degree. We start each search at
#     mdeg rather than at 0, and every vertex below it is known to be gone.
#   - a vertex whose degree changes must be pulled out of the middle of its old
#     bucket, so buckets need O(1) removal. Here that is a Python set; the
#     vendored codes use doubly linked lists (MMD's fwd/bwd, AMD's Next/Last)
#     because Fortran and C of that era had no such container.
#
# Keeping mdeg correct is the whole of the difficulty, and it is a lower bound
# rather than the true minimum on purpose: it may lag, and the walk fixes it.
# What it must never do is overshoot, since a bucket below mdeg is never
# examined and a vertex sitting there would never be chosen.

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
    """Unchanged from md5: absorb, prune, merge."""
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


def minimum_degree(graph):
    """Quotient graph, supervariables, maintained degrees, and degree buckets."""
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
    buckets = [set() for _ in range(n + 1)]      # buckets[d] holds degree-d vertices
    for v in range(n):
        buckets[degrees[v]].add(v)
    mdeg = min(degrees) if n else 0              # lower bound on the minimum
    probes = 0                                   # bucket slots examined, the metric

    def refile(i, new_degree):
        """Move i from its old bucket to the one for new_degree."""
        nonlocal mdeg
        buckets[degrees[i]].discard(i)
        degrees[i] = new_degree
        buckets[new_degree].add(i)
        if new_degree < mdeg:                    # the bound may only ever fall
            mdeg = new_degree

    def show_state(neighbors=None, absorbed=None, pruned=None, merged=None,
                   refreshed=None, walked=None):
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
        print(f"         A       = {a}")
        print(f"         C       = {c}")
        print(f"         cliques = {cl}")
        print(f"         weights = {w}   storage {storage(A, cliques)}")
        print(f"         degrees = {dg}")
        print(f"         buckets = {bk}   mdeg = {mdeg}")
        if neighbors is None:
            print("         neighbors = none, absorbed = none, pruned = none, "
                  "merged = none, refreshed = all (initial)")
        else:
            ab = ", ".join(f"c{x}" for x in sorted(absorbed)) if absorbed else "none"
            pr = ", ".join(f"{u}-{w2}" for u, w2 in sorted(pruned)) if pruned else "none"
            mg = ", ".join(map(str, merged)) if merged else "none"
            rf = ", ".join(map(str, sorted(refreshed))) if refreshed else "none"
            print(f"         neighbors = {sorted(neighbors)}, absorbed = {ab}, "
                  f"pruned = {pr}, merged = {mg}")
            print(f"         refreshed = {rf}   (walked {walked} bucket(s) to find the pivot)")

    print("start: no cliques yet, degrees from A, vertices filed by degree")
    show_state()
    step = 0
    while any(not eliminated[v] for v in range(n)):
        # PICK: walk up from mdeg to the first non-empty bucket. No scan.
        walked = 0
        while not buckets[mdeg]:
            mdeg += 1
            walked += 1
        probes += walked + 1
        pivot = min(buckets[mdeg])               # lowest index, to match md5's tie-break
        d = degrees[pivot]
        w = weight[pivot]
        buckets[d].discard(pivot)

        neighbors, absorbed, pruned, merged = eliminate(
            A, C, cliques, weight, eliminated, pivot)
        for i in merged:
            members[pivot] += members[i]
            buckets[degrees[i]].discard(i)       # merged vertices leave the buckets
        pivots.append(pivot)

        refreshed = [i for i in sorted(neighbors) if not eliminated[i]]
        for i in refreshed:
            refile(i, weighted_degree(A, C, cliques, weight, i))

        wp = weight[pivot]
        ext = sum(weight[j] for j in cliques[pivot])
        nnz_L += wp * ext + wp * (wp - 1) // 2 + wp

        print(f"step {step}: eliminate {pivot} (degree {d}, weight {w} -> {wp})")
        show_state(neighbors, absorbed, pruned, merged, refreshed, walked)
        step += 1

    order = [v for p in pivots for v in members[p]]
    print(f"supervariable pivots = {pivots}")
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"bucket probes: {probes}   (md5 would scan every live vertex per step)")
    return order


# %%
# The same three graphs as md1, md3, md4 and md5.
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
    order = minimum_degree(g)
    print("order:", order)
    print()
