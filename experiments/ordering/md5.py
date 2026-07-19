# %% [markdown]
# # Minimum degree, step 5: maintained degrees
#
# Every version so far has recomputed a reachable set for EVERY live vertex at
# EVERY step, just to find the smallest, then thrown all but one away. On a 3D
# grid that is roughly ten times the necessary work, and the ratio grows with n.
# Section 5.7 of archive/sparse_factorization.md.
#
# The waste is easy to see once stated: eliminating a pivot can only change the
# degrees of the vertices it REACHED. Every other vertex has the same A, the same
# cliques, and the same weights as before, so its degree is still whatever it was.
#
# So we keep a degrees[] array and refresh only the reached set. The picker then
# scans cached integers instead of building set unions. This is the first half of
# what both MMD and AMD do before their ideas diverge; md6 adds the second half,
# degree buckets, which removes the scan itself.
#
# Why refreshing the clique alone is enough, in three parts:
#
#   - PRUNING and clique membership change only for members of the new clique
#   - ABSORPTION deletes cliques the pivot belonged to, and every member of such
#     a clique is reachable from the pivot, hence in the new clique
#   - MERGING removes a vertex i, but i merged only because everything it could
#     see lay inside the new clique, so nobody outside sees i disappear
#
# Nothing else in the graph can tell that an elimination happened.

# %%
def reachable(A, C, cliques, i):
    """The true neighbors of live variable i."""
    neighbors = set(A[i])
    for c in C[i]:
        neighbors |= cliques[c]
    neighbors.discard(i)
    return neighbors


def weighted_degree(A, C, cliques, weight, i):
    """External degree, weighted. Called only on refresh now, not on every pick."""
    return sum(weight[j] for j in reachable(A, C, cliques, i))


def storage(A, cliques):
    return sum(len(a) for a in A) + sum(len(m) for m in cliques.values())


def eliminate(A, C, cliques, weight, eliminated, pivot):
    """Turn the pivot into a clique, absorb, prune, and merge.

    Returns (neighbors, absorbed, pruned, merged). Same as md4; the difference in
    this file is entirely in who recomputes a degree afterwards.
    """
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
    """Quotient graph, supervariables, and degrees that are maintained."""
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

    # The cache. Built once here, then touched only where it can be wrong.
    degrees = [len(A[v]) for v in range(n)]
    refreshes = n                      # count them: this is the whole point

    def show_state(neighbors=None, absorbed=None, pruned=None, merged=None,
                   refreshed=None):
        a = "{" + ", ".join(f"{v}: {sorted(A[v])}" for v in range(n)) + "}"
        c = "{" + ", ".join(
            f"{v}: {[f'c{x}' for x in sorted(C[v])]}" for v in range(n)) + "}"
        cl = "{" + ", ".join(
            f"c{x}: {sorted(m)}" for x, m in sorted(cliques.items())) + "}"
        w = "{" + ", ".join(f"{v}: {weight[v]}" for v in range(n)) + "}"
        dg = "{" + ", ".join(
            f"{v}: {degrees[v] if not eliminated[v] else '-'}" for v in range(n)) + "}"
        print(f"         A       = {a}")
        print(f"         C       = {c}")
        print(f"         cliques = {cl}")
        print(f"         weights = {w}   storage {storage(A, cliques)}")
        print(f"         degrees = {dg}")
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
            print(f"         refreshed = {rf}")

    print("start: no cliques yet, every neighbor explicit, degrees from A")
    show_state()
    step = 0
    while any(not eliminated[v] for v in range(n)):
        # The pick is now a scan over cached integers, not over set unions.
        pivot = min((v for v in range(n) if not eliminated[v]),
                    key=lambda v: degrees[v])
        d = degrees[pivot]
        w = weight[pivot]
        neighbors, absorbed, pruned, merged = eliminate(
            A, C, cliques, weight, eliminated, pivot)
        for i in merged:
            members[pivot] += members[i]
        pivots.append(pivot)

        # REFRESH, and only here. The surviving members of the new clique are
        # exactly the vertices whose degree can have changed.
        refreshed = [i for i in sorted(neighbors) if not eliminated[i]]
        for i in refreshed:
            degrees[i] = weighted_degree(A, C, cliques, weight, i)
        refreshes += len(refreshed)

        wp = weight[pivot]
        ext = sum(weight[j] for j in cliques[pivot])
        nnz_L += wp * ext + wp * (wp - 1) // 2 + wp

        print(f"step {step}: eliminate {pivot} (degree {d}, weight {w} -> {wp})")
        show_state(neighbors, absorbed, pruned, merged, refreshed)
        step += 1

    order = [v for p in pivots for v in members[p]]
    print(f"supervariable pivots = {pivots}")
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"degree computations: {refreshes}   "
          f"(md4 would do one per live vertex per step)")
    return order


# %%
# The same three graphs as md1, md3 and md4.
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
