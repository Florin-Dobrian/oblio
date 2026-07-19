# %% [markdown]
# # Minimum degree, step 4: supervariables and mass elimination
#
# md3 stopped the graph from growing. This step stops it from being as wide, by
# noticing that eliminating a pivot often makes some of its neighbors
# INDISTINGUISHABLE from it: their whole remaining neighborhood lies inside the
# new clique, so eliminating them next costs no fill at all. Rather than let the
# picker rediscover them one at a time, we merge them into the pivot on the spot
# and give the group a WEIGHT. Section 5.5 of archive/sparse_factorization.md.
#
# What changes from md3, and it is only these three things:
#
#   - weight[p] counts how many original vertices the supervariable p stands for
#   - the degree is WEIGHTED: a neighbor of weight 3 counts 3, not 1
#   - after forming a clique, any member i with A[i] empty and C[i] == {p} is
#     merged into p (MASS ELIMINATION) and stops being a vertex
#
# The elimination order comes out over supervariables, so a final pass expands
# each one into consecutive numbers. That expansion is why merged vertices must
# be numbered together, and it is the one place this layer changes the answer
# rather than just the cost.

# %%
def reachable(A, C, cliques, i):
    """The true neighbors of live variable i: its explicit neighbors A[i]
    together with the members of every clique it belongs to, minus itself."""
    neighbors = set(A[i])
    for c in C[i]:
        neighbors |= cliques[c]
    neighbors.discard(i)
    return neighbors


def degree(A, C, cliques, weight, i):
    """External degree, WEIGHTED: each neighbor counts for the number of
    original vertices it stands for. This is the quantity md3 counted with
    len(), and the only reason it differs is that supervariables now exist."""
    return sum(weight[j] for j in reachable(A, C, cliques, i))


def storage(A, cliques):
    """Explicit endpoints plus clique members."""
    return sum(len(a) for a in A) + sum(len(m) for m in cliques.values())


def eliminate(A, C, cliques, weight, eliminated, pivot):
    """Turn the pivot into a clique, then absorb every member it makes
    indistinguishable. Returns (neighbors, absorbed, pruned, merged)."""
    neighbors = reachable(A, C, cliques, pivot)
    absorbed = set(C[pivot])
    for c in absorbed:
        del cliques[c]
    cliques[pivot] = set(neighbors)

    pruned = []
    for i in neighbors:
        redundant = A[i] & neighbors     # both ends inside the new clique
        for j in redundant:
            if i < j:
                pruned.append((i, j))
        A[i] -= redundant
        A[i].discard(pivot)
        C[i] -= absorbed
        C[i].add(pivot)

    # MASS ELIMINATION. i is indistinguishable from the pivot when everything it
    # still sees lies inside the new clique: nothing explicit left, and no other
    # clique to reach through. Eliminating it next would create no fill, so we
    # merge it into the pivot now and let the weight remember how many there are.
    merged = []
    for i in sorted(neighbors):
        if not A[i] and C[i] == {pivot}:
            weight[pivot] += weight[i]
            weight[i] = 0
            C[i].clear()
            eliminated[i] = True
            merged.append(i)
    for i in merged:                     # a merged vertex is no longer a vertex
        cliques[pivot].discard(i)
        for c in cliques.values():
            c.discard(i)

    A[pivot].clear()
    C[pivot].clear()
    eliminated[pivot] = True
    return neighbors, absorbed, pruned, merged


def minimum_degree(graph):
    """Quotient-graph minimum degree with supervariables."""
    n = len(graph)
    nnz_tril_A = sum(len(neigh) for neigh in graph) // 2 + n
    A = [set(neigh) for neigh in graph]
    C = [set() for _ in range(n)]
    cliques = {}
    weight = [1] * n                  # original vertices per supervariable
    members = [[v] for v in range(n)]  # which ones, for the final expansion
    eliminated = [False] * n
    pivots = []                       # supervariable order
    nnz_L = 0

    def show_state(neighbors=None, absorbed=None, pruned=None, merged=None):
        a = "{" + ", ".join(f"{v}: {sorted(A[v])}" for v in range(n)) + "}"
        c = "{" + ", ".join(
            f"{v}: {[f'c{x}' for x in sorted(C[v])]}" for v in range(n)) + "}"
        cl = "{" + ", ".join(
            f"c{x}: {sorted(m)}" for x, m in sorted(cliques.items())) + "}"
        w = "{" + ", ".join(f"{v}: {weight[v]}" for v in range(n)) + "}"
        print(f"         A       = {a}")
        print(f"         C       = {c}")
        print(f"         cliques = {cl}")
        print(f"         weights = {w}   storage {storage(A, cliques)}")
        if neighbors is None:
            print("         neighbors = none, absorbed = none, pruned = none, "
                  "merged = none   (nothing eliminated yet)")
        else:
            ab = ", ".join(f"c{x}" for x in sorted(absorbed)) if absorbed else "none"
            pr = ", ".join(f"{u}-{w2}" for u, w2 in sorted(pruned)) if pruned else "none"
            mg = ", ".join(map(str, merged)) if merged else "none"
            print(f"         neighbors = {sorted(neighbors)}, absorbed = {ab}, "
                  f"pruned = {pr}, merged = {mg}")

    print("start: no cliques yet, every neighbor explicit, every weight 1")
    show_state()
    step = 0
    while any(not eliminated[v] for v in range(n)):
        pivot = min((v for v in range(n) if not eliminated[v]),
                    key=lambda v: degree(A, C, cliques, weight, v))
        d = degree(A, C, cliques, weight, pivot)
        w = weight[pivot]
        neighbors, absorbed, pruned, merged = eliminate(
            A, C, cliques, weight, eliminated, pivot)
        for i in merged:
            members[pivot] += members[i]
        pivots.append(pivot)

        # A supervariable of weight w is w consecutive columns of L. Its
        # EXTERNAL degree is what remains of the clique after the merges, since
        # a merged vertex joins the supervariable instead of neighboring it.
        # The first column then has ext + w - 1 entries below its diagonal, the
        # next ext + w - 2, down to ext, and each contributes its own diagonal.
        wp = weight[pivot]
        ext = sum(weight[j] for j in cliques[pivot])
        nnz_L += wp * ext + wp * (wp - 1) // 2 + wp

        print(f"step {step}: eliminate {pivot} (degree {d}, weight {w} -> {wp})")
        show_state(neighbors, absorbed, pruned, merged)
        step += 1

    order = [v for p in pivots for v in members[p]]
    print(f"supervariable pivots = {pivots}")
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    return order


# %%
# The same three graphs as md1 and md3.
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

for name, g in [("graph1", graph1), ("graph2", graph2),
                ("graph3", graph3), ("graph4", graph4)]:
    print(f"=== {name} ===")
    order = minimum_degree(g)
    print("order:", order)
    print()
