# %% [markdown]
# # Minimum degree, step 3: the quotient graph
#
# Same ordering as md1, computed WITHOUT ever storing fill. When a vertex is
# eliminated it becomes a CLIQUE on the vertices it would have joined. A clique is
# fully described by its vertex list, so every edge inside it is implicit, and that
# cuts twice:
#
#   - the fill edges are never added, and
#   - the edges ALREADY present between two members are now redundant, so they
#     are pruned from the explicit adjacency.
#
# So an elimination adds nothing and removes something. Each A[i] only ever
# shrinks, which is why this representation never needs more room than the
# original graph. Section 5.2 of archive/sparse_factorization.md.
#
# A live variable i is stored as A[i], its remaining explicit variable neighbors,
# and C[i], the cliques it belongs to. Its true neighborhood is the union of the
# two, formed only when asked.
#
# Naming: the literature calls these objects ELEMENTS and writes E_i for C[i].
# They are cliques; we name them for what they are.
#
# The order and the per-step degrees match md1 exactly: same algorithm, cheaper
# storage. What this layer does NOT yet fix is that the degree is still a full
# union every time it is asked; a cheap degree is a later layer.

# %%
def reachable(A, C, cliques, i):
    """The true neighbors of live variable i: its explicit neighbors A[i]
    together with the members of every clique it belongs to, minus itself."""
    neighbors = set(A[i])
    for c in C[i]:
        neighbors |= cliques[c]
    neighbors.discard(i)
    return neighbors


def storage(A, cliques):
    """What the quotient graph currently costs: explicit endpoints plus clique
    members. Watch it fall. The naive graph's edge count only rises."""
    return sum(len(a) for a in A) + sum(len(m) for m in cliques.values())


def eliminate(A, C, cliques, eliminated, pivot):
    """Turn the pivot into a clique. Returns (neighbors, absorbed, pruned).

    neighbors the pivot's reachable set. Stored once as a clique instead of as
              fill edges; in that role the doc calls it the pattern L_pivot,
              which is also the nonzero pattern of column pivot of L.
    absorbed  cliques the pivot belonged to, swallowed by the new one
    pruned    explicit edges deleted because the new clique now implies them
    """
    neighbors = reachable(A, C, cliques, pivot)
    absorbed = set(C[pivot])
    for c in absorbed:
        del cliques[c]
    cliques[pivot] = set(neighbors)    # becomes L_pivot: the clique's pattern

    pruned = []
    for i in neighbors:
        redundant = A[i] & neighbors    # both ends inside the new clique
        for j in redundant:
            if i < j:
                pruned.append((i, j))
        A[i] -= redundant               # implicit now: delete the explicit copy
        A[i].discard(pivot)             # the pivot is no longer a variable
        C[i] -= absorbed                # its absorbed cliques are gone
        C[i].add(pivot)                 # i belongs to the new clique instead

    A[pivot].clear()
    C[pivot].clear()
    eliminated[pivot] = True
    return neighbors, absorbed, pruned


def minimum_degree(graph):
    """Quotient-graph minimum degree. Same result as md1, no fill stored."""
    n = len(graph)
    nnz_tril_A = sum(len(neigh) for neigh in graph) // 2 + n     # before we mutate
    A = [set(neigh) for neigh in graph]   # explicit variable neighbors
    C = [set() for _ in range(n)]         # cliques each variable belongs to
    cliques = {}                           # clique id -> its live members
    eliminated = [False] * n
    order = []
    degree_sum = 0                # sum of pivot degrees == sum of column counts of L

    def show_state(neighbors=None, absorbed=None, pruned=None):
        """A, C and cliques on their own lines; then the three results of the
        eliminate call, in the order it returns them."""
        a = "{" + ", ".join(f"{v}: {sorted(A[v])}" for v in range(n)) + "}"
        c = "{" + ", ".join(
            f"{v}: {[f'c{x}' for x in sorted(C[v])]}" for v in range(n)) + "}"
        cl = "{" + ", ".join(
            f"c{x}: {sorted(m)}" for x, m in sorted(cliques.items())) + "}"
        print(f"         A       = {a}")
        print(f"         C       = {c}")
        print(f"         cliques = {cl}   storage {storage(A, cliques)}")
        if neighbors is None:
            print("         neighbors = none, absorbed = none, pruned = none"
                  "   (nothing eliminated yet)")
        else:
            nb = sorted(neighbors)
            ab = ", ".join(f"c{x}" for x in sorted(absorbed)) if absorbed else "none"
            pr = ", ".join(f"{u}-{w}" for u, w in sorted(pruned)) if pruned else "none"
            print(f"         neighbors = {nb}, absorbed = {ab}, pruned = {pr}")

    print("start: no cliques yet, every neighbor explicit")
    show_state()
    for step in range(n):
        pivot = min((v for v in range(n) if not eliminated[v]),
                    key=lambda v: len(reachable(A, C, cliques, v)))
        degree = len(reachable(A, C, cliques, pivot))
        neighbors, absorbed, pruned = eliminate(A, C, cliques, eliminated, pivot)
        order.append(pivot)
        degree_sum += degree

        print(f"step {step}: eliminate {pivot} (degree {degree})")
        show_state(neighbors, absorbed, pruned)

    # The degree of a pivot at elimination is the count of its column of L, so
    # the degrees already computed give nnz(L) with no extra work (Section 5.1).
    nnz_L = degree_sum + n
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    return order


# %%
# Two examples.
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
