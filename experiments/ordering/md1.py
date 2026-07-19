# %% [markdown]
# # Minimum degree, step 1: the smallest version
#
# Naive minimum degree, nothing else. Eliminate the vertex of least degree,
# make its neighbors a clique, repeat. The new edges are FILL: the whole point
# of the ordering is to keep them few. Section 5.1 of
# archive/sparse_factorization.md as code. We build on it later.
#
# It names each fill edge as it is created, so the ordering can be seen earning
# (or wasting) its keep, step by step.

# %%
def eliminate(graph, pivot):
    """Make the pivot's neighbors a clique, then remove the pivot.

    Returns the fill edges created: pairs of neighbors that were not already
    adjacent. Those new edges are the fill.
    """
    neighbors = set(graph[pivot])
    fill = []
    for u in neighbors:
        for w in neighbors:
            if u < w and w not in graph[u]:   # a genuinely new edge
                graph[u].add(w)
                graph[w].add(u)
                fill.append((u, w))
    for u in neighbors:
        graph[u].discard(pivot)
    graph[pivot].clear()
    return fill


def storage(graph, eliminated):
    """What the graph currently costs: one entry per edge endpoint. Compare with
    md3, where the same number falls monotonically. Here fill pushes it back up."""
    return sum(len(graph[v]) for v in range(len(graph)) if not eliminated[v])


def minimum_degree(graph):
    """Eliminate the least-degree vertex each step, naming the fill it makes."""
    n = len(graph)
    nnz_tril_A = sum(len(graph[v]) for v in range(n)) // 2 + n   # before we mutate it
    eliminated = [False] * n
    order = []
    total_fill = 0
    degree_sum = 0                # sum of pivot degrees == sum of column counts of L

    def show_state():
        live = {v: sorted(graph[v]) for v in range(n) if not eliminated[v]}
        print(f"         graph:   {live}   storage {storage(graph, eliminated)}")

    print("start: every edge explicit, no fill yet")
    show_state()
    for step in range(n):
        pivot = min((v for v in range(n) if not eliminated[v]),
                    key=lambda v: len(graph[v]))
        degree = len(graph[pivot])
        fill = eliminate(graph, pivot)
        eliminated[pivot] = True
        order.append(pivot)
        total_fill += len(fill)
        degree_sum += degree

        fill_edges = ", ".join(f"{u}-{w}" for u, w in fill) if fill else "none"
        print(f"step {step}: eliminate {pivot} (degree {degree})")
        print(f"         fill: {fill_edges}   (fill so far {total_fill})")
        show_state()

    # The degree of a pivot at elimination is the count of its column of L, so
    # the degrees already computed give nnz(L) with no extra work (Section 5.1).
    nnz_L = degree_sum + n
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}   (edges counted: {total_fill})")
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
    {1, 3},   # 0
    {0, 2},   # 1
    {1, 3},   # 2
    {0, 2},   # 3
]
graph2 = [
    {1, 2},      # 0
    {0, 3},      # 1
    {0, 4},      # 2
    {1, 4, 5},   # 3
    {2, 3, 5},   # 4
    {3, 4},      # 5
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
