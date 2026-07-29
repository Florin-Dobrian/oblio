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
import sys

def md1_show(A, title=None, eliminated=None):
    """Print a graph: adjacency sets."""
    n = len(A)
    width = len(str(max(n - 1, 0)))
    alive_vertices = [u for u in range(n) if eliminated is None or not eliminated[u]]
    num_alive_edges = sum(len(A[u]) for u in alive_vertices) // 2
    if title:
        print(title)
    alive_vertices_text = f"{n}" if eliminated is None else f"{len(alive_vertices)} of {n}"
    print(f"num alive vertices = {alive_vertices_text}, "
          f"num alive edges = {num_alive_edges}, "
          f"storage = {2 * num_alive_edges}")
    for u in alive_vertices:
        adjacency_text = " ".join(f"{v:>{width}}" for v in sorted(A[u]))
        print(f"  {u:>{width}}: {{{adjacency_text}}} degree {len(A[u])}")
    print()

def md1_storage(A):
    """What the graph currently costs: one entry per edge endpoint. Compare with
    md2, where the same number falls monotonically. Here fill pushes it back up."""
    return sum(len(adjacency) for adjacency in A)

def md1_eliminate(A, eliminated, pivot):
    """Make the pivot's neighbors a clique, then remove the pivot.

    Returns (neighbors, fill_edges): the pivot's adjacency at elimination, which is the
    pattern of its column of L, and the fill edges created among those neighbors,
    pairs that were not already adjacent.
    """
    neighbors = set(A[pivot])
    fill_edges = []
    for u in neighbors:
        for v in neighbors:
            if u < v and v not in A[u]:   # a genuinely new edge
                A[u].add(v)
                A[v].add(u)
                fill_edges.append((u, v))
    for u in neighbors:
        A[u].discard(pivot)
    A[pivot].clear()
    eliminated[pivot] = True
    return neighbors, fill_edges

def md1_minimum_degree(G):
    """Eliminate the least-degree vertex each step, naming the fill it makes."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n   # before we mutate it
    A = [set(adjacency) for adjacency in G]
    eliminated = [False] * n
    order = []
    total_fill = 0
    degree_sum = 0                # sum of pivot degrees == sum of column counts of L

    md1_show(A, "start: every edge explicit, no fill yet", eliminated=eliminated)
    for step in range(n):
        pivot = min((u for u in range(n) if not eliminated[u]), key=lambda u: len(A[u]))
        neighbors, fill_edges = md1_eliminate(A, eliminated, pivot)
        degree = len(neighbors)
        order.append(pivot)
        total_fill += len(fill_edges)
        degree_sum += degree

        fill_edges_text = ", ".join(f"{u}-{v}" for u, v in fill_edges) if fill_edges else "none"
        md1_show(A,
             (f"step {step}: eliminate {pivot} (degree {degree}), "
              f"fill edges: {fill_edges_text}, fill so far: {total_fill}"),
             eliminated=eliminated)

    # The degree of a pivot at elimination is the count of its column of L, so
    # the degrees already computed give nnz(L) with no extra work (Section 5.1).
    nnz_L = degree_sum + n
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A} (fill edges counted: {total_fill})")
    print(f"order: {order}")
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
# weight that merged. And super_members ends with a hole in the middle, slot 4
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

# All of them by default. To run just one, pass its number: python3 md1.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    md1_minimum_degree(g)
    print()
