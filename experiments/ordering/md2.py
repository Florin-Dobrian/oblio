# %% [markdown]
# # Minimum degree, step 2: the quotient graph
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
# So an elimination adds nothing and removes something. Each A[u] only ever
# shrinks, which is why this representation never needs more room than the
# original graph. Section 5.3 of archive/sparse_factorization.md.
#
# A live vertex u is stored as A[u], its remaining explicit vertex neighbors, and
# I[u], the ids of the cliques that contain u; C[c] holds the members of clique c,
# so an incidence is stored twice, once from each side, just as an edge is. The
# true neighborhood of u is the union of the two, formed only when asked.
#
# Naming: the literature calls the cliques ELEMENTS and writes A_i and E_i for
# what we call A[u] and I[u]. They are cliques; we name them for what they are.
#
# The order and the per-step degrees match md1 exactly: same algorithm, cheaper
# storage. What this layer does NOT yet fix is that the degree is still a full
# union every time it is asked; a cheap degree is a later layer.

# %%
import sys

# I[u] cliques that contain u
# C[c] vertices that c contains

def md2_show(A, I, C, title=None, eliminated=None):
    """Print a quotient graph: adjacency sets, incidence sets, cliques."""
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
        degree = len(md2_neighbors(A, I, C, u))
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} degree {degree}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in sorted(C[c]))
        print(f"  c{c}: {{{clique_members_text}}}")
    print()

def md2_storage(A, I, C):
    """Entries actually stored. Each edge costs two, one per endpoint in A. Each
    incidence costs two as well, the clique id in I and the member in C. Watch
    the total fall monotonically; the naive graph's only rises."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def md2_neighbors(A, I, C, u):
    """The neighbors of live vertex u: its explicit adjacency A[u] together with
    the members of every clique that contains u, minus u itself, which the
    cliques always carry. This is George and Liu's reachable set, and it is what
    the elimination graph would hold explicitly."""
    neighbors = set(A[u])
    for c in I[u]:
        neighbors |= C[c]
    neighbors.discard(u)
    return neighbors

def md2_eliminate(A, I, C, eliminated, pivot):
    """Turn the pivot into a clique.

    Returns (neighbors, absorbed_cliques, pruned_edges): the pivot's neighbor set,
    which becomes the clique and the pattern of its column of L; the cliques that
    the new one swallows; and the explicit edges the new clique makes redundant.
    The last two are reported for display; only neighbors is used by the caller.
    """
    neighbors = md2_neighbors(A, I, C, pivot)
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

    A[pivot].clear()
    I[pivot].clear()
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges

def md2_minimum_degree(G):
    """Same heuristic as md1, on the quotient graph. No fill is ever stored."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n
    A = [set(adjacency) for adjacency in G]    # explicit variable neighbors
    I = [set() for _ in range(n)]          # cliques each variable belongs to
    C = {}                           # clique id -> member set
    eliminated = [False] * n
    order = []
    degree_sum = 0

    md2_show(A, I, C, "start: every edge explicit, no clique yet", eliminated=eliminated)
    for step in range(n):
        pivot = min((u for u in range(n) if not eliminated[u]), key=lambda u: len(md2_neighbors(A, I, C, u)))
        neighbors, absorbed_cliques, pruned_edges = md2_eliminate(A, I, C, eliminated, pivot)
        degree = len(neighbors)
        order.append(pivot)
        degree_sum += degree

        absorbed_cliques_text = ", ".join(f"c{c}" for c in sorted(absorbed_cliques)) if absorbed_cliques else "none"
        pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
        md2_show(A, I, C,
                 (f"step {step}: eliminate {pivot} (degree {degree}), "
                  f"absorbed cliques: {absorbed_cliques_text}, pruned edges: {pruned_edges_text}"),
                 eliminated=eliminated)

    nnz_L = degree_sum + n
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
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

# All of them by default. To run just one, pass its number: python3 md2.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    md2_minimum_degree(g)
    print()
