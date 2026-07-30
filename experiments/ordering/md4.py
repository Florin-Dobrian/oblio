# %% [markdown]
# # Minimum degree, step 4: maintained degrees
#
# Every version so far has recomputed a reachable set for EVERY live vertex at
# EVERY step, just to find the smallest, then thrown all but one away. On a 3D
# grid that is roughly ten times the necessary work, and the ratio grows with n.
# Section 5.7 of archive/sparse_factorization.md.
#
# The waste is easy to see once stated: eliminating a pivot can only change the
# degrees of the vertices it REACHED. Every other vertex has the same A, the same
# cliques and the same live neighbors as before, so its degree is still whatever
# it was.
#
# So we keep a degrees[] array and refresh only the reached set. The picker then
# scans cached integers instead of building set unions. This is the first half of
# what both MMD and AMD do before their ideas diverge; md5 adds the second half,
# degree buckets, which removes the scan itself.
#
# Why refreshing the clique alone is enough, in three parts:
#
#   - PRUNING and clique membership change only for members of the new clique
#   - ABSORPTION deletes cliques the pivot belonged to, and every member of such
#     a clique is reachable from the pivot, hence in the new clique
#   - MERGING removes a vertex u, but u merged only because everything it could
#     see lay inside the new clique, so nobody outside sees u disappear
#
# Nothing else in the graph can tell that an elimination happened.

# %%
import sys

# I[u] cliques that contain u
# C[c] vertices that c contains

def md4_show(A, I, C, degrees, title=None, eliminated=None):
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

def md4_show_state(degrees, super_members, eliminated, pivots, title=None):
    """Print the state arrays: degrees, members, eliminated, and the order so far."""
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
    print(f"  degrees: [{degrees_text}]")
    print(f"  members: [{super_members_text}]")
    print(f"  eliminated: [{eliminated_text}]")
    print(f"  pivots: {pivots}")
    print(f"  order: {[u for pivot in pivots for u in super_members[pivot]]}")
    print()

def md4_storage(A, I, C):
    """Entries actually stored, as in md3. The degree cache costs one number per
    vertex, in a slot that already exists, and holds no graph structure."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def md4_neighbors(A, I, C, u):
    """The neighbors of live vertex u, exactly as in md3: its explicit adjacency
    A[u] together with the members of every clique that contains u, minus u."""
    neighbors = set(A[u])
    for c in I[u]:
        neighbors |= C[c]
    neighbors.discard(u)
    return neighbors

def md4_eliminate(A, I, C, eliminated, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md3_eliminate: this layer changes who
    recomputes a degree afterwards, not what an elimination does.

    Returns (neighbors, absorbed_cliques, pruned_edges, merged_vertices).
    """
    neighbors = md4_neighbors(A, I, C, pivot)
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
    # the same closed neighborhood, md4_neighbors(u) | {u} == md4_neighbors(pivot)
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

def md4_minimum_degree(G):
    """Same as md3, with the degrees kept in an array instead of recomputed.
    The picker reads cached integers; only the new clique's members are refreshed."""
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

    # The cache, and the count of degree computations, which is what this layer
    # exists to reduce. Built once, then touched only where it can be wrong.
    degrees = [len(A[u]) for u in range(n)]
    num_degree_computations = n

    md4_show(A, I, C, degrees, "start: every edge explicit, no clique yet",
             eliminated=eliminated)
    md4_show_state(degrees, super_members, eliminated, pivots)
    step = 0
    while num_eliminated < n:
        pivot = min((u for u in range(n) if not eliminated[u]),
                    key=lambda u: degrees[u])
        neighbors, absorbed_cliques, pruned_edges, merged_vertices = md4_eliminate(
            A, I, C, eliminated, pivot)
        degree = len(neighbors)
        pivots.append(pivot)
        num_eliminated += 1 + len(merged_vertices)
        for u in merged_vertices:              # the pivot now stands for them too
            super_members[pivot] += super_members[u]
            super_members[u] = []

        # Only the new clique's surviving members can have a different degree.
        # Everything else has the same A, the same cliques and the same live
        # neighbors as before, so its cached value is still correct.
        refreshed_vertices = sorted(C[pivot])
        for u in refreshed_vertices:
            degrees[u] = len(md4_neighbors(A, I, C, u))
        num_degree_computations += len(refreshed_vertices)
        degrees[pivot] = 0
        for u in merged_vertices:
            degrees[u] = 0

        # A supervariable of size w is w consecutive columns of L. Its external
        # degree is what remains of the clique after the merges, since a merged
        # vertex joins the supervariable instead of neighboring it, and every
        # member left there is a live vertex standing for itself alone. The first
        # column then holds ext + w - 1 entries below its diagonal, the next
        # ext + w - 2, down to ext, and each column contributes its own diagonal.
        super_size = len(super_members[pivot])
        external_degree = len(C[pivot])
        nnz_L += (super_size * external_degree
                  + super_size * (super_size - 1) // 2
                  + super_size)

        absorbed_cliques_text = ", ".join(f"c{c}" for c in sorted(absorbed_cliques)) if absorbed_cliques else "none"
        pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
        merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
        refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
        md4_show(A, I, C, degrees,
                 (f"step {step}: eliminate {pivot} (degree {degree}, size {super_size}, "
                  f"external degree {external_degree}), "
                  f"absorbed cliques: {absorbed_cliques_text}, pruned edges: {pruned_edges_text}, "
                  f"merged vertices: {merged_vertices_text}, refreshed: {refreshed_vertices_text}"),
                 eliminated=eliminated)
        md4_show_state(degrees, super_members, eliminated, pivots)
        step += 1

    order = [u for pivot in pivots for u in super_members[pivot]]
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"degree computations: {num_degree_computations}")
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

# All of them by default. To run just one, pass its number: python3 md4.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    md4_minimum_degree(g)
    print()
