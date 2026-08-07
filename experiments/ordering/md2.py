# %% [markdown]
# # Minimum degree, iteration 2: the quotient graph
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
# The order and the per-iteration degrees match md1 exactly: same algorithm, cheaper
# storage. What this layer does NOT yet fix is that the degree is still a full
# union every time it is asked; a cheap degree is a later layer.

# %%
import sys

# The mark array is a set and the tag names it, so a tag must never repeat: a
# repeat makes a stale stamp read as a match, which is wrong silently. The tag
# only ever climbs, so the ceiling is where it has to be swept back. Half the
# positive range of the C++ twin's int32_t, which is a pragmatic choice and not
# a derived one: nothing here stores anything but a tag, so the true ceiling is
# the type's own maximum, and the room left over is against a later layer
# wanting some of it.
TAG_CEILING = 2**30 - 1

# I[u] cliques that contain u
# C[c] vertices that c contains

def md2_show(A, I, C, mark, tag, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, in the order the
    structure holds them. Returns the advanced tag, since the degree it prints is
    recomputed through md2_neighbors."""
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
        adjacency_text = " ".join(f"{v:>{width}}" for v in A[u])
        incidence_text = " ".join(f"c{c}" for c in I[u])
        neighbors, tag = md2_neighbors(A, I, C, mark, tag, u)
        degree = len(neighbors)
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} degree {degree}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in C[c])
        print(f"  c{c}: {{{clique_members_text}}}")
    print()
    return tag

def md2_storage(A, I, C):
    """Entries actually stored. Each edge costs two, one per endpoint in A. Each
    incidence costs two as well, the clique id in I and the member in C. Watch
    the total fall monotonically; the naive graph's only rises."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def md2_neighbors(A, I, C, mark, tag, u):
    """The neighbors of live vertex u: its explicit adjacency A[u] together with
    the members of every clique that contains u, minus u itself, which the
    cliques always carry. This is George and Liu's reachable set, and it is what
    the elimination graph would hold explicitly.

    In set terms this is one line, and it is worth keeping in view because the code
    below is that line with the set taken away:

        reach(u) = ( A[u] | C[c] for every c in I[u] ) - {u}

    The mark array IS the set. mark[v] == tag is the membership test, one
    comparison; mark[v] = tag is the insertion, one store. So the union costs one
    pass per source rather than a hash per member, and nothing is allocated.
    Python could write the union directly, and earlier drafts of this file did, but
    the C++ twin cannot afford it and the two must agree line for line.

    One pass per source, with the mark array doing the deduplication, so the cost
    is linear in what is touched. Returns (neighbors, tag): nothing is sorted, and
    the order is the order the sources were walked in.
    """
    tag += 1
    neighbors = []
    mark[u] = tag                      # never its own neighbor
    for v in A[u]:
        if mark[v] != tag:
            mark[v] = tag
            neighbors.append(v)
    for c in I[u]:
        for v in C[c]:
            if mark[v] != tag:
                mark[v] = tag
                neighbors.append(v)
    return neighbors, tag

def md2_eliminate(A, I, C, mark, tag, eliminated, pivot):
    """Turn the pivot into a clique.

    Returns (neighbors, absorbed_cliques, pruned_edges, tag): the pivot's neighbor
    set, which becomes the clique and the pattern of its column of L; the cliques
    that the new one swallows; the explicit edges the new clique makes redundant;
    and the advanced tag. The middle two are reported for display; only neighbors
    is used by the caller.

    Set view of the whole function, in the order the code does it:

        C[pivot] = reach(pivot)                    absorb into C[pivot]
        C        = C - I[pivot]                    reclaim I[pivot]
        for u in C[pivot]:
            A[u] = A[u] - C[pivot] - {pivot}       prune
            I[u] = ( I[u] - I[pivot] ) | {pivot}   absorb into C[pivot], reclaim I[pivot]

    The new clique is C[pivot] and gets no name of its own, so the first line reads
    as what an elimination IS: the pivot stops being a vertex with a reachable set
    and becomes a clique holding that same set. The last line is the first two
    written on the I side, since u is in C[c] exactly when c is in I[u].

    Three set differences, and not one of them builds a set. Each is a single stamp
    of the subtrahend followed by one compaction pass over the minuend, which turns
    |A[u]| * |C[pivot]| comparisons into |A[u]| + |C[pivot]|.
    """
    neighbors, tag = md2_neighbors(A, I, C, mark, tag, pivot)
    absorbed_cliques = list(I[pivot])
    for c in absorbed_cliques:
        del C[c]
    C[pivot] = list(neighbors)      # becomes the column pattern of the pivot

    # Stamp the new clique once, and the absorbed cliques once. Membership is then
    # a comparison, and both loops below are compactions in place. clique_tag is
    # the set C[pivot] and absorbed_tag is the set I[pivot], each built in one pass
    # and then queried for free.
    tag += 1
    clique_tag = tag
    for v in neighbors:
        mark[v] = clique_tag
    tag += 1
    absorbed_tag = tag
    for c in absorbed_cliques:
        mark[c] = absorbed_tag

    pruned_edges = []
    for u in neighbors:
        kept = []
        for v in A[u]:
            if v == pivot:              # the pivot is no longer a variable
                continue
            if mark[v] == clique_tag:   # both ends inside the new clique
                if u < v:
                    pruned_edges.append((u, v))
                continue                # implicit now: delete the explicit copy
            kept.append(v)
        A[u] = kept                     # what survives is A[u] - C[pivot] - {pivot}

        kept = [c for c in I[u] if mark[c] != absorbed_tag]   # I[u] - I[pivot]
        kept.append(pivot)              # u joins the new clique, whose id is the pivot
        I[u] = kept

    A[pivot] = []
    I[pivot] = []
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges, tag

def md2_minimum_degree(G):
    """Same heuristic as md1, on the quotient graph. No fill is ever stored."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n
    # The input is given as sets, so sort once here to match the C++ literals.
    # After this nothing is sorted: the order is whatever the structure produces.
    A = [sorted(adjacency) for adjacency in G]    # explicit vertex neighbors
    I = [[] for _ in range(n)]             # cliques each vertex belongs to
    C = {}                                 # clique id -> member list
    mark = [-1] * n                        # scratch for membership, stamped with tag
    tag = 0
    # Calls to the eliminate procedure, one per pivot. Not the count of vertices
    # removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    # an iteration batches several eliminations before one degree update pass. The three
    # counts coincide only where both of those are absent.
    num_eliminations = 0
    # Summed over the eliminations, |C[p]| being the new clique AFTER the trim, so
    # in supernodal terms the update rather than the front. It is the raw reach of
    # the eliminations, undeduplicated: where a layer deduplicates, the degree
    # update count comes out below this, and the gap is what the batching saved.
    # In md2 it is nnz(L) - n, there being no mass elimination to shrink a clique.
    num_clique_entries = 0
    # Passes of the outer loop, each one a batch of eliminations followed by one
    # degree update pass. Here the batch is always a single elimination, so this
    # equals num_eliminations; from mmd1 up the two come apart.
    num_iterations = 0
    num_degree_computations = 0
    # Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    # is here as the witness that the guard is inert rather than as a statistic.
    num_tag_sweeps = 0
    eliminated = [False] * n
    order = []
    degree_sum = 0

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    tag = md2_show(A, I, C, mark, tag, "start: every edge explicit, no clique yet",
                   eliminated=eliminated)
    for iteration in range(n):
        num_iterations += 1
        # Sweep the tag back before it can wrap. Two sites in this layer, one before
        # each region that advances the tag, and each placed where nothing in mark is
        # live. Here the region is the pivot search, which calls md2_neighbors once
        # per alive vertex; every call stamps what it reads in the same call, so
        # there is nothing to erase. Never observed to fire.
        if tag >= TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        pivot, best = -1, 0
        for u in range(n):                 # the key advances the tag, so no lambda
            if eliminated[u]:
                continue
            num_degree_computations += 1
            candidate, tag = md2_neighbors(A, I, C, mark, tag, u)
            if pivot == -1 or len(candidate) < best:
                pivot, best = u, len(candidate)
        # The second site. Not inside md2_eliminate: it holds clique_tag and
        # absorbed_tag live across the whole prune loop, so a sweep in there would
        # erase marks about to be read.
        if tag >= TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        neighbors, absorbed_cliques, pruned_edges, tag = md2_eliminate(
            A, I, C, mark, tag, eliminated, pivot)
        num_eliminations += 1
        num_clique_entries += len(C[pivot])
        degree = len(neighbors)
        order.append(pivot)
        degree_sum += degree

        absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
        pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
        # NOT PRODUCTION: display only. The trace is what makes these files teachable and
        # is the whole reason they exist; nothing downstream reads it.
        tag = md2_show(A, I, C, mark, tag,
                       (f"iteration {iteration}: eliminate {pivot} (degree {degree}), "
                        f"absorbed cliques: {absorbed_cliques_text}, "
                        f"pruned edges: {pruned_edges_text}"),
                       eliminated=eliminated)

    nnz_L = degree_sum + n
    print(f"n = {n}, nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"iterations: {num_iterations}")
    print(f"eliminations: {num_eliminations}")
    print(f"sum of |C[p]|: {num_clique_entries}")
    print(f"degree computations: {num_degree_computations}")
    print(f"tag sweeps: {num_tag_sweeps}")
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
# merge test declines a genuine supervariable. At the iteration whose pivot is 0 and
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
# already has 1 as a child when 0 merges into it. The merge happens at iteration 2 of
# 5, so the run continues afterwards and the selection degree, 3 over {2, 3, 4},
# differs from the external degree, 2 over {2, 3}, with the difference being the
# weight that merged. And super_members ends with a hole in the middle, slot 4
# empty between two used ones, while no pivot equals its own iteration number. See the
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

# graph7, five vertices and six edges. The pairwise case: at the iteration whose pivot
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
