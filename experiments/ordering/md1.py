# %% [markdown]
# # Minimum degree, iteration 1: the smallest version
#
# Naive minimum degree, nothing else. Eliminate the vertex of least degree,
# make its neighbors a clique, repeat. The new edges are FILL: the whole point
# of the ordering is to keep them few. Section 5.1 of
# notes/SPARSE_FACTORIZATION.md as code. We build on it later.
#
# It names each fill edge as it is created, so the ordering can be seen earning
# (or wasting) its keep, iteration by iteration.
#
# The adjacency of a vertex is a plain list, UNSORTED, and membership comes from a
# mark array stamped with a tag, one comparison per query. That is what the C++
# twin does and what the vendored codes do; see the README section on complexity.
# The one sort is at construction, since the input is given as sets here and as
# ascending literals there.

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

# Above this n, nothing is printed from inside the run: no initial state, no
# per-iteration trace. That output is for reading a small example by eye, and at
# any size worth calling large it is O(n) lines of O(n) each, so it is unreadable
# and slow to produce. What still prints at every size is the end of the run, the
# counters and the order, since each is O(1) lines and that is what the twin
# comparison comes down to. To watch a larger run, raise this.
SHOW_THRESHOLD = 32

def md1_show(A, title=None, eliminated=None):
    """Print a graph: adjacency lists, in the order the structure holds them."""
    n = len(A)
    width = len(str(max(n - 1, 0)))
    live_vertices = [u for u in range(n) if eliminated is None or not eliminated[u]]
    num_live_edges = sum(len(A[u]) for u in live_vertices) // 2
    if title:
        print(title)
    live_vertices_text = f"{n}" if eliminated is None else f"{len(live_vertices)} of {n}"
    print(f"num live vertices = {live_vertices_text}, "
          f"num live edges = {num_live_edges}, "
          f"storage = {2 * num_live_edges}")
    for u in live_vertices:
        adjacency_text = " ".join(f"{v:>{width}}" for v in A[u])
        print(f"  {u:>{width}}: {{{adjacency_text}}} degree {len(A[u])}")
    print()

def md1_storage(A):
    """What the graph currently costs: one entry per edge endpoint. Compare with
    md2, where the same number falls monotonically. Here fill pushes it back up."""
    return sum(len(adjacency) for adjacency in A)

def md1_eliminate(A, mark, tag, eliminated, pivot):
    """Make the pivot's neighbors a clique, then remove the pivot.

    Returns (neighbors, fill_edges, tag): the pivot's adjacency at elimination,
    which is the pattern of its column of L; the fill edges created among those
    neighbors, pairs that were not already adjacent; and the advanced tag.

    In set terms this is the elimination game itself, three lines:

        for u in A[pivot]:
            fill(u) = A[pivot] - A[u] - {u}       what was not already there
            A[u]    = ( A[u] | fill(u) ) - {pivot}
        A[pivot] = {}

    The loop below is that difference without a set: stamp A[u], then walk
    A[pivot] and keep whatever is unstamped. Two passes over lists rather than a
    hash per element, and A[pivot] is captured before the loop because u is inside
    it.

    Nothing is sorted. Membership comes from the mark array, one stamp per query,
    which is what the vendored codes do and what keeps every pass linear in what
    it touches.

    One tag per neighbor, each about one neighbor u and labeling what u already
    sees: A[u], with u and the pivot stamped alongside so they fail the test below
    and never become fill. A tag is about a vertex here, never about a set of
    vertices shared by several of them, because there is no such set in this layer:
    the fill is pairwise and each neighbor's missing edges are its own. So the set a
    stamp belongs to changes every time round the loop, and the eliminator's advance
    is len(A[pivot]), at most n - 1, which is the whole of this layer's cost in the
    tag-overflow table.

    md2 is where that collapses, and the reason is what its tags are about rather
    than how many there are. The quotient graph gives the fill a name, the clique,
    so its eliminator stamps sets that belong to the pivot rather than to each
    neighbor: the pivot's reach, the members of the clique that reach becomes, and
    the ids of the cliques it absorbs. Three, whatever the pivot's degree, against
    one per neighbor here.
    """
    neighbors = list(A[pivot])
    fill_edges = []
    for u in neighbors:
        tag += 1                          # stamp what u already sees
        for v in A[u]:
            mark[v] = tag
        mark[u] = tag                     # never adjacent to itself
        mark[pivot] = tag                 # the pivot is leaving anyway
        for v in neighbors:               # A[pivot] - A[u] - {u, pivot}, one pass
            if mark[v] != tag:            # a genuinely new edge
                mark[v] = tag
                A[u].append(v)
                if u < v:
                    fill_edges.append((u, v))
    for u in neighbors:                   # A[u] - {pivot}, compacting in place
        A[u] = [v for v in A[u] if v != pivot]
    A[pivot] = []
    eliminated[pivot] = True
    return neighbors, fill_edges, tag

def md1_minimum_degree(G):
    """Eliminate the least-degree vertex each iteration, naming the fill it makes."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n   # before we mutate it
    # The input is given as sets, so sort once here to match the C++ literals.
    # After this nothing is sorted: the order is whatever the structure produces.
    A = [sorted(adjacency) for adjacency in G]
    mark = [-1] * n                # scratch for membership, stamped with tag
    tag = 0
    # Calls to the eliminate procedure, one per pivot. Not the count of vertices
    # removed: a pivot can carry mass-merged vertices out with it, and from mmd1 up
    # an iteration batches several eliminations before one degree update pass. The three
    # counts coincide only where both of those are absent.
    num_eliminations = 0
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
    total_fill = 0
    degree_sum = 0                # sum of pivot degrees == sum of column counts of L

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    if n <= SHOW_THRESHOLD:
        md1_show(A, "start: every edge explicit, no fill yet", eliminated=eliminated)
    for iteration in range(n):
        num_iterations += 1
        # The scan asks every live vertex for its degree, so the count is the live_vertices
        # count summed over iterations, n(n+1)/2 here since exactly one vertex leaves
        # per iteration. Two things keep it from being comparable with the layers
        # above. There is no initial build to charge for, degrees being computed here
        # and nowhere else, so this starts at 0 where md4 and md5 start at n. And a
        # degree computation is len(A[u]) rather than a union over A[u] and the
        # cliques in I[u], so it is the same count of a much cheaper operation.
        live_vertices = [u for u in range(n) if not eliminated[u]]
        num_degree_computations += len(live_vertices)
        pivot = min(live_vertices, key=lambda u: len(A[u]))
        # Sweep the tag back before it can wrap. Here because nothing in mark is
        # live between eliminations: every pass inside md1_eliminate stamps what
        # it reads in the same pass, so there is nothing to erase. One elimination
        # advances the tag once per neighbor of the pivot, at most n, which is the
        # room the ceiling has to leave and does. Never observed to fire.
        if tag > TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        neighbors, fill_edges, tag = md1_eliminate(A, mark, tag, eliminated, pivot)
        num_eliminations += 1
        degree = len(neighbors)
        order.append(pivot)
        total_fill += len(fill_edges)
        degree_sum += degree

        # NOT PRODUCTION: display only. The trace is what makes these files teachable and
        # is the whole reason they exist; nothing downstream reads it.
        if n <= SHOW_THRESHOLD:
            fill_edges_text = ", ".join(f"{u}-{v}" for u, v in fill_edges) if fill_edges else "none"
            md1_show(A,
                     (f"iteration {iteration}: eliminate {pivot} (degree {degree}), "
                      f"fill edges: {fill_edges_text}, fill so far: {total_fill}"),
                     eliminated=eliminated)

    # The degree of a pivot at elimination is the count of its column of L, so
    # the degrees already computed give nnz(L) with no extra work (Section 5.1).
    nnz_L = degree_sum + n
    print(f"n = {n}, nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A} (fill edges counted: {total_fill})")
    print(f"iterations: {num_iterations}")
    print(f"eliminations: {num_eliminations}")
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



# A square grid graph, four-neighbor, for running the counters at a size the seven examples
# cannot reach. It is here rather than among them because it is not an example: nothing about it
# illustrates a mechanism, and above SHOW_THRESHOLD its trace is not printed at all.
#
# It must match the C++ twin's gridGraph exactly, vertex for vertex, or `make test` would be
# diffing two different problems.
def grid_graph(side):
    n = side * side
    graph = [[] for _ in range(n)]
    for r in range(side):
        for c in range(side):
            u = r * side + c
            if r > 0:
                graph[u].append(u - side)
            if c > 0:
                graph[u].append(u - 1)
            if c + 1 < side:
                graph[u].append(u + 1)
            if r + 1 < side:
                graph[u].append(u + side)
    return graph

examples = [("graph1", graph1), ("graph2", graph2),
            ("graph3", graph3), ("graph4", graph4),
            ("graph5", graph5), ("graph6", graph6),
            ("graph7", graph7)]

# Grid mode, spelled the same way in every layer and in both twins, so `make test` can diff at a
# size the seven examples cannot reach. That matters: two defects in amd2 left every example byte
# for byte identical while the ordering was wrong on any grid of 10 a side or more. Nothing is
# filtered: the run is silent above SHOW_THRESHOLD and prints its closing lines as always.
#
#   python3 md1.py grid 22
if len(sys.argv) > 2 and sys.argv[1] == "grid":
    grid_side = int(sys.argv[2])
    print(f"=== grid {grid_side}x{grid_side} (n = {grid_side * grid_side}) ===")
    md1_minimum_degree(grid_graph(grid_side))
    sys.exit(0)

# All of them by default. To run just one, pass its number: python3 md1.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    md1_minimum_degree(g)
    print()
