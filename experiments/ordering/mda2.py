# %% [markdown]
# # Minimum degree, iteration 2a: the quotient graph, with the degree bounded
#
# md2 with one line changed: the picker asks for a BOUND on the degree instead of
# the degree. Everything else, the quotient graph, the elimination, the pruning,
# the absorption, is md2's and untouched.
#
# This is a fork, not a rung. md1 has no cliques, so there is nothing to
# approximate and no choice to make; md2 is the FIRST layer where the choice
# exists, and it exists because the degree has become a union over a structure
# rather than a length. From here the ladder runs twice, once computing the degree
# exactly and once bounding it, and mda2 is the earliest point where the two can be
# compared side by side.
#
# The exact degree unites:
#
#     reach(u)  = ( A[u] | C[c] for every c in I[u] ) - {u}
#     degree(u) = |reach(u)|
#
# The bound adds instead:
#
#     bound(u)  = |A[u]| + sum over c in I[u] of ( |C[c]| - 1 )
#
# **Why the minus one, and why it is exact.** u is in C[c] for every c in I[u],
# that being what an incidence means, and reach excludes u. So |C[c] - {u}| is
# |C[c]| - 1 with no test needed.
#
# **Why A[u] contributes its whole length.** It cannot overlap any of the cliques.
# When u joined c the prune took C[c] out of A[u], and A[u] only ever shrinks
# afterwards, so A[u] and every C[c] in I[u] are disjoint for the rest of the run.
# That is md2's pruning invariant doing work it was not written for, and it means
# the only source of overcounting is one vertex lying in two cliques of I[u].
#
# **And the whole point is in the shape of the sum.** |C[c]| - 1 depends on the
# clique alone, not on u, so it is one number per clique read by every vertex that
# names it. A union cannot be decomposed that way: whether two cliques overlap is a
# different question for each u, so the union has to be redone per vertex. The
# bound turns a walk over the cliques' members into one addition per clique.
#
# So bound(u) >= degree(u) always, with equality exactly when no vertex belongs to
# two cliques of I[u]. Both are printed at every iteration, side by side, so the
# looseness is visible rather than argued. Computing the exact degree defeats the
# purpose and is here only to be looked at.
#
# **This is the PIVOT-FREE bound, and it is not the tighter one.** The other form,
# the one used from mdam2 onward and by amd1, is stated against the new clique:
#
#     bound(u) = |A[u] - C[pivot]|
#              + |C[pivot] - {u}|
#              + sum over c in I[u], c != pivot, of |C[c] - C[pivot]|
#
# The two differ in their third term and nowhere else:
#
#     PIVOT-FREE            sum over c in I[u]              of ( |C[c]| - 1 )
#     AGAINST C[pivot]      sum over c in I[u], c != pivot, of |C[c] - C[pivot]|
#
# and the second is strictly better. u is one member of C[pivot], so |C[c]| - 1 is
# |C[c] - {u}|, and {u} is a subset of C[pivot], so removing C[pivot] removes at
# least as much:
#
#     |C[c] - C[pivot]|   <=   |C[c] - {u}|   =   |C[c]| - 1
#
# The gap is exactly the double counting the new clique creates: a vertex in
# C[c] & C[pivot] other than u is already counted by |C[pivot] - {u}|, and the
# pivot-free form counts it a second time.
#
# **mda2 cannot use the tighter form, and the reason is structural.** It is stated
# against C[pivot], so it applies only to a vertex IN C[pivot]. mda2's picker
# produces a number for every live vertex at every iteration, touched or not, and for an
# untouched u there is no C[pivot] in I[u] and no group to state the bound against.
# What makes the tighter form available is MAINTAINED degrees, which narrow the
# refresh set to exactly C[pivot]'s members: that is mdam2, and it is the earliest
# file in the ladder that computes what amd1 computes.
#
# So the pivot-free bound is not a stage on the way to the other one. It is what a
# recomputing picker is left with, and mda2 is the only file in the ladder that
# uses it. See the section "Zooming in on md2" in README.md, and 5.13 of
# archive/sparse_factorization.md.

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

def mda2_show(A, I, C, mark, tag, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, in the order the
    structure holds them, with the bound and the exact degree side by side.

    The degree IS recomputed here, through mda2_neighbors, which costs the very
    union the bound exists to avoid. It is here to be looked at and appears nowhere
    in the picker, which is why the tag advances and why this function returns
    it."""
    n = len(A)
    width = len(str(max(n - 1, 0)))
    alive_vertices = [u for u in range(n) if eliminated is None or not eliminated[u]]
    num_alive_edges = sum(len(A[u]) for u in alive_vertices) // 2
    num_alive_incidences = sum(len(I[u]) for u in alive_vertices)
    num_alive_cliques = len(C)
    loose = 0
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
        neighbors, tag = mda2_neighbors(A, I, C, mark, tag, u)
        degree = len(neighbors)
        bound = mda2_bound(A, I, C, u)
        if bound > degree:
            loose += 1
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} "
              f"bound {bound} degree {degree}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in C[c])
        print(f"  c{c}: {{{clique_members_text}}}")
    if loose:
        print(f"  bound loose on {loose} of {len(alive_vertices)} live vertices")
    print()
    return tag

def mda2_storage(A, I, C):
    """Entries actually stored. Each edge costs two, one per endpoint in A. Each
    incidence costs two as well, the clique id in I and the member in C. Identical
    to md2's: bounding the degree changes what is computed, never what is kept."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def mda2_neighbors(A, I, C, mark, tag, u):
    """The reachable set of u: its explicit neighbors together with the members of
    every clique it belongs to, u itself removed, each vertex once.

    Set view:

        reach(u) = ( A[u] | C[c] for every c in I[u] ) - {u}

    One mark stamp deduplicates the union in one pass per source, so the cost is
    the total size of the sources and not a search per element.

    md2 calls this to pick a pivot AND to eliminate. Here only the eliminator calls
    it: the picker uses mda2_bound, which never opens a clique. That is the whole
    difference between the two files, and it is one call site."""
    tag += 1
    mark[u] = tag                      # never its own neighbor
    neighbors = []
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

def mda2_bound(A, I, C, u):
    """An upper bound on |reach(u)|, computed without forming reach(u).

    Term by term, as the return line computes them:

        bound(u) = |A[u]|                               -> len(A[u])
                 + sum over c in I[u] of ( |C[c]| - 1 )  -> sum of len(C[c]) - 1

    TWO terms, where the bound against C[pivot] has three. There is no designated
    clique to separate out here, so where u belongs to the new clique it is simply
    one more c in the sum. And nothing is precomputed: |C[c]| - 1 is a length read
    straight off C[c], so this form needs no array where the other needs one.

    Two facts make the terms what they are, and both are md2's doing rather than
    this layer's. u belongs to every C[c] with c in I[u], so |C[c] - {u}| is
    |C[c]| - 1 and needs no test. And A[u] is disjoint from every one of those
    cliques, because joining c pruned C[c] out of A[u] and A[u] only shrinks after,
    so the explicit term needs no correction either.

    So the ONLY overcount is a vertex lying in two cliques of I[u], counted once
    per clique. bound(u) >= |reach(u)|, with equality exactly when the cliques of
    I[u] are pairwise disjoint outside u.

    Cost is |I[u]| additions and no clique is opened, against a walk over every
    member of every clique. The saving is entirely in |C[c]| - 1 depending on c and
    not on u: one number per clique, read by every vertex that names it. A union
    cannot be shared that way, since whether two cliques overlap is a different
    question for every u.

    **This is the PIVOT-FREE form**, and the bound AGAINST C[pivot] used from mdam2
    onward is strictly tighter: it subtracts the whole new clique from each other
    clique where this subtracts only {u}, and u is one member of C[pivot]. mda2
    cannot use it, because it is stated against C[pivot] and this picker must bound
    untouched vertices too. See the header."""
    return len(A[u]) + sum(len(C[c]) - 1 for c in I[u])

def mda2_eliminate(A, I, C, mark, tag, eliminated, pivot):
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

    This file bounds the degree and this function is untouched by that, which is
    worth saying here rather than only in the header. The new clique IS the
    reachable set, so forming it needs the members and not a count and there is
    nothing to approximate. The approximation buys nothing at this call and
    everything at the picker's, which is the asymmetry the whole bounded branch
    rests on: one union per pivot either way, and the exact branch pays another per
    vertex it refreshes.
    """
    neighbors, tag = mda2_neighbors(A, I, C, mark, tag, pivot)
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

def mda2_minimum_degree(G):
    """md2's algorithm with the picker reading a bound instead of a degree.

    The order can differ from md2's, and that is the result rather than a defect:
    where the bound is loose the picker may prefer a vertex the exact degree would
    not have chosen. Whether that costs fill is what the two files are for."""
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
    # Every bound this layer computes. No split into a build and updates, because
    # nothing is maintained: the picker recomputes each candidate's bound from
    # scratch on every iteration, as md2 and md3 do with exact degrees.
    num_bound_computations = 0
    # Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    # is here as the witness that the guard is inert rather than as a statistic.
    num_tag_sweeps = 0
    eliminated = [False] * n
    order = []
    degree_sum = 0
    # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    loose_picks = 0

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    tag = mda2_show(A, I, C, mark, tag, "start: every edge explicit, no clique yet",
                    eliminated=eliminated)
    for iteration in range(n):
        num_iterations += 1
        pivot, best = -1, 0
        for u in range(n):                 # no tag advances here: the bound reads lengths
            if eliminated[u]:
                continue
            num_bound_computations += 1
            candidate_bound = mda2_bound(A, I, C, u)
            if pivot == -1 or candidate_bound < best:
                pivot, best = u, candidate_bound
        # Sweep the tag back before it can wrap. One site in this layer, unlike its
        # three neighbors in the square: mda2_bound reads lengths and takes no mark
        # or tag, and the recomputing column has no refresh, so the elimination is
        # the only region that spends a tag. Not inside mda2_eliminate, which holds
        # clique_tag and absorbed_tag live across the whole prune loop. Never
        # observed to fire.
        if tag >= TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        neighbors, absorbed_cliques, pruned_edges, tag = mda2_eliminate(
            A, I, C, mark, tag, eliminated, pivot)
        num_eliminations += 1
        num_clique_entries += len(C[pivot])
        degree = len(neighbors)
        # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
        # No union is needed for it: the eliminator has just formed reach(pivot) as
        # the new clique, so degree IS the pivot's exact degree, free of charge.
        if best > degree:
            loose_picks += 1
        order.append(pivot)
        degree_sum += degree

        absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
        pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
        # NOT PRODUCTION: display only. The trace is what makes these files teachable and
        # is the whole reason they exist; nothing downstream reads it.
        tag = mda2_show(A, I, C, mark, tag,
                        (f"iteration {iteration}: eliminate {pivot} (bound {best}, degree {degree}), "
                         f"absorbed cliques: {absorbed_cliques_text}, "
                         f"pruned edges: {pruned_edges_text}"),
                        eliminated=eliminated)

    nnz_L = degree_sum + n
    print(f"n = {n}, nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"iterations: {num_iterations}")
    print(f"eliminations: {num_eliminations}")
    print(f"sum of |C[p]|: {num_clique_entries}")
    # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    print(f"loose picks = {loose_picks} of {n}")
    print(f"bound computations: {num_bound_computations}")
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

# All of them by default. To run just one, pass its number: python3 mda2.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    mda2_minimum_degree(g)
    print()
