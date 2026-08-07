# %% [markdown]
# # Approximate minimum degree, iteration 2: the two extras
#
# amd1 has the idea, the degree bound. This file adds the two mechanisms the
# vendored amd_2 carries beyond it, and nothing else. Both ride along with the
# bound rather than being about the degree at all, and both are cheap only because
# the bound's work has already been done. Section 5.13 of
# archive/sparse_factorization.md.
#
#   1. AGGRESSIVE ABSORPTION, which kills a clique the moment the bound work shows
#      it lies inside the new one. |C[c] - C[pivot]| has just been computed for
#      every touched clique; if it is zero, C[c] lies entirely inside C[pivot], so
#      that clique is dead. Ordinary absorption kills only the cliques the PIVOT
#      touched; this kills cliques that ANY reached vertex touched.
#
#   2. HASH SUPERVARIABLE DETECTION, which finds vertices indistinguishable from
#      EACH OTHER rather than from the pivot. Mass elimination merges into the
#      pivot, and two vertices can match each other while neither matches the
#      pivot. The hash is a filter and never the decision: a collision costs a
#      comparison, and no merge is ever missed, since equal sets give equal keys.
#
# **This file is production's Amd2, layer for layer**, and `make test` checks it by
# PERMUTATION rather than by fill. What amd3 adds beyond here is not ordering ideas
# at all: dense row detection, the pattern of A + A', the postorder and the
# Control/Info interface, none of which production has or wants.
#
# **One thing comes from amd3 and not from amd1**, and it is not optional. A hash
# merge leaves the merged vertex in place with weight zero rather than removing it
# from every list, so the walks must skip eliminated vertices. amd1 has no such
# vertices and its core takes no `eliminated` argument; this file's does. That is
# the prototype's version of what production calls live merges.
#
# The other fork from md5. Section 5.13 of archive/sparse_factorization.md.
#
# md5 has the quotient graph, supervariables, maintained degrees and buckets, and
# returns exactly md1's ordering. What is left costing anything is the refresh
# itself, which for each reached vertex u unites the members of every clique in
# I[u] and counts the result. That union is the expensive object.
#
# MMD made the refresh RARE. AMD makes each one CHEAP, and the two are the same
# answer reached from opposite ends: do the expensive thing less.
#
# THE BOUND. Rather than uniting the cliques, sum their separate contributions:
#
#   degree(u) <= min( n - k - weight(u),                nothing exceeds what remains
#                     degree_old[u] + |C[p] - {u}|,     it can only grow by the new clique
#                     |A[u] - C[p]| + |C[p] - {u}|
#                                   + sum |C[c] - C[p]| )   over c in I[u] - {p}
#
# where p is the pivot, so C[p] is the new clique, k is the count of original vertices
# eliminated so far, and weight(u) is the size of u's supervariable. The third line
# OVERCOUNTS, because two cliques may overlap outside C[p] and the overlap is counted
# twice. So it is an upper bound, not the degree.
#
# WHY THAT IS FAST, which is the entire point and is easy to miss. The quantity
# |C[c] - C[p]| depends only on the clique c, not on the vertex u, so it is
# computed ONCE PER CLIQUE and then read by every vertex whose incidence list
# holds c. The exact degree costs, per vertex, a walk over the members of all its
# cliques. The bound costs, per vertex, one addition per clique. Both are counted
# below, and the gap widens with the size of the cliques, which is to say with the
# amount of fill, which is to say exactly where it matters.
#
# WHAT IS GIVEN UP, and it is a different kind of loss from mmd1's. Every layer up
# to here picks a true minimum-degree vertex and differs only in how it finds one
# or how ties fall. This one can pick the WRONG vertex outright, because an
# overcounted bound can hide the true minimum. It is the first layer whose
# heuristic changes rather than its implementation, and the first whose pivot is
# not guaranteed to be minimal at all.
#
# The trace prints the exact degree beside the bound, so the gap is visible at
# every iteration, and the closing lines count how often the bound was loose.
#
# THIS FILE IS THE IDEA ALONE. Aggressive absorption, hash supervariable
# detection, the two-pass update, dense row handling and the rest of amd_1 and
# amd_2 are amd2's business, exactly as mmd1 held only the batching and mmd2 took
# the rest of genmmd.
#
# Run: python3 amd2.py       every example
#      python3 amd2.py 3     just the third

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

# I[u] cliques that contain u
# C[c] vertices that c contains

def amd2_show(A, I, C, degrees, exact, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, and both degrees, in
    the order the structure holds them. The stored value is the BOUND; the exact
    degree is printed beside it so the gap is visible at every iteration."""
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
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} "
              f"bound {degrees[u]} exact {exact[u]}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in C[c])
        print(f"  c{c}: {{{clique_members_text}}}")
    print()

def amd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots,
                   title=None):
    """Print the state arrays: degrees, buckets, min degree, members, eliminated,
    and the order so far."""
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
    buckets_text = "  ".join(f"{d}: [{' '.join(str(v) for v in buckets[d])}]"
                             for d in range(len(buckets)) if buckets[d])
    print(f"  degrees: [{degrees_text}]")
    print(f"  buckets: {buckets_text if buckets_text else 'all empty'}")
    print(f"  min degree: {min_degree}")
    print(f"  members: [{super_members_text}]")
    print(f"  eliminated: [{eliminated_text}]")
    print(f"  pivots: {pivots}")
    print(f"  order: {[u for pivot in pivots for u in super_members[pivot]]}")
    print()

def amd2_storage(A, I, C):
    """Entries actually stored, as in md5. The bound changes what goes into the
    degree cache, not what the quotient graph holds."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def amd2_neighbors(A, I, C, eliminated, mark, tag, u):
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
        if mark[v] != tag and not eliminated[v]:
            mark[v] = tag
            neighbors.append(v)
    for c in I[u]:
        for v in C[c]:
            if mark[v] != tag and not eliminated[v]:
                mark[v] = tag
                neighbors.append(v)
    return neighbors, tag

def amd2_exact_degree(A, I, C, eliminated, super_members, mark, tag, u):
    """The degree md5 would have computed: the union of A[u] with the members of
    every clique in I[u], counted in original vertices. Kept only so the trace can
    show the bound beside the truth and count how often the bound is loose.

    Set view: sum of |super_members[v]| over v in reach(u). It is the union the
    bound exists to avoid, so this function is instrumentation and nothing more."""
    neighbors, tag = amd2_neighbors(A, I, C, eliminated, mark, tag, u)
    return sum(len(super_members[v]) for v in neighbors), tag

def amd2_eliminate(A, I, C, mark, tag, eliminated, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md5_eliminate: this layer changes how a degree
    is estimated afterwards, not what an elimination does.

    Returns (neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag): as
    in md2, plus the vertices folded into the pivot by mass elimination. The middle
    three are reported for display; only neighbors is used by the caller.

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

    Mass elimination adds two more lines, and breaks the identity in the first one:
    from here C[pivot] is reach(pivot) minus what the pivot absorbed.

        merged   = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
        C[pivot] = C[pivot] - merged
    """
    neighbors, tag = amd2_neighbors(A, I, C, eliminated, mark, tag, pivot)
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
            # A hash merge leaves the merged vertex in place with weight zero rather
            # than removing it from every list, so a dead vertex can still sit in
            # A[u]. Dropping it here matters beyond tidiness: the mass elimination
            # test below asks whether A[u] is EMPTY, and a stale entry makes it
            # answer no. Without this line the pivot count is higher and the
            # permutation the same, since a vertex missed by mass elimination is
            # simply eliminated on its own later.
            if eliminated[v]:
                continue
            kept.append(v)
        A[u] = kept                     # what survives is A[u] - C[pivot] - {pivot}

        kept = [c for c in I[u] if mark[c] != absorbed_tag]   # I[u] - I[pivot]
        kept.append(pivot)              # u joins the new clique, whose id is the pivot
        I[u] = kept

    # Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    # the same closed neighborhood, amd2_neighbors(u) | {u} == amd2_neighbors(pivot)
    # | {pivot}, as it stood before the iteration. Equivalently, now that the clique is
    # formed, when everything u can still reach lies inside it. The test below is
    # a cheap sufficient condition for that: nothing explicit left and no clique
    # but the new one means u sees exactly what the pivot sees, so eliminating it
    # next would cost no fill. Fold it into the pivot now and strip it from the
    # cliques, since it is no longer a vertex.
    # merged = { u in C[pivot] : A[u] == {} and I[u] == {pivot} }
    merged_vertices = []
    for u in neighbors:
        if not A[u] and len(I[u]) == 1 and I[u][0] == pivot:
            I[u] = []
            eliminated[u] = True
            merged_vertices.append(u)
    if merged_vertices:                 # C[pivot] - merged, in one compaction pass
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]

    A[pivot] = []
    I[pivot] = []
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag

def amd2_file(buckets, filed, d, u):
    """Push u at the head of bucket d, which is the O(1) end of the C++ twin's
    linked list.

    Set view: buckets[d] is the set of live vertices whose current degree is d, and
    the three functions here are add, discard and move between two of them. A
    linked list gives all three in O(1) and gives the head in O(1) too, which is
    everything the picker asks of it. What it does not give is a minimum, which is
    why min_degree walks. A sorted container would hand over the minimum directly
    and charge a log on every file, and files outnumber picks."""
    buckets[d].insert(0, u)
    filed[u] = True

def amd2_unfile(buckets, filed, d, u):
    """Take u out of bucket d, if it is there. Set view: buckets[d].discard(u), and
    discard rather than remove, since a caller may unfile a vertex twice."""
    if not filed[u]:
        return
    buckets[d].remove(u)
    filed[u] = False

def amd2_refile(buckets, filed, degrees, u, new_degree):
    """Move u from the bucket for its old degree to the one for new_degree. Set
    view: buckets[old].discard(u) then buckets[new].add(u)."""
    amd2_unfile(buckets, filed, degrees[u], u)
    degrees[u] = new_degree
    amd2_file(buckets, filed, new_degree, u)

def amd2_minimum_degree(G):
    """Same as md5, with the exact refresh replaced by the approximate bound.
    Everything else, the quotient graph, mass elimination, the buckets, is md5's."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n
    # The input is given as sets, so sort once here to match the C++ literals.
    # After this nothing is sorted: the order is whatever the structure produces.
    A = [sorted(adjacency) for adjacency in G]    # explicit vertex neighbors
    I = [[] for _ in range(n)]                 # cliques that contain each vertex
    C = {}                                     # clique id -> member list
    # Twice n: vertices are stamped below n and cliques at c + n, so the exact
    # comparison in the hash pass cannot confuse a vertex with a clique.
    mark = [-1] * (2 * n)                      # scratch for membership, with tag
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
    super_members = [[u] for u in range(n)]    # the vertices each pivot stands for
    eliminated = [False] * n
    pivots = []                                # the order over supervariables
    num_eliminated_vertices = 0                         # a counter, not a scan of eliminated
    # Live ORIGINAL vertices, which is not n - num_eliminated_vertices. num_eliminated_vertices
    # counts what has left the SELECTION, and a hash merge folds v into a LIVE u,
    # so v stops being selectable while the vertices it stands for are still live
    # inside u. The first cap of the bound needs the second reading, so it gets its
    # own counter: only an elimination reduces it.
    #
    # amd1 has no such counter and does not need one, since it has no hash merges
    # and every increment of num_eliminated_vertices really is an original leaving. This
    # file inherited that line along with the driver, which made the cap too tight
    # and drove the bound BELOW the true degree, 22 times on a 10 by 10 grid. The
    # vendored amd_2 is the oracle here: its nel is advanced only at the pivot and
    # at mass elimination, never at the hash merge, which moves the weight with
    # Nv[i] += Nv[j] and leaves nel alone.
    num_live = n
    nnz_L = 0

    # The cache, as in md5, except that from the first elimination it holds a
    # BOUND rather than a degree. exact[] is carried alongside for the trace only.
    degrees = [len(A[u]) for u in range(n)]
    exact = list(degrees)
    # Only the updates are counted. The total, including the initial pass over all
    # n vertices, is that plus n, so the report derives it. That first pass finds
    # |A[u]| with no clique yet formed, which is the bound formula on an empty
    # clique set and so is exact; the bound becomes a bound from the first
    # elimination on.
    num_bound_updates = 0
    # Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    # is here as the witness that the guard is inert rather than as a statistic.
    num_tag_sweeps = 0
    num_member_visits = 0                      # what an exact refresh would cost
    num_clique_reads = 0                       # what the bound costs instead
    # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    num_absorbed = 0                           # cliques killed aggressively
    num_hash_merges = 0                        # pairs found by the hash
    hash_bucket = [[] for _ in range(n + 1)]   # Amd.cpp's Head[hval]
    num_bound_checks = 0
    num_loose_bounds = 0
    num_bounds_below_exact = 0             # an invariant, not a measurement

    # The buckets, and min_degree, a LOWER BOUND on the current minimum degree.
    # The search starts at min_degree rather than at 0, so it never looks at
    # buckets known to be empty. The bound may lag, and the walk corrects it; what
    # it must never do is overshoot, since a vertex below it would never be seen.
    #
    # n buckets is exactly right. A live vertex counts only live neighbors, so its
    # degree is at most n - 1, and the walk stops at the first non-empty bucket,
    # which exists while anything is live.
    buckets = [[] for _ in range(n)]           # buckets[d] holds the live degree-d
    filed = [False] * n                        # whether u is in a bucket at all
    for u in range(n):
        amd2_file(buckets, filed, degrees[u], u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    if n <= SHOW_THRESHOLD:
        amd2_show(A, I, C, degrees, exact,
                  "start: every edge explicit, no clique yet, degrees exact",
                  eliminated=eliminated)
        amd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
    iteration = 0
    while num_eliminated_vertices < n:
        num_iterations += 1
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1
        pivot = buckets[min_degree][0]         # the head, whatever was filed last

        # Sweep the tag back before it can wrap. THREE sites in this layer, unlike
        # the two everywhere else, and each placed where nothing in mark is live.
        # Note the sweep fills 2n here: the hash pass stamps cliques at c + n, so
        # this is the one layer whose mark array is not n long. Not inside
        # amd2_eliminate, which holds three stamps live in turn: clique_tag and
        # absorbed_tag across the prune loop, then the merged set across the
        # C[pivot] compaction. Never observed to fire.
        if tag > TAG_CEILING:
            mark = [-1] * (2 * n)
            tag = 0
            num_tag_sweeps += 1
        neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag = amd2_eliminate(
            A, I, C, mark, tag, eliminated, pivot)
        num_eliminations += 1
        num_clique_entries += len(C[pivot])
        degree = len(neighbors)
        pivots.append(pivot)
        num_eliminated_vertices += 1 + len(merged_vertices)
        for u in merged_vertices:              # the pivot now stands for them too
            super_members[pivot] += super_members[u]
            super_members[u] = []
        num_live -= len(super_members[pivot])  # every original the pivot stands for

        amd2_unfile(buckets, filed, degrees[pivot], pivot)   # the pivot has left
        degrees[pivot] = 0
        for u in merged_vertices:               # and so have the merged vertices
            amd2_unfile(buckets, filed, degrees[u], u)
            degrees[u] = 0

        # ---- the BOUND, in place of md5's exact refresh --------------------
        # C[pivot] is the new clique. Everything it reached needs a new degree, and
        # the bound replaces the union with a sum of separate contributions.
        #
        # Set view of the three quantities, none of which is built as a set:
        #
        #     pivot_clique = C[pivot]
        #     degme        = |C[pivot]|             weighted, original vertices
        #     outside[c]   = |C[c] - C[pivot]|      ONE value per CLIQUE
        #
        # The last line is the whole idea. |C[c] - C[pivot]| depends on c and not
        # on the vertex asking, so it is computed once and read many times, where
        # the exact degree recomputes a union per vertex. mark[v] == in_clique is
        # the membership test for C[pivot]; a second tag makes touched_cliques a
        # set too, so a clique is listed once however many vertices reach it.
        # The second site, guarding one contiguous region: in_clique, seen_clique,
        # the outside[c] loop, aggressive absorption's dead_tag, and the bound
        # loop's per-vertex amd2_exact_degree calls. in_clique is stamped here and
        # still read inside the outside[c] loop, so no sweep may land between them.
        # The whole region advances the tag by about n.
        if tag > TAG_CEILING:
            mark = [-1] * (2 * n)
            tag = 0
            num_tag_sweeps += 1
        pivot_clique = list(C[pivot])
        tag += 1
        in_clique = tag                         # membership of C[pivot], one test
        for v in pivot_clique:
            mark[v] = in_clique
        degme = sum(len(super_members[v]) for v in pivot_clique)

        # |C[c] - C[pivot]| ONCE PER CLIQUE. This is the whole reason the bound is
        # cheap: the quantity depends on c alone, so every vertex whose incidence
        # list holds c reads it rather than recomputing it.
        touched_cliques = []
        tag += 1
        seen_clique = tag
        for u in pivot_clique:
            for c in I[u]:
                if c != pivot and mark[c] != seen_clique:
                    mark[c] = seen_clique
                    touched_cliques.append(c)
        outside = {}
        for c in touched_cliques:
            total = 0
            for v in C[c]:
                if mark[v] != in_clique and not eliminated[v]:
                    total += len(super_members[v])
            outside[c] = total
            num_member_visits += len(C[c])      # what an exact degree pays PER VERTEX


        # AGGRESSIVE ABSORPTION, the first of amd2's two extras. Set view:
        # dead = { c : C[c] <= C[pivot] }, the containment decided by the count
        # already computed for the bound, since |C[c] - C[pivot]| == 0 IS
        # C[c] <= C[pivot]. Then I[u] = I[u] - dead for every u in C[pivot], one
        # stamp and one compaction pass, the same shape as the absorption in the
        # eliminator. Ordinary absorption kills only what the PIVOT touched; this
        # kills what ANY reached vertex touched, and it is free because the
        # quantity was computed for the bound anyway.
        dead_cliques = [c for c in touched_cliques if outside[c] == 0]
        for c in dead_cliques:
            del C[c]
        if dead_cliques:
            tag += 1
            dead_tag = tag
            for c in dead_cliques:
                mark[c] = dead_tag
            for u in pivot_clique:
                I[u] = [c for c in I[u] if mark[c] != dead_tag]   # I[u] - dead
            num_absorbed += len(dead_cliques)

        num_left = num_live
        refreshed_vertices = pivot_clique
        for u in refreshed_vertices:
            # bound = |A[u]| + |C[pivot] - {u}| + sum |C[c] - C[pivot]| over the
            # cliques in I[u] - {pivot}, against the exact
            # |( A[u] | C[c] for c in I[u] ) - {u}|. The bound replaces the union
            # by a sum, so an overlap outside C[pivot] is counted once per clique
            # that holds it, which is exactly where it overcounts.
            explicit = sum(len(super_members[v]) for v in A[u])
            cliques = [c for c in I[u] if c != pivot]
            num_clique_reads += len(cliques)    # what the bound pays instead
            bound = explicit + degme - len(super_members[u])
            for c in cliques:
                bound += outside[c]
            bound = min(bound,
                        num_left - len(super_members[u]),
                        degrees[u] + degme - len(super_members[u]))
            # NOT PRODUCTION: instrumentation. This computes the very union the bound exists to
            # avoid, and its only purpose is to show the truth beside the estimate.
            exact_u, tag = amd2_exact_degree(A, I, C, eliminated, super_members, mark, tag, u)
            exact[u] = exact_u
            num_bound_checks += 1
            if bound > exact_u:
                num_loose_bounds += 1
            # The other direction, which is not a quality signal but an INVARIANT. A bound
            # may exceed the degree by any amount and still be a bound; falling below it is
            # the one thing it must never do, since the picker would then be told a vertex is
            # cheaper than it is. Counted because it was not: the layer measured looseness
            # only, and a cap taken from the wrong counter drove this negative 22 times on a
            # 10 by 10 grid while every example stayed green. Anything but zero here is a
            # defect. See the amd2 subsection of README.md.
            if bound < exact_u:
                num_bounds_below_exact += 1
            amd2_refile(buckets, filed, degrees, u, bound)

        # HASH SUPERVARIABLE DETECTION, the second extra. Vertices indistinguishable
        # from EACH OTHER, which the pivot test cannot see: mass elimination merges
        # a vertex into the pivot, and two vertices can match each other while
        # neither matches the pivot. Hash first so the exact comparison runs only
        # within a bucket; the hash is a filter, never the decision.
        #
        # The buckets are an array indexed by the hash value, allocated once and
        # cleared only where it was used, which is Amd.cpp's Head[hval]. A map keyed
        # by the hash would cost a log per insertion and a node per group, for a
        # quantity that is already an index into 0 .. n.
        survivors = [u for u in pivot_clique if not eliminated[u]]
        used_keys = []
        for u in survivors:
            # The hash stands for the PAIR of sets (A[u], I[u]), so equal sets
            # always collide and unequal ones usually do not. A collision costs a
            # comparison, not a wrong merge, and no merge is ever missed, since
            # equal sets give equal keys by construction.
            #
            # A SUM, because addition has no order and neither do the sets. Sorting
            # to build a key would be a log factor for nothing; Amd.cpp sums the
            # indices and reduces modulo n for the same reason. The + 1 makes index
            # zero contribute, and the stride keeps cliques and vertices apart.
            key = 0
            for v in A[u]:
                if not eliminated[v]:
                    key += v + 1
            for c in I[u]:
                key += (c + 1) * (n + 1)       # cliques and vertices must not mix
            k = key % (n + 1)
            if not hash_bucket[k]:
                used_keys.append(k)
            hash_bucket[k].append(u)
        hash_pairs = []
        for k in used_keys:
            group = hash_bucket[k]
            if len(group) < 2:
                continue
            for x in range(len(group)):
                u = group[x]
                if eliminated[u]:
                    continue
                for y in range(x + 1, len(group)):
                    v = group[y]
                    if eliminated[v]:
                        continue
                    # The exact test, which the hash only filters for:
                    #     A[u] - {v} == A[v] - {u}  and  I[u] == I[v]
                    # Decided by stamping one side and counting matches on the
                    # other, as every other membership test in this file is, so it
                    # costs one pass and no sort. Each vertex is removed from the
                    # other's adjacency first: indistinguishable vertices are
                    # adjacent to each other, so without that no pair would match.
                    # The third site, and the reason this layer has one more than
                    # the others. `other` advances once per PAIR TESTED rather than
                    # once per pass, and the pair count is quadratic in the bucket
                    # sizes with no clean bound, so a check before the pass would
                    # leave the gap between checks unbounded. Safe at the top of a
                    # pair because the previous pair's stamps are spent and
                    # hash_bucket, used_keys and eliminated are separate structures.
                    if tag > TAG_CEILING:
                        mark = [-1] * (2 * n)
                        tag = 0
                        num_tag_sweeps += 1
                    tag += 1
                    other = tag
                    size_v = 0
                    for w in A[v]:
                        if w != u and not eliminated[w]:
                            mark[w] = other
                            size_v += 1
                    for c in I[v]:
                        mark[c + n] = other    # cliques stamped past the vertices
                        size_v += 1
                    size_u = 0
                    same = True
                    for w in A[u]:
                        if w == v or eliminated[w]:
                            continue
                        size_u += 1
                        if mark[w] != other:
                            same = False
                            break
                    if same:
                        for c in I[u]:
                            size_u += 1
                            if mark[c + n] != other:
                                same = False
                                break
                    if same and size_u == size_v:
                        # v is folded into u and left exactly where it lies, with a
                        # weight of zero, which is Amd.cpp's Nv[v] = 0. Removing it
                        # from every clique and every adjacency would cost a pass
                        # over the whole structure per merge, against O(1) for this.
                        #
                        # Nothing is lost by leaving it. The merge required
                        # A[u] - {v} == A[v] - {u}, so every list holding v holds u
                        # as well: v is redundant wherever it appears, never the
                        # only way to reach anything. The walks below skip it.
                        #
                        # And no degree is recomputed. u keeps the bound just
                        # written: the two were adjacent, v leaves the graph, and an
                        # external degree excludes u's own supervariable, so u's
                        # reachable set is unchanged. What changes is u's WEIGHT.
                        super_members[u] += super_members[v]
                        super_members[v] = []
                        amd2_unfile(buckets, filed, degrees[v], v)
                        A[v] = []
                        I[v] = []
                        eliminated[v] = True
                        num_eliminated_vertices += 1
                        hash_pairs.append((u, v))
                        num_hash_merges += 1

        for k in used_keys:                    # only what was touched
            hash_bucket[k].clear()

        num_bound_updates += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])

        # A supervariable of size w is w consecutive columns of L. Its external
        # degree is what remains of the clique after the merges, since a merged
        # vertex joins the supervariable instead of neighboring it, and every
        # member left there is a live vertex standing for itself alone. The first
        # column then holds ext + w - 1 entries below its diagonal, the next
        # ext + w - 2, down to ext, and each column contributes its own diagonal.
        super_size = len(super_members[pivot])
        # Weighted, because hash detection folds a vertex into a LIVE one, so a
        # member of the new clique can stand for several original vertices, and a
        # merged one stands for none and is still lying there. amd1 can use the
        # plain count and this file cannot, which is the same inheritance the
        # num_live counter above records.
        external_degree = sum(len(super_members[v]) for v in C[pivot]
                              if not eliminated[v])
        nnz_L += (super_size * external_degree
                  + super_size * (super_size - 1) // 2
                  + super_size)

        absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
        pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
        merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
        refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
        # NOT PRODUCTION: display only. The trace is what makes these files teachable and
        # is the whole reason they exist; nothing downstream reads it.
        if n <= SHOW_THRESHOLD:
            amd2_show(A, I, C, degrees, exact,
                      (f"iteration {iteration}: eliminate {pivot} (degree {degree}, size {super_size}, "
                      f"external degree {external_degree}), "
                      f"absorbed cliques: {absorbed_cliques_text}, "
                      f"pruned edges: {pruned_edges_text}, "
                      f"merged vertices: {merged_vertices_text}, "
                      f"refreshed vertices: {refreshed_vertices_text}"),
                     eliminated=eliminated)
            amd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
        iteration += 1

    order = [u for pivot in pivots for u in super_members[pivot]]
    print(f"n = {n}, nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"iterations: {num_iterations}")
    print(f"eliminations: {num_eliminations}")
    print(f"sum of |C[p]|: {num_clique_entries}")
    print(f"bound computations: {num_bound_updates + n}, "
          f"bound updates: {num_bound_updates}, "
          f"bucket probes: {num_bucket_probes}")
    print(f"clique-member visits an exact degree would need: {num_member_visits}")
    print(f"clique reads the bound needed:                    {num_clique_reads}")
    # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    print(f"aggressively absorbed: {num_absorbed}, hash merges: {num_hash_merges}")
    # NOT PRODUCTION: instrumentation, counting how often the bound was loose.
    print(f"bound below exact {num_bounds_below_exact} times, "
          f"which must be zero")
    print(f"bound was loose {num_loose_bounds} times out of {num_bound_checks}")
    print(f"tag sweeps: {num_tag_sweeps}")
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
# two cliques that overlap outside the new one, which needs enough eliminations
# to have made several cliques and enough fill for them to intersect. Every
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
# size of what merged. And super_members ends with a hole in the middle, slot 4
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
#   python3 amd2.py grid 22
if len(sys.argv) > 2 and sys.argv[1] == "grid":
    grid_side = int(sys.argv[2])
    print(f"=== grid {grid_side}x{grid_side} (n = {grid_side * grid_side}) ===")
    amd2_minimum_degree(grid_graph(grid_side))
    sys.exit(0)

# All of them by default. To run just one, pass its number: python3 amd2.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    amd2_minimum_degree(g)
    print()
