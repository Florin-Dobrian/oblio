# %% [markdown]
# # Minimum degree, iteration 5: degree buckets
#
# md4 stopped recomputing degrees that could not have changed. What it left in
# place is the scan: the picker still walks every live vertex to find the
# smallest cached degree, O(n) per iteration, now over integers rather than set
# unions. Cheap, but still the only remaining O(n) per pivot.
#
# The fix is to file each supervariable in a bucket indexed by its degree, so
# the minimum can be found by walking UP from the last known minimum rather than
# looking at everything. Section 5.9 of notes/SPARSE_FACTORIZATION.md describes
# this as common ground: both MMD and AMD do it, neither invented it.
#
# Two things make the walk cheap:
#
#   - min_degree, a LOWER BOUND on the current minimum degree. Each search starts
#     at min_degree rather than at 0, and every vertex below it is known to be gone.
#   - a vertex whose degree changes must be pulled out of the middle of its old
#     bucket, so buckets need O(1) removal. That is a doubly linked list, which is
#     what MMD's fwd/bwd and AMD's Next/Last are, and what the C++ twin uses:
#     head[d], next[u], prev[u], all O(1) to push, pop and splice. The Python
#     mirrors it with a list whose position 0 is the head, so both hold the same
#     sequence at every iteration and pick the same pivot.
#
# Keeping min_degree correct is the whole of the difficulty, and it is a lower
# bound rather than the true minimum on purpose: it may lag, and the walk fixes
# it. What it must never do is overshoot, since a bucket below min_degree is
# never examined and a vertex sitting there would never be chosen.

#
# COMPLEXITY. The goal is the same asymptotic cost as the vendored routines,
# without their coding style. The containers are flat: A, I and the clique member
# lists are unsorted lists, membership comes from a mark array with a tag, and a
# bucket is a linked list. Every pass is linear in what it touches.
#
# THE SET VIEW, kept in comments. Flat containers cost the Python its set algebra,
# which is what made the earlier layers readable. So every place a set operation is
# being computed by stamping and compacting says which set operation it is, in a
# comment above the loop. The reading stays "union, difference, containment" and
# the code stays flat, which is the only arrangement that serves both the teaching
# and the twin.
#
# THE TIE-BREAK CHANGES HERE, and it is the price of the buckets. md1 through md4
# scan ascending and keep the first strict minimum, so ties go to the lowest index.
# A bucket is a linked list pushed at the head and popped at the head, so the
# winner is whatever was filed last. That is what the vendored codes do, it is O(1)
# where an index-ordered pop is not, and it means md5 and mmd1 may return a
# different permutation from md1 through md4. Different, not worse: the pivots are
# still exact minima, and only the choice among equals moves.
#
# The Python mirrors the C++ structure rather than using a set: buckets[d] is a
# list whose position 0 is the head, so push is insert(0, u), pop is buckets[d][0],
# and unfile is remove(u). Same sequence, same pivot, same trace. The Python pays
# O(bucket) for insert and remove where the C++ splices in O(1), which is the one
# remaining place it is asymptotically worse than its twin.

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

def md5_show(A, I, C, degrees, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, degrees, in the order
    the structure holds them. The degree shown is the stored one, never
    recomputed, which is the point of this layer."""
    n = len(A)
    width = len(str(max(n - 1, 0)))
    live_vertices = [u for u in range(n) if eliminated is None or not eliminated[u]]
    num_live_edges = sum(len(A[u]) for u in live_vertices) // 2
    num_live_incidences = sum(len(I[u]) for u in live_vertices)
    num_live_cliques = len(C)
    if title:
        print(title)
    live_vertices_text = f"{n}" if eliminated is None else f"{len(live_vertices)} of {n}"
    print(f"num live vertices = {live_vertices_text}, "
          f"num live edges = {num_live_edges}, "
          f"num live cliques = {num_live_cliques}, "
          f"storage = {2 * num_live_edges} + {2 * num_live_incidences} = {2 * (num_live_edges + num_live_incidences)}")
    for u in live_vertices:
        adjacency_text = " ".join(f"{v:>{width}}" for v in A[u])
        incidence_text = " ".join(f"c{c}" for c in I[u])
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} degree {degrees[u]}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in C[c])
        print(f"  c{c}: {{{clique_members_text}}}")
    print()

def md5_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots,
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

def md5_storage(A, I, C):
    """Entries actually stored, as in md4. The degree cache and the buckets hold
    no graph structure: every live vertex appears in exactly one bucket."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def md5_neighbors(A, I, C, mark, tag, u):
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

    A reach tag, about vertex u, labeling reach(u) together with u. Not about any
    clique: the cliques in I[u] are read here as sources of members, never stamped
    as ids. Consumed before this function returns, unlike the eliminator's two,
    which stay live across its whole prune loop; that is why the sweep guard may
    sit before this call and not before those.
    """
    tag += 1
    neighbors = []
    mark[u] = tag                      # never its own neighbor
    for v in A[u]:
        mark[v] = tag
        neighbors.append(v)
    for c in I[u]:
        for v in C[c]:
            if mark[v] != tag:
                mark[v] = tag
                neighbors.append(v)
    return neighbors, tag

def md5_eliminate(A, I, C, mark, tag, eliminated, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md4_eliminate: this layer changes how the
    minimum is found, not what an elimination does.

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
    neighbors, tag = md5_neighbors(A, I, C, mark, tag, pivot)
    absorbed_cliques = list(I[pivot])
    for c in absorbed_cliques:
        del C[c]
    C[pivot] = list(neighbors)      # becomes the column pattern of the pivot

    # Stamp the new clique once, and the absorbed cliques once. Membership is then
    # a comparison, and both loops below are compactions in place.
    #
    # Three tags in this eliminator. Two are md2's, about cliques, differing in
    # which SIDE of a clique they name:
    #
    #     pivot_clique_tag      about the PIVOT'S clique, labels its MEMBERS,
    #                           so stamps VERTICES
    #     absorbed_cliques_tag  about the ABSORBED cliques, labels their IDS,
    #                           so stamps CLIQUE IDS
    #
    # Each side is what one loop below needs: pruning A[u] asks whether a VERTEX is
    # in C[pivot], pruning I[u] asks whether a CLIQUE ID is one of the absorbed.
    #
    # The third is mass elimination's, further down, and it is the first tag in the
    # family on the MEMBER side alongside pivot_clique_tag: it labels the merged
    # vertices so C[pivot] can be compacted against them. Two consequences. It is
    # conditional, fired only when something merged, so this layer's eliminate
    # advance is a bound, at most 4, where md2's is exactly 3. And md2's two tags
    # could have shared one value, their sides being disjoint; here they could not,
    # because two of the three now stamp vertices.
    tag += 1
    pivot_clique_tag = tag
    for v in neighbors:
        mark[v] = pivot_clique_tag
    tag += 1
    absorbed_cliques_tag = tag
    for c in absorbed_cliques:
        mark[c] = absorbed_cliques_tag

    pruned_edges = []
    # One name for both compactions, as in the C++ twin, so the two read as the same
    # code. What the name means changes four lines apart, so each use is labeled.
    # The twin has one buffer it clears and refills, the swap handing it the list it
    # just replaced; here the name is simply rebound and there is nothing to reuse.
    for u in neighbors:
        kept = []                       # KEPT IS ADJACENCY here: A[u] - C[pivot] - {pivot}
        for v in A[u]:
            if v == pivot:              # the pivot is no longer a variable
                continue
            if mark[v] == pivot_clique_tag:   # both ends inside the new clique
                if u < v:
                    pruned_edges.append((u, v))
                continue                # implicit now: delete the explicit copy
            kept.append(v)
        A[u] = kept                     # what survives is A[u] - C[pivot] - {pivot}

        # KEPT IS INCIDENCE here: I[u] - I[pivot], plus the pivot
        kept = [c for c in I[u] if mark[c] != absorbed_cliques_tag]
        kept.append(pivot)              # u joins the new clique, whose id is the pivot
        I[u] = kept

    # Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    # the same closed neighborhood, md5_neighbors(u) | {u} == md5_neighbors(pivot)
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
        # The merged tag: about the vertices this step absorbed into
        # the pivot, labeling them directly. On the MEMBER side, like
        # pivot_clique_tag, and the only conditional advance in the file.
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]

    A[pivot] = []
    I[pivot] = []
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag

def md5_file(buckets, filed, d, u):
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

def md5_unfile(buckets, filed, d, u):
    """Take u out of bucket d, if it is there. Set view: buckets[d].discard(u), and
    discard rather than remove, since a caller may unfile a vertex twice."""
    if not filed[u]:
        return
    buckets[d].remove(u)
    filed[u] = False

def md5_refile(buckets, filed, degrees, u, new_degree):
    """Move u from the bucket for its old degree to the one for new_degree. Set
    view: buckets[old].discard(u) then buckets[new].add(u)."""
    md5_unfile(buckets, filed, degrees[u], u)
    degrees[u] = new_degree
    md5_file(buckets, filed, new_degree, u)

def md5_minimum_degree(G):
    """Same as md4, with the live vertices filed in buckets by degree. The picker
    walks up from min_degree to the first non-empty bucket instead of scanning."""
    n = len(G)
    nnz_tril_A = sum(len(G[u]) for u in range(n)) // 2 + n
    # The input is given as sets, so sort once here to match the C++ literals.
    # After this nothing is sorted: the order is whatever the structure produces.
    A = [sorted(adjacency) for adjacency in G]    # explicit vertex neighbors
    I = [[] for _ in range(n)]                 # cliques that contain each vertex
    C = {}                                     # clique id -> member list
    mark = [-1] * n                            # scratch for membership, with tag
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
    nnz_L = 0

    # The cache, as in md4, unchanged by this layer.
    degrees = [len(A[u]) for u in range(n)]
    # Only the updates are counted. The total, including the initial pass over all
    # n vertices, is that plus n, so the report derives it rather than keeping a
    # second counter that could drift from this one.
    num_degree_updates = 0
    # Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    # is here as the witness that the guard is inert rather than as a statistic.
    num_tag_sweeps = 0

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
        md5_file(buckets, filed, degrees[u], u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    if n <= SHOW_THRESHOLD:
        md5_show(A, I, C, degrees, "start: every edge explicit, no clique yet",
                 eliminated=eliminated)
        md5_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
    iteration = 0
    while num_eliminated_vertices < n:
        num_iterations += 1
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1
        pivot = buckets[min_degree][0]         # the head, whatever was filed last

        # Sweep the tag back before it can wrap. Two sites in this layer, one before
        # each region that advances the tag, and each placed where nothing in mark is
        # live. The bucket walk above spends no tag, so the first region is the
        # elimination. Not inside md5_eliminate, which holds three stamps live in
        # turn: pivot_clique_tag and absorbed_cliques_tag across the prune loop, then the merged
        # set across the C[pivot] compaction. Never observed to fire.
        if tag > TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag = md5_eliminate(
            A, I, C, mark, tag, eliminated, pivot)
        num_eliminations += 1
        num_clique_entries += len(C[pivot])
        pivots.append(pivot)
        num_eliminated_vertices += 1 + len(merged_vertices)
        for u in merged_vertices:              # the pivot now stands for them too
            super_members[pivot] += super_members[u]
            super_members[u] = []

        md5_unfile(buckets, filed, degrees[pivot], pivot)   # the pivot has left
        degrees[pivot] = 0
        for u in merged_vertices:               # and so have the merged vertices
            md5_unfile(buckets, filed, degrees[u], u)
            degrees[u] = 0

        # Only the new clique's surviving members can have a different degree.
        # Everything else has the same A, the same cliques and the same live
        # neighbors as before, so its cached value is still correct.
        # Set view: the refresh set is exactly C[pivot], because reach(u) can only
        # change when a source of it changed, and the iteration touched no source
        # outside C[pivot].
        refreshed_vertices = list(C[pivot])
        # The second site, before the degree update pass. Safe here because
        # md5_eliminate's stamps are spent and the bucket work between touches no
        # mark, and because every md5_neighbors call stamps what it reads in the
        # same call.
        if tag > TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        for u in refreshed_vertices:
            neighbors_u, tag = md5_neighbors(A, I, C, mark, tag, u)
            md5_refile(buckets, filed, degrees, u, len(neighbors_u))
        num_degree_updates += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])

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

        # NOT PRODUCTION: display only, and silent above the threshold. The trace is
        # what makes these files teachable and is the whole reason they exist; nothing
        # downstream reads it. Everything it needs is built INSIDE the guard, so a run
        # above the threshold formats nothing: the four joins and the pivot's degree are
        # per elimination, and on a grid that is work for a line nobody prints.
        if n <= SHOW_THRESHOLD:
            degree = len(neighbors)            # the reach the eliminator found
            absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
            pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
            merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
            refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
            md5_show(A, I, C, degrees,
                     (f"iteration {iteration}: eliminate {pivot} (degree {degree}, size {super_size}, "
                      f"external degree {external_degree}), "
                      f"absorbed cliques: {absorbed_cliques_text}, "
                      f"pruned edges: {pruned_edges_text}, "
                      f"merged vertices: {merged_vertices_text}, "
                      f"refreshed vertices: {refreshed_vertices_text}"),
                     eliminated=eliminated)
            md5_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
        iteration += 1

    order = [u for pivot in pivots for u in super_members[pivot]]
    print(f"n = {n}, nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"iterations: {num_iterations}")
    print(f"eliminations: {num_eliminations}")
    print(f"sum of |C[p]|: {num_clique_entries}")
    print(f"degree computations: {num_degree_updates + n}, "
          f"degree updates: {num_degree_updates}, "
          f"bucket probes: {num_bucket_probes}")
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
#   python3 md5.py grid 22
if len(sys.argv) > 2 and sys.argv[1] == "grid":
    grid_side = int(sys.argv[2])
    print(f"=== grid {grid_side}x{grid_side} (n = {grid_side * grid_side}) ===")
    md5_minimum_degree(grid_graph(grid_side))
    sys.exit(0)

# All of them by default. To run just one, pass its number: python3 md5.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    md5_minimum_degree(g)
    print()
