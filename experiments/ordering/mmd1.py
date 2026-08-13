# %% [markdown]
# # Multiple minimum degree
#
# md5 finished the cheap wins. It has the quotient graph, supervariables,
# maintained degrees and degree buckets, and it returns exactly the ordering md1
# returns, only far faster. Everything left costs something.
#
# This is the first layer that changes the ANSWER. Section 5.11 of
# archive/sparse_factorization.md.
#
# The idea, from Liu (1985), is the M in MMD. Refreshing degrees is the expensive
# iteration, so do it less often: eliminate a whole INDEPENDENT SET of least-degree
# vertices before refreshing anything. Non-adjacent pivots cannot disturb each
# other's degrees, so every pivot in a batch is still a true minimum-degree
# vertex when it is taken.
#
# We never search for the independent set. It falls out of the bookkeeping:
# eliminating a pivot EVICTS every vertex it reached from the degree buckets, so
# whatever is still sitting in the bucket was not reached, hence is non-adjacent
# to everything already taken this iteration. Draining the bucket drains an
# independent set.
#
# WHAT THIS GIVES UP, and it is not what one would guess. The pivots are exact,
# but the vertices the batch evicted are invisible for the rest of the iteration, so
# the choice is made among the untouched remainder rather than among all
# candidates. The batch does not pick a worse vertex, it picks a different vertex
# OF THE SAME DEGREE. Minimum degree is famously sensitive to tie-breaks, so the
# fill moves by a fraction of a percent, in either direction.
#
# Batching across connected components is free (5.4); batching within one is the
# wager. This code does not distinguish them, exactly as the vendored MMD does
# not.
#
# WHAT IS HERE, AND WHAT MMD2 ADDS. This file is the idea alone. Everything else
# genmmd does is deliberately left to mmd1's successor, which completes it:
#
#   - the PREPASS that numbers degree 0 and 1 vertices before the main loop,
#     leaving their neighbors' degrees stale (genmmd, the loop over head[1])
#   - mmdupd's q2h path. mmdelm stashes each reached vertex's pruned adjacency
#     count as fwd[rn] = nq+1, and mmdupd routes the nq==1 cases into a separate
#     list where it merges indistinguishable PAIRS. The merge test here catches
#     only vertices indistinguishable from the pivot, so MMD's supervariables are
#     at least as coarse as ours and sometimes coarser.
#   - OUTMATCHED marking, bwd[nd] = -maxint, which takes a vertex out of the
#     degree lists without merging it.
#   - the filing convention: MMD files at `dg - qsize[en] + 1` floored at 1, so
#     its least bucket is 1 where ours is 0, and it never uses bucket 0. Plus the
#     ncsub subscript statistic.
#
# NO WEIGHT ARRAY, for the same reason md3 through md5 have none: mass elimination
# merges only into the PIVOT, which is eliminated in the same call, so no live
# vertex ever stands for more than one original vertex, and a supervariable's size
# is len(super_members[pivot]) whenever it is wanted. mmd2 needs one, because its
# q2h merge folds a vertex into a LIVE one.
#
# TIE-BREAKS. The buckets are linked lists pushed and popped at the head, exactly
# as MMD's are, so the winner among equal degrees is whatever was filed last. That
# is the O(1) structure rather than a quality choice, and it is why md5 and mmd1
# may order differently from md1 through md4.
#
# The tag/marker machinery with its maxint overflow reset is not modeled at all.
# It exists because the marks live in reusable integer arrays; the C++ twin's mark
# and touched_iteration arrays do the same job with an explicit tag, and the Python
# uses sets where the trace does not depend on the order.

#
# COMPLEXITY, AND ONE PLACE THE PYTHON PAYS MORE THAN THE C++. The goal is the
# same asymptotic cost as the vendored routines, without their coding style. Two
# things were wrong and are fixed: the driver loop counts eliminations rather than
# scanning `eliminated` (O(n) per iteration before, O(1) now), and the mass elimination
# block strips a merged vertex from C[pivot] alone rather than from every clique,
# which is sound because I[u] was {pivot}. On a 20 by 20 grid those two cost 14800
# and 4247 elementary iterations before, against 34 and 47 after, with the real
# neighbor work at 26408.
#
# The containers are flat in BOTH twins: A, I and the clique member lists are plain
# unsorted vectors, C is indexed by clique id, membership comes from a mark array
# with a monotone tag, and a bucket is a linked list, head[d] with next and prev
# over n, so filing, unfiling and popping are O(1). With that the C++ performs the
# same operations at the same cost as the vendored genmmd; what it does not yet
# have is genmmd's remaining features, which are mmd2's business.
#
# The Python gave up set algebra to stay a twin of that, since the two must agree
# line for line and a set union costs a hash per element the C++ cannot pay. What
# it did not give up is the set VIEW: every place a set operation is being computed
# by stamping and compacting carries the set expression it stands for in a comment,
# so the reading stays "union, difference, containment" while the code stays flat.
#
# The tie-break follows the structure: the winner is whatever was filed last, not
# the lowest index, so md5 and mmd1 may return a different permutation from md1
# through md4. Different, not worse: the pivots are still exact minima and only the
# choice among equals moves.

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

def mmd1_show(A, I, C, degrees, title=None, eliminated=None):
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

def mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots,
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

def mmd1_storage(A, I, C):
    """Entries actually stored, as in md5. Batching changes when degrees are
    refreshed, not what the quotient graph holds."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def mmd1_neighbors(A, I, C, mark, tag, u):
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

def mmd1_eliminate(A, I, C, mark, tag, eliminated, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md5_eliminate: this layer changes how often
    degrees are refreshed, not what an elimination does.

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
    neighbors, tag = mmd1_neighbors(A, I, C, mark, tag, pivot)
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
    # the same closed neighborhood, mmd1_neighbors(u) | {u} == mmd1_neighbors(pivot)
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

def mmd1_file(buckets, filed, d, u):
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

def mmd1_unfile(buckets, filed, d, u):
    """Take u out of bucket d, if it is there. Set view: buckets[d].discard(u), and
    discard rather than remove, since a caller may unfile a vertex twice."""
    if not filed[u]:
        return
    buckets[d].remove(u)
    filed[u] = False

def mmd1_refile(buckets, filed, degrees, u, new_degree):
    """Move u from the bucket for its old degree to the one for new_degree. Set
    view: buckets[old].discard(u) then buckets[new].add(u)."""
    mmd1_unfile(buckets, filed, degrees[u], u)
    degrees[u] = new_degree
    mmd1_file(buckets, filed, new_degree, u)

def mmd1_minimum_degree(G, delta=0):
    """Multiple elimination: a batch of independent pivots per degree refresh.

    delta is signed: negative means one pivot per iteration. It is compared against a
    degree and its useful range stops at n - 1, so in the C++ twin it is a
    std::int32_t, an index-like quantity rather than a count.

    delta widens the batch to vertices within delta of the minimum degree, which
    buys still fewer refreshes for a real concession, since those vertices are not
    minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
    pivot per iteration, which is md5's behavior reached through this code path.
    """
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
    super_members = [[u] for u in range(n)]    # the vertices each pivot stands for
    eliminated = [False] * n
    pivots = []                                # the order over supervariables
    num_eliminated_vertices = 0                         # a counter, not a scan of eliminated
    nnz_L = 0

    degrees = [len(A[u]) for u in range(n)]
    # Only the updates are counted. The total, including the initial pass over all
    # n vertices, is that plus n, so the report derives it rather than keeping a
    # second counter that could drift from this one.
    num_degree_updates = 0
    # Sweeps of the tag back to zero. Expected to be 0 at every size we run, so it
    # is here as the witness that the guard is inert rather than as a statistic.
    num_tag_sweeps = 0

    buckets = [[] for _ in range(n)]           # buckets[d] holds the live degree-d
    filed = [False] * n                        # whether u is in a bucket at all
    for u in range(n):
        mmd1_file(buckets, filed, degrees[u], u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0
    num_iterations = 0                             # batches, the metric this layer adds
    touched_iteration = [-1] * n                   # the iteration in which u was last evicted

    # NOT PRODUCTION: display only. The trace is what makes these files teachable and
    # is the whole reason they exist; nothing downstream reads it.
    if n <= SHOW_THRESHOLD:
        mmd1_show(A, I, C, degrees, "start: every edge explicit, no clique yet",
                  eliminated=eliminated)
        mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
    while num_eliminated_vertices < n:
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1

        # ---- one BATCH, no degree refreshed inside it ----------------------
        # Take pivots from buckets [min_degree, min_degree + delta]. Eviction is
        # what keeps them independent: eliminating a pivot pulls every vertex it
        # reached out of the buckets, so whatever is still filed was not reached,
        # hence is not adjacent to anything taken this iteration.
        #
        # Set view of the invariant the eviction maintains, where reached is the
        # union of C[p] over the pivots taken so far:
        #
        #     filed = live - reached,  so  batch & reached == {}
        #
        # No set is built for either side. Membership in filed is the filed[] flag,
        # and touched_iteration[] is the same idea one level up: it stamps the iteration a
        # vertex was evicted in, so the refresh set is accumulated without a set
        # and without a sort.
        # Clamped: a degree is at most n - 1, so a wider window would walk the
        # bucket array off its end.
        batch_limit = min(min_degree + delta, n - 1) if delta >= 0 else min_degree
        batch = []
        touched = []                           # first-touch order, no set and no sort
        while True:
            if not buckets[min_degree]:        # this degree is drained
                if min_degree >= batch_limit:
                    break
                min_degree += 1
                num_bucket_probes += 1
                continue
            pivot = buckets[min_degree][0]     # the head, whatever was filed last
            degree = degrees[pivot]
            mmd1_unfile(buckets, filed, degree, pivot)

            # Sweep the tag back before it can wrap. Two sites in this layer, one
            # before each region that advances the tag, and each placed where nothing
            # in mark is live. This one is INSIDE the batch loop rather than before
            # it, since a batch takes several pivots and each calls the eliminator.
            # Safe between eliminations because the eviction that follows stamps
            # touched_iteration and filed, which are separate arrays. Not inside
            # mmd1_eliminate, which holds three stamps live in turn: pivot_clique_tag and
            # absorbed_cliques_tag across the prune loop, then the merged set across the
            # C[pivot] compaction. Never observed to fire.
            if tag > TAG_CEILING:
                mark = [-1] * n
                tag = 0
                num_tag_sweeps += 1
            neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag = mmd1_eliminate(
                A, I, C, mark, tag, eliminated, pivot)
            num_eliminations += 1
            num_clique_entries += len(C[pivot])
            batch.append(pivot)
            pivots.append(pivot)
            num_eliminated_vertices += 1 + len(merged_vertices)
            for u in merged_vertices:          # the pivot now stands for them too
                super_members[pivot] += super_members[u]
                super_members[u] = []
                mmd1_unfile(buckets, filed, degrees[u], u)
                degrees[u] = 0
            degrees[pivot] = 0

            for u in C[pivot]:                 # EVICT, with a stale degree
                mmd1_unfile(buckets, filed, degrees[u], u)
                if touched_iteration[u] != num_iterations:   # a marker, so O(1) per eviction
                    touched_iteration[u] = num_iterations
                    touched.append(u)

            super_size = len(super_members[pivot])
            external_degree = len(C[pivot])
            nnz_L += (super_size * external_degree
                      + super_size * (super_size - 1) // 2
                      + super_size)

            # NOT PRODUCTION: display only, and silent above the threshold.
            if n <= SHOW_THRESHOLD:
                absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
                pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
                merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
                evicted_text = ", ".join(str(u) for u in C[pivot]) if C[pivot] else "none"
                print(f"iteration {num_iterations}: eliminate {pivot} (degree {degree}, size {super_size}, "
                      f"external degree {external_degree}), "
                      f"absorbed cliques: {absorbed_cliques_text}, pruned edges: {pruned_edges_text}, "
                      f"merged vertices: {merged_vertices_text}, evicted: {evicted_text}")
            if delta < 0:                      # one pivot per iteration, as md5 does
                break

        # ---- one REFRESH, for everything the batch reached -----------------
        refreshed_vertices = [u for u in touched if not eliminated[u]]
        # The second site, before the degree update pass. Safe here because the
        # batch's stamps are all spent, and because every mmd1_neighbors call stamps
        # what it reads in the same call.
        if tag > TAG_CEILING:
            mark = [-1] * n
            tag = 0
            num_tag_sweeps += 1
        for u in refreshed_vertices:
            neighbors_u, tag = mmd1_neighbors(A, I, C, mark, tag, u)
            degrees[u] = len(neighbors_u)
            mmd1_file(buckets, filed, degrees[u], u)
        num_degree_updates += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
        num_iterations += 1

        batch_text = ", ".join(str(u) for u in batch)
        refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
        # NOT PRODUCTION: display only. The trace is what makes these files teachable and
        # is the whole reason they exist; nothing downstream reads it.
        if n <= SHOW_THRESHOLD:
            mmd1_show(A, I, C, degrees,
                      (f"iteration {num_iterations - 1} done: batch of {len(batch)}: {batch_text}, "
                       f"refreshed vertices: {refreshed_vertices_text}"),
                      eliminated=eliminated)
            mmd1_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)

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
#   python3 mmd1.py grid 22
if len(sys.argv) > 2 and sys.argv[1] == "grid":
    grid_side = int(sys.argv[2])
    print(f"=== grid {grid_side}x{grid_side} (n = {grid_side * grid_side}) ===")
    mmd1_minimum_degree(grid_graph(grid_side))
    sys.exit(0)

# All of them by default. To run just one, pass its number: python3 mmd1.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    mmd1_minimum_degree(g)
    print()
