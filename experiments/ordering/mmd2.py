# %% [markdown]
# # Multiple minimum degree, complete
#
# mmd1 has the idea: a batch of independent pivots per degree refresh. What it
# does not have is the rest of what genmmd does, and this file adds it, one pass
# at a time. Section 5.11 of archive/sparse_factorization.md, plus the vendored
# routine itself in vendored/vendored_mmd.cpp.
#
# The list, against the vendored routine that carries each:
#
#   1. the PREPASS over head[1], numbering degree 0 and 1 vertices before the main
#      loop and leaving their neighbors' degrees stale                    [done]
#   2. mmdupd's q2h path, split on mmdelm's fwd[rn] = nq + 1 stash        [done]
#   3. the pairwise merge inside the q2h walk, which folds a vertex into a
#      LIVE one, so a candidate can stand for several vertices             [done]
#   4. OUTMATCHED marking, bwd[nd] = -maxint                              [done]
#   5. the filing convention, dg - qsize[en] + 1 floored at 1             [done]
#   6. the counters: ncsub, and the early termination on num + qsize      [done]
#
# PASS 1, THE PREPASS. genmmd numbers every vertex in the degree-1 list before the
# main loop starts, marks each marker[mn] = maxint, and never refreshes a
# neighbor. Two things travel with it. mmdint files a degree-0 vertex under degree
# 1, `if(dg==0)dg=1`, so isolated and degree-1 vertices are numbered together, and
# the bucket a vertex sits in stops being its true degree. And head[1] = 0
# afterwards, with the main loop starting at mdeg = 2.
#
# A prepass vertex is NOT eliminated in the quotient-graph sense: no clique is
# formed, nothing is pruned, and its neighbors keep degrees that still count it.
# It is simply numbered and then skipped, which marker[mn] = maxint does there and
# `eliminated[u]` does here. The neighbor query skips such vertices, and the prune
# loop drops them when it compacts, which is what mmdelm's `marker[nb] < tag` test
# does.
#
# PASS 2, THE q2h SPLIT. mmdupd does not walk a flat list of reached vertices. It
# walks the ELEMENTS created this round, `el = list[el]`, and for each one it
# computes dg0 once, the weighted size of that element, then visits the element's
# members. A member is classified by what it has left BESIDES the new element:
# mmdelm stashes fwd[rn] = nq + 1 where nq counts the survivors of the compaction,
# which here is len(A[u]) + len(I[u]) - 1. nq == 1 puts the vertex on the q2h list,
# anything else on qxh.
#
# The q2h case is answered without a union. Everything the vertex reaches is
# either inside the element, already counted in dg0, or comes from that one other
# source, so the walk adds only what the other source contributes. The marker
# stops the element's own members being counted twice, which is what mt does in
# mmdupd. The qxh case does the full union, as before.
#
# The degrees come out identical either way, which is the check on this pass. What
# does move is the ORDER of the filing, since the refresh is now element by element
# with q2h before qxh, and filing order decides what a bucket holds.
#
# A vertex reached by two pivots in the same round is refreshed once: mmdupd skips
# it on the second visit with `if(bwd[en]!=0) goto n2200`, since a refiled vertex
# has a bucket again. Here that is the `filed` flag.
#
# PASS 3, THE PAIRWISE MERGE AND OUTMATCHED MARKING. Both live in the same branch
# of the q2h walk, reached when a member of the one other source is ALSO a member
# of the new element:
#
#   else if(bwd[nd]==0){
#       if(fwd[nd]==2){qsize[en]+=qsize[nd];qsize[nd]=0;marker[nd]=maxint;
#                      fwd[nd]=-en;bwd[nd]=-maxint;}
#       else if(bwd[nd]==0)bwd[nd]=-maxint;}
#
# MERGE. If nd is q2h too, its only other source is that same element, so en and nd
# reach exactly the same vertices and are indistinguishable. en absorbs nd. This is
# the first merge in the whole sequence that folds a vertex into a LIVE one, which
# matters here: from here a candidate can stand for several original vertices, and
# every degree has to count them rather than count entries. It is also what makes
# MMD's supervariables coarser than md3's, whose test only ever compares a vertex
# against the pivot.
#
# The count comes from len(super_members[v]), which is O(1), so there is still no
# weight array. MMD keeps qsize because its members are a chain, not a list.
#
# OUTMATCHED. If nd is not q2h it has other sources besides these two, so its reach
# contains en's. It can never be the minimum before en, and MMD withdraws it from
# the degree lists rather than refiling it: bwd[nd] = -maxint. It is not merged and
# not eliminated, just held out until something reaches it again, at which point
# mmdelm restores it with bwd[rn] = 0. Here that is the `outmatched` flag, cleared
# in mmd2_eliminate for every vertex the new clique reaches.
#
# PASS 4, THE FILING CONVENTION. mmdupd does not file a vertex under its degree.
# It files under `dg = dg - qsize[en] + 1`, floored at 1, where dg was the weighted
# reach INCLUDING en's own members. So the bucket index is the external degree plus
# one, and the floor catches the case where a vertex reaches nothing outside itself.
#
# mmdint, meanwhile, files at the plain degree, `dg = xadj[nd+1] - xadj[nd]`, with
# only the zero case lifted to 1. So MMD runs on two scales: the initial buckets
# hold degrees, every refiled bucket holds degree + 1. That is genuine, not a
# misreading, and it tilts the pivot choice slightly against refreshed vertices,
# which sit one bucket higher than an untouched vertex of the same reach.
#
# From here degrees[] holds the FILED value, which is what the picker compares and
# what min_degree tracks. The nnz(L) accounting does not use it: that sums weights
# over the live members of C[pivot] and is unaffected.
#
# PASS 5, THE COUNTERS. Two small things in genmmd's main loop.
#
# ncsub, `*ncsub += mdeg + qsize[mn] - 2`, accumulated per pivot. It is the
# statistic genmmd returns alongside the permutation, an estimate of the subscript
# storage the factor will need, and it is computed from the values the loop already
# has rather than from anything extra.
#
# The early termination, `if((num+qsize[mn])>neqns)goto n1000`, checked after the
# pivot is numbered and before it is eliminated. When the last supervariable is
# reached there is nothing left to update, so genmmd skips the elimination and goes
# straight to the numbering. Ours is the same test on num_eliminated: take the
# pivot, account for it, and stop rather than eliminate into an empty graph.
#
# Run: python3 mmd2.py       every example
#      python3 mmd2.py 3     just the third

# %%
import sys

# I[u] cliques that contain u
# C[c] vertices that c contains

def mmd2_show(A, I, C, degrees, title=None, eliminated=None):
    """Print a quotient graph: adjacency, incidence, cliques, degrees, in the order
    the structure holds them. The degree shown is the stored one, never
    recomputed, which is the point of this layer."""
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
        print(f"  {u:>{width}}: {{{adjacency_text}}} {{{incidence_text}}} degree {degrees[u]}")
    for c in sorted(C):
        clique_members_text = " ".join(f"{u:>{width}}" for u in C[c])
        print(f"  c{c}: {{{clique_members_text}}}")
    print()

def mmd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots,
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

def mmd2_storage(A, I, C):
    """Entries actually stored, as in md5. Batching changes when degrees are
    refreshed, not what the quotient graph holds."""
    return (sum(len(adjacency) for adjacency in A)
            + sum(len(incidence) for incidence in I)
            + sum(len(clique_members) for clique_members in C.values()))

def mmd2_neighbors(A, I, C, eliminated, mark, tag, u):
    """The neighbors of live vertex u: its explicit adjacency A[u] together with
    the members of every clique that contains u, minus u itself, which the
    cliques always carry. This is George and Liu's reachable set, and it is what
    the elimination graph would hold explicitly.

    A vertex numbered by the prepass is skipped, which is what marker[mn] = maxint
    does in genmmd: it is ordered but never eliminated, so it is still sitting in
    the adjacency lists of its neighbors.

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

def mmd2_degree(A, I, C, eliminated, super_members, mark, tag, u):
    """The weighted degree of u: its neighbors counted in ORIGINAL vertices, since
    a neighbor may stand for several. Returns (degree, tag).

    The count of a supervariable is len(super_members[v]), which is O(1), so no
    weight array is kept: md3 through mmd1 have none for the same reason. One
    becomes necessary only when the members are chains over a flat array, where a
    size stops being free."""
    neighbors, tag = mmd2_neighbors(A, I, C, eliminated, mark, tag, u)
    return sum(len(super_members[v]) for v in neighbors), tag

def mmd2_eliminate(A, I, C, mark, tag, eliminated, outmatched, pivot):
    """Turn the pivot into a clique, then merge in every member it makes
    indistinguishable. Identical to md5_eliminate: this layer changes how often
    degrees are refreshed, not what an elimination does.

    Returns (neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag): as
    in md2, plus the vertices folded into the pivot by mass elimination. The middle
    three are reported for display; only neighbors is used by the caller.
    """
    neighbors, tag = mmd2_neighbors(A, I, C, eliminated, mark, tag, pivot)
    absorbed_cliques = list(I[pivot])
    for c in absorbed_cliques:
        del C[c]
    C[pivot] = list(neighbors)      # becomes L_pivot, the column pattern

    # Stamp the new clique once, and the absorbed cliques once. Membership is then
    # a comparison, and both loops below are compactions in place.
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
        outmatched[u] = False       # mmdelm's bwd[rn] = 0: back in the running
        kept = []
        for v in A[u]:
            if v == pivot:              # the pivot is no longer a variable
                continue
            if eliminated[v]:           # numbered by the prepass, gone for good
                continue
            if mark[v] == clique_tag:   # both ends inside the new clique
                if u < v:
                    pruned_edges.append((u, v))
                continue                # implicit now: delete the explicit copy
            kept.append(v)
        A[u] = kept

        kept = [c for c in I[u] if mark[c] != absorbed_tag]   # absorbed are gone
        kept.append(pivot)              # u joins the new clique, whose id is the pivot
        I[u] = kept

    # Mass elimination. u is INDISTINGUISHABLE from the pivot when the two have
    # the same closed neighborhood, mmd2_neighbors(u) | {u} == mmd2_neighbors(pivot)
    # | {pivot}, as it stood before the step. Equivalently, now that the clique is
    # formed, when everything u can still reach lies inside it. The test below is
    # a cheap sufficient condition for that: nothing explicit left and no clique
    # but the new one means u sees exactly what the pivot sees, so eliminating it
    # next would cost no fill. Fold it into the pivot now and strip it from the
    # cliques, since it is no longer a vertex.
    merged_vertices = []
    for u in neighbors:
        if not A[u] and len(I[u]) == 1 and I[u][0] == pivot:
            I[u] = []
            eliminated[u] = True
            merged_vertices.append(u)
    if merged_vertices:                 # one compaction pass, not a removal each
        tag += 1
        for u in merged_vertices:
            mark[u] = tag
        C[pivot] = [v for v in C[pivot] if mark[v] != tag]

    A[pivot] = []
    I[pivot] = []
    eliminated[pivot] = True
    return neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag

def mmd2_file(buckets, filed, d, u):
    """Push u at the head of bucket d, which is the O(1) end of the C++ twin's
    linked list."""
    buckets[d].insert(0, u)
    filed[u] = True

def mmd2_unfile(buckets, filed, d, u):
    """Take u out of bucket d, if it is there."""
    if not filed[u]:
        return
    buckets[d].remove(u)
    filed[u] = False

def mmd2_refile(buckets, filed, degrees, u, new_degree):
    """Move u from the bucket for its old degree to the one for new_degree."""
    mmd2_unfile(buckets, filed, degrees[u], u)
    degrees[u] = new_degree
    mmd2_file(buckets, filed, new_degree, u)

def mmd2_minimum_degree(G, delta=0):
    """Multiple elimination: a batch of independent pivots per degree refresh.

    delta is signed: negative means one pivot per round. It is compared against a
    degree and its useful range stops at n - 1, so in the C++ twin it is a
    std::int32_t, an index-like quantity rather than a count.

    delta widens the batch to vertices within delta of the minimum degree, which
    buys still fewer refreshes for a real concession, since those vertices are not
    minimal. delta = 0 keeps the batch to true minima. A negative delta takes one
    pivot per round, which is md5's behavior reached through this code path.
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
    super_members = [[u] for u in range(n)]    # the vertices each pivot stands for
    eliminated = [False] * n
    outmatched = [False] * n                   # withheld from the buckets, not merged
    pivots = []                                # the order over supervariables
    num_eliminated = 0                         # a counter, not a scan of eliminated
    nnz_L = 0

    degrees = [len(A[u]) for u in range(n)]     # weighted from here on
    num_degree_computations = n

    # mmdint files a degree-0 vertex under degree 1, `if(dg==0)dg=1`, so the
    # bucket a vertex sits in is max(degree, 1) rather than its degree. From here
    # degrees[] holds that filed value, which is what MMD compares and files by.
    buckets = [[] for _ in range(n)]           # buckets[d] holds the live degree-d
    filed = [False] * n                        # whether u is in a bucket at all
    for u in range(n):
        degrees[u] = max(degrees[u], 1)
        mmd2_file(buckets, filed, degrees[u], u)
    min_degree = min(degrees) if n else 0
    num_bucket_probes = 0
    num_rounds = 0                             # batches, the metric this layer adds
    ncsub = 0                                  # genmmd's subscript estimate
    pair_merges = 0                            # q2h merges, the coarser supervariables
    outmatched_count = 0                       # vertices withheld rather than refiled
    touched_round = [-1] * n                   # the round in which u was last evicted

    mmd2_show(A, I, C, degrees, "start: every edge explicit, no clique yet",
              eliminated=eliminated)
    mmd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)

    # ---- the PREPASS -------------------------------------------------------
    # Number everything in bucket 1, which after the floor above holds the
    # isolated and the degree-1 vertices together, and then leave the bucket
    # empty. Nothing is eliminated in the quotient-graph sense: no clique is
    # formed, nothing is pruned, and the neighbors keep degrees that still count
    # these vertices. That staleness is the point, and it is what genmmd does.
    prepass_vertices = list(buckets[1]) if n > 1 else list(buckets[1])
    for u in prepass_vertices:
        mmd2_unfile(buckets, filed, degrees[u], u)
        external_degree = sum(1 for v in A[u] if not eliminated[v])
        nnz_L += external_degree + 1
        eliminated[u] = True
        pivots.append(u)
        num_eliminated += 1
    if prepass_vertices:
        prepass_text = ", ".join(str(u) for u in prepass_vertices)
        mmd2_show(A, I, C, degrees,
                  f"prepass: numbered {len(prepass_vertices)}: {prepass_text}",
                  eliminated=eliminated)
        mmd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)
    min_degree = 2 if n > 2 else min_degree     # head[1] = 0, and mdeg starts at 2

    while num_eliminated < n:
        while not buckets[min_degree]:         # walk up to the first live bucket
            min_degree += 1
            num_bucket_probes += 1
        num_bucket_probes += 1

        # ---- one BATCH, no degree refreshed inside it ----------------------
        # Take pivots from buckets [min_degree, min_degree + delta]. Eviction is
        # what keeps them independent: eliminating a pivot pulls every vertex it
        # reached out of the buckets, so whatever is still filed was not reached,
        # hence is not adjacent to anything taken this round.
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
            mmd2_unfile(buckets, filed, degree, pivot)

            neighbors, absorbed_cliques, pruned_edges, merged_vertices, tag = mmd2_eliminate(
                A, I, C, mark, tag, eliminated, outmatched, pivot)
            batch.append(pivot)
            pivots.append(pivot)
            num_eliminated += 1 + len(merged_vertices)
            for u in merged_vertices:          # the pivot now stands for them too
                super_members[pivot] += super_members[u]
                super_members[u] = []
                mmd2_unfile(buckets, filed, degrees[u], u)
                degrees[u] = 0
            degrees[pivot] = 0

            for u in C[pivot]:                 # EVICT, with a stale degree
                mmd2_unfile(buckets, filed, degrees[u], u)
                if touched_round[u] != num_rounds:   # a marker, so O(1) per eviction
                    touched_round[u] = num_rounds
                    touched.append(u)

            ncsub += degree + len(super_members[pivot]) - 2   # genmmd's *ncsub
            super_size = len(super_members[pivot])
            external_degree = sum(len(super_members[v]) for v in C[pivot]
                                  if not eliminated[v])
            nnz_L += (super_size * external_degree
                      + super_size * (super_size - 1) // 2
                      + super_size)

            absorbed_cliques_text = ", ".join(f"c{c}" for c in absorbed_cliques) if absorbed_cliques else "none"
            pruned_edges_text = ", ".join(f"{u}-{v}" for u, v in pruned_edges) if pruned_edges else "none"
            merged_vertices_text = ", ".join(str(u) for u in merged_vertices) if merged_vertices else "none"
            evicted_text = ", ".join(str(u) for u in C[pivot]) if C[pivot] else "none"
            print(f"round {num_rounds}: eliminate {pivot} (degree {degree}, size {super_size}, "
                  f"external degree {external_degree}), "
                  f"absorbed cliques: {absorbed_cliques_text}, pruned edges: {pruned_edges_text}, "
                  f"merged vertices: {merged_vertices_text}, evicted: {evicted_text}")
            if num_eliminated >= n:            # genmmd's num + qsize[mn] > neqns:
                break                          # nothing left to update
            if delta < 0:                      # one pivot per round, as md5 does
                break

        # ---- one REFRESH, walked ELEMENT BY ELEMENT ------------------------
        # mmdupd walks the elements this round created, not the vertices it
        # reached, and computes dg0 once per element: the size of that element,
        # which every member of it reaches in full. A member with exactly one other
        # source goes on the q2h list and is answered from dg0 plus that source;
        # everything else goes on qxh and pays for the full union. Same degrees,
        # different work, and a different filing order.
        refreshed_vertices = []
        for element in batch:
            element_members = [u for u in C[element] if not eliminated[u]]
            tag += 1                            # dg0's members, marked once
            element_tag = tag
            for v in element_members:
                mark[v] = element_tag
            dg0 = sum(len(super_members[v]) for v in element_members)

            q2h, qxh = [], []
            for u in element_members:
                if filed[u] or outmatched[u]:   # already refreshed this round, or
                    continue                    # withheld as outmatched
                other_sources = len(A[u]) + len(I[u]) - 1
                (q2h if other_sources == 1 else qxh).append(u)

            for u in q2h:
                if eliminated[u] or outmatched[u]:   # merged or withheld by an
                    continue                         # earlier q2h vertex
                # Everything u reaches is in the element or comes from its one
                # other source. dg0 counts the element, minus u itself. Two mark
                # levels, as mmdupd has: element_tag says "already in dg0" and
                # survives the whole element, while vertex_tag is fresh per vertex,
                # so one q2h vertex cannot hide a neighbor from the next.
                tag += 1
                vertex_tag = tag
                degree = dg0 - len(super_members[u])
                for v in A[u]:
                    if eliminated[v] or mark[v] == vertex_tag:
                        continue
                    if mark[v] == element_tag:
                        continue                # already in dg0
                    mark[v] = vertex_tag
                    degree += len(super_members[v])
                for c in I[u]:
                    if c == element:
                        continue
                    for v in C[c]:
                        if v == u or eliminated[v] or mark[v] == vertex_tag:
                            continue
                        if mark[v] == element_tag:
                            # v is in the new element AND in this same other
                            # source, so it sees at least what u sees.
                            if filed[v] or outmatched[v]:
                                continue
                            if len(A[v]) + len(I[v]) - 1 == 1:
                                # v is q2h too, so its only other source is this
                                # one: identical reach, and u absorbs it.
                                super_members[u] += super_members[v]
                                super_members[v] = []
                                eliminated[v] = True
                                num_eliminated += 1
                                pair_merges += 1
                            else:
                                # v reaches more than u does, so it can never be
                                # the minimum first. Withhold it from the buckets.
                                outmatched[v] = True
                                outmatched_count += 1
                            continue
                        mark[v] = vertex_tag
                        degree += len(super_members[v])
                degrees[u] = max(degree + 1, 1)   # dg - qsize[en] + 1, floored
                mmd2_file(buckets, filed, degrees[u], u)
                refreshed_vertices.append(u)

            for u in qxh:
                if eliminated[u] or outmatched[u]:
                    continue
                degree, tag = mmd2_degree(A, I, C, eliminated, super_members, mark, tag, u)
                degrees[u] = max(degree + 1, 1)   # dg - qsize[en] + 1, floored
                mmd2_file(buckets, filed, degrees[u], u)
                refreshed_vertices.append(u)

        num_degree_computations += len(refreshed_vertices)
        min_degree = min([min_degree] + [degrees[u] for u in refreshed_vertices])
        num_rounds += 1

        batch_text = ", ".join(str(u) for u in batch)
        refreshed_vertices_text = ", ".join(str(u) for u in refreshed_vertices) if refreshed_vertices else "none"
        mmd2_show(A, I, C, degrees,
                  (f"round {num_rounds - 1} done: batch of {len(batch)}: {batch_text}, "
                   f"refreshed: {refreshed_vertices_text}"),
                  eliminated=eliminated)
        mmd2_show_state(degrees, buckets, min_degree, super_members, eliminated, pivots)

    order = [u for pivot in pivots for u in super_members[pivot]]
    print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
          f"fill = {nnz_L - nnz_tril_A}")
    print(f"degree computations: {num_degree_computations}, "
          f"bucket probes: {num_bucket_probes}, rounds: {num_rounds}, "
          f"pair merges: {pair_merges}, outmatched: {outmatched_count}, "
          f"ncsub: {ncsub}")
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

# All of them by default. To run just one, pass its number: python3 mmd2.py 3
selected = int(sys.argv[1]) if len(sys.argv) > 1 else 0
for number, (name, g) in enumerate(examples, start=1):
    if selected and number != selected:
        continue
    print(f"=== {name} ===")
    mmd2_minimum_degree(g)
    print()
