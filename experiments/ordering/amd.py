# %% [markdown]
# # Approximate minimum degree
#
# The other fork from md6. Section 5.13 of archive/sparse_factorization.md.
#
# md6 has the quotient graph, supervariables, maintained degrees and buckets,
# and returns exactly md1's ordering. What is left costing anything is the
# refresh itself, which for each reached vertex i unites the members of every
# element in E_i and counts the result. That union is the expensive object.
#
# MMD made the refresh RARE (mmd.py). AMD makes each one CHEAP, and the two are
# the same answer reached from opposite ends: do the expensive thing less.
#
# THE BOUND. Rather than uniting the elements, sum their separate contributions:
#
#   degree[i] <= min( n - k,                                 nothing exceeds what remains
#                     degree_old[i] + |L \ i|,               it can only grow by the new element
#                     |A_i \ L| + |L \ i| + sum |L_e \ L| )  over e in E_i, e != the new element
#
# The third line overcounts, because two elements may overlap outside L and the
# overlap is counted twice. So it is an UPPER BOUND, not the degree.
#
# WHY THAT IS FAST, which is the entire point and is easy to miss. The quantity
# |L_e \ L| depends only on the element e, not on the vertex i, so it is computed
# ONCE PER ELEMENT and then read by every vertex whose element list contains e.
# The exact degree costs, per vertex, a walk over the members of all its
# elements. The bound costs, per vertex, one addition per element. Counting both
# below shows the gap widening with the size of the elements, which is to say
# with the amount of fill, which is to say exactly where it matters.
#
# TWO MORE MECHANISMS, both beyond md6 and neither about the degree:
#
#   - AGGRESSIVE ABSORPTION. If |L_e \ L| == 0, element e lies entirely inside
#     the new element, so it is dead and can be absorbed at once. Ordinary
#     absorption only kills the elements the PIVOT touched; this kills elements
#     that any reached vertex touched. Cheap, since |L_e \ L| has just been
#     computed anyway.
#   - HASH SUPERVARIABLE DETECTION. Our mass-elimination test only merges
#     vertices indistinguishable from the PIVOT. Two vertices can be
#     indistinguishable from EACH OTHER without either being absorbable into the
#     pivot. AMD hashes (A_i, E_i), compares only within a hash bucket, and
#     merges on an exact match. MMD reaches the same vertices through mmdupd's
#     q2h list; the goal is shared and the mechanism is not.
#
# WHAT IS GIVEN UP. Unlike MMD, whose pivots are exact and whose ordering moves
# only through tie-breaks, AMD can genuinely pick the wrong vertex: an
# overcounted bound can hide the true minimum. So the quality loss here is of a
# different kind, and this is the first prototype whose pivot is not guaranteed
# to be a minimum-degree vertex at all.

#
# HOW THIS DIFFERS FROM THE VENDORED Amd.cpp. One algorithmic difference, and it
# shows in the output; see section 5.13.
#
#   - THE UPDATE RUNS IN ONE PASS HERE, IN TWO THERE. Amd.cpp computes a bound
#     that excludes the new element, mass-eliminating and shrinking degme as it
#     goes, and then a SECOND loop adds the final degme to every survivor. We fold
#     both into one loop, so a vertex handled early sees a larger degme than one
#     handled late, where the vendored code gives them all the same value. Since
#     degme only shrinks, our early vertices get looser bounds. On eleven test
#     graphs this changes the ordering on four, all grids, with fill moving a few
#     tenths of a percent either way. One pass is easier to read; two is faithful.
#
# And several missing features, none of them about the degree: dense-row handling
# (amd_preprocess, AMD_DENSE), which changes orderings on matrices that have such
# rows; the postorder (amd_postorder), which leaves fill alone but changes the
# permutation; amd_aat, which forms A + A' so unsymmetric input can be ordered on
# its symmetric pattern; amd_valid; the Control and Info arrays, where we hardcode
# the choices; and the workspace compression, which is a memory strategy for flat
# arrays with no counterpart in a prototype built on sets.

# %%
def reachable(A, C, cliques, i):
    neighbors = set(A[i])
    for c in C[i]:
        neighbors |= cliques[c]
    neighbors.discard(i)
    return neighbors


def exact_degree(A, C, cliques, weight, i):
    """Kept only to measure the bound against the truth."""
    return sum(weight[j] for j in reachable(A, C, cliques, i))


def storage(A, cliques):
    return sum(len(a) for a in A) + sum(len(m) for m in cliques.values())


def amd(graph, aggressive=True, verbose=True):
    n = len(graph)
    nnz_tril_A = sum(len(neigh) for neigh in graph) // 2 + n
    A = [set(neigh) for neigh in graph]
    C = [set() for _ in range(n)]
    cliques = {}
    weight = [1] * n
    members = [[v] for v in range(n)]
    eliminated = [False] * n
    pivots = []
    nnz_L = 0

    degrees = [len(A[v]) for v in range(n)]
    buckets = [set() for _ in range(n + 1)]
    for v in range(n):
        buckets[degrees[v]].add(v)
    mdeg = min(degrees) if n else 0

    # work counters: the whole argument for the bound is the gap between these
    member_visits = 0   # what the EXACT degree would cost (element members walked)
    element_reads = 0   # what the BOUND costs (one addition per element)
    overcounts = 0      # times the bound exceeded the true degree
    bound_checks = 0
    eliminated_count = 0

    def show(label, indent="  "):
        if not verbose:
            return
        a = "{" + ", ".join(f"{v}: {sorted(A[v])}" for v in range(n)) + "}"
        c = "{" + ", ".join(
            f"{v}: {[f'c{x}' for x in sorted(C[v])]}" for v in range(n)) + "}"
        cl = "{" + ", ".join(
            f"c{x}: {sorted(m)}" for x, m in sorted(cliques.items())) + "}"
        w = "{" + ", ".join(f"{v}: {weight[v]}" for v in range(n)) + "}"
        dg = "{" + ", ".join(
            f"{v}: {degrees[v] if not eliminated[v] else '-'}" for v in range(n)) + "}"
        ex = "{" + ", ".join(
            f"{v}: {exact_degree(A, C, cliques, weight, v) if not eliminated[v] else '-'}"
            for v in range(n)) + "}"
        bk = "{" + ", ".join(f"{d}: {sorted(b)}"
                             for d, b in enumerate(buckets) if b) + "}"
        print(f"{indent}A       = {a}")
        print(f"{indent}C       = {c}")
        print(f"{indent}cliques = {cl}")
        print(f"{indent}weights = {w}   storage {storage(A, cliques)}")
        print(f"{indent}bounds  = {dg}   {label}")
        print(f"{indent}exact   = {ex}")
        print(f"{indent}buckets = {bk}   mdeg = {mdeg}")

    if verbose:
        print("start: no cliques yet, degrees exact (no elements to approximate over)")
    show("(exact at the start)")

    while eliminated_count < n:
        while mdeg <= n and not buckets[mdeg]:
            mdeg += 1
        if mdeg > n:
            break
        pivot = min(buckets[mdeg])
        d_pivot = degrees[pivot]
        buckets[d_pivot].discard(pivot)

        # ---- form the new element -----------------------------------------
        L = reachable(A, C, cliques, pivot)
        absorbed = set(C[pivot])
        for c in absorbed:
            del cliques[c]
        cliques[pivot] = set(L)
        degme = sum(weight[j] for j in L)

        for i in L:
            A[i] -= (A[i] & L)
            A[i].discard(pivot)
            C[i] -= absorbed
            C[i].add(pivot)
        A[pivot].clear()
        C[pivot].clear()
        eliminated[pivot] = True
        eliminated_count += weight[pivot]   # merged vertices are counted as they merge

        # ---- |L_e \ L| ONCE PER ELEMENT, and aggressive absorption --------
        touched_elements = set()
        for i in L:
            touched_elements |= (C[i] - {pivot})
        we = {}
        dead = []
        for e in touched_elements:
            outside = sum(weight[j] for j in cliques[e] if j not in L)
            member_visits += len(cliques[e])   # the exact degree pays this PER VERTEX
            we[e] = outside
            if aggressive and outside == 0:
                dead.append(e)                 # e lies inside L: absorb it now
        for e in dead:
            del cliques[e]
            for i in range(n):
                C[i].discard(e)

        # ---- the bound, plus mass elimination -----------------------------
        nleft = n - eliminated_count
        merged = []
        for i in sorted(L):
            if eliminated[i]:
                continue
            outside_A = sum(weight[j] for j in A[i])
            elems = C[i] - {pivot}
            element_reads += len(elems)        # the bound pays only this
            bound = outside_A + sum(we[e] for e in elems if e in we)

            if not A[i] and not elems:
                # MASS ELIMINATION: nothing outside the new clique
                weight[pivot] += weight[i]
                eliminated_count += weight[i]   # counted here, not via the pivot
                degme -= weight[i]              # i joins the pivot, so it is no
                                                # longer part of the new element
                weight[i] = 0
                C[i].clear()
                eliminated[i] = True
                merged.append(i)
                continue

            bound += degme - weight[i]
            bound = min(bound, nleft - weight[i], degrees[i] + degme - weight[i])
            true_degree = exact_degree(A, C, cliques, weight, i)
            bound_checks += 1
            if bound > true_degree:
                overcounts += 1
            buckets[degrees[i]].discard(i)
            degrees[i] = bound
            buckets[bound].add(i)
            if bound < mdeg:
                mdeg = bound

        for i in merged:
            members[pivot] += members[i]
            buckets[degrees[i]].discard(i)
            cliques[pivot].discard(i)
            for c in cliques.values():
                c.discard(i)

        # ---- HASH SUPERVARIABLE DETECTION ---------------------------------
        # Vertices indistinguishable from EACH OTHER, which the pivot test above
        # cannot see. Hash first so the exact comparison runs only within a
        # bucket; the hash is a filter, never the decision.
        survivors = [i for i in sorted(L) if not eliminated[i]]
        by_hash = {}
        for i in survivors:
            h = (hash(frozenset(A[i])), hash(frozenset(C[i])))
            by_hash.setdefault(h, []).append(i)
        pairs = []
        for group in by_hash.values():
            if len(group) < 2:
                continue
            for x in range(len(group)):
                i = group[x]
                if eliminated[i]:
                    continue
                for y in range(x + 1, len(group)):
                    j = group[y]
                    if eliminated[j]:
                        continue
                    # exact test: same explicit neighbors and same elements,
                    # once each is removed from the other's neighbor set
                    if (A[i] - {j}) == (A[j] - {i}) and C[i] == C[j]:
                        weight[i] += weight[j]
                        weight[j] = 0
                        members[i] += members[j]
                        buckets[degrees[j]].discard(j)
                        A[j].clear()
                        C[j].clear()
                        eliminated[j] = True
                        # j is absorbed into the live supervariable i, not
                        # eliminated: its weight is counted when i is chosen
                        for c in cliques.values():
                            c.discard(j)
                        for v in range(n):
                            A[v].discard(j)
                        pairs.append((i, j))

        pivots.append(pivot)
        wp = weight[pivot]
        ext = sum(weight[j] for j in cliques[pivot])
        nnz_L += wp * ext + wp * (wp - 1) // 2 + wp

        if verbose:
            mg = ", ".join(map(str, merged)) if merged else "none"
            ab = ", ".join(f"c{e}" for e in sorted(dead)) if dead else "none"
            pr = ", ".join(f"{i}<-{j}" for i, j in pairs) if pairs else "none"
            print(f"eliminate {pivot} (bound {d_pivot}, weight -> {wp})"
                  f"   merged = {mg}   absorbed = {ab}   hash-merged = {pr}")
            show("(bounds after)")

    order = [v for p in pivots for v in members[p]]
    if verbose:
        print(f"nnz(L) = {nnz_L} against nnz(tril A) = {nnz_tril_A}, "
              f"fill = {nnz_L - nnz_tril_A}")
        print(f"element-member visits an exact degree would need: {member_visits}")
        print(f"element reads the bound needed:                   {element_reads}")
        print(f"bound was loose {overcounts} times out of {bound_checks}")
    return order


# %%
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

for name, g in [("graph1", graph1), ("graph2", graph2),
                ("graph3", graph3), ("graph4", graph4)]:
    print(f"=== {name} ===")
    order = amd(g)
    print("order:", order)
    print()
