# %% [markdown]
# # Minimum degree, layer 1: naive, materializing fill
#
# The learning version of the ordering, in Python so we can step through it.
# This is the NAIVE minimum degree of Section 5.1 (archive/sparse_factorization.md):
# it materializes fill, the one thing a real code refuses to do. That is the
# point. It is the vertex elimination game transcribed line for line, the ground
# truth the quotient graph (layer 2) will optimize without changing the result.
#
# Scope: correct, not equivalent. It produces a valid minimum degree ordering,
# not the vendored genmmd's exact one; matching genmmd's tie-breaks is a later,
# separate goal. The vendored src/Mmd.cpp is the oracle for that.

# %%
from copy import deepcopy

# The elimination graph: one neighbor set per vertex. A set keeps each
# neighborhood sorted-enough and duplicate-free, and reads plainly.
Graph = list  # list[set[int]], spelled out for the reader


def graph_from_edges(n, edges):
    """Build an undirected graph on n vertices from a list of (u, w) edges.

    Pythonic input for the notebook. The real orderer takes CSC (colPtr/rowIdx);
    we bridge to that only when we compare against the oracle, later.
    """
    graph = [set() for _ in range(n)]
    for u, w in edges:
        if u != w:
            graph[u].add(w)
            graph[w].add(u)  # keep it symmetric
    return graph


# %%
def eliminate(graph, pivot):
    """Eliminate one vertex: make its neighbors a clique, then detach it.

    Making the neighbors mutually adjacent is where fill is created, and here we
    actually store it. This function *is* the definition of a pivot's effect on
    the graph, and the clique loop is the one line layer 2 will remove.
    """
    neighbors = set(graph[pivot])  # copy: we mutate the graph below
    for u in neighbors:
        for w in neighbors:
            if u != w:
                graph[u].add(w)  # the fill edge, materialized
    for u in neighbors:
        graph[u].discard(pivot)
    graph[pivot].clear()


def minimum_degree(graph):
    """Repeatedly eliminate the live vertex of least degree.

    Returns the elimination order; order[k] is the vertex eliminated k-th, which
    read as a sequence is the new-to-old permutation. Does not mutate the input.
    """
    graph = deepcopy(graph)
    n = len(graph)
    eliminated = [False] * n
    order = []

    for _ in range(n):
        # PICK the live vertex of least degree. Linear scan, first index wins a
        # tie (strict <), matching the C++ layer 1 so the two agree.
        pivot = -1
        least_degree = 0
        for vertex in range(n):
            if eliminated[vertex]:
                continue
            degree = len(graph[vertex])
            if pivot == -1 or degree < least_degree:
                pivot = vertex
                least_degree = degree

        eliminate(graph, pivot)
        eliminated[pivot] = True
        order.append(pivot)

    return order


def count_fill(graph, order):
    """Count the fill an ordering produces, replaying it on a fresh copy."""
    graph = deepcopy(graph)
    fill = 0
    for pivot in order:
        neighbors = set(graph[pivot])
        for u in neighbors:
            for w in neighbors:
                if u < w and w not in graph[u]:
                    graph[u].add(w)
                    graph[w].add(u)
                    fill += 1
        for u in neighbors:
            graph[u].discard(pivot)
        graph[pivot].clear()
    return fill


# %%
def minimum_degree_steps(graph):
    """Same as minimum_degree, but yields state after each elimination.

    This is why we are in a notebook. Each yield carries the step number, the
    pivot chosen, its degree at the moment of choice, and a snapshot of the live
    graph, so we can watch degrees change and fill appear as the cell runs.
    """
    graph = deepcopy(graph)
    n = len(graph)
    eliminated = [False] * n

    for step in range(n):
        pivot = -1
        least_degree = 0
        for vertex in range(n):
            if eliminated[vertex]:
                continue
            degree = len(graph[vertex])
            if pivot == -1 or degree < least_degree:
                pivot = vertex
                least_degree = degree

        eliminate(graph, pivot)
        eliminated[pivot] = True
        live = {v: sorted(graph[v]) for v in range(n) if not eliminated[v]}
        yield step, pivot, least_degree, live


# %%
# Arrowhead: vertex 0 is a hub joined to 1..n-1; every other vertex touches only
# the hub. Example 1 of the doc: eliminate the hub first and its neighbors become
# a full clique (much fill); eliminate it last and there is none. Minimum degree
# finds the good order on its own.
n = 6
arrowhead = graph_from_edges(n, [(0, leaf) for leaf in range(1, n)])

md_order = minimum_degree(arrowhead)
natural_order = list(range(n))  # hub first

print("minimum degree order:     ", md_order,
      " fill =", count_fill(arrowhead, md_order))
print("natural order (hub first):", natural_order,
      " fill =", count_fill(arrowhead, natural_order))

# %%
# The interactive payoff: watch it eliminate. Run this cell and read the pivots
# and degrees. Notice the hub's degree falling by one each time a leaf goes.
for step, pivot, degree, live in minimum_degree_steps(arrowhead):
    print(f"step {step}: eliminate {pivot} (degree {degree}), live now {live}")

# %%
# A slightly richer graph to experiment with: a 3x3 grid. Change it, re-run, see
# how the ordering and the fill move. This is the cell to play in.
def grid_edges(rows, cols):
    edges = []
    for r in range(rows):
        for c in range(cols):
            v = r * cols + c
            if c + 1 < cols:
                edges.append((v, v + 1))
            if r + 1 < rows:
                edges.append((v, v + cols))
    return edges


grid = graph_from_edges(9, grid_edges(3, 3))
grid_order = minimum_degree(grid)
print("3x3 grid, minimum degree:", grid_order, " fill =", count_fill(grid, grid_order))
print("3x3 grid, natural order: ", list(range(9)),
      " fill =", count_fill(grid, list(range(9))))
