"""
Bond enumeration via Ganter's Next-Closure algorithm over the cycle-matroid
flat lattice, returning bonds as vertex partitions (list of lists of lists).
"""

def preprocess_graph(G):
    edges = list(G.edges())
    nodes = list(G.nodes())
    vid = {v: i for i, v in enumerate(nodes)}
    edge_u, edge_v = [], []
    for u, v in edges:
        edge_u.append(vid[u])
        edge_v.append(vid[v])
    return edges, nodes, vid, edge_u, edge_v


class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _components(mask, edge_u, edge_v, n_vertices):
    """DSU over all vertices using the selected edges; returns comp-root array."""
    dsu = DSU(n_vertices)
    m, i = mask, 0
    while m:
        if m & 1:
            dsu.union(edge_u[i], edge_v[i])
        m >>= 1
        i += 1
    return [dsu.find(v) for v in range(n_vertices)]


def bond_closure(mask, edge_u, edge_v, n_vertices):
    comp = _components(mask, edge_u, edge_v, n_vertices)
    out, n_edges = 0, len(edge_u)
    for i in range(n_edges):
        if comp[edge_u[i]] == comp[edge_v[i]]:
            out |= (1 << i)
    return out


def mask_to_partition(mask, edge_u, edge_v, nodes):
    """Convert a (closed) edge bitmask to a vertex partition: list of blocks,
    each block a sorted list of original node labels. Isolated vertices are
    their own singleton blocks."""
    n_vertices = len(nodes)
    comp = _components(mask, edge_u, edge_v, n_vertices)
    groups = {}
    for v in range(n_vertices):
        groups.setdefault(comp[v], []).append(nodes[v])
    blocks = [sorted(b, key=lambda x: (str(type(x)), x)) for b in groups.values()]
    # sort blocks by their smallest element for a canonical ordering
    blocks.sort(key=lambda b: (len(b), [str(x) for x in b]))
    return blocks


def A_oplus_i(A, i, edge_u, edge_v, n_vertices):
    B = A & ((1 << i) - 1)
    return bond_closure(B | (1 << i), edge_u, edge_v, n_vertices)


def icover_check(A, B, i):
    if not (B & (1 << i)):
        return False
    diff = A ^ B
    if diff == 0:
        return False
    lowest = (diff & -diff).bit_length() - 1
    return lowest == i


def one_step(A, edge_u, edge_v, n_vertices, n_edges):
    for i in reversed(range(n_edges)):
        if A & (1 << i):
            continue
        candidate = A_oplus_i(A, i, edge_u, edge_v, n_vertices)
        if icover_check(A, candidate, i):
            return candidate
    return None


def next_closure_algo(G, as_partitions=True):
    """Enumerate all bonds of G.

    Returns
    -------
    If as_partitions (default): a list of bonds, each bond a list of blocks,
      each block a sorted list of vertices  ->  list[list[list]].
    Else: the raw list of edge-set bitmasks (original behaviour).
    """
    edges, nodes, vid, edge_u, edge_v = preprocess_graph(G)
    n_vertices = len(nodes)
    n_edges = len(edges)

    A = 0
    masks = [A]
    while True:
        A = one_step(A, edge_u, edge_v, n_vertices, n_edges)
        if A is None:
            break
        masks.append(A)

    if not as_partitions:
        return masks
    return [mask_to_partition(m, edge_u, edge_v, nodes) for m in masks]