def preprocess_graph(G):
    edges = list(G.edges())
    nodes = list(G.nodes())

    vid = {v: i for i, v in enumerate(nodes)}

    edge_u = []
    edge_v = []

    for u, v in edges:
        edge_u.append(vid[u])
        edge_v.append(vid[v])

    return edges, vid, edge_u, edge_v

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

def bond_closure(mask, edge_u, edge_v, n_vertices):
    """
    Closure of an edge set represented as a bitmask.
    """

    n_edges = len(edge_u)

    dsu = DSU(n_vertices)

    # Step 1: build connectivity from selected edges
    m = mask
    i = 0
    while m:
        if m & 1:
            dsu.union(edge_u[i], edge_v[i])
        m >>= 1
        i += 1

    # Step 2: compute vertex components
    comp = [dsu.find(v) for v in range(n_vertices)]

    # Step 3: build closure: all edges inside components
    out = 0
    for i in range(n_edges):
        if comp[edge_u[i]] == comp[edge_v[i]]:
            out |= (1 << i)

    return out

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

def next_closure_algo(G):
    edges, vid, edge_u, edge_v = preprocess_graph(G)

    n_vertices = len(G.nodes())
    n_edges = len(edges)

    A = 0
    bonds = [A]

    while True:
        A = one_step(A, edge_u, edge_v, n_vertices, n_edges)

        if A is None:
            break

        bonds.append(A)

    return bonds
