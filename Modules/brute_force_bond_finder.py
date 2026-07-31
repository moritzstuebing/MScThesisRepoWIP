import networkx as nx
import more_itertools as mit

def bond_checker(partition, G):

    """
        Checks whether a given partition of vertices of a graph G defines a bond
    """

    for block in partition:
        subgraph = nx.subgraph(G, block)
        if nx.is_connected(subgraph) == True:
            continue
        else:
            return False

    return True

def bond_sorter(bonds, n):  # try and work out if possible to do this without needing n input, purely working out number of nodes n from the bonds
    
    """
        Sorts bonds by rank (number of blocks) order
    """
    
    return sorted(bonds, key=lambda pi: n - len(pi))

def brute_force_bond_finder(G):

    """
        Finds all vertex partitions that define bonds of the graph G. Does this by brute force,
        computing every possible partition of the set and checking whether this defines a bond
    """

    # All partitions of the vertices
    partitions = mit.set_partitions(G.nodes)

    # Loop through all set partitions, check they are bonds, add to list of bonds if so
    bonds = []
    for p in partitions:
        verdict = bond_checker(p, G)
        if verdict == True:
            bonds.append(frozenset(frozenset(block) for block in p)) # turn blocks and partitions into sets and sort
    
    return bond_sorter(bonds, G.number_of_nodes())



