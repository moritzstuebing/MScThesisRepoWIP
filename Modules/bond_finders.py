import networkx as nx
import itertools
import more_itertools as mit
from collections import deque

def bond_checker(G, partition):

    """
        Checks whether a given partition of vertices of a graph defines a bond.

        Inputs:
        partition - candidate partition of graph vertices to be checked, frozenset of frozensets
        G - graph to check the candidate partition in (NetworkX graph object).

        Outputs:
        True/False depending on whether partition is a bond or not.
    """

    # Loop over blocks in the partition
    for block in partition:

        # Check if the induced subgraph of the block is connected
        subgraph = nx.subgraph(G, block)
        if nx.is_connected(subgraph) == True:
            continue
        else:
            return False

    # If all blocks induce connected subgraphs, return True
    return True

def bond_viewer(bond):

    """
        Takes the bond as a frozenset of frozensets and outputs it in a nice format.

        Inputs:
        bond - bond partition as frozenset of frozensets.

        Outputs:
        prints the bond in the nice format, just numbers with midlines separating the blocks.
    """

    # Sort the bond's blocks in order of the smallest element in each block
    blocks = sorted(bond, key=min)

    # Correctly formatted and print
    formatted = " | ".join(
        " ".join(str(v + 1) for v in sorted(block))
        for block in blocks
    )
    print(formatted)

def sort_key(bond):

    """
        Key function for sorting bonds

        Input:
        bond - frozenset of frozensets

        Output:
        (-len(blocks, blocks)) - tuple with the negative number of blocks (for sorting from 
        finest to coarsest), and the bond itself, sorted lexicographically
    """

    blocks = sorted(sorted(block) for block in bond)
    return (-len(blocks), blocks)


def bond_sorter(bonds):  # try and work out if possible to do this without needing n input, purely working out number of nodes n from the bonds
    
    """
        Sorts bonds by rank (number of blocks) order

        Inputs:
        bonds - the bonds (list of frozensets of frozensets)

        Outputs:
        Sorted bonds list by rank order
    """
    
    return sorted(bonds, key= sort_key)

def brute_force_bond_finder(G, nice_format_view=False):

    """
        Finds all vertex partitions that define bonds of the graph G. Does this by brute force,
        computing every possible partition of the set and checking whether this defines a bond

        Inputs:
        G - the graph
        nice_format_view - for presentation purposes, will output the bonds in order from
                            finest to coarsest, formatted nicely

        Outputs:
        if nice_format_view == True - prints the bonds in rank order formatted nicely
        else - returns the bonds as a list of frozensets of frozensets
    """

    # All partitions of the vertices
    partitions = mit.set_partitions(G.nodes)

    # Loop through all set partitions, check they are bonds, add to list of bonds if so
    bonds = []
    for p in partitions:
        verdict = bond_checker(G, p)
        if verdict == True:
            bonds.append(frozenset(frozenset(block) for block in p)) # turn blocks and partitions into frozenset of frozensets

    # If we just want to present the bonds we format them nicely and print here
    if nice_format_view == True:
        for bond in sorted(bonds, key=sort_key):
            bond_viewer(bond)

    # If we actually want to use the bonds for further operation we use this
    else:
        return bond_sorter(bonds)

def edge_bond_finder(G):

    """
        Finds bonds of a graph by finding all edge subsets and their corresponding vertex partitions.

        Inputs:
        G - the graph (NetworkX object)

        Outputs:
        bond_sorter(bonds) - the sorted bonds of the graph

    """

    # Edge enumeration and initialising bonds list
    dict = {i: e for i, e in enumerate(G.edges)} # edge enumeration
    inv = {v: k for k, v in dict.items()} # inverse list to easily go from edge to edge number
    bonds = set()

    # Loop over all possible sizes of edge subsets
    for i in range(len(G.edges) + 1):

        # Loop over all possible edge sets of size i
        for edge_set in itertools.combinations(G.edges, i):

            # Find edges not present in the edge set and remove them from the graph to from subgraph
            present_edge_numbers = [inv[v] for v in edge_set]
            not_present_numbers = set(range(len(G.edges))) - set(present_edge_numbers)
            removing_edges = [dict[i] for i in not_present_numbers]
            subgraph = G.copy()
            subgraph.remove_edges_from(removing_edges)

            # Get vertices of connected components in subgraph in blocks and add to bonds list
            nodes = frozenset(frozenset(component) for component in nx.connected_components(subgraph))
            if nodes not in bonds: # deduplicating
                bonds.add(nodes)

    # Sort the bonds
    return bond_sorter(bonds)

def merger(G, edges):

    """
        Given a graph and an edge set, produces all edge sets that merge two components.

        Inputs:
        G - the graph (NetworkX object)
        edges - the edge set to find resulting merged edge subsets from (list of tuples)

        Outputs:
        tuples - list of tuples of each resulting edge set and corresponding bond
    """

    # Form spanning subgraph given edges, retrieve nodes of the connected components, and initialise
    # list
    subgraph = G.copy()
    subgraph.remove_edges_from(list(set(G.edges()) - set(edges)))
    component_nodes = [nodes for nodes in nx.connected_components(subgraph)]
    bonds = set()
    tuples = []

    # Loop over all edges in G
    for edge in G.edges():

        # Check if the edge connects separate components of the spanning subgraph (if one endpoint is
        # in a given component and the other isn't, this edge must connect separate components)
        verdict = (
            any(edge[0] in nodes and edge[1] not in nodes for nodes in component_nodes)
        )

        # If edge connects separate components, form another subgraph from spanning subgraph with
        # just this edge added, and retrieve the given bond and edge subset, adding to lists
        if verdict is True:
            edge_subgraph = subgraph.copy()
            edge_subgraph.add_edge(edge[0], edge[1])
            bond = frozenset(frozenset(nodes) for nodes in nx.connected_components(edge_subgraph))

            # Ensure we dont add any duplicated bonds or corresponding edge subsets
            if bond not in bonds:
                tuples.append((bond, [e for e in edge_subgraph.edges()]))
                bonds.add(bond)


    return tuples

def walk_bond_finder(G):

    """
        Finds all bonds of a graph using a walk/branch scheme from the bottom of the lattice to top.

        Inputs:
        G - the graph (NetworkX object)

        Outputs:
        bond_sorter(bonds) - sorted list of bonds, where each bond is a frozenset of frozensets
    """

    # Bonds set (needs a set to check for previously visited bonds later) and wait_list for edge sets
    # waiting to be expanded
    bonds = {frozenset(frozenset([i]) for i in range(G.number_of_nodes()))}
    wait_list = deque([[]])

    while wait_list:  # while non-empty

        # Extract left-most edge subset and find all edge subsets resulting from merging connected
        # components of graph given by this edge subset
        edges = wait_list.popleft()
        tuples = merger(G, edges)

        # Add edge subsets to wait list and bonds to bonds set only if the bond has not been seen
        wait_list.extend(edge_set for (bond, edge_set) in tuples if bond not in bonds)
        bonds.update(bond for (bond, _) in tuples if bond not in bonds)
    
    # bond_sorter returns a list
    return bond_sorter(bonds)
