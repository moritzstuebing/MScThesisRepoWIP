import networkx as nx
import more_itertools as mit

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
        " ".join(str(v) for v in sorted(block))
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



