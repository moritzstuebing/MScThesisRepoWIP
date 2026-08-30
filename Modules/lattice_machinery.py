import numpy as np
import networkx as nx
from math import log2
from Modules.bond_finders import brute_force_bond_finder, bond_sorter

def partition_entropy(partition):

    """
        Computes the entropy of a partition. Used for creating weighted Hasse diagrams used in
        BVI computation.

        Inputs:
        partition - the partition of variables (frozenset of frozensets)

        Outputs:
        entropy - the partition entropy
    """
    
    # List for block sizes
    block_sizes = [len(i) for i in partition]

    # Random variables
    rv = np.asarray(block_sizes) / sum(block_sizes)

    # Entropy
    entropy = - sum(rv * np.log(rv) / np.log(2))

    return entropy

def cover_check(pi, sigma):

    """
        Checks whether sigma covers pi by checking that sigma is obtained from pi by merging 
        exactly two blocks.

        Inputs:
        pi - the bond partition to check for being covered (frozenset of frozensets)
        sigma - the bond partition to check for covering pi (frozenset of frozensets)

        Outputs:
        True/False depending on if sigma covers pi or not
    """

    # Checks if number of blocks match the cover relation
    if len(sigma) != len(pi) - 1:
        return False

    verdicts = []
    seen = set()

    # Loop over blocks in pi
    for block_1 in pi:
        seen.add(block_1) # Keep a list of blocks we have done in outer loop to avoid duplicates

        for block_2 in pi:
            if block_2 not in seen: # Preventing duplicates

                # Check whether this union of blocks of pi gives a merged block in sigma
                block_3 = block_1 | block_2
                check_1 = block_3 in sigma

                # If two blocks merge to give block in sigma, check that all other blocks of sigma
                # are blocks of pi.
                if check_1 is True:
                    check_2 = [block_4 in pi for block_4 in sigma if block_4 != block_3]
                    verdict = check_1 and all(check_2)
                    verdicts.append(verdict)

    # Only returns true if strictly one such pair of blocks exists
    if sum(verdicts) == 1:
        return True

    else:
        return False

def hasse_creator(G, algo, weight=False):

    """
        Constructs the Hasse diagram of the bond lattice of G. Used for BVI calculation.

        Inputs:
        G - the graph (NetworkX object)
        algo - the chosen bond-finding algorithm
        weight - whether or not to include edge weights given by the entropy differences
    """

    # Retrieve bonds, initialise Hasse diagram object, and add bond nodes
    bonds = algo(G)
    hasse = nx.Graph()
    for bond in bonds:
        hasse.add_node(bond)

    # Check cover relations and add edges if cover relations satisfied
    for bond in bonds:
        for bondd in bonds:
            verdict = cover_check(bond, bondd)
            if verdict == True and weight == False:
                hasse.add_edge(bond, bondd)

            # Add weighted edge if required
            elif verdict == True and weight == True:
                            hasse.add_edge(
                                bond, bondd,
                                weight = abs(partition_entropy(bondd) - partition_entropy(bond)))

    return hasse


def partition_meet(p_1, p_2):

    """
        Takes the partition meet of two partitions.

        Inputs:
        p_1, p_2 - partitions (both frozensets of frozensets)

        Outputs:
        meet - the partition meet of p_1 and p_2 (frozenset of frozensets)
    """

    meet_list = []

    # Check all intersections of blocks in p_1 and p_2, keeping only those that are non-empty
    for b_1 in p_1:
        for b_2 in p_2:
            inter = b_1.intersection(b_2)
            if len(inter) != 0: 
                meet_list.append(inter)

    # Ensuring output as frozenset of frozensets
    meet = frozenset(block for block in meet_list)

    return meet

def bond_meet(p_1, p_2, G, algo):

    """
        Finds the bond meet of two bonds p_1 and p_2. This is the greatest bond partition that 
        refines both p_1 and p_2. 

        Inputs:
        p_1, p_2 - partitions (both frozensets of frozensets).
        G - the graph (NetworkX object).
        algo - the chosen bond-finding algorithm.

        Outputs:
        meet - the bond meet of p_1 and p_2 (frozenset of frozensets).
    """
    
    # Get bonds and initialise refinement candidate list
    bonds = algo(G)
    candidates = []
    
    # Loop over bonds and blocks in bonds
    for bond in bonds:
        block_verdicts = []

        for block in bond:

            # Check whether this block is contained in some block of p_1 and p_2
            in_p_1 = any([block <= block_1 for block_1 in p_1])
            in_p_2 = any([block <= block_2 for block_2 in p_2])

            # Check contained in some block of both p_1 and p_2
            verdict = in_p_1 and in_p_2
            block_verdicts.append(verdict)

        # If all blocks of the bond are contained in some block of p_1 and p_2, append to refinment
        # candidate list
        if all(block_verdicts) == True:
            candidates.append(bond)

   # The greatest lower bound will always be at the end of the list because bonds is sorted rank order
   # and this must be unique
    meet = candidates[-1]

    return meet

def partition_sublattice_checker(G, algo):

    """
        Checks whether the bond lattice of a graph G is a sublattice of the partition lattice, by
        checking whether the partition meet of all bond pairs is in the bond lattice (closed under
        partition meet).

        Inputs:
        G - the graph (NetworkX object).
        algo - the chosen bond-finding algorithm.

        Outputs:
        False if any partition meet of bonds does not define a bond. Otherwise true
    """

    bonds = algo(G)
    
    for p_1 in bonds:
        for p_2 in bonds:
            part_meet = partition_meet(p_1, p_2)

            verdict = part_meet in bonds

            if verdict is False:
                return False

    return True

def is_block_graph(G):

    """
        Checks if a graph is a block graph by checking that all biconnected components are cliques.

        Inputs:
        G - the graph (NetworkX object)

        Outputs:
        True/False - verdict of whether or not the graph is a block graph
    """

    if G.number_of_nodes() <= 1:
        return True

    for bi_comp in nx.biconnected_components(G):
        subgraph = G.subgraph(bi_comp)
        k = len(bi_comp)
        if subgraph.number_of_edges() != k * (k - 1) // 2:  # is it a clique?
            return False
    return True