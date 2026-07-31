import numpy as np
import networkx as nx
from Modules.brute_force_bond_finder import brute_force_bond_finder, bond_sorter

def partition_meet(p_1, p_2):

    """
        Takes the partition meet of two partitions p_1 and p_2, each defined as a frozen set of
        frozen sets.
    """

    meet_list = []

    for b_1 in p_1:
        for b_2 in p_2:
            inter = set(b_1).intersection(set(b_2))
            if len(inter) != 0: 
                meet_list.append(list(inter))

    return frozenset(frozenset(block) for block in meet_list)

def bond_meet(p_1, p_2, G, algo):

    """
        Finds the bond meet of two bonds p_1 and p_2. This is the greatest bond partition that refines
        both p_1 and p_2. algo is the chosen bond-finding algorithm.
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
    return candidates[-1]

def partition_sublattice_checker(G, algo):

    """
        Checks whether the bond lattice of a graph G is a sublattice of the partition lattice, by
        checking whether the partition meet of all bond pairs is in the bond lattice (closed under
        partition meet).
    """

    bonds = algo(G)
    
    for p_1 in bonds:
        for p_2 in bonds:
            part_meet = partition_meet(p_1, p_2)

            verdict = part_meet in bonds

            if verdict is False:
                return False

    return True