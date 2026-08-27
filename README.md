# Bond Lattice, Bond Information and Bond Variation of Information

This repository provides various programs for the Bond Lattice, it's Moebius coefficients, and the two measures derived
from the bond lattice that have been used in my thesis. All graph-theoretic objects have been coded using the NetworkX package.
Modules contains scripts with useful functions for these purposes, which can be imported into notebooks to when required.
Any figure or result in the write-up has it's own Notebook. Let us now go over each file and it's contents in detail. 
For details on specific functions, check the docstrings and comments in the Modules scripts.

## Modules

#### bond_finders
Contains three different algorithms to find the bonds of a given graph, and functions required for them. brute_force_bond_finder finds all possible partitions of the vertices of the graph, and evaluates
whether each block induces a connected subgraph. edge_bond_finder finds all possible edge subsets and the corresponding bonds of these subsets, taking care with uniqueness (various edge subsets can map to the same bond). 
walk_bond_finder uses a sliding-window type algorithm. Given some edge subset, it finds all edge subsets that can be produced
by connecting separate components of the edge subset. It then applies the same process to those new edge subsets, and hence iteratively
produces all bonds from the bottom-up, taking care to not consider duplicate edge subsets and bonds. All of these algorithms require
only the graph itself (as a NetworkX Graph object) as input.

#### gaussian_bond_information
Contains various functions used for generation of Gaussian data respecting some bond structure, as well as computation of the
analytic formula for the Bond Information with Gaussian data.

#### general_bond_information
Contains a function to calculate the bond information on for a general graph-structured dataset, requiring only the observed dataset itself
and the graph structure.

#### lattice_machinery
Contains various functions that allow for lattice-theoretic operations, including the creation of the Hasse diagram of the
bond lattice of any graph (unweighted or weighted if required for BVI calculation), partition and bond meets of partitions/bonds,
and a function that checks if the bond lattice of a given graph is a sublattice of the partition lattice.

#### moebius_machinery
Contains functions required for calculating the Moebius coefficients of the bond lattice of a given graph, with separate 
functions for using either a full matrix inversion, or just a matrix-vector solve if one only requires the Moebius coefficients
from each bond partition to the top of the lattice.

