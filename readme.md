# The Bond Lattice: Structure and Information Measures

This repository develops various codes to analyse the Bond Lattice, as well as information-theoretic measures derived from it. All graph-theoretic requirements use NetworkX. The Modules folder provides various scripts that include various useful functions, and the Notebooks folder contains Jupyter notebook experiments and calculations to analyse the behaviours of the lattice and associated information measures. We now give a brief overview of each Modules file. For deeper explanations, check the docstrings and comments in corresponding scripts.

## Modules

#### bond_finders
Contains three algorithms of varying computational complexity, to find all bonds of any given graph.

#### gaussian_bond_information
Contains functions required to compute the Bond Information for Gaussian data, which has an analytic form.

#### general_bond_information
Contains a Bond Information function that can be used for any dataset with a graph structure.

#### lattice_machinery
Contains various lattice-theoretic functions. The Hasse creation function with weights is used to compute the Bond Variation of Information.

#### moebius_machinery
Functions to compute the Moebius coefficients of any given graph's bond lattice.

## Example Usage

The following is quick example illustrating how to use the main functions in this repository. We use the 5-vertex path graph $P_5$.

```python
import networkx as nx
from Modules.bond_finders import *
from Modules.gaussian_bond_information import *
from Modules.lattice_machinery import *
from Modules.moebius_machinery import *

# 5-vertex path graph P5
P_5 = nx.path_graph(5)
nx.draw(P_5)

# Retrieving the bonds using the brute force bond-finder, good for small graphs. For larger graphs use walk_bond_finder.
bonds = brute_force_bond_finder(G)
print(f"{len(bonds)} bond partitions")

# Calculating the analytic Bond Information for multivariate Gaussian data. The dependence/factorisation structure of the data is encoded through the covariance matrix. Note we cannot encode pure correlation since otherwise the covariance matrix becomes singular, so we use a correlation strength parameter rho < 1.

rho = 0.99
sigma_full = np.full((5, 5), 0.99)
bi = analytic_gaussian_bi(P_5, sigma_full, rho)
print(f"Bond Information of fully-dependent data: {bi:.4f}")

# Now we calculate the Bond Information when the Gaussian data factorises according to a bond partition (can also be calculated for any partition of the variables).

bonds = bonds[4] # choosing arbitrary bond
print(bond)
sigma_bond = covariance_matrix_bond(P_5, bond, rho) # bond-induced covariance
bi_fact = analytic_gaussian_bi(P_5, sigma_bond, rho)
print(f"Bond Information of factorised data: {bi_fact:.4f}")

# Now we show how to create the Hasse diagram of the bond lattice of a graph, which is itself a graph. We also show how to find the Bond Variation of Information between two bonds, and the path in the lattice producing it.

hasse = hasse_creator(G, brute_force_bond_finder, weight=True) # any of the algorithms can be used in the second argument
pi = frozenset(frozenset([0, 1]), frozenset([2, 3, 4]))
sigma = frozenset(frozenset([0]), frozenset([1, 2]), frozenset([3, 4]))
bvi = nx.shortest_path_length(hasse, source=pi, target=sigma, weight="weight")
print(f"Bond Variation of Information: {bvi}")
path = nx.shortest_path(hasse, source=pi, target=sigma, weight="weight")
print(f"Shortest path: {path}")
```