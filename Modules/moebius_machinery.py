from Modules.brute_force_bond_finder import *
import numpy as np

def zeta_matrix(G, algo):

    """
        Calculates the zeta matrix of bond lattice for any graph.

        Inputs:
        G - the graph (NetworkX object).
        algo - the chosen bond-finding algorithm.

        Outputs:
        zeta - matrix (2D np.array) where i,j entry = 1 if i'th bond refines j, 0 otherwise.
    """

    # Automatically non-singular since upper triangular

    # Retrieve bonds
    bonds = algo(G)
    # Get number of bonds and initialise zeta matrix
    zeta = np.zeros((len(bonds), len(bonds)))

    # Refinement check / filling in zeta matrix
    for i, pi in enumerate(bonds):
        for j, sigma in enumerate(bonds):
            verdict = all(any(b <= B for B in sigma) for b in pi)
            if verdict == True:
                zeta[i, j] = 1

    return zeta

def full_inv_mu(G, algo):

    """
        For a given graph G, calculates Moebius coefficients.

        Inputs:
        G - the graph (NetworkX object).
        algo - the chosen bond-finding algorithm.

        bonds - the bonds of G (list of frozensets of frozensets).
        zeta - the zeta matrix (2D np.array) of L_G.
        mu_mat - the mu matrix (2D np.array) of L_G, inverse of zeta.
        mu_values - the Moebius coefficients we are interested in, correspond to intervals from
                    each bond to the top of the lattice.
    """

    # Retrieve bonds and form zeta matrix
    bonds = algo(G)
    zeta = zeta_matrix(G, algo)

    # Form mu matrix via inversion
    mu_mat = np.linalg.inv(zeta)

    # Extract required Moebius coefficients for our purposes
    mu_values = mu_mat[:, -1]

    return bonds, zeta, mu_mat, mu_values

def mat_vec_solve_mu(G, algo):
    
    """
        For some graph, returns required Moebius coefficients via a matrix-vector system solve.

        Inputs:
        G - the graph (NetworkX object)
        algo - the chosen bond-finding algorithm

        Outputs:
        bonds - the bonds of G (list of frozensets of frozensets)
        zeta - the zeta matrix (2D np.array)
        mu_values - the required Moebius coefficients for our purposes
    """

    # Retrieve bonds and form zeta matrix
    bonds = algo(G)
    zeta = zeta_matrix(G, algo)

    # Form canonical basis vector and solve matrix-vector system
    e = np.zeros(len(bonds)); e[len(bonds) - 1] = 1
    mu_values = np.linalg.solve(zeta, e)

    return bonds, zeta, mu_values