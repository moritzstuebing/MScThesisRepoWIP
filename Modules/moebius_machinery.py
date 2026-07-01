from Modules.next_closure_bond_finder import next_closure_algo
import numpy as np

# Needed for well-defined set inclusion for partitions
def set_maker(pi):

    """
        Turns a list of lists into a frozen set of a frozen sets.
    """

    return frozenset(frozenset(block) for block in pi)

def zeta_matrix(partitions):

    """
        Calculates the zeta matrix for any partition lattice or sublattice, ordered by refinement
    """

    # Automatically non-singular since upper triangular

    # Turn partitions and blocks within partitions into frozen sets so that set inclusion operations work
    bonds = [set_maker(pi) for pi in partitions]

    # Get number of bonds and initialise zeta matrix
    n = len(partitions)
    zeta = np.zeros((n, n))

    # Refinement check / filling in zeta matrix
    for i, pi in enumerate(bonds):
        for j, sigma in enumerate(bonds): # can possibly simplify further by working out how many bonds at each rank, then checking less bonds
            verdict = all(any(b <= B for B in sigma) for b in pi)
            #print(pi, sigma, verdict)
            if verdict == True:
                zeta[i, j] = 1

    return zeta

def full_inv_mu(partitions):

    """
        For a given set of partitions in refinement order, calculates the zeta and moebius matrices, 
        and returns the required Moebius coefficients.
    """

    zeta = zeta_matrix(partitions)
    mu_mat = np.linalg.inv(zeta)
    mu_values = mu_mat[:, -1]

    return partitions, zeta, mu_mat, mu_values

def full_inv_mu_graph(G):

    """
        For a specified graph G, finds the bonds ordered by refinement, calculates the zeta matrix, calculates the moebius
        matrix by inversion, extracts the required moebius coefficients
    """
    bonds = next_closure_algo(G)

    return full_inv_mu(bonds)

def mat_vec_solve_mu(partitions):
    
    """
        For some set of partitions in refinement order, calculates the zeta matrix, and returns the required Moebius
        coefficients by a matrix-vector solve
    """

    zeta = zeta_matrix(partitions)
    e = np.zeros(len(partitions)); e[len(partitions) - 1] = 1
    mu_values = np.linalg.solve(zeta, e)

    return partitions, zeta, mu_values

def mat_vec_mu_graph(G, algo):

    
    """
        For a specified graph G, finds the bonds ordered by refinement, calculates the zeta matrix, extracts the required 
        moebius coefficients by a matrix-vector solve.
    """

    bonds = algo(G)

    return mat_vec_solve_mu(bonds)