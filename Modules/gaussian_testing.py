import numpy as np
import itertools
from Modules.brute_force_bond_finder import brute_force_bond_finder
from Modules.moebius_machinery import mat_vec_solve_mu

### Covariance and gaussian generator never actually used

def covariance_matrix_graph(G, rho):

    """
        Generate covariance matrix that respects the independence structure introduced by the graph G. The correlation strength between
        variables is rho.
    """

    cov = np.zeros((len(list(G.nodes)), len(list(G.nodes))))
    for e in G.edges:
        cov[e[0], e[1]] = rho
        cov[e[1], e[0]] = rho

    for i in range(len(list(G.nodes))):
        cov[i, i] = 1

    return cov

def covariance_matrix_bond(d, bond, rho):

    """
        Generates the covariance matrix respecting the factorisation structure of a given bond. d is the number of variables,
        rho is the correlation strength between dependent variables
    """

    cov = np.zeros((d, d))

    # Loop over all blocks
    for block in bond:

        # Generate all pairs in the block
        if len(block) > 1:
            pairs = itertools.combinations(block, 2)

            # Fill in covariance matrix entries
            for pair in pairs:
                cov[pair[0], pair[1]] = rho
                cov[pair[1], pair[0]] = rho

    for i in range(d):
        cov[i, i] = 1

    return cov

def generate_gaussian(cov, n_samples):

    """
        Generates multivariate Gaussian data that respects the independence/factorisation structure induced by the covariance matrix cov
    """

    # Zero mean for each variable
    mean = np.zeros(cov.shape[0])
    
    # Return data
    return np.random.multivariate_normal(mean, cov, n_samples)

def restrict(sigma, bond):

    """
        Restricts a known covariance matrix giving factorisation structure of a joint (Gaussian) distribution to give the covariance for
        factorisation restricted by some given bond
    """

    # Copy to prevent overwriting
    sigma_bond = sigma.copy()

    # Number of variables
    d = sigma.shape[0]

    # Update covariance to match bond factorisation structure
    for block in bond:
        for i in block:
            for j in range(d):
                if j not in block:
                    sigma_bond[i, j] = 0

    return sigma_bond

def analytic_gaussian_bi(G, algo, sigma, alpha=0.5):

    """
        G the graph, algo the bond-finding algorithm, sigma gives the factorisation structure of the joint, alpha the specified parameter for 
        the TA divergence (normally 0.5). This function calculates the analytic bond information for some graph G when the dataset is
        multivariate Gaussian that factorises according to the covariance matrix sigma
    """

    # Get bonds and moebius coefficients
    bonds = algo(G)
    d = len(bonds[0])
    _, _, mu = mat_vec_solve_mu(bonds)

    # Initialise running sum
    sum = 0.0

    # Main calculation loop
    for i, bond in enumerate(bonds):
        cov_pi = restrict(sigma, bond)
        div = ((np.linalg.det(np.linalg.inv(cov_pi))) ** (alpha / 2)) / (np.linalg.det(alpha * np.linalg.inv(cov_pi) + (1 - alpha) * np.eye(cov_pi.shape[0]))) ** (1/2)
        sum += mu[i] * ((div - 1.0) / (alpha - 1.0))

    return sum