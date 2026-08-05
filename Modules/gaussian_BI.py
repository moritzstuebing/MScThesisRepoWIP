import numpy as np
import networkx as nx
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
    _, _, mu = mat_vec_solve_mu(bonds)

    # Initialise running sum
    sum = 0.0

    # Main calculation loop
    for i, bond in enumerate(bonds):
        cov_pi = restrict(sigma, bond)
        div = ((np.linalg.det(np.linalg.inv(cov_pi))) ** (alpha / 2)) / (np.linalg.det(alpha * np.linalg.inv(cov_pi) + (1 - alpha) * np.eye(cov_pi.shape[0]))) ** (1/2)
        sum += mu[i] * ((div - 1.0) / (alpha - 1.0))

    return sum

def bond_nonbond_list_maker(G, algo):

    """
        This is a function to sort the list of partitions for a d-vertex graph into bonds and non-bonds
        of the graph G.
    """

    d = len(G.nodes)
    partitions = algo(nx.complete_graph(d))
    bonds = algo(G)

    # Initiate separated set with all the bonds of C_4
    sep_list = [bond for bond in bonds]

    # Now append non-bond partitions
    for partition in partitions:
        if partition not in bonds:
            sep_list.append(partition)

    return bonds, partitions, sep_list

def factorised_bi(G, algo, rho):

    """
        Calculates the Bond Information for all factorisations of the joint distribution. First
        calculates for bond factorisations, then non-bond factorisations.
    """

    bonds, partitions, sep_list = bond_nonbond_list_maker(G, algo)

    bond_facts = {}
    nonbond_facts = {}

    # Calculate bond informations
    for bond in sep_list[:len(bonds)]:

        # Covariance to induce joint factorisation according to current bond. Correlation strength between 
        # dependent variables = 0.99
        sigma = covariance_matrix_bond(4, bond, rho)

        # Calculate and append
        bond_facts[bond] = (analytic_gaussian_bi(G, brute_force_bond_finder, sigma))

    for part in sep_list[len(bonds):]:

        # Covariance to induce joint factorisation according to current bond. Correlation strength between 
        # dependent variables = 0.99
        sigma = covariance_matrix_bond(4, bond, rho)
        
        # Calculate and append
        nonbond_facts[part] = (analytic_gaussian_bi(G, brute_force_bond_finder, sigma))

    return bond_facts, nonbond_facts


def bond_to_label(bond, offset=1):

    """
        Generates LaTeX label from bond given as a frozenset of frozensets
    """
    blocks = sorted(
        (sorted(v + offset for v in block) for block in bond),
        key=lambda b: b[0],
    )
    body = "".join("P_{" + "".join(map(str, b)) + "}" for b in blocks)
    return "$" + body + "$"