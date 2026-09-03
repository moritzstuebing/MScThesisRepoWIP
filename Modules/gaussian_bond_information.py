import numpy as np
import networkx as nx
import itertools
from Modules.bond_finders import *
from Modules.moebius_machinery import mat_vec_solve_mu

def covariance_matrix_graph(G, rho):

    """
        Generate covariance matrix that respects the depoendence structure introduced by the
        graph G.

        Inputs:
        G - the graph.
        rho - correlation strength between dependent variables.

        cov - covariance matrix (2D np.array) that respects the dependence structure of the graph.
    """

    # Initialise array for covariance
    cov = np.zeros((len(list(G.nodes)), len(list(G.nodes))))

    # Fill in correlation values
    for e in G.edges:
        cov[e[0], e[1]] = rho
        cov[e[1], e[0]] = rho

    # Fill in diagonals (a variable is always perfectly correlated with itself)
    for i in range(len(list(G.nodes))):
        cov[i, i] = 1

    return cov


def covariance_matrix_bond(bond, rho):

    """
        Generates the covariance matrix respecting the factorisation structure of a given bond.

        Inputs:
        bond - the given bond that the covariance matrix's dependence structure should respect.
        rho - correlation strength between dependent variables.

        Outputs:
        sigma - covariance matrix (2D np.array) that respects the dependence structure of the bond.
    """

    # Initialise array for covariance
    d = sum([1 for block in bond for b in block])
    sigma = np.zeros((d, d))

    # Loop over all blocks
    for block in bond:

        # Generate all pairs in the block
        if len(block) > 1:
            pairs = itertools.combinations(block, 2)

            # Fill in covariance matrix entries
            for pair in pairs:
                sigma[pair[0], pair[1]] = rho
                sigma[pair[1], pair[0]] = rho

    for i in range(d):
        sigma[i, i] = 1

    return sigma

def generate_gaussian(sigma, n_samples):

    """
        Generates multivariate Gaussian data that respects the independence/factorisation
        structure induced by the covariance matrix cov.

        Inputs:
        sigma - the covariance matrix with expected dependence structure encoded.
        n_samples - the number of samples we want to return.

        Outputs:
        data - the generated multivariate Gaussian data (2D np.array), means zero and covariance
                specified by cov.
    """

    # Zero mean for each variable
    mean = np.zeros(sigma.shape[0])

    # Form data
    data = np.random.multivariate_normal(mean, sigma, n_samples)
    
    return data

def restrict(sigma, bond):

    """
        Further restricts a known covariance matrix giving factorisation structure of a joint
        (Gaussian) distribution to give the covariance for factorisation restriction given by some 
        given bond.

        Inputs:
        sigma - covariance matrix (2D np.array) encoding overall dependence structure of the 
                joint distribution to be further factorised according to bond. 
        bond - the bond according to which the joint distribution should be further factorised
                according to.

        Outputs:
        sigma_bond - the covariance matrix that has been further factorised according to the bond.
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
        Finds the analytic Gaussian bond information. Works because TA divergence has a closed
        form on Gaussian data.

        Inputs:
        G - the graph
        algo - the desired bond finding algorithm
        sigma - the covariance matrix encoding the overall factorisation structure of the
                gaussian data (this gets further restricted for each bond)
        alpha - parameter in TA divergence
        
        Output:
        bi - the analytic bond information found for this theoretical Gaussian distribution
    """

    # Get bonds and moebius coefficients
    bonds = algo(G)
    _, _, mu = mat_vec_solve_mu(G, algo)

    # Initialise running sum
    bi = 0.0

    # Main calculation loop
    for i, bond in enumerate(bonds):

        # Further restrict the covariance according to the current bond
        cov_pi = restrict(sigma, bond)

        # Calculate the analytic TA divergence term for this bond
        div = (
            (np.linalg.det(np.linalg.inv(cov_pi))) ** (alpha / 2)) / (
                np.linalg.det(alpha * np.linalg.inv(cov_pi) + (1 - alpha) * np.eye(
                    cov_pi.shape[0]))) ** (1/2)

        # Add the divergence term multiplied by the Moebius coefficient to the running total
        bi += mu[i] * ((div - 1.0) / (alpha - 1.0))

    return bi

def bond_nonbond_list_maker(G, algo):

    """
        This is a function to sort the list of partitions for a d-vertex graph into bonds and 
        non-bonds of the graph G.

        Inputs:
        G - the graph (NetworkX object).
        algo - the desired bond-finding algorithm/function.

        Outputs:
        bonds - the bonds of G, list of frozensets of frozensets.
        partitions - all possible partitions of the vertices of G, bond or not. list of frozensets
                        of frozensets.
        sep_list - organised list of partitions of vertices of G, bonds first and non-bonds second.
    """

    # Retrieve all bonds
    d = len(G.nodes)
    partitions = algo(nx.complete_graph(d))
    bonds = algo(G)

    # Initiate separated list with all the bonds of the graph
    sep_list = [bond for bond in bonds]

    # Now append non-bond partitions
    for partition in partitions:
        if partition not in bonds:
            sep_list.append(partition)

    return bonds, partitions, sep_list

def factorised_bi(G, algo, rho):

    """
        Calculates the Bond Information for all factorisations of the joint distribution, given
        by all partitions of the vertices of G. First calculates for bond factorisations, then 
        non-bond factorisations.

        Inputs:
        G - the graph (NetworkX object).
        algo - the desired bond-finding algorithm/function.
        rho - correlation strength of dependent variables, scalar.

        Outputs:
        bond_facts - dictionary giving the bond and the value of analytic gaussian bond
                        information when respecting the dependence structure of this bond
        nonbond_facts - dictionary giving the non-bond partition and the value of analytic 
                        gaussian bond information when respecting the dependence structure of this 
                        partition

    """

    # Get bonds and separated list of bonds / non-bonds
    bonds, _, sep_list = bond_nonbond_list_maker(G, algo)

    # Initialise dictinaries for value storing
    bond_facts = {}
    nonbond_facts = {}

    # Calculate bond informations
    for bond in sep_list[:len(bonds)]:

        # Covariance to induce joint factorisation according to current bond
        sigma = covariance_matrix_bond(G, bond, rho)

        # Calculate and add to dictionary
        bond_facts[bond] = (analytic_gaussian_bi(G, brute_force_bond_finder, sigma))

    for part in sep_list[len(bonds):]:

        # Covariance to induce joint factorisation according to current bond. Correlation strength between 
        # dependent variables = 0.99
        sigma = covariance_matrix_bond(G, part, rho)
        
        # Calculate and add to dictionary
        nonbond_facts[part] = (analytic_gaussian_bi(G, brute_force_bond_finder, sigma))

    return bond_facts, nonbond_facts


def bond_to_label(bond):

    """
        Generates LaTeX label from bond given as a frozenset of frozensets.

        Inputs:
        bond - the bond to be displayed (frozenset of frozensets).

        Outputs:
        a nice string displaying the factorisation of the joint distribution.

    """

    # Sort the blocks in order of the first element in each block
    blocks = sorted(
        (sorted(v + 1 for v in block) for block in bond),
        key= min,
    )
    body = "".join("P_{" + "".join(map(str, b)) + "}" for b in blocks)
    return "$" + body + "$"