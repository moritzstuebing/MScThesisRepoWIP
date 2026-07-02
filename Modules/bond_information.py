import numpy as np
from Modules.moebius_machinery import mat_vec_mu_graph

# Taken from Tim's repo

def shuffle_data(data, shuffle_indices):
    # Create a copy of the original data to avoid modifying the input in place.
    shuffled_data = data.copy()

    for indices in shuffle_indices:
        # Generate a random permutation for the indices.
        permutation = np.random.permutation(shuffled_data.shape[0])

        # Apply the permutation to the specified indices.
        for idx in indices:
            shuffled_data[:, idx] = shuffled_data[permutation, idx]

    return shuffled_data

def div_simplifier(bond):

    """
        Simplifies the argument of the divergence by cancelling singleton indices appearing in both numerator and denominator.
    """

    numerator = []

    for block in bond:
        if len(block) != 1:
            numerator.append(block)

    # Turn the correct indices into singletons for fully factorised denominator
    denominator = [[i] for block in numerator for i in block]

    return numerator, denominator

def bond_information(G, X, div, algo):
    
    """
        Calculates the bond information for a given graph G, dataset X, divergence measure div, and bond finding algorithm algo.
    """

    # Retrieve all bonds of the graph and initialise array for divergence storing
    #bonds = algo(G)

    # Calculate the mu coefficients of the given graph
    bonds, _, mu = mat_vec_mu_graph(G, algo)
    divs = np.zeros(len(bonds))

    # Check that shapes match
    assert mu.shape[0] == len(bonds)

    # Data shuffling for denominator
    indices = [[i] for i in range(X.shape[1])]
    X_shuffled_denom = shuffle_data(X, indices)

    # Loop over all bonds
    for i, bond in enumerate(bonds):

        # Simplify the argument of the divergence (remove singleton indices)
        num, denom = div_simplifier(bond)

        # Introduce required numerator independence
        X_shuffled_num = shuffle_data(X, num)

        # Turn denom into one list for indexing
        denom_one_list = [i for block in denom for i in block]

        # Guard against pure singleton case
        if len(denom_one_list) == 0:
            divs[i] = 0.0

        else:
            # Calculate divergence
            divs[i] = div(X_shuffled_num[:, denom_one_list], X_shuffled_denom[:, denom_one_list])
        
    # Return sum of divergences weighted by mu's
    return np.sum(mu * divs)