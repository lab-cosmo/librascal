"""
A collection of utility functions to manipulate the radial basis.
"""

from scipy.special import legendre, gamma
from copy import deepcopy
import numpy as np
import rascal.representations as representations


def radial_basis_functions_dvr(
    radial_grid, max_radial, interaction_cutoff, gaussian_sigma
):
    """
    Evaluate DVR radial basis function for R_n on a grid.
    These are equivalent to the basis used to expand the density in the
    librascal DVR implementation. Implementation follows reference

    radial_grid : array_like
        R_n is evaluated on the grid (0, r_c + 3 * σ]

    max_radial : int
        max_radial in rascal

    interaction_cutoff : float
        interaction_cutoff in rascal

    gaussian_sigma : float
        gaussian_sigma_constan in rascal

    Returns:
    -------
    R_n : array_like
        R_n is evaluated on the `radial_grid`

    Reference:
    ----------
    Light, J. C., & Carrington Jr, T. (2000). Discrete-variable
    representations and their utilization. Advances in Chemical Physics,
    114, 263-310. http://light-group.uchicago.edu/dvr-rev.pdf
    """

    legendre_points, legendre_weights = np.polynomial.legendre.leggauss(max_radial)
    grid_end = interaction_cutoff + 3 * gaussian_sigma
    assert np.min(radial_grid) > 0
    assert np.max(radial_grid) <= grid_end
    legendre_grid = (radial_grid - grid_end / 2) / (grid_end / 2)
    legendre_f = [
        legendre(n_radial) * np.sqrt((2 * n_radial + 1) / 2)
        for n_radial in range(max_radial)
    ]
    # Transformation matrix in the reference
    T = np.array(
        [
            [
                np.sqrt(legendre_weights[alpha]) * legendre_f[j](legendre_points[alpha])
                for alpha in range(max_radial)
            ]
            for j in range(max_radial)
        ]
    )
    legendre_f = np.array(
        [legendre_f[n_radial](legendre_grid) for n_radial in range(max_radial)]
    )
    rR_n = T.T @ legendre_f
    # rascal does the expansion for r*R_n to obtain the right density on the
    # power spectrum level, therefore we divide here this r factor out
    R_n = rR_n / radial_grid
    return R_n


def gto_sigma(cutoff, n, n_max):
    """
    Compute GTO sigma

    ---Arguments---
    cutoff: environment cutoff
    n: order of the GTO
    n_max: maximum order of the GTO

    ---Returns---
    sigma: GTO sigma parameter
    """
    return np.maximum(np.sqrt(n), 1) * cutoff / n_max


def gto_width(sigma):
    """
    Compute GTO width

    ---Arguments---
    sigma: GTO sigma parameter

    ---Returns---
    b: GTO (Gaussian) width
    """
    return 1.0 / (2 * sigma ** 2)


def gto_prefactor(n, sigma):
    """
    Compute GTO prefactor

    ---Arguments---
    n: order of the GTO
    sigma: GTO sigma parameter

    ---Returns---
    N: GTO prefactor (normalization factor)
    """
    return np.sqrt(2 / (sigma ** (2 * n + 3) * gamma(n + 1.5)))


def gto(r, n, sigma):
    """
    Compute GTO

    ---Arguments---
    r: grid on which to evaluate the GTO
    n: order of the GTO
    sigma: GTO sigma parameter

    ---Returns---
    R_n: GTO of order n evaluated on the provided grid
    """
    b = gto_width(sigma)
    N = gto_prefactor(n, sigma)
    return N * r ** (n + 1) * np.exp(-b * r ** 2)  # why n+1?


def gto_overlap(n, m, sigma_n, sigma_m):
    """
    Compute overlap of two GTOs

    ---Arguments---
    n: order of the first GTO
    m: order of the second GTO
    sigma_n: sigma parameter of the first GTO
    sigma_m: sigma parameter of the second GTO

    ---Returns---
    S: overlap of the two GTOs
    """
    b_n = gto_width(sigma_n)
    b_m = gto_width(sigma_m)
    N_n = gto_prefactor(n, sigma_n)
    N_m = gto_prefactor(m, sigma_m)
    nm = 0.5 * (3 + n + m)
    return 0.5 * N_n * N_m * (b_n + b_m) ** (-nm) * gamma(nm)


def gto_S(hypers):
    """
    Compute the GTO overlap

    ---Arguments---
    hypers: dictionary of SOAP hyperparameters
    r_grid: grid of radial distances on which to evaluate the GTOs

    ---Returns---
    R_n: orthogonalized GTOs evaluated on `r_grid`
    """

    # Setup grids of the expansion orders
    n_grid = np.arange(0, hypers["max_radial"])
    sigma_grid = gto_sigma(hypers["interaction_cutoff"], n_grid, hypers["max_radial"])

    # Compute radial normalization factor based on the GTO overlap
    S = gto_overlap(
        n_grid[:, np.newaxis],
        n_grid[np.newaxis, :],
        sigma_grid[:, np.newaxis],
        sigma_grid[np.newaxis, :],
    )
    return S


def radial_basis_functions_gto(radial_grid, max_radial, interaction_cutoff):
    """
    Evaluate GTO radial basis function for R_n on a grid.
    These are equivalent to the basis used to expand the density in the
    librascal GTO implementation.

    radial_grid : array_like
        R_n is evaluated on the grid (0, r_c + 3 * σ]

    max_radial : int
        max_radial in rascal

    interaction_cutoff : float
        interaction_cutoff in rascal

    Returns:
    -------
    R_n : array_like
        R_n is evaluated on the `radial_grid`
    """

    # Setup grids of the expansion orders
    n_grid = np.arange(0, max_radial)
    sigma_grid = gto_sigma(interaction_cutoff, n_grid, max_radial)

    # Compute radial normalization factor based on the GTO overlap
    S = gto_overlap(
        n_grid[:, np.newaxis],
        n_grid[np.newaxis, :],
        sigma_grid[:, np.newaxis],
        sigma_grid[np.newaxis, :],
    )
    # S = fractional_matrix_power(S, -0.5)
    eva, eve = np.linalg.eigh(S)
    S = eve @ np.diag(1 / np.sqrt(eva)) @ eve.T

    # Compute GTOs, shape (n_max, len(r_grid))
    R_n = np.matmul(
        S,
        gto(
            radial_grid[np.newaxis, :], n_grid[:, np.newaxis], sigma_grid[:, np.newaxis]
        ),
    )

    return R_n


def get_radial_basis_covariance(spherical_expansion, spherical_expansion_by_species):
    """
    spherical_expansion : array_like
        an instance of a spherical_expansion calculator that defines
        how the expansion coefficients should be computed

    spherical_expansion_by_species : dictionary
        a dictionary containing spherical expansion coefficients computed for each
        neighbor species, typically obtained by applying
        SphericalExpansion.transform().get_features_by_species to a set of atomic structures

    Returns:
    -------
    invariant_covariance_matrices : array_like
        covariance matrix for each species and angular channel
        of shape (n_species, max_angular+1, max_radial, max_radial)

    """

    # infers shape of the expansion coefficients
    n_environments = len(
        spherical_expansion_by_species[list(spherical_expansion_by_species)[0]]
    )
    n_species = len(spherical_expansion_by_species.keys())
    max_radial = spherical_expansion.hypers["max_radial"]
    max_angular = spherical_expansion.hypers["max_angular"]

    invariant_covariance_matrices = {}

    for species in spherical_expansion_by_species.keys():
        # reshape features of shape (n_environments, n_features) so we can access the species, radial and angular channels separately
        spherical_expansion_coefficients = spherical_expansion_by_species[
            species
        ].reshape((n_environments, max_radial, (max_angular + 1) ** 2))
        invariant_covariance_matrices[species] = np.zeros(
            (max_angular + 1, max_radial, max_radial)
        )

        # 1/n_envs \sum_i sum_m <anlm|ρ> <aklm|ρ> for environment i, species a, radial n k, angular l, magnetic m

        # awkard indexing handles the compact storage of the (l,m) terms in the spex coefficients
        for l in range(max_angular + 1):
            invariant_covariance_matrices[species][l] = np.tensordot(
                spherical_expansion_coefficients[:, :, l * l : (l + 1) * (l + 1)],
                spherical_expansion_coefficients[:, :, l * l : (l + 1) * (l + 1)],
                ([0, 2], [0, 2]),
            )
        invariant_covariance_matrices[species] /= n_environments

    return invariant_covariance_matrices


def get_radial_basis_pca(covariance_matrices):
    """
    Computes the PCA of the invariant covariance matrices, as obtained
    by get_radial_basis_covariance

    covariance_matrices : array_like
        list of invariant covariance matrices,
        shape (n_species, max_angular+1, max_radial, max_radial)

    Returns:
    -------
    principal_values : dictionary of array_like
        dictionary of principal values, shape shape (max_angular+1, max_radial), arranged
        by species
    principal_vectors  : dictionary of array_like
        dictionary of principal vectors (stored as columns),
        shape (n_species, max_angular+1, max_radial, max_radial)
        arranged by species

    """

    max_angular_plusone, max_radial, _ = list(covariance_matrices.values())[0].shape

    principal_values = {}
    principal_vectors = {}
    for species in covariance_matrices.keys():
        principal_values[species] = np.zeros((max_angular_plusone, max_radial))
        principal_vectors[species] = np.zeros(
            (max_angular_plusone, max_radial, max_radial)
        )
        for l in range(max_angular_plusone):
            covariance_matrix = covariance_matrices[species][l]
            eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
            principal_values[species][l] = eigenvalues[::-1]
            principal_vectors[species][l] = eigenvectors[:, ::-1]
    return principal_values, principal_vectors


def get_radial_basis_projections(principal_vectors, n_optimal):
    """
    Extracts the top n_optimal projectors from the principal vectors, as obtained
    by get_radial_basis_pca, and format them in a way that is suitable to be used
    as the projection_matrices hyperparameter of SphericalExpansion.

    principal_vectors  : array_like
        list of principal vectors (stored as columns)
        shape shape (n_species, max_angular+1, max_radial, max_radial)
    n_optimal : int
        the number of principal components that should be retained in the
        projection matrices

    Returns:
    -------
    projection matrics : dictionary of array_like
        dictionary of projection matrices retaining n_optimal basis functions,
        arranged by species, and formatted so that they can be passed to SphericalExpansion

    """
    max_angular_plusone, max_radial, _ = list(principal_vectors.values())[0].shape

    projection_matrices = {}
    for species in principal_vectors.keys():
        projection_matrices[species[0]] = []
        for l in range(max_angular_plusone):
            projection_matrices[species[0]].append(
                principal_vectors[species][l][:, :n_optimal].T.tolist()
            )
    return projection_matrices


def get_optimal_radial_basis_hypers(hypers, frames, expanded_max_radial=-1):
    """
    Helper function to compute an optimal radial basis following
    Goscinski et al, arxiv:2105.08717.

    hypers: dictionary
        hyperparameters for the desired representation. "max_radial" indicates
        the desired size of the optimal basis
    frames: ase.Atoms
        a list of structures used to estimate the optimal radial basis
    expanded_max_radial: int
        number of intermediate basis to be used to estimate the optimal basis.
        defaults to -1, in which case it is taken to be 2*max_radial

    Returns:
    -------
    optimal_hypers: dictionary
        hyperparameters including the optimal basis projectors
    """

    spherical_expansion_hypers = deepcopy(hypers)

    # removes parameters that don't make sense for a spherical expansion
    spherical_expansion_hypers.pop("normalize", None)
    spherical_expansion_hypers.pop("soap_type", None)
    spherical_expansion_hypers.pop("compute_gradients", None)

    if "optimization" in spherical_expansion_hypers:
        spherical_expansion_hypers["optimization"].pop("RadialDimReduction", None)

    if expanded_max_radial == -1:
        expanded_max_radial = 2 * hypers["max_radial"]
    spherical_expansion_hypers["max_radial"] = expanded_max_radial

    spex = representations.SphericalExpansion(**spherical_expansion_hypers)
    # computes density expansion coefficients
    feats = spex.transform(frames).get_features_by_species(spex)

    # computes covariance and principal components
    cov = get_radial_basis_covariance(spex, feats)
    p_val, p_vec = get_radial_basis_pca(cov)

    p_mat = get_radial_basis_projections(p_vec, hypers["max_radial"])

    # assemble the updated hypers
    optimal_hypers = deepcopy(hypers)
    if not "optimization" in optimal_hypers:
        optimal_hypers["optimization"] = {}
    optimal_hypers["optimization"] = {
        "RadialDimReduction": {"projection_matrices": p_mat},
    }

    if not "Spline" in optimal_hypers["optimization"]:
        optimal_hypers["optimization"]["Spline"] = {"accuracy": 1e-8}

    return optimal_hypers
