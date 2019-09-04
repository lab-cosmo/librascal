import json

from .base import CalculatorFactory
from ..neighbourlist import AtomsList
import numpy as np


class SphericalExpansion(object):

    """
    Computes the spherical expansion of the neighbour density [soap]

    Hyperparameters
    ----------
    interaction_cutoff : float
        Maximum pairwise distance for atoms to be considered in
        expansion

    cutoff_smooth_width : float
        The distance over which the the interaction is smoothed to zero

    max_radial : int
        Number of radial basis functions

    max_angular : int
        Highest angular momentum number (l) in the expansion

    n_species : int
        Number of species to be considered separately

    gaussian_sigma_type : str
        How the Gaussian atom sigmas (smearing widths) are allowed to
        vary -- fixed ('Constant'), by species ('PerSpecies'), or by
        distance from the central atom ('Radial').

    gaussian_sigma_constant : float
        Specifies the atomic Gaussian widths, in the case where they're
        fixed.

    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [soap] Bartók, Kondor, and Csányi, "On representing chemical
        environments", Phys. Rev. B. 87(18), p. 184115
        http://link.aps.org/doi/10.1103/PhysRevB.87.184115

    """

    def __init__(self, interaction_cutoff, cutoff_smooth_width,
                 max_radial, max_angular, gaussian_sigma_type,
                 gaussian_sigma_constant=0., cutoff_function_type="Cosine",
                 n_species=1, radial_basis="GTO",
                 method='thread', n_workers=1, disable_pbar=False):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'sphericalexpansion'
        self.hypers = dict()
        self.update_hyperparameters(
            max_radial=max_radial, max_angular=max_angular,
            n_species=n_species
        )

        cutoff_function = dict(
            type=cutoff_function_type,
            cutoff=dict(
                value=interaction_cutoff,
                unit='A'
            ),
            smooth_width=dict(
                value=cutoff_smooth_width,
                unit='A'
            ),
        )
        gaussian_density = dict(
            type=gaussian_sigma_type,
            gaussian_sigma=dict(
                value=gaussian_sigma_constant,
                unit='A'
            ),
        )
        radial_contribution = dict(
            type=radial_basis,
        )
        self.update_hyperparameters(cutoff_function=cutoff_function,
                                    gaussian_density=gaussian_density,
                                    radial_contribution=radial_contribution)

        self.nl_options = [
            dict(name='centers', args=dict()),
            dict(name='neighbourlist', args=dict(cutoff=interaction_cutoff)),
            dict(name="centercontribution", args=dict()),
            dict(name='strict', args=dict(cutoff=interaction_cutoff))
        ]

        hypers_str = json.dumps(self.hypers)
        self.rep_options = dict(name=self.name, args=[hypers_str])

        n_features = self.get_num_coefficients()
        self.feature_options = dict(name='blocksparse_double', args=[
                                    n_features, hypers_str])

        self.misc = dict(method=method, n_workers=n_workers,
                         disable_pbar=disable_pbar)

        self._representation = CalculatorFactory(self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like _representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_constant', 'n_species', 'gaussian_density', 'cutoff_function',
                        'radial_contribution'}
        hypers_clean = {key: hypers[key] for key in hypers
                        if key in allowed_keys}
        self.hypers.update(hypers_clean)

    def transform(self, frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------

            Object containing the representation

        """
        if not isinstance(frames, AtomsList):
            frames = AtomsList(frames, self.nl_options)

        self._representation.compute(frames.managers)
        return frames

    def get_num_coefficients(self):
        """Return the number of coefficients in the spherical expansion

        (this is the descriptor size per atomic centre)

        """
        return (self.hypers['n_species'] * self.hypers['max_radial']
                * (self.hypers['max_angular'] + 1)**2)
