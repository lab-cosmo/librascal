import json

from ..neighbourlist import get_neighbourlist
from .base import RepresentationFactory, FeatureFactory
from ..utils import get_full_name
from ..neighbourlist.structure_manager import convert_to_structure
from ..neighbourlist.base import NeighbourListFactory
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
                 gaussian_sigma_constant=0., n_species=1,
                 method='thread', n_workers=1, disable_pbar=False):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'sphericalexpansion'
        self.hypers = dict()
        self.update_hyperparameters(
            interaction_cutoff=interaction_cutoff,
            cutoff_smooth_width=cutoff_smooth_width,
            max_radial=max_radial, max_angular=max_angular,
            gaussian_sigma_type=gaussian_sigma_type,
            gaussian_sigma_constant=gaussian_sigma_constant,
            n_species=n_species
        )

        self.nl_options = [
            dict(name='centers', args=[]),
            dict(name='neighbourlist', args=[interaction_cutoff]),
            dict(name='strict', args=[interaction_cutoff])
        ]

        neighbourlist_full_name = get_full_name(self.nl_options)
        self.name = self.name + '_' + neighbourlist_full_name
        hypers_str = json.dumps(self.hypers)
        self.rep_options = dict(name=self.name, args=[hypers_str])

        n_features = self.get_num_coefficients()
        self.feature_options = dict(name='blocksparse_double', args=[
                                    n_features, hypers_str])

        self.misc = dict(method=method, n_workers=n_workers,
                         disable_pbar=disable_pbar)

        self.manager = NeighbourListFactory(self.nl_options)
        self.representation = RepresentationFactory(
            self.manager, self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_constant', 'n_species'}
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
        FeatureManager.blocksparse_double
            Object containing the representation

        """
        structures = [convert_to_structure(frame) for frame in frames]

        n_atoms = [0]+[len(structure['atom_types'])
                       for structure in structures]
        structure_ids = np.cumsum(n_atoms)[:-1]
        n_centers = np.sum(n_atoms)

        features = FeatureFactory(self.feature_options)

        for structure in structures:
            self.manager.update(**structure)
            self.representation.compute()
            features.append(self.representation)

        return features

    def get_num_coefficients(self):
        """Return the number of coefficients in the spherical expansion

        (this is the descriptor size per atomic centre)

        """
        return (self.hypers['n_species'] * self.hypers['max_radial']
                * (self.hypers['max_angular'] + 1)**2)
