import numpy as np
import json

from ..utils import get_strict_neighbourlist
from ..lib import RepresentationManager, FeatureManager
from .base import RepresentationFactory

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
                 gaussian_sigma_constant=-1.)
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'sphericalexpansion'
        #TODO is this really necessary if the underlying C++ object
        #     stores all the hypers in dict-like format anyway?
        self.hypers_dict = dict()
        self.update_hyperparameters(
            interaction_cutoff=interaction_cutoff,
            cutoff_smooth_width=cutoff_smooth_width,
            max_radial=max_radial, max_angular=max_angular,
            gaussian_sigma_type=gaussian_sigma_type,
            gaussian_sigma_constant=gaussian_sigma_constant)

    def update_hyperparameters(self, **hypers)
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_width'}
        hypers_clean = {key: hypers[key] if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self,frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------
        FeatureManager.Dense_double
            Object containing the representation

        """
        n_frames = len(frames)

        managers = list(map(
            get_strict_neighbourlist,
            frames, [self.hypers['interaction_cutoff'], ]*n_frames))

        hypers_str = json.dumps(self.hypers)

        n_features = self.get_num_components()
        features = FeatureManager.Dense_double(n_features, hypers_str)

        cms = map(RepresentationFactory(self.name,self.options),
                     managers,[inp]*Nframe)

        for cm in cms:
            cm.compute()
            features.append(cm)

        return features

    #TODO deprecate in favour of get_feature_size() and get_center_size()
    #     already exposed through the bindings?
    def get_num_components(self):
        return 

