import json

from ..neighbourlist import get_neighbourlist
from .base import RepresentationFactory, FeatureFactory
from ..utils import get_full_name
from ..neighbourlist.structure_manager import convert_to_structure
from ..neighbourlist.base import NeighbourListFactory
import numpy as np


class SOAPCovariant(object):

    """
    Computes a SOAPCovariant representation, e.g. lambda spectrum.

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

    soap_type : string
        Specifies the type of representation to be computed.

    inversion_symmetry : Boolean
        Specifies whether inversion invariance should be enforced.

    lam : int
        Order of the lambda spectrum.

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
                 cutoff_function_type="Cosine", normalize=True,
                 radial_basis="GTO",
                 soap_type="LambdaSpectrum", inversion_symmetry=True,
                 lam=0):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = 'soapcovariant'
        self.hypers = dict()
        self.update_hyperparameters(
            max_radial=max_radial, max_angular=max_angular,
            n_species=n_species,
            soap_type=soap_type,
            normalize=normalize,
            inversion_symmetry=inversion_symmetry,
            lam=lam)

        cutoff_function = dict(
            type=cutoff_function_type,
            cutoff=dict(
                value=interaction_cutoff,
                unit='AA'
            ),
            smooth_width=dict(
                value=cutoff_smooth_width,
                unit='AA'
            ),
        )
        gaussian_density = dict(
            type=gaussian_sigma_type,
            gaussian_sigma=dict(
                value=gaussian_sigma_constant,
                unit='AA'
            ),
        )
        radial_contribution = dict(
            type=radial_basis,
        )

        self.update_hyperparameters(cutoff_function=cutoff_function,
                                    gaussian_density=gaussian_density,
                                    radial_contribution=radial_contribution)

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

        self.manager = NeighbourListFactory(self.nl_options)
        self.representation = RepresentationFactory(
            self.manager, self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_constant', 'n_species', 'soap_type',
                        'inversion_symmetry', 'lam', 'cutoff_function', 'normalize', 'gaussian_density', 'radial_contribution'}
        hypers_clean = {key: hypers[key] for key in hypers
                        if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self, frames, features=None):
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
        if features is None:
            features = FeatureFactory(self.feature_options)

        structures = [convert_to_structure(frame) for frame in frames]

        n_atoms = [0]+[len(structure['atom_types'])
                       for structure in structures]
        structure_ids = np.cumsum(n_atoms)[:-1]
        n_centers = np.sum(n_atoms)

        ii = 0
        for structure in structures:
            try:
                self.manager.update(**structure)
            except:
                print("Structure NL {} failed".format(ii))

            try:
                self.representation.compute()
            except:
                print("Structure Rep computation {} failed".format(ii))

            try:
                features.append(self.representation)
            except:
                print("Structure data gather {} failed".format(ii))

            ii += 1

        return features

    def get_num_coefficients(self):
        """Return the number of coefficients in the representation

        (this is the descriptor size per atomic centre)

        """
        if self.hypers['soap_type'] == 'LambdaSpectrum':
            if self.hypers['inversion_symmetry'] == True:
                n_col = np.ceil((self.hypers['max_angular'] + 1)**2/2.0) - \
                    (1.0 + np.floor((self.hypers['lam'] - 1)/2.0))**2 - \
                    np.floor((self.hypers['max_angular'] + 1 -
                              self.hypers['lam'])**2/2.0)*(self.hypers['lam'] % 2) - \
                    (np.ceil((self.hypers['max_angular'] + 1 -
                              self.hypers['lam'])**2/2.0) -
                     (self.hypers['max_angular'] -
                      self.hypers['lam'] + 1)) * \
                    (1.0 - self.hypers['lam'] % 2)
                if (self.hypers['lam'] % 2 == 1):
                    n_col = -n_col + 0.5*(2.0 + self.hypers['lam'] -
                                          3*self.hypers['lam']**2 +
                                          2*self.hypers['max_angular'] +
                                          4*self.hypers['lam']*self.hypers['max_angular'])
                n_col *= (2*self.hypers['lam'] + 1)
                return int(n_col*self.hypers['n_species']**2*self.hypers['max_radial']**2)
            else:
                return (self.hypers['n_species']**2*self.hypers['max_radial']**2 *
                        int((2 + self.hypers['lam'] - 3*self.hypers['lam']**2 +
                             2*self.hypers['max_angular'] +
                             4*self.hypers['lam']*self.hypers['max_angular'])/2) *
                        (2*self.hypers['lam'] + 1))
        else:
            raise RuntimeError('Only soap_type = LambdaSpectrum '
                               'implemented for now')
