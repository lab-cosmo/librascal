import json

from .base import CalculatorFactory
from ..neighbourlist import AtomsList
import numpy as np


class SphericalCovariants(object):

    """
    Computes a SphericalCovariants representation, i.e. lambda spectrum.

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


    .. [lambda-soap] Grisafi, A., Wilkins, D. M., Csányi, G., & Ceriotti, M.
    (2018). Symmetry-Adapted Machine Learning for Tensorial Properties of
    Atomistic Systems. Physical Review Letters, 120(3), 036002.
    https://doi.org/10.1103/PhysRevLett.120.036002

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
        self.name = 'sphericalcovariants'
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
            dict(name='centers', args=dict()),
            dict(name='neighbourlist', args=dict(cutoff=interaction_cutoff)),
            dict(name="centercontribution", args=dict()),
            dict(name='strict', args=dict(cutoff=interaction_cutoff))
        ]

        hypers_str = json.dumps(self.hypers)
        self.rep_options = dict(name=self.name, args=[hypers_str])

        self._representation = CalculatorFactory(self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {'interaction_cutoff', 'cutoff_smooth_width',
                        'max_radial', 'max_angular', 'gaussian_sigma_type',
                        'gaussian_sigma_constant', 'n_species', 'soap_type',
                        'inversion_symmetry', 'lam', 'cutoff_function',
                        'normalize', 'gaussian_density', 'radial_contribution'}
        hypers_clean = {key: hypers[key] for key in hypers
                        if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self, frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms)
            List of atomic structures.

        Returns
        -------


        """
        if not isinstance(frames, AtomsList):
            frames = AtomsList(frames, self.nl_options)

        self._representation.compute(frames.managers)
        return frames

    def get_num_coefficients(self):
        """Return the number of coefficients in the representation

        (this is the descriptor size per atomic centre)

        """
        if self.hypers['soap_type'] == 'LambdaSpectrum':
            if self.hypers['inversion_symmetry'] == True:
                n_col = (np.ceil((self.hypers['max_angular'] + 1)**2/2.0) -
                         (1.0 + np.floor((self.hypers['lam'] - 1)/2.0))**2 -
                         np.floor((self.hypers['max_angular'] + 1 -
                                   self.hypers['lam'])**2/2.0)
                         * (self.hypers['lam'] % 2) -
                         (np.ceil((self.hypers['max_angular'] + 1 -
                                   self.hypers['lam'])**2/2.0) -
                          (self.hypers['max_angular'] -
                             self.hypers['lam'] + 1)) *
                         (1.0 - self.hypers['lam'] % 2))
                if (self.hypers['lam'] % 2 == 1):
                    n_col = -n_col + 0.5*(2.0 + self.hypers['lam'] -
                                          3 * self.hypers['lam']**2 +
                                          2 * self.hypers['max_angular'] +
                                          4 * self.hypers['lam']
                                            * self.hypers['max_angular'])
                n_col *= (2*self.hypers['lam'] + 1)
                return int(n_col * self.hypers['n_species']**2
                           * self.hypers['max_radial']**2)
            else:
                return (self.hypers['n_species']**2 *
                        self.hypers['max_radial']**2 *
                        int((2 + self.hypers['lam'] - 3*self.hypers['lam']**2 +
                             2 * self.hypers['max_angular'] +
                             4 * self.hypers['lam']
                             * self.hypers['max_angular'])/2) *
                        (2*self.hypers['lam'] + 1))
        else:
            raise ValueError('Only soap_type = LambdaSpectrum '
                             'implemented for now')
