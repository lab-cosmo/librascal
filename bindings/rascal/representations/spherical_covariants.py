import json

from .base import (
    CalculatorFactory,
    cutoff_function_dict_switch,
    check_optimization_for_spherical_representations,
)
from ..neighbourlist import AtomsList
import numpy as np
from ..utils import BaseIO
from copy import deepcopy


class SphericalCovariants(BaseIO):
    """
    Computes a SphericalCovariants representation, i.e. lambda spectrum.

    Attributes
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

    cutoff_function_type : string
        Choose the type of smooth cutoff function used to define the local
        environment. Can be either 'ShiftedCosine' or 'RadialScaling'.

        If 'ShiftedCosine', the functional form of the switching function is:

        .. math::

            sc(r) = \\begin{cases}
            1 &r < r_c - sw,\\\\
            0.5 + 0.5 \cos(\pi * (r - r_c + sw) / sw) &r_c - sw < r <= r_c, \\\\
            0 &r_c < r,
            \\end{cases}

        where :math:`r_c` is the interaction_cutoff and :math:`sw` is the
        cutoff_smooth_width.

        If 'RadialScaling', the functional form of the switching function is
        as expressed in equation 21 of https://doi.org/10.1039/c8cp05921g:

        .. math::

            rs(r) = sc(r) u(r),

        where

        .. math::

            u(r) = \\begin{cases}
            \\frac{1}{(r/r_0)^m} &\\text{if c=0,}\\\\
            1 &\\text{if m=0,} \\\\
            \\frac{c}{c+(r/r_0)^m} &\\text{else},
            \\end{cases}

        where :math:`c` is the rate, :math:`r_0` is the scale, :math:`m` is the
        exponent.

    normalize : boolean
        Whether to normalize so that the kernel between identical environments
        is 1.  Default and highly recommended: True.

    radial_basis :  string
        Specifies the type of radial basis R_n to be computed
        ("GTO" for Gaussian typed orbitals and "DVR" discrete variable
        representation using Gauss-Legendre quadrature rule)

    soap_type : string
        Specifies the type of representation to be computed.

    inversion_symmetry : boolean
        Specifies whether inversion invariance should be enforced.

    covariant_lambda : int
        Order of the lambda spectrum.

    cutoff_function_parameters : dict
        Additional parameters for the cutoff function.
        if cutoff_function_type == 'RadialScaling' then it should have the form

        .. code:: python

            dict(rate=...,
                 scale=...,
                 exponent=...)

        where :code:`...` should be replaced by the desired positive float.

    optimization : dict, default None
        Optional arguments for optimization of the computation of spherical
        expansion coefficients. "Spline" and "RadialDimReduction" are available.

        Spline: Enables cubic splining for the radial basis functions.

            accuracy : float
                accuracy of the cubic spline

        RadialDimReduction: Projection matrices to optimize radial basis,
                            requires Spline to be set

            projection_matrices : dict
                Contains or each species a list of projection matrices for each
                angular channel. A projection matrix for an angular channel has
                the shape (max_radial, expanded_max_radial). A number of
                `expanded_max_radial` radial basis are computed
                to be then projected to `max_radial` radial basis. The projected
                radial basis is then splined for each species and angular channel

        Default settings is using spline

        .. code: python

            dict(Spline=dict(accuracy=1e-8))


    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [lambda-soap] Grisafi, A., Wilkins, D. M., Csányi, G., & Ceriotti, M.
    (2018). Symmetry-Adapted Machine Learning for Tensorial Properties of
    Atomistic Systems. Physical Review Letters, 120(3), 036002.
    https://doi.org/10.1103/PhysRevLett.120.036002

    """

    def __init__(
        self,
        interaction_cutoff,
        cutoff_smooth_width,
        max_radial,
        max_angular,
        gaussian_sigma_type,
        gaussian_sigma_constant=0.3,
        cutoff_function_type="ShiftedCosine",
        normalize=True,
        radial_basis="GTO",
        optimization=None,
        optimization_args=None,
        soap_type="LambdaSpectrum",
        inversion_symmetry=True,
        covariant_lambda=0,
        cutoff_function_parameters=dict(),
    ):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """
        self.name = "sphericalcovariants"
        self.hypers = dict()
        self.update_hyperparameters(
            max_radial=max_radial,
            max_angular=max_angular,
            soap_type=soap_type,
            normalize=normalize,
            inversion_symmetry=inversion_symmetry,
            covariant_lambda=covariant_lambda,
        )

        self.cutoff_function_parameters = deepcopy(cutoff_function_parameters)

        cutoff_function_parameters.update(
            interaction_cutoff=interaction_cutoff,
            cutoff_smooth_width=cutoff_smooth_width,
        )
        cutoff_function = cutoff_function_dict_switch(
            cutoff_function_type, **cutoff_function_parameters
        )

        gaussian_density = dict(
            type=gaussian_sigma_type,
            gaussian_sigma=dict(value=gaussian_sigma_constant, unit="AA"),
        )
        optimization = check_optimization_for_spherical_representations(
            optimization, optimization_args
        )

        radial_contribution = dict(type=radial_basis, optimization=optimization)
        radial_contribution = dict(
            type=radial_basis,
        )

        self.update_hyperparameters(
            cutoff_function=cutoff_function,
            gaussian_density=gaussian_density,
            radial_contribution=radial_contribution,
        )

        self.nl_options = [
            dict(name="centers", args=dict()),
            dict(name="neighbourlist", args=dict(cutoff=interaction_cutoff)),
            dict(name="centercontribution", args=dict()),
            dict(name="strict", args=dict(cutoff=interaction_cutoff)),
        ]

        self.rep_options = dict(name=self.name, args=[self.hypers])

        self._representation = CalculatorFactory(self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like representation

        """
        allowed_keys = {
            "interaction_cutoff",
            "cutoff_smooth_width",
            "max_radial",
            "max_angular",
            "gaussian_sigma_type",
            "gaussian_sigma_constant",
            "soap_type",
            "inversion_symmetry",
            "covariant_lambda",
            "cutoff_function",
            "normalize",
            "gaussian_density",
            "radial_contribution",
            "cutoff_function_parameters",
        }
        hypers_clean = {key: hypers[key] for key in hypers if key in allowed_keys}
        self.hypers.update(hypers_clean)
        return

    def transform(self, frames):
        """Compute the representation.

        Parameters
        ----------
        frames : list(ase.Atoms) or AtomsList
            List of atomic structures.

        Returns
        -------
           AtomsList : Object containing the representation

        """
        if not isinstance(frames, AtomsList):
            frames = AtomsList(frames, self.nl_options)

        self._representation.compute(frames.managers)
        return frames

    def get_num_coefficients(self, n_species=1):
        """Return the number of coefficients in the representation

        (this is the descriptor size per atomic centre)

        """
        if self.hypers["soap_type"] == "LambdaSpectrum":
            if self.hypers["inversion_symmetry"]:
                n_col = (
                    np.ceil((self.hypers["max_angular"] + 1) ** 2 / 2.0)
                    - (1.0 + np.floor((self.hypers["covariant_lambda"] - 1) / 2.0)) ** 2
                    - np.floor(
                        (
                            self.hypers["max_angular"]
                            + 1
                            - self.hypers["covariant_lambda"]
                        )
                        ** 2
                        / 2.0
                    )
                    * (self.hypers["covariant_lambda"] % 2)
                    - (
                        np.ceil(
                            (
                                self.hypers["max_angular"]
                                + 1
                                - self.hypers["covariant_lambda"]
                            )
                            ** 2
                            / 2.0
                        )
                        - (
                            self.hypers["max_angular"]
                            - self.hypers["covariant_lambda"]
                            + 1
                        )
                    )
                    * (1.0 - self.hypers["covariant_lambda"] % 2)
                )
                if self.hypers["covariant_lambda"] % 2 == 1:
                    n_col = -n_col + 0.5 * (
                        2.0
                        + self.hypers["covariant_lambda"]
                        - 3 * self.hypers["covariant_lambda"] ** 2
                        + 2 * self.hypers["max_angular"]
                        + 4
                        * self.hypers["covariant_lambda"]
                        * self.hypers["max_angular"]
                    )
                n_col *= 2 * self.hypers["covariant_lambda"] + 1
                return int(n_col * n_species ** 2 * self.hypers["max_radial"] ** 2)
            else:
                return (
                    n_species ** 2
                    * self.hypers["max_radial"] ** 2
                    * int(
                        (
                            2
                            + self.hypers["covariant_lambda"]
                            - 3 * self.hypers["covariant_lambda"] ** 2
                            + 2 * self.hypers["max_angular"]
                            + 4
                            * self.hypers["covariant_lambda"]
                            * self.hypers["max_angular"]
                        )
                        / 2
                    )
                    * (2 * self.hypers["covariant_lambda"] + 1)
                )
        else:
            raise ValueError("Only soap_type = LambdaSpectrum " "implemented for now")

    def get_keys(self, species):
        """
        return the proper list of keys used to build the representation
        """
        keys = []
        if self.hypers["soap_type"] == "LambdaSpectrum":
            for sp1 in species:
                for sp2 in species:
                    if sp1 > sp2:
                        continue
                    keys.append([sp1, sp2])
        else:
            raise ValueError("Only soap_type = LambdaSpectrum " "implemented for now")

        return keys

    def _get_init_params(self):
        gaussian_density = self.hypers["gaussian_density"]
        cutoff_function = self.hypers["cutoff_function"]
        radial_contribution = self.hypers["radial_contribution"]

        init_params = dict(
            interaction_cutoff=cutoff_function["cutoff"]["value"],
            cutoff_smooth_width=cutoff_function["smooth_width"]["value"],
            max_radial=self.hypers["max_radial"],
            max_angular=self.hypers["max_angular"],
            soap_type=self.hypers["soap_type"],
            inversion_symmetry=self.hypers["inversion_symmetry"],
            normalize=self.hypers["normalize"],
            gaussian_sigma_type=gaussian_density["type"],
            gaussian_sigma_constant=gaussian_density["gaussian_sigma"]["value"],
            lam=self.hypers["covariant_lambda"],
            cutoff_function_type=cutoff_function["type"],
            radial_basis=radial_contribution["type"],
            optimization=radial_contribution["optimization"],
            cutoff_function_parameters=self.cutoff_function_parameters,
        )
        return init_params

    def _set_data(self, data):
        super()._set_data(data)

    def _get_data(self):
        return super()._get_data()
