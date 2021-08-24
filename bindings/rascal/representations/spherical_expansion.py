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


class SphericalExpansion(BaseIO):
    """
    Computes the spherical expansion of the neighbour density [soap]

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
        vary. Only fixed smearing width ('Constant') are implemented.

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

    radial_basis : string
        Specifies the type of radial basis R_n to be computed
        ("GTO" for Gaussian typed orbitals and "DVR" discrete variable representation using Gaussian quadrature rule)

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

    expansion_by_species_method : string
        Specifies the how the species key of the invariant are set-up.
        Possible values: 'environment wise', 'user defined', 'structure wise'.
        The descriptor is computed for each atomic enviroment and it is indexed
        using tuples of atomic species that are present within the environment.
        This index is by definition sparse since a species tuple will be non
        zero only if the atomic species are present inside the environment.
        'environment wise' means that each environmental representation
        will only contain the minimal set of species tuples needed by each
        atomic environment.
        'structure wise' means that within a structure the species tuples
        will be the same for each environment coefficients.
        'user defined' uses global_species to set-up the species tuples.

        These different settings correspond to different trade-off between
        the memory efficiency of the invariants and the computational
        efficiency of the kernel computation.
        When computing a kernel using 'environment wise' setting does not allow
        for efficent matrix matrix multiplications which is ensured when
        'user defined' is used. 'structure wise' is a balance between the
        memory footprint and the use of matrix matrix products.

        Note that the sparsity of the gradient coefficients and their use to
        build kernels does not allow for clear efficiency gains so their
        sparsity is kept irrespective of expansion_by_species_method.

    global_species : list
        list of species to use to set-up the species key of the invariant. It
        should contain all the species present in the structure for which
        invariants will be computed

    compute_gradients : bool
        control the computation of the representation's gradients w.r.t. atomic
        positions.

    cutoff_function_parameters : dict
        Additional parameters for the cutoff function.
        if cutoff_function_type == 'RadialScaling' then it should have the form

        .. code:: python

            dict(rate=...,
                 scale=...,
                 exponent=...)

        where :code:`...` should be replaced by the desired positive float.

    Methods
    -------
    transform(frames)
        Compute the representation for a list of ase.Atoms object.


    .. [soap] Bartók, Kondor, and Csányi, "On representing chemical
        environments", Phys. Rev. B. 87(18), p. 184115
        http://link.aps.org/doi/10.1103/PhysRevB.87.184115

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
        radial_basis="GTO",
        optimization=None,
        optimization_args=None,
        expansion_by_species_method="environment wise",
        global_species=None,
        compute_gradients=False,
        cutoff_function_parameters=dict(),
    ):
        """Construct a SphericalExpansion representation

        Required arguments are all the hyperparameters named in the
        class documentation
        """

        self.name = "sphericalexpansion"
        self.hypers = dict()

        if global_species is None:
            global_species = []
        elif not isinstance(global_species, list):
            global_species = list(global_species)

        self.update_hyperparameters(
            max_radial=max_radial,
            max_angular=max_angular,
            expansion_by_species_method=expansion_by_species_method,
            global_species=global_species,
            compute_gradients=compute_gradients,
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
            gaussian_sigma=dict(value=gaussian_sigma_constant, unit="A"),
        )
        optimization = check_optimization_for_spherical_representations(
            optimization, optimization_args
        )

        radial_contribution = dict(type=radial_basis, optimization=optimization)

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

        n_features = self.get_num_coefficients()

        self._representation = CalculatorFactory(self.rep_options)

    def update_hyperparameters(self, **hypers):
        """Store the given dict of hyperparameters

        Also updates the internal json-like _representation

        """
        allowed_keys = {
            "interaction_cutoff",
            "cutoff_smooth_width",
            "max_radial",
            "max_angular",
            "gaussian_sigma_type",
            "gaussian_sigma_constant",
            "gaussian_density",
            "cutoff_function",
            "radial_contribution",
            "compute_gradients",
            "cutoff_function_parameters",
            "expansion_by_species_method",
            "global_species",
        }
        hypers_clean = {key: hypers[key] for key in hypers if key in allowed_keys}
        self.hypers.update(hypers_clean)

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
        """Return the number of coefficients in the spherical expansion

        (this is the descriptor size per atomic centre)

        """
        return (
            n_species
            * self.hypers["max_radial"]
            * (self.hypers["max_angular"] + 1) ** 2
        )

    def get_keys(self, species):
        """
        return the proper list of keys used to build the representation
        """
        keys = []
        for sp in species:
            keys.append([sp])
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
            expansion_by_species_method=self.hypers["expansion_by_species_method"],
            global_species=self.hypers["global_species"],
            compute_gradients=self.hypers["compute_gradients"],
            gaussian_sigma_type=gaussian_density["type"],
            gaussian_sigma_constant=gaussian_density["gaussian_sigma"]["value"],
            cutoff_function_type=cutoff_function["type"],
            radial_basis=radial_contribution["type"],
            optimization=radial_contribution["optimization"],
            cutoff_function_parameters=self.cutoff_function_parameters,
        )
        return init_params

    def _set_data(self, data):
        super()._set_data(data)
        # allows to load deprecated models
        if "cpp_representation" in data.keys():
            self._representation = self._representation.from_dict(
                data["cpp_representation"]
            )
        else:
            print(
                "WARNING: a deprecated model was loaded. Key "
                "'cpp_representation' was not found in model. Please dump and "
                "reload the model to update it. The model parameters will not "
                "change, only the format."
            )

    def _get_data(self):
        data = super()._get_data()
        data.update(cpp_representation=self._representation.to_dict())
        return data
