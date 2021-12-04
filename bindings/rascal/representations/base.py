from ..lib import representation_calculators
from ..neighbourlist.base import NeighbourListFactory
from ..neighbourlist.structure_manager import convert_to_structure_list

import numpy as np
import queue

from copy import deepcopy


_supported_optimization = ["Spline", "RadialDimReduction"]
# Register Calculators
_representations_list = [
    "sortedcoulomb",
    "sphericalexpansion",
    "sphericalinvariants",
    "sphericalcovariants",
    "kspacesphericalexpansion",
    "kspacesphericalinvariants",
]
_representations = {}
for k, v in representation_calculators.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _representations_list:
            if name in kl:
                _representations[kl] = v


def CalculatorFactory(rep_options):
    name = rep_options["name"]
    if name not in _representations:
        raise NameError(
            (
                "The representations factory {} has not been registered. "
                + "The available combinations are: {}"
            ).format(name, list(_representations.keys()))
        )
    return _representations[name](*rep_options["args"])


def cutoff_function_dict_switch(cutoff_function_type, **kwargs):
    """
    return appropriate dict for the cutoff function parameters
    """
    if cutoff_function_type == "ShiftedCosine":
        cutoff_function_dict = dict(
            type=cutoff_function_type,
            cutoff=dict(value=kwargs["interaction_cutoff"], unit="AA"),
            smooth_width=dict(value=kwargs["cutoff_smooth_width"], unit="AA"),
        )
    elif cutoff_function_type == "RadialScaling":
        cutoff_function_dict = dict(
            type=cutoff_function_type,
            cutoff=dict(value=kwargs["interaction_cutoff"], unit="AA"),
            smooth_width=dict(value=kwargs["cutoff_smooth_width"], unit="AA"),
            rate=dict(value=kwargs["rate"], unit="AA"),
            scale=dict(value=kwargs["scale"], unit="AA"),
            exponent=dict(value=kwargs["exponent"], unit="AA"),
        )
    else:
        raise NotImplementedError(
            "cutoff_function: " + cutoff_function_type + " has not been implemented."
        )

    return cutoff_function_dict


def check_optimization_for_spherical_representations(optimization, optimization_args):
    """
    Checks if the arguments in optimzation have been set correctly
    for representation based on spherical expansion.
    """
    # Soft backwards compatibility (remove this whole if-statement after 01.11.2021)
    if optimization_args is not None:
        # TODO(veit) TODO(alex) replace with logger as soon as merged
        print(
            "Warning: The 'optimization_args' parameter is deprecated "
            "(see 'optimization' parameter instead).\n"
            "This message will become an error after 2021-11-01."
        )
        optimization = dict(Spline=dict(accuracy=optimization_args["accuracy"]))

    if optimization is None:
        optimization = dict(Spline=dict(accuracy=1e-8))
    optimization = deepcopy(optimization)
    # check supported optimization keys
    for key in optimization:
        if key not in _supported_optimization:
            print(
                f"Warning: Optimization argument {key} is not supported and therefore ignored."
            )
    if "RadialDimReduction" in optimization and "Spline" not in optimization:
        raise ValueError(
            "Optimization key `RadialDimReduction` can only be used with `Spline`, please set the spline arguments in `optimization`"
        )
    return optimization
