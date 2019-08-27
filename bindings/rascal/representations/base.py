from ..lib import RepresentationManager
from ..utils.pool_worker import FactoryPool
from ..neighbourlist.base import NeighbourListFactory
from ..neighbourlist.structure_manager import convert_to_structure_list

import numpy as np
import queue

_representations_list = ["sortedcoulomb", "sphericalexpansion",
                         "sphericalinvariants", "sphericalcovariants"]
_representations = {}
for k, v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _representations_list:
            if name in kl:
                _representations[kl] = v


def CalculatorFactory(rep_options):
    name = rep_options['name']
    if name not in _representations:
        raise NameError('The representations factory {} has not been registered. The available combinations are: {}'.format(
            name, list(_representations.keys())))
    return _representations[name](*rep_options['args'])

