import numpy as np
from .base import NeighbourListFactory, unpack_ase, is_valid_structure


def get_neighbourlist(frame, options):
    names = []
    args = []
    full_name = []
    for opt in options:
        full_name.insert(0, opt['name'])
        name = '_'.join(full_name)
        names.append(name)
        args.append(opt['args'])

    if not is_valid_structure(frame):
        structure = unpack_ase(frame)
    else:
        structure = frame

    managers = [NeighbourListFactory(names[0], *args[0])]
    for name, arg in zip(names[1:], args[1:]):
        manager = NeighbourListFactory(name, managers[-1], *arg)
        managers.append(manager)
    manager = managers[-1]
    manager.update(**structure)
    return managers[-1]


def get_neighbourlist_full_name(options):
    full_name = []
    for opt in options:
        full_name.insert(0, opt['name'])
    return '_'.join(full_name)
