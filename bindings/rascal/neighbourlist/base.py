from ..lib import NeighbourList
import numpy as np

_neighbourlist_list = ["centers", "neighbourlist",
                       "strict", "maxorder", "halflist", "fulllist"]

_neighbourlists = {}
for k, v in NeighbourList.__dict__.items():
    if "make_adapted_manager" in k or "make_structure_manager" in k:
        name = k.lower().replace('make_adapted_manager_',
                                 '').replace('make_structure_manager_', '')
        _neighbourlists[name] = v


def NeighbourListFactory(nl_options):
    names = []
    args = []
    full_name = []
    for opt in nl_options:
        full_name.insert(0, opt['name'])
        name = '_'.join(full_name)
        names.append(name)
        args.append(opt['args'])

        if name not in _neighbourlists:
            raise NameError('The neighbourlist factory {} has not been registered. The available combinations are: {}'.format(
                name, list(_neighbourlists.keys())))

    managers = [_neighbourlists[names[0]](*args[0])]
    for name, arg in zip(names[1:], args[1:]):
        manager = _neighbourlists[name](managers[-1], *arg)
        managers.append(manager)

    return managers[-1]


def is_valid_structure(structure):
    keys = ['cell', 'positions', 'atom_types', 'pbc']
    is_valid = True
    if isinstance(structure, dict):
        for k in keys:
            if k not in structure:
                is_valid = False
    else:
        is_valid = False

    return is_valid


def adapt_structure(cell, positions, atom_types, pbc):
    cell = np.array(cell.T, order='F')
    positions = np.array(positions.T, order='F')
    atom_types = atom_types.reshape(-1, 1)
    pbc = pbc.reshape(3, 1)
    return dict(cell=cell, positions=positions, atom_types=atom_types, pbc=pbc)
