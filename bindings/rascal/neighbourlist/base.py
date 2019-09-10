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

_structure_collections = {}
for k, v in NeighbourList.__dict__.items():
    if "ManagerCollection" in k:
        name = k.lower().replace('managercollection_', '')
        _structure_collections[name] = v

def NeighbourListFactory(nl_options):
    names = []
    kargs = []
    full_name = []
    for opt in nl_options:
        full_name.insert(0, opt['name'])
        name = '_'.join(full_name)
        names.append(name)
        kargs.append(opt['args'])

        if name not in _neighbourlists:
            raise NameError('The neighbourlist factory {} has not been registered. The available combinations are: {}'.format(
                name, list(_neighbourlists.keys())))

    managers = [_neighbourlists[names[0]](**kargs[0])]
    for name, karg in zip(names[1:], kargs[1:]):
        manager = _neighbourlists[name](managers[-1], **karg)
        managers.append(manager)

    return managers[-1]

def StructureCollectionFactory(nl_options):
    import json

    args = []
    full_name = []
    for opt in nl_options:
        full_name.insert(0, opt['name'])
        name = '_'.join(full_name)
        args.append(dict(name=opt['name'], initialization_arguments=opt['args']))

    if name not in _structure_collections:
        raise NameError('The StructureCollection factory {} has not been registered. The available combinations are: {}'.format(
            name, list(_structure_collections.keys())))
    # remove the arguments relative to stucture manager centers
    args.pop(0)
    agrs_str = json.dumps(args)
    structure_collection = _structure_collections[name](agrs_str)

    return structure_collection


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
