from ..lib import NeighbourList
import numpy as np

_neighbourlist_list = ["centers","neighbourlist","strict","maxorder","halflist","fulllist"]

_neighbourlists = {}
for k,v in NeighbourList.__dict__.items():
    if "make_adapted_manager" in k or "make_structure_manager" in k:
        name = k.lower().replace('make_adapted_manager_','').replace('make_structure_manager_','')
        _neighbourlists[name] = v

def NeighbourListFactory(name,*args):
    if name not in _neighbourlists:
        raise NameError('The neighbourlist factory {} has not been registered. The available combinations are: {}'.format(name,list(_neighbourlists.keys())))
    return _neighbourlists[name](*args)

def is_valid_structure(structure):
    keys = ['cell','positions','atom_types','pbc']
    is_valid = True
    if isinstance(structure, dict):
        for k in keys:
            if k not in structure:
                is_valid = False
    else:
        is_valid = False

    return is_valid

def adapt_structure(cell, positions, atom_types, pbc):
    cell = np.array(cell.T,order='F')
    positions = np.array(positions.T,order='F')
    atom_types = atom_types.reshape(-1,1)
    pbc = pbc.reshape(3,1)
    return dict(cell=cell,positions=positions,atom_types=atom_types,pbc=pbc)

def unpack_ase(frame):
    """
    Convert ASE Atoms object to rascal's equivalent

    Parameters
    ----------
    frame : ase.Atoms
        Atomic structure

    Returns
    -------
    StructureManagerCenters
        base structure manager.
    """
    cell = frame.get_cell()
    positions = frame.get_positions()
    numbers = frame.get_atomic_numbers()
    pbc = frame.get_pbc().astype(int)

    return adapt_structure(cell=cell,positions=positions,atom_types=numbers,pbc=pbc)
