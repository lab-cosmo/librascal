import numpy as np
from .base import NeighbourListFactory

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

    cell = np.array(cell.T,order='F')
    positions = np.array(positions.T,order='F')
    numbers = numbers.reshape(-1,1)
    pbc = pbc.reshape(3,1)
    return dict(cell=cell,positions=positions,atom_types=numbers,pbc=pbc)

def get_neighbourlist(frame,options):
    names = []
    args = []
    full_name = []
    for opt in options:
        full_name.insert(0,opt['name'])
        name = '_'.join(full_name)
        names.append(name)
        args.append(opt['args'])

    structure = unpack_ase(frame)

    managers = [NeighbourListFactory(names[0],*args[0])]
    for name,arg in zip(names[1:],args[1:]):
        manager = NeighbourListFactory(name,managers[-1],*arg)
        managers.append(manager)
    manager = managers[-1]
    manager.update(**structure)
    return managers[-1]
