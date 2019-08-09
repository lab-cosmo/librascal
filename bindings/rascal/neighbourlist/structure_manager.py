import numpy as np
from .base import NeighbourListFactory, is_valid_structure, adapt_structure


def get_neighbourlist(structure, options):
    manager = NeighbourListFactory(options)
    manager.update(**structure)
    return manager


def convert_to_structure(frame):
    if is_valid_structure(frame):
        structure = frame
    else:
        if is_ase_Atoms(frame):
            structure = unpack_ase(frame)
        else:
            raise RuntimeError(
                'Cannot convert structure of type {}'.format(type(frame)))
    return structure


def is_ase_Atoms(frame):
    is_ase = True
    if not hasattr(frame, 'get_cell'):
        is_ase = False
    if not hasattr(frame, 'get_positions'):
        is_ase = False
    if not hasattr(frame, 'get_atomic_numbers'):
        is_ase = False
    if not hasattr(frame, 'get_pbc'):
        is_ase = False
    return is_ase


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

    if "center_atoms_mask" in frame.arrays.keys():
        center_atoms_mask = frame.get_array("center_atoms_mask")
    else:
        center_atoms_mask = np.ones_like(numbers, dtype=bool)

    return adapt_structure(cell=cell, positions=positions, atom_types=numbers, pbc=pbc, center_atoms_mask=center_atoms_mask)
