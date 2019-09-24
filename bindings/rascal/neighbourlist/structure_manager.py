from .base import (NeighbourList, NeighbourListFactory, is_valid_structure,
 adapt_structure, StructureCollectionFactory)
from collections.abc import Iterable


class AtomsList(object):
    """
    A wrapper class for a stack of managers precompiled on the C++ side of the form Strict->NeighbourList->Center.
    A container for atoms/centers/atomic environments.

    Attributes
    ----------
    nl_options : dict
        Parameters for each layer of the wrapped structure manager. Parameters can be specified for these layers: center, neighbourlist and strict.

    Methods
    -------
    """
    def __init__(self, frames, nl_options, start=None, length=None):
        self.nl_options = nl_options
        self._frames = frames

        if isinstance(frames, str):
            # if filename
            managers = StructureCollectionFactory(nl_options)
            if start is None and length is None:
                managers.add_structures(frames)
            elif start is not None and length is None:
                managers.add_structures(frames, start=start)
            elif start is not None and not length is None:
                managers.add_structures(frames, start=start, length=length)
            elif start is None and not length is None:
                managers.add_structures(frames, length=length)
        else:
            # if python structure
            structures = convert_to_structure_list(frames)
            managers = StructureCollectionFactory(nl_options)
            try:
                managers.add_structures(structures)
            except:
                print("""Neighbourlist of structures failed. trying
                one at a time.""")
                ii = 0
                for structure, manager in zip(structures, managers):
                    try:
                        manager.update(structure)
                    except:
                        print("Structure Rep computation {} failed".format(ii))

        self.managers = managers

    def __iter__(self):
        return self.managers
    def __getitem__(self, key):
        return self.managers[key]
    def get_dense_feature_matrix(self, calculator):
        """
        Parameters
        -------
        calculator : Calculator (an object owning a _representation object)

        Returns
        -------
        represenation_matrix : ndarray
            returns the representation bound to the calculator as dense matrix.
        """
        return self.managers.get_dense_feature_matrix(calculator._representation)


def get_neighbourlist(structure, options):
    manager = NeighbourListFactory(options)
    manager.update(**structure)
    return manager


def convert_to_structure_list(frames):
    if not isinstance(frames, Iterable):
        frames = [frames]
    structure_list = NeighbourList.AtomicStructureList()
    for frame in frames:
        if is_valid_structure(frame):
            structure = frame
        else:
            if is_ase_Atoms(frame):
                structure = unpack_ase(frame)
            else:
                raise RuntimeError(
                    'Cannot convert structure of type {}'.format(type(frame)))
        structure_list.append(**structure)
    return structure_list


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

    return adapt_structure(cell=cell, positions=positions, atom_types=numbers, pbc=pbc)
