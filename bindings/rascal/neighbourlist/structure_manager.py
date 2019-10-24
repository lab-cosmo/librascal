from collections.abc import Iterable
from functools import reduce
from operator import and_

import numpy as np

from ..lib import neighbour_list
from .base import (NeighbourListFactory, is_valid_structure,
                   adapt_structure, StructureCollectionFactory)

class AtomsList(object):
    """
    A wrapper class for a stack of managers precompiled on the C++ side of the
    form Strict->NeighbourList->Center.  A container for atoms/centers/atomic
    environments.

    Attributes
    ----------
    nl_options : dict
        Parameters for each layer of the wrapped structure manager. Parameters
        can be specified for these layers: center, neighbourlist and strict.

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
            elif start is not None and length is not None:
                managers.add_structures(frames, start=start, length=length)
            elif start is None and length is not None:
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
        return self.managers.get_dense_feature_matrix(
                calculator._representation)


def get_neighbourlist(structure, options):
    manager = NeighbourListFactory(options)
    manager.update(**structure)
    return manager


def convert_to_structure_list(frames):
    if not isinstance(frames, Iterable):
        frames = [frames]
    structure_list = neighbour_list.AtomicStructureList()
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

    If the frame has an ase.atoms.arrays entry called
    'center_atoms_mask' then it will be used as the center mask
    (surprise) for any representations computed on this
    StructureManager.
    """
    cell = frame.get_cell()
    positions = frame.get_positions()
    numbers = frame.get_atomic_numbers()
    pbc = frame.get_pbc().astype(int)

    if "center_atoms_mask" in frame.arrays.keys():
        center_atoms_mask = frame.get_array("center_atoms_mask")
    else:
        center_atoms_mask = np.ones_like(numbers, dtype=bool)

    return adapt_structure(cell=cell, positions=positions,
                           atom_types=numbers, pbc=pbc,
                           center_atoms_mask=center_atoms_mask)


def mask_center_atoms_by_id(frame, id_select=None, id_blacklist=None):
    """Mask the centers (center-select) of an ASE atoms object, by index

    Parameters
    ----------
    frame: ase.Atoms
        Atomic structure to mask

    id_select: list of int
        List of atom IDs to select

    id_blacklist: list of int
        List of atom IDs to exclude

    Returns
    -------
    None (the Atoms object is modified directly)

    Notes
    -----
    The default is to select all atoms.  If `id_select` is provided,
    select only those atoms.  If only `id_blacklist` is provided, select
    all atoms *except* those in the blacklist.  If both are provided,
    atoms are first selected based on `id_select` and then excluded based
    on `id_blacklist`.  If the atoms object already has a mask, then
    `id_select` is applied first using the `or` operation, then
    `id_blacklist` is applied using the `and not` operation (so the order
    of precedence is: blacklist, selection, previous mask).

    This logic allows this function to be combined with
    `mask_center_atoms_by_species` to allow mixed species/id-based
    masking.
    """
    if 'center_atoms_mask' not in frame.arrays:
        # add a default mask
        if id_select is not None:
            mask = np.zeros((frame.get_number_of_atoms(),), dtype='bool')
        else:
            mask = np.ones((frame.get_number_of_atoms(),), dtype='bool')
    else:
        mask = frame.arrays['center_atoms_mask']
    if id_select is not None:
        mask[id_select] = True
    if id_blacklist is not None:
        mask[id_blacklist] = False
    frame.arrays['center_atoms_mask'] = mask


def mask_center_atoms_by_species(frame, species_select=[],
                                 species_blacklist=[]):
    """Mask the centers of an ASE atoms object, by atomic species

    Parameters
    ----------
    frame: ase.Atoms
        Atomic structure to mask

    species_select: list of int or str
        List of atomic numbers, or species symbols, to select.
        Should be of consistent type across list.

    species_blacklist: list of int or str
        List of atomic numbers, or species symbols, to exclude.
        Should be of consistent type across list.

    Returns
    -------
    None (the Atoms object is modified directly)

    Notes
    -----
    The default is to select all atoms.  If `species_select` is
    provided, select only those atoms whose species is in the list.  If
    only `species_blacklist` is provided, select all atoms *except*
    those whose species is in the blacklist.  If both are provided,
    atoms are first selected based on `species_select` and then excluded
    based on `species_blacklist`.  If the atoms object already has a
    mask, then `species_select` is applied first using the `or`
    operation, then `species_blacklist` is applied using the `and not`
    operation (so the order of precedence is: blacklist, selection,
    previous mask).

    This logic allows this function to be combined with
    `mask_center_atoms_by_id` to allow mixed species/id-based masking.
    """
    select_is_str = reduce(and_,
                           map(lambda x: isinstance(x, str), species_select),
                           True)
    select_is_int = reduce(and_,
                           map(lambda x: isinstance(x, int), species_select),
                           True)
    blacklist_is_str = reduce(
            and_, map(lambda x: isinstance(x, str), species_blacklist), True)
    blacklist_is_int = reduce(
            and_, map(lambda x: isinstance(x, int), species_blacklist), True)
    if select_is_str:
        id_select = np.isin(frame.get_chemical_symbols(), species_select)
    elif select_is_int:
        id_select = np.isin(frame.get_atomic_numbers(), species_select)
    else:
        raise ValueError("Species select must be either all string or all int")
    if blacklist_is_str:
        id_blacklist = np.isin(frame.get_chemical_symbols(), species_blacklist)
    elif blacklist_is_int:
        id_blacklist = np.isin(frame.get_atomic_numbers(), species_blacklist)
    else:
        raise ValueError(
                "Species blacklist must be either all string or all int")
    if 'center_atoms_mask' not in frame.arrays:
        # add a default mask
        if species_select:
            old_mask = np.zeros((frame.get_number_of_atoms(),), dtype='bool')
        else:
            old_mask = np.ones((frame.get_number_of_atoms(),), dtype='bool')
    else:
        old_mask = frame.arrays['center_atoms_mask']
    # Python's "bitwise" operators do per-element logical operations in NumPy
    # see for instance https://docs.scipy.org/doc/numpy/reference/ufuncs.html#comparison-functions
    mask = (old_mask | id_select) & ~id_blacklist
    frame.arrays['center_atoms_mask'] = mask
