from collections.abc import Iterable
from functools import reduce
from operator import and_

import numpy as np

from ase.geometry import wrap_positions
from ase import Atoms

from ..lib import neighbour_list
from .base import (
    NeighbourListFactory,
    is_valid_structure,
    adapt_structure,
    StructureCollectionFactory,
)


class AtomsList(object):
    """
    A container for the neighbourlist and representation data associated with a list of atomic structures.

    This is a wrapper class for the `StructureManagerCollection` that have between precompiled on the C++ side.

    Attributes
    ----------
    nl_options : dict
        Parameters for each layer of the wrapped structure manager. Parameters
        can be specified for these layers: center, neighbourlist and strict.
    managers : StructureManagerCollection
        C++ object from rascal that holds the neighbourlist and the data associated with representations.

    """

    def __init__(self, frames, nl_options, start=None, length=None, managers=None):
        """Build a new AtomsList with only the selected atomic structures and
        corresponding neighborlist and representations (if present).

        Parameters
        -------
        frames :
            list of atomic structures.
        nl_options : dict
            Parameters for each layer of the wrapped structure manager. For
            example to initialize a neighbourlist for computing `SphericalInvariants` representation using a linked cell algorithm
        .. code:: python
            nl_options = [
            dict(name="centers", args=dict()),
            dict(name="neighbourlist", args=dict(cutoff=interaction_cutoff)),
            dict(name="centercontribution", args=dict()),
            dict(name="strict", args=dict(cutoff=interaction_cutoff)),
            ]

        managers : `StructureManagerCollection`
            Take directly a `StructureManagerCollection` without recomputing
            anything.
        """
        self.nl_options = nl_options
        self._frames = frames

        if managers is not None:
            self.managers = managers
        elif isinstance(frames, str):
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
            self.managers = managers
        else:
            # if python structure
            structures = convert_to_structure_list(frames)
            managers = StructureCollectionFactory(nl_options)
            try:
                managers.add_structures(structures)
            except Exception as e:
                raise RuntimeError(
                    "Neighbourlist of structures failed " + "because: " + str(e)
                )
            self.managers = managers

    def __iter__(self):
        return self.managers.__iter__()

    def __len__(self):
        return len(self.managers)

    def __getitem__(self, key):
        return self.managers[key]

    def get_subset(self, selected_ids):
        """Build a new AtomsList with only the selected atomic structures and
        corresponding neighborlist and representations (if present).

        Parameters
        -------
        selected_ids : list/array of indices

        Returns
        -------
        new_atom_list : AtomsList
        """
        selected_ids = list(map(int, selected_ids))
        new_managers = self.managers.get_subset(selected_ids)
        new_atom_list = AtomsList(
            [self._frames[idx] for idx in selected_ids],
            self.nl_options,
            managers=new_managers,
        )
        return new_atom_list

    def get_features(self, calculator, species=None):
        """
        Parameters
        -------
        calculator : Calculator (an object owning a _representation object)

        species :  list of atomic number to use for building the dense feature
        matrix computed with calculators of name Spherical*

        Returns
        -------
        represenation_matrix : ndarray
            returns the representation bound to the calculator as dense matrix.
        """

        if species is None:
            X = self.managers.get_features(calculator._representation)
        else:
            keys_list = calculator.get_keys(species)
            X = self.managers.get_features(calculator._representation, keys_list)

        return X

    def get_features_gradient(self, calculator, species=None):
        """
        Parameters
        -------
        calculator : Calculator (an object owning a _representation object)

        species :  list of atomic number to use for building the dense feature
        matrix computed with calculators of name Spherical*

        Returns
        -------
        dX_dr : ndarray of size (3*(n_neighbor+n_atom), n_features)
            returns the gradient of representation with respect to the atomic
            positions that have been computed with the calculator as a dense
            matrix. The method `get_gradients_info` provides the
            necessary information for operating on the dX_dr matrix.
        """

        if species is None:
            X = self.managers.get_features_gradient(calculator._representation, [])
        else:
            keys_list = calculator.get_keys(species)
            X = self.managers.get_features_gradient(
                calculator._representation, keys_list
            )

        return X

    def get_features_by_species(self, calculator):
        """
        Parameters
        -------
        calculator : one of the representation calculators named Spherical*

        Returns
        -------
        representation_matrix : dict of ndarray
            returns a dictionary associating tuples of atomic numbers sorted
            alphabetically to the corresponding feature matrices
        """
        return self.managers.get_features_by_species(calculator._representation)

    def get_gradients_info(self):
        """
        Returns
        -------
        ij : np.array of size (n_neighbor+n_atom, 5)
            Get informations necessary to the computation of gradients returned
            by `get_features_gradient`. It has as many rows as as the number
            gradients and each columns correspond to the index of the atomic
            structure, central atom, the neighbor atom and their atomic species.
        """
        return self.managers.get_gradients_info()

    def get_representation_info(self):
        """
        Returns
        -------
        ij : np.array of size (n_atoms, 3)
            Get informations necessary to the computation of predictions using
            the representation from `get_features`. It has as many rows as the
            number representations and they correspond to the index of the
            structure, the central atom and its atomic species.
        """
        return self.managers.get_representation_info()

    def get_direction_vectors(self):
        """
        Returns
        -------
        direction_vector : np.array of size (n_neighbor+n_atom, 3)
            Get the direction vectors from the atoms to their neighbors.
        """
        return self.managers.get_direction_vectors()

    def get_distances(self):
        """
        Returns
        -------
        distance : np.array of size (n_neighbor+n_atom, 4)
            Get the distances from the atoms to their neighbors.
        """
        return self.managers.get_distances()


def store_features_ase_atoms(atoms_list, feature_matrix, key="features"):
    """Store a feature matrix directly into a list of ASE Atoms

    This helps with compatibility and portability outside librascal

    Parameters
    ----------
    atoms_list: list(ase.Atoms)
        List of ASE Atoms to store the features into

    feature_matrix: np.ndarray, NxQ
        Array of (atom-wise) features calculated for the atoms list.
        First index is atoms, concatenated across structures.
        Second index is features.

    key : str
        Dictionary key under which to store the features.
        WARNING: If features are already present under this key,
        they will be overwritten!

    Returns
    -------
    Returns the modified atoms list, but note that this also modifies
    the atoms list directly!
    """
    offset = 0
    for structure in atoms_list:
        natoms = len(structure)
        # TODO maybe this should be done with Atoms.add_array() or similar
        structure.arrays[key] = feature_matrix[offset : offset + natoms]
        offset += natoms
    return atoms_list


def retrieve_features_ase_atoms(atoms_list, key="features"):
    """Retrieve a feature matrix stored in a list of ASE Atoms

    Format is expected as stored by store_features_ase_atoms()

    Parameters
    ----------
    atoms_list: list(ase.Atoms)
        List of ASE Atoms to store the features into

    key : str
        Dictionary key under which the features are stored.

    Returns
    -------
    feature_matrix: np.ndarray, NxQ
        Array of (atom-wise) features calculated for the atoms list.
        First index is atoms, concatenated across structures.
        Second index is features.

    Raises
    ------
    IndexError
        If the given key does not exist in one or more of the structures
        (the first structure with the missing key is indicated)

    """
    features_matrix = []
    for idx, structure in enumerate(atoms_list):
        try:
            features_matrix.append(structure.arrays[key])
        except IndexError as ide:
            raise IndexError(f"Key not found in structure number {idx:d}") from ide
    return np.concatenate(features_matrix)

def get_neighbourlist(structure, options):
    manager = NeighbourListFactory(options)
    manager.update(**structure)
    return manager


def convert_to_structure_list(frames):
    """Convert an atomic structure to format used internaly in rascal.
    Parameters
    ----------
    frames : ase.Atoms or list(ase.Atoms) or list(dict)
        atomic structure(s) in various formats.
        Note that the dictionary must be valid structure as per is_valid_structure
    Returns
    -------
    structure_list : AtomicStructureList

    """
    if isinstance(frames, Atoms):
        frames = [frames]
    elif not isinstance(frames, Iterable):
        raise ValueError("Must pass either an ase.Atoms object or an iterable")

    structure_list = neighbour_list.AtomicStructureList()
    for frame in frames:
        if is_valid_structure(frame):
            structure = frame
        else:
            if is_ase_Atoms(frame):
                structure = unpack_ase(frame)
            else:
                raise RuntimeError(
                    "Cannot convert structure of type {}".format(type(frame))
                )

        structure = sanitize_non_periodic_structure(structure)
        structure_list.append(**structure)
    return structure_list


def sanitize_non_periodic_structure(structure):
    """
    Rascal expects a unit cell that contains all the atoms even if the
    structure is not periodic.
    If the cell is set to 0 and the structure is not periodic then it
    is adapted to contain the atoms and the atoms are shifted inside the unit
    cell.

    Parameters
    ----------
    structure : a valid structure as per is_valid_structure


    Returns
    -------
    a valid structure as per is_valid_structure
        cell and positions have been modified if structure is not periodic
    """

    if np.all(structure["pbc"] == 0):
        cell = structure["cell"]
        if np.allclose(cell, np.zeros((3, 3))):
            pos = structure["positions"]
            bounds = np.array([pos.min(axis=1), pos.max(axis=1)])
            bounding_box_lengths = (bounds[1] - bounds[0]) * 1.05
            new_cell = np.diag(bounding_box_lengths)
            CoM = pos.mean(axis=1)
            disp = 0.5 * bounding_box_lengths - CoM
            new_pos = pos + disp[:, None]
            structure["positions"] = new_pos
    return structure


def is_ase_Atoms(frame):
    is_ase = True
    if not hasattr(frame, "get_cell"):
        is_ase = False
    if not hasattr(frame, "get_positions"):
        is_ase = False
    if not hasattr(frame, "get_atomic_numbers"):
        is_ase = False
    if not hasattr(frame, "get_pbc"):
        is_ase = False
    return is_ase


def unpack_ase(frame, wrap_pos=False):
    """
    Convert ASE Atoms object to rascal's equivalent

    Parameters
    ----------
    frame : ase.Atoms
        Atomic structure

    Returns
    -------
    dict containing the atomic structure informations readily readable by librascal.

    If the frame has an ase.atoms.arrays entry called
    'center_atoms_mask' then it will be used as the center mask
    (surprise) for any representations computed on this
    StructureManager.
    """
    cell = frame.get_cell()
    positions = frame.get_positions()
    numbers = frame.get_atomic_numbers()
    pbc = frame.get_pbc().astype(int)

    if wrap_pos:
        positions = wrap_positions(positions, cell, frame.get_pbc(), eps=1e-11)

    if "center_atoms_mask" in frame.arrays.keys():
        center_atoms_mask = frame.get_array("center_atoms_mask")
    else:
        center_atoms_mask = np.ones_like(numbers, dtype=bool)

    return adapt_structure(
        cell=cell,
        positions=positions,
        atom_types=numbers,
        pbc=pbc,
        center_atoms_mask=center_atoms_mask,
    )


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
    if "center_atoms_mask" not in frame.arrays:
        # add a default mask
        if id_select is not None:
            mask = np.zeros((len(frame),), dtype="bool")
        else:
            mask = np.ones((len(frame),), dtype="bool")
    else:
        mask = frame.arrays["center_atoms_mask"]
    if id_select is not None:
        mask[id_select] = True
    if id_blacklist is not None:
        mask[id_blacklist] = False
    frame.arrays["center_atoms_mask"] = mask


def mask_center_atoms_by_species(frame, species_select=[], species_blacklist=[]):
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
    select_is_str = reduce(
        and_, map(lambda x: isinstance(x, str), species_select), True
    )
    select_is_int = reduce(
        and_, map(lambda x: isinstance(x, int), species_select), True
    )
    blacklist_is_str = reduce(
        and_, map(lambda x: isinstance(x, str), species_blacklist), True
    )
    blacklist_is_int = reduce(
        and_, map(lambda x: isinstance(x, int), species_blacklist), True
    )
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
        raise ValueError("Species blacklist must be either all string or all int")
    if "center_atoms_mask" not in frame.arrays:
        # add a default mask
        if species_select:
            old_mask = np.zeros((len(frame),), dtype="bool")
        else:
            old_mask = np.ones((len(frame),), dtype="bool")
    else:
        old_mask = frame.arrays["center_atoms_mask"]
    # Python's "bitwise" operators do per-element logical operations in NumPy
    # see for instance https://docs.scipy.org/doc/numpy/reference/ufuncs.html#comparison-functions
    mask = (old_mask | id_select) & ~id_blacklist
    frame.arrays["center_atoms_mask"] = mask
