import logging
from ase import Atoms
import numpy as np
from .io import BaseIO
from ..models.sparse_points import SparsePoints
from ..representations.spherical_invariants import SphericalInvariants
from ..neighbourlist import retrieve_features_ase_atoms

LOGGER = logging.getLogger(__name__)

try:
    from skcosmo._selection import _FPS, _CUR
except ImportError as ie:
    LOGGER.warn(
        "Warning: skcosmo module not found. CUR and FPS filters will be unavailable."
    )
    LOGGER.warn("Original error:\n" + str(ie))
    _FPS = _CUR = None

# Index conversion utilities


def _indices_manager_to_perstructure(managers, selected_ids_global):
    """Convert manager-global center indexing to per-structure format

    That is, change them to a list of lists of structure-local atom indices
    that is accepted as input to SparsePoints

    Parameters
    ----------
    managers : AtomsList or list(ase.Atoms)
        list of atomic structures
    selected_ids_global : list of int
        global indices to convert

    Returns
    -------
    selected_ids : list of list(int)
        list the atom indices (within their structure) that have been selected
    """
    selected_ids = []
    selected_ids_global = np.array(selected_ids_global)
    natoms_list = [len(manager) for manager in managers]
    split_idces = np.cumsum(natoms_list)
    structure_start_idx = 0
    # Handle any out-of-range indices
    if np.any(selected_ids_global >= split_idces[-1]):
        bad_indices = selected_ids_global[selected_ids_global >= split_idces[-1]]
        bad_indices_str = np.array2string(
            bad_indices, threshold=5, edgeitems=2, separator=", "
        )
        raise ValueError(f"Selected index(es): {bad_indices_str} out of range")
    # Do the actual conversion if OK
    for structure_end_idx in split_idces:
        this_structure_idces = selected_ids_global[
            (selected_ids_global >= structure_start_idx)
            & (selected_ids_global < structure_end_idx)
        ]
        selected_structure_idces = this_structure_idces - structure_start_idx
        selected_ids.append(list(selected_structure_idces))
        structure_start_idx = structure_end_idx
    assert sum(len(ids) for ids in selected_ids) == len(selected_ids_global)
    return selected_ids


def _indices_perspecies_manager_to_perstructure(managers, selected_ids_by_sp, sps):
    """Convert per-species, manager-global center indexing to per-structure

    That is, change them to a list of lists of structure-local atom indices
    that is accepted as input to SparsePoints

    See _get_index_mappings_sample_per_species() to make the intput

    Parameters
    ----------
    managers : AtomsList
        list of atomic structures
    selected_ids_by_sp : dict
        indices to convert; each dictionary entry (keyed by species)
        is a list of indices into the array of all atoms in the manager
        of only that species
    sps : list(int) or set(int)
        unique center atom species present in managers

    Returns
    -------
    selected_ids : list of lists
        list the atom indices (within their structure) that have been selected
        They are ordered, first by species (sorted, irrespective of the order
        passed in), then by selection order.

    Notes
    -----
    Asking for indices for a species not present in the AtomsList will
    result in an out-of-range error (since the slice of atoms of that species
    is of size zero).

    This function does not yet support lists of ASE Atoms (instead of a list
    of managers) as input, but it should be fairly easy to add support in
    the future if required.
    """
    if len(set(sps)) != len(sps):
        raise ValueError(f"List of species contains duplicated entries: {sps}")
    selected_ids = [[] for ii in range(len(managers))]
    structure_sp_start_idx = {sp: 0 for sp in sps}
    for sp in sps:
        selected_ids_by_sp[sp] = np.array(selected_ids_by_sp[sp])
    for structure, structure_selected_ids in zip(managers, selected_ids):
        perspecies_counter = {sp: 0 for sp in sps}
        structure_index_mapping = {sp: [] for sp in sps}
        for perstructure_idx, atom in enumerate(structure):
            atom_sp = atom.atom_type
            if atom_sp not in sps:
                raise ValueError(
                    f"Atom of type {atom_sp} found but was not listed in sps: {sps}"
                )
            structure_index_mapping[atom_sp].append(perstructure_idx)
            perspecies_counter[atom_sp] += 1
        for sp in sorted(list(sps)):
            selected_ids_sp = selected_ids_by_sp[sp]
            structure_sp_end_idx = structure_sp_start_idx[sp] + perspecies_counter[sp]
            this_structure_sp_idces = (
                selected_ids_sp[
                    (selected_ids_sp >= structure_sp_start_idx[sp])
                    & (selected_ids_sp < structure_sp_end_idx)
                ]
                - structure_sp_start_idx[sp]
            )
            structure_selected_ids.extend(
                [
                    structure_index_mapping[sp][sp_idx]
                    for sp_idx in this_structure_sp_idces
                ]
            )
            structure_sp_start_idx[sp] += perspecies_counter[sp]
    for sp in sps:
        selected_out_of_range = (
            np.array(selected_ids_by_sp[sp]) >= structure_sp_start_idx[sp]
        )
        if np.any(selected_out_of_range):
            bad_indices = np.array(selected_ids_by_sp[sp])[selected_out_of_range]
            bad_indices_str = np.array2string(
                bad_indices, threshold=5, edgeitems=2, separator=", "
            )
            error_str = (
                f"Selected index(es): {bad_indices_str} for species {sp} out of range"
            )
            if 0 in bad_indices:
                error_str += " (species does not appear to be present)"
            raise ValueError(error_str)
    # Check that we haven't missed anything
    assert sum(len(ids) for ids in selected_ids) == sum(
        len(selected_ids_by_sp[sp]) for sp in sps
    )
    return selected_ids


def _split_feature_matrix_by_species(managers, X, sps):
    """Does exactly what it says on the tin

    Parameters
    ----------
    managers : AtomsList or list(ase.Atoms)
        list of atomic structures
    X : np.ndarray (2-D)
        feature matrix computed from managers; rows must correspond to atoms
    sps : list(int)
        list of unique center atom species present in managers

    Returns
    -------
    dict(np.ndarray)
        The feature matrix split into matrices each corresponding to one
        of the atomic species requested

    Warnings
    --------
    This function does not check that the list of species provided actually
    corresponds to those present in managers; it only performs the selection
    (which would be empty for nonexistent species).
    """
    X_per_species = {}
    global_species_list = []
    for structure in managers:
        if isinstance(structure, Atoms):
            global_species_list.extend(structure.get_atomic_numbers())
        else:
            # Assuming rascal StructureManagerCollection
            global_species_list.extend([atom.atom_type for atom in structure])
    global_species_list = np.array(global_species_list)
    for sp in sps:
        X_per_species[sp] = X[global_species_list == sp]
    return X_per_species


class AtomsListWrapper:

    """Thin wrapper around a list of ASE Atoms objects to make it
    compatible with the Filter utilities
    """

    def __init__(self, atoms_list):
        self.atoms_list = atoms_list

    def __iter__(self):
        return iter(self.atoms_list)

    def get_features(self, representation_key):
        """Get the features stored in the list of atoms

        See librascal.neighbourlist.store_features_ase_atoms() for
        information on how to store the features in the first place

        Instead of a representation manager, takes a simple string key
        that is used to retrieve the features from the ASE Atoms'
        arrays dictionaries
        """
        return retrieve_features_ase_atoms(self.atoms_list, representation_key)


class Filter(BaseIO):
    """
    A super class for filtering representations based upon a standard
    sample or feature selection class.

    This is mainly a wrapper around selectors (implemented e.g. in
    scikit-cosmo) that handles the semantic-index transformations
    required after selection.

    Parameters
    ----------

    representation : Calculator
        Representation calculator associated with the kernel

    Nselect: int
        number of points to select. If act_on='sample per species' then it should
        be a dictionary mapping atom type to the number of samples, e.g.
        Nselect = {1:200,6:100,8:50}.

    selector: selector to use for filtering. The selector should
            have a `fit` function, which when called will select from the input
            matrix the desired features / samples and a `get_support` function
            which takes parameters `indices` and `ordered`, and returns a list
            of selection indices, in the order that they were selected,
            when `indices=True` and `ordered=True`.

    act_on: string
        Select how to apply the selection. Can be either of 'sample',
        'sample per species','feature'.  Default 'sample per species'.
        Note that for 'feature' mode only the SphericalInvariants
        representation is supported.

    """

    def __init__(
        self,
        representation,
        Nselect,
        selector,
        act_on="sample per species",
    ):
        self._representation = representation
        self.Nselect = Nselect
        if self.act_on is None:
            self._check_set_mode(act_on)
        # effectively selected list of indices at the filter step
        # the indices have been reordered for effiency and compatibility with
        # the c++ routines
        self.selected_ids = None
        # for 'sample' selection
        self.selected_sample_ids = None
        # for 'sample per species' selection
        self.selected_sample_ids_by_sp = None
        # for feature selection
        self.selected_feature_ids_global = None
        self._selector = selector

    def _check_set_mode(
        self, act_on, modes=["sample", "sample per species", "feature"]
    ):
        """Check that the supplied act_on is one of the supported modes

        Set the mode if it is valid, aise a ValueError with a helpful
        message otherwise

        A list of valid modes can be supplied in case it differs from the
        superclass default.
        """
        if act_on in modes:
            self.act_on = act_on
        else:
            valid_modes = ['"{}"'.format(mode) for mode in modes]
            if len(valid_modes) > 1:
                valid_modes[-1] = "or " + valid_modes[-1]
            valid_modes_str = ", ".join(valid_modes)
            raise ValueError('"act_on" should be one of: ' + valid_modes_str)

    def select(self, managers):
        """Perform selection of samples/features.

        Parameters
        ----------
        managers : AtomsList
            list of structures containing features computed with representation

        Returns
        -------
        Filter (self)
            Returns self; use `filter()` to perform the actual filtering
            operation

        """
        X = managers.get_features(self._representation)
        if self.act_on == "sample per species":
            self.selected_sample_ids_by_sp = {}
            sps = list(self.Nselect.keys())
            LOGGER.info(
                f"The number of pseudo points selected by central atom species is: {self.Nselect}"
            )
            X_by_sp = _split_feature_matrix_by_species(managers, X, sps)
            for sp in sps:
                LOGGER.info(f"Selecting species: {sp}")
                if self._selector[sp] is not None:
                    self._selector[sp].fit(X_by_sp[sp])

                    in_sample_indices = np.array(
                        self._selector[sp].get_support(indices=True, ordered=True),
                        dtype=int,
                    )
                    self.selected_sample_ids_by_sp[sp] = in_sample_indices
                else:
                    self.selected_sample_ids_by_sp[sp] = []
            return self
        else:
            self._selector.fit(X)
            if self.act_on == "sample":
                self.selected_sample_ids = self._selector.get_support(
                    indices=True, ordered=True
                )
            else:
                self.selected_feature_ids_global = self._selector.get_support(
                    indices=True, ordered=True
                )
            return self

    def filter(self, managers, n_select=None):
        """Apply the fitted selection to a new set of managers

        Parameters
        ----------
        managers : AtomsList
            list of structures containing features computed with representation

        n_select : int
            number of selections to return, must be less than self.Nselect

        Returns
        -------
        SparsePoints or 2-d array or dict
            Selected samples.  The format depends on self.act_on - if it is
            "sample" or "sample per species", then a SparsePoints instance
            is directly returned (in the case that 'managers' is a rascal
            StructureManagerCollection -- if it is a AtomsListWrapper, the
            sparse feature matrix itself is returned).
            If it is "feature", then it is a dictionary
            containing the "coefficient_subselection" key that can be directly
            passed to the SphericalInvariants constructor.

        Warnings
        --------
        Note that the selected points are sorted in order of selection,
        _except_ if self.act_on=="sample", in which case the sparse points
        are afterwards sorted by species.

        Raises
        ------
        ValueError
            if requesting more selected samples or features than were used
            to initialize the representation

        """
        if n_select is None:
            n_select = self.Nselect
        else:
            if n_select > self.Nselect:
                raise ValueError(
                    f"It is only possible to filter {self.Nselect} {self.act_on}(s), "
                    f"you have requested {n_select}"
                )
        if self.act_on == "sample per species":
            sps = list(n_select.keys())
            selected_ids_by_sp = {
                key: val[: n_select[key]]
                for key, val in self.selected_sample_ids_by_sp.items()
            }
            if isinstance(managers, AtomsListWrapper):
                X = managers.get_features(self._representation)
                X_by_sp = _split_feature_matrix_by_species(managers, X, sps)
                # This splits the sparse points by species, which is what the
                # internal rascal version does anyway
                sparse_points = dict()
                for sp in sps:
                    sparse_points[sp] = X_by_sp[sp][selected_ids_by_sp[sp]]
            else:
                # Assuming rascal StructureManagerCollection
                self.selected_ids = _indices_perspecies_manager_to_perstructure(
                    managers, selected_ids_by_sp, sps
                )
                sparse_points = SparsePoints(self._representation)
                sparse_points.extend(managers, self.selected_ids)
            return sparse_points
        elif self.act_on == "sample":
            selected_ids_global = self.selected_sample_ids[:n_select]
            if isinstance(managers, AtomsListWrapper):
                X = managers.get_features(self._representation)
                sparse_points = X[selected_ids_global]
            else:
                # Assuming rascal StructureManagerCollection
                self.selected_ids = _indices_manager_to_perstructure(
                    managers, selected_ids_global
                )
                # The sparse points will be reordered since they're not per-species
                # but the resulting object is still usable
                sparse_points = SparsePoints(self._representation)
                sparse_points.extend(managers, self.selected_ids)
            return sparse_points
        elif self.act_on == "feature":
            if not isinstance(self._representation, SphericalInvariants):
                raise ValueError(
                    "Feature filtering currently only supported for SphericalInvariants"
                )
            feat_idx2coeff_idx = self._representation.get_feature_index_mapping(
                managers
            )
            self.selected_ids = {key: [] for key in feat_idx2coeff_idx[0].keys()}
            selected_ids_sorting = np.argsort(
                self.selected_feature_ids_global[:n_select]
            )
            selected_feature_ids = self.selected_feature_ids_global[
                selected_ids_sorting
            ]
            for idx in selected_feature_ids:
                coef_idx = feat_idx2coeff_idx[idx]
                for key in self.selected_ids.keys():
                    self.selected_ids[key].append(int(coef_idx[key]))
            self.selected_ids = dict(coefficient_subselection=self.selected_ids)
            # keep the global indices and ordering for ease of use
            self.selected_ids[
                "selected_feature_ids_global"
            ] = selected_feature_ids.tolist()
            self.selected_ids[
                "selected_feature_ids_global_selection_ordering"
            ] = selected_ids_sorting.tolist()
            return self.selected_ids

    def select_and_filter(self, managers):
        return self.select(managers).filter(managers)

    def _get_data(self):
        return dict(
            selected_ids=self.selected_ids,
            selected_sample_ids=self.selected_sample_ids,
            selected_sample_ids_by_sp=self.selected_sample_ids_by_sp,
            selected_feature_ids_global=self.selected_feature_ids_global,
        )

    def _set_data(self, data):
        self.selected_ids = data["selected_ids"]
        self.selected_sample_ids = data["selected_sample_ids"]
        self.selected_sample_ids_by_sp = data["selected_sample_ids_by_sp"]
        self.selected_feature_ids_global = data["selected_feature_ids_global"]

    def _get_init_params(self):
        return dict(
            representation=self._representation,
            Nselect=self.Nselect,
            act_on=self.act_on,
        )


class CURFilter(Filter):
    def __init__(
        self,
        representation,
        Nselect,
        act_on="sample per species",
        selector_args={},
        **kwargs,
    ):
        modes = ["sample", "sample per species", "feature"]
        self._check_set_mode(act_on, modes)
        if act_on == "sample":
            selector = _CUR(
                selection_type="sample", n_to_select=Nselect, **selector_args
            )
        elif act_on == "sample per species":
            selector = {
                n: _CUR(
                    selection_type="sample", n_to_select=Nselect[n], **selector_args
                )
                if Nselect[n] > 0
                else None
                for n in Nselect
            }

        else:
            assert act_on == "feature"
            selector = _CUR(
                selection_type="feature", n_to_select=Nselect, **selector_args
            )

        super().__init__(
            representation=representation,
            Nselect=Nselect,
            act_on=act_on,
            selector=selector,
            **kwargs,
        )


class FPSFilter(Filter):
    def __init__(
        self,
        representation,
        Nselect,
        act_on="sample per species",
        selector_args={},
        **kwargs,
    ):
        modes = ["sample", "sample per species", "feature"]
        self._check_set_mode(act_on, modes)
        if act_on == "sample":
            selector = _FPS(
                selection_type="sample", n_to_select=Nselect, **selector_args
            )
        elif act_on == "sample per species":
            selector = {
                n: _FPS(
                    selection_type="sample", n_to_select=Nselect[n], **selector_args
                )
                if Nselect[n] > 0
                else None
                for n in Nselect
            }

        else:
            assert act_on == "feature"
            selector = _FPS(
                selection_type="feature", n_to_select=Nselect, **selector_args
            )

        super().__init__(
            representation=representation,
            Nselect=Nselect,
            act_on=act_on,
            selector=selector,
            **kwargs,
        )

    def get_fps_distances(self):
        """Return the Hausdorff distances over the course of selection

        This may be a useful (rough) indicator for choosing how many points to
        select, as a small distance generally indicates that the selected point
        is close to the existing set of selected points and therefore probably
        does not add much additional information.

        Returns either an array of Hausdorff distances, or a species-indexed
        dict of arrays (for the "sample per species" mode).
        """
        if self.act_on == "sample per species":
            return {
                sp: self._selector[sp].get_select_distance() for sp in self._selector
            }
        else:
            return self._selector.get_select_distance()
