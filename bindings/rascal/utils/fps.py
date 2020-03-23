from ..lib import sparsification
import numpy as np
from ..models import PseudoPoints


def fps(feature_matrix, n_select, starting_index=None,
        method='simple', restart=None):
    """
    Farthest Point Sampling [1] routine using librascal.

    Parameters
    ----------
    feature_matrix : numpy.ndarray[float64[n, m], flags.c_contiguous]
        Feature matrix with n samples and m features.
    n_select : int
        Number of sample to select
    starting_index : int, optional
        Index of the first sample to select (the default is None,
                           which corresponds to starting_index == 0)
    restart : dictionary, optional (only valid for method="simple")
        the return value of a previous call to FPS. the selection
        will continue where it has been left off. the dictionary
        should contain at least "fps_indices"; if "fps_minmax_d2" and
        "fps_hausdorff_d2" are present, they will be used to restart,
        otherwise they'll be recomputed

    method : str, optional
        Select which kind of FPS selection to perform:
        + 'simple' is the basic fps selection [1].
        + 'voronoi' uses voronoi polyhedra to avoid some
                    distance computations.

    Returns
    -------
    A dictionary containing the following quantities
    "fps_indices" : numpy.ndarray[int[n_select]]
        Selected indices refering to the feature_matrix order.
    "fps_minmax_d2": numpy.ndarray[float64[n_select]]
        MIN-MAX distance^2 at each step in the selection
    "fps_hausforff_d2": numpy.ndarray[float64[n]]
        array of Hausdorff distances between the n points and the
        n_select FPS points
    In addition, when method="voronoi", the following arrays are returned
    "fps_voronoi_indices": numpy.ndarray[int[n]]
        the indices that assign each of the n inputs to the closest FPS point
    "fps_voronoi_r2": numpy.ndarray[float64[n_select]]
        the squared "Voronoi radius" of each FPS point (the largest distance to
        one of the points associated to its Voronoi cell


    .. [1] Ceriotti, M., Tribello, G. A., & Parrinello, M. (2013).
        Demonstrating the Transferability and the Descriptive Power of Sketch-Map.
        Journal of Chemical Theory and Computation, 9(3), 1521–1532.
        https://doi.org/10.1021/ct3010563
    """

    if starting_index is None:
        starting_index = 0

    return_dict = {}
    if method == 'simple':
        if restart is None:
            sparse_indices, sparse_minmax_d2, lmin_d2 = \
                sparsification.fps(feature_matrix, n_select, starting_index)
        else:
            res_tuple = (restart["fps_indices"],
                         (restart["fps_minmax_d2"] if "fps_minmax_d2" in restart
                          else np.zeros(0, float)),
                         (restart["fps_hausdorff_d2"] if "fps_hausdorff_d2" in restart
                          else np.zeros(0, float)))
            sparse_indices, sparse_minmax_d2, lmin_d2 = \
                sparsification.fps(feature_matrix, n_select,
                                   starting_index, res_tuple)

    elif method == 'voronoi':
        sparse_indices, sparse_minmax_d2, lmin_d2, \
            voronoi_indices, voronoi_r2 = \
            sparsification.fps_voronoi(feature_matrix,
                                       n_select, starting_index)
        return_dict["fps_voronoi_indices"] = voronoi_indices
        return_dict["fps_voronoi_r2"] = voronoi_r2

    else:
        raise Exception('Unknown FPS algorithm {}'.format(method))

    return_dict["fps_indices"] = sparse_indices
    return_dict["fps_minmax_d2"] = sparse_minmax_d2
    return_dict["fps_hausdorff_d2"] = lmin_d2

    return return_dict


class FPSFilter(object):
    """FPS selection of samples in a given feature matrix.

    Parameters
    ----------
    representation : Calculator
        Representation calculator associated with the kernel

    Nselect: int
        number of points to select. 

    act_on: string
        Select how to apply the selection. Can be either of 'sample',
        'sample per species','feature'.
        For the moment only 'sample_per_species' is implemented.

    first: int
        index of first sample to select

    """

    def __init__(self, representation, Nselect, act_on='sample per specie', first=0):
        super(FPSFilter, self).__init__()
        
        self._representation = representation
        self.Nselect = Nselect
        if act_on in ['sample', 'sample per specie', 'feature']:
            self.act_on = act_on
        else:
            raise 'Wrong input: {}'.format(act_on)
        self.first = first
        self.selected_ids = None

    def fit_transform(self, managers):
        """Perform CUR selection of samples/features.

        Parameters
        ----------
        managers : AtomsList
            list of structures containing features computed with representation

        Returns
        -------
        PseudoPoints
            Selected samples

        Raises
        ------
        ValueError
            [description]
        NotImplementedError
            [description]
        """

        if self.act_on in ['sample per specie']:
            # get the dense feature matrix
            X = managers.get_features(self._representation)

            sps = list(self.Nselect.keys())

            # get various info from the structures about the center atom species and indexing
            (strides_by_sp, global_counter, map_by_manager,
             indices_by_sp) = self.get_index_mappings_sample_per_species(managers)

            print('The number of pseudo points selected by central atom species is: {}'.format(
                self.Nselect))

            # organize features w.r.t. central atom type
            X_by_sp = {}
            for sp in sps:
                X_by_sp[sp] = X[indices_by_sp[sp]]
            self._XX = X_by_sp

            # split the dense feature matrix by center species and apply CUR decomposition
            selected_ids_by_sp = {}
            for sp in sps:
                print('Selecting species: {}'.format(sp))
                fps_sel = fps(X_by_sp[sp], self.Nselect[sp], self.first)
                selected_ids_by_sp[sp] = np.sort(fps_sel["fps_indices"])

            self.selected_ids = self.convert_selected_global_index2rascal_sample_per_species(
                managers, selected_ids_by_sp, strides_by_sp, map_by_manager)

            # build the pseudo points
            pseudo_points = PseudoPoints(self._representation)
            pseudo_points.extend(managers, self.selected_ids)

            return pseudo_points
        else:
            raise NotImplementedError("method: {}".format(self.act_on))

    def get_index_mappings_sample_per_species(self, managers):
        # get various info from the structures about the center atom species and indexing
        sps = list(self.Nselect.keys())
        types = []
        strides_by_sp = {sp: [0] for sp in sps}
        global_counter = {sp: 0 for sp in sps}
        indices_by_sp = {sp: [] for sp in sps}
        map_by_manager = [{} for ii in range(len(managers))]
        for i_man, man in enumerate(managers):
            counter = {sp: 0 for sp in sps}
            for i_at, at in enumerate(man):
                types.append(at.atom_type)
                if at.atom_type in sps:
                    map_by_manager[i_man][global_counter[at.atom_type]] = i_at
                    counter[at.atom_type] += 1
                    global_counter[at.atom_type] += 1
                else:
                    raise ValueError('Atom type {} has not been specified in fselect: {}'.format(
                        at.atom_type, self.Nselect))
            for sp in sps:
                strides_by_sp[sp].append(counter[sp])

        for sp in sps:
            strides_by_sp[sp] = np.cumsum(strides_by_sp[sp])

        for ii, sp in enumerate(types):
            indices_by_sp[sp].append(ii)

        return strides_by_sp, global_counter, map_by_manager, indices_by_sp

    def convert_selected_global_index2rascal_sample_per_species(self, managers, selected_ids_by_sp, strides_by_sp, map_by_manager):
        # convert selected center indexing into the rascal format
        selected_ids = [[] for ii in range(len(managers))]
        sps = list(self.Nselect.keys())
        i_manager = {sp: 0 for sp in sps}
        for sp in sps:
            for idx in selected_ids_by_sp[sp]:
                carry_on = True
                while carry_on:
                    if idx >= strides_by_sp[sp][i_manager[sp]] and idx < strides_by_sp[sp][i_manager[sp] + 1]:
                        selected_ids[i_manager[sp]].append(
                            map_by_manager[i_manager[sp]][idx])
                        carry_on = False
                    else:
                        i_manager[sp] += 1
        for ii in range(len(selected_ids)):
            selected_ids[ii] = list(np.sort(selected_ids[ii]))
        return selected_ids
