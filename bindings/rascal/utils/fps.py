from .io import BaseIO
from ..lib import sparsification

from .filter_utils import (get_index_mappings_sample_per_species,
                           convert_selected_global_index2perstructure_index_per_species,
                           get_index_mappings_sample,
                           convert_selected_global_index2perstructure_index)
import numpy as np


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


class FPSFilter(BaseIO):
    """Farthest Point Sampling (FPS) to select samples or features in a given feature matrix.
    Wrapper around the fps function for convenience.
    Parameters
    ----------
    representation : Calculator
        Representation calculator associated with the kernel
    Nselect: int
        number of points to select. if act_on='sample per specie' then it should
        be a dictionary mapping atom type to the number of samples, e.g.
        Nselect = {1:200,6:100,8:50}.
    act_on: string
        Select how to apply the selection. Can be either of 'sample',
        'sample per species','feature'.
    starting_index: int
        Index used to start the FPS selection with.

    """

    def __init__(self, representation, Nselect, act_on='sample per species',
                 starting_index=0):
        self._representation = representation
        self.Nselect = Nselect
        self.starting_index = starting_index
        modes = ['sample', 'sample per species', 'feature']
        if act_on in modes:
            self.act_on = act_on
        else:
            raise ValueError(
                '"act_on" should be either of: "{}", "{}", "{}"'.format(*modes))

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

        self.fps_minmax_d2_by_sp = None
        self.fps_minmax_d2 = None

    def select(self, managers):
        """Perform FPS selection of samples/features.
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

        # get the dense feature matrix
        X = managers.get_features(self._representation)

        if self.act_on == 'sample per species':
            sps = list(self.Nselect.keys())

            # get various info from the structures about the center atom species and indexing
            (strides_by_sp, global_counter, map_by_manager,
             indices_by_sp) = get_index_mappings_sample_per_species(managers, sps)

            print('The number of pseudo points selected by central atom species is: {}'.format(
                self.Nselect))

            # organize features w.r.t. central atom type
            X_by_sp = {}
            for sp in sps:
                X_by_sp[sp] = X[indices_by_sp[sp]]
            self._XX = X_by_sp

            # split the dense feature matrix by center species and apply CUR decomposition
            self.selected_sample_ids_by_sp = {}
            self.fps_minmax_d2_by_sp = {}
            self.fps_hausforff_d2_by_sp = {}
            for sp in sps:
                print('Selecting species: {}'.format(sp))
                fps_out = fps(X_by_sp[sp], self.Nselect[sp],
                              starting_index=self.starting_index)
                self.selected_sample_ids_by_sp[sp] = fps_out['fps_indices']
                self.fps_minmax_d2_by_sp[sp] = fps_out['fps_minmax_d2']

            return self
        elif self.act_on == 'feature':
            fps_out = fps(X.T, self.Nselect, starting_index=self.starting_index)
            self.selected_feature_ids_global = fps_out['fps_indices']
            self.fps_minmax_d2 = fps_out['fps_minmax_d2']
        elif self.act_on == 'sample':
            fps_out = fps(X, self.Nselect, starting_index=self.starting_index)
            self.selected_sample_ids = fps_out['fps_indices']
            self.fps_minmax_d2 = fps_out['fps_minmax_d2']
        else:
            raise NotImplementedError("method: {}".format(self.act_on))

    def filter(self, managers, n_select=None):
        from ..models import SparsePoints
        if n_select is None:
            n_select = self.Nselect

        if self.act_on == 'sample per species':
            sps = list(n_select.keys())
            # get various info from the structures about the center atom species and indexing
            (strides_by_sp, global_counter, map_by_manager,
             indices_by_sp) = get_index_mappings_sample_per_species(managers, sps)
            selected_ids_by_sp = {key: np.sort(val[:n_select[key]])
                                  for key, val in self.selected_sample_ids_by_sp.items()}
            self.selected_ids = convert_selected_global_index2perstructure_index_per_species(
                managers, selected_ids_by_sp, strides_by_sp, map_by_manager, sps)
            # return self.selected_ids
            # build the pseudo points
            pseudo_points = SparsePoints(self._representation)
            pseudo_points.extend(managers, self.selected_ids)
            return pseudo_points

        elif self.act_on == 'sample':
            selected_ids_global = np.sort(self.selected_sample_ids[:n_select])
            strides, _, map_by_manager = get_index_mappings_sample(managers)
            self.selected_ids = convert_selected_global_index2perstructure_index(managers,
                                                                                 selected_ids_global, strides, map_by_manager)
            return self.selected_ids
            # SparsePoints is not compatible with a non center atom species
            # dependant sparse points
            # # build the pseudo points
            # pseudo_points = SparsePoints(self._representation)
            # pseudo_points.extend(managers, self.selected_ids)
            # return pseudo_points

        elif self.act_on == 'feature':
            feat_idx2coeff_idx = self._representation.get_feature_index_mapping(
                managers)
            self.selected_ids = {key: []
                                 for key in feat_idx2coeff_idx[0].keys()}
            selected_ids_sorting = np.argsort(
                self.selected_feature_ids_global[:n_select])
            selected_feature_ids = self.selected_feature_ids_global[selected_ids_sorting]
            for idx in selected_feature_ids:
                coef_idx = feat_idx2coeff_idx[idx]
                for key in self.selected_ids.keys():
                    self.selected_ids[key].append(int(coef_idx[key]))
            # keep the global indices and ordering for ease of use
            self.selected_ids['selected_features_global_ids'] = selected_feature_ids.tolist(
            )
            self.selected_ids['selected_features_global_ids_fps_ordering'] = selected_ids_sorting.tolist(
            )
            self.selected_ids = dict(coefficient_subselection=self.selected_ids)
            return self.selected_ids

    def select_and_filter(self, managers):
        return self.select(managers).filter(managers)

    def plot_fps_error(self):
        import matplotlib.pyplot as plt
        if self.fps_minmax_d2_by_sp is None:
            plt.semilogy(self.fps_minmax_d2, label=self.act_on)

        else:
            for sp in self.fps_minmax_d2_by_sp:
                plt.semilogy(self.fps_minmax_d2_by_sp[sp],
                             label='{} species {}'.format(self.act_on, sp))
            plt.legend()
        plt.title('FPSFilter')
        plt.ylabel('fps minmax d^2')

    def _get_data(self):
        return dict(selected_ids=self.selected_ids,
                    selected_sample_ids=self.selected_sample_ids,
                    selected_sample_ids_by_sp=self.selected_sample_ids_by_sp,
                    selected_feature_ids_global=self.selected_feature_ids_global)

    def _set_data(self, data):
        self.selected_ids = data['selected_ids']
        self.selected_sample_ids = data['selected_sample_ids']
        self.selected_sample_ids_by_sp = data['selected_sample_ids_by_sp']
        self.selected_feature_ids_global = data['selected_feature_ids_global']

    def _get_init_params(self):
        return dict(representation=self._representation,
                    Nselect=self.Nselect,
                    act_on=self.act_on,
                    starting_index=self.starting_index,)
