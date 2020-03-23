from ..lib import sparsification
# from ..models import PseudoPoints
from .filter_utils import (get_index_mappings_sample_per_species,
convert_selected_global_index2rascal_sample_per_species,
get_index_mappings_sample,
convert_selected_global_index2rascal_sample)
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


class FPSFilter(object):
    """Farther Point Sampling (FPS) to select samples or features in a given feature matrix.
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

    def __init__(self, representation, Nselect, act_on='sample per specie', starting_index=0):
        super(FPSFilter, self).__init__()
        self._representation = representation
        self.Nselect = Nselect
        self.starting_index = starting_index
        if act_on in ['sample', 'sample per specie', 'feature']:
            self.act_on = act_on
        else:
            raise 'Wrong input: {}'.format(act_on)

        self.selected_ids = None
        self.fps_minmax_d2_by_sp = None
        self.fps_minmax_d2 = None

    def fit(self, managers):
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

        # get the dense feature matrix
        X = managers.get_features(self._representation)

        if self.act_on in ['sample per specie']:
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
            self.selected_ids_by_sp = {}
            self.fps_minmax_d2_by_sp = {}
            self.fps_hausforff_d2_by_sp = {}
            for sp in sps:
                print('Selecting species: {}'.format(sp))
                fps_out = fps(X_by_sp[sp], self.Nselect[sp], starting_index=self.starting_index)
                self.selected_ids_by_sp[sp] = fps_out['fps_indices']
                self.fps_minmax_d2_by_sp[sp] = fps_out['fps_minmax_d2']

            return self
        elif self.act_on in ['feature']:
            fps_out = fps(X.T, self.Nselect, starting_index=self.starting_index)
            self.selected_ids = fps_out['fps_indices']
            self.fps_minmax_d2 = fps_out['fps_minmax_d2']
        elif self.act_on in ['sample']:
            fps_out = fps(X, self.Nselect, starting_index=self.starting_index)
            self.selected_ids_global = fps_out['fps_indices']
            self.fps_minmax_d2 = fps_out['fps_minmax_d2']
        else:
            raise NotImplementedError("method: {}".format(self.act_on))

    def transform(self, managers):
        if self.act_on in ['sample per specie']:
            sps = list(self.Nselect.keys())
            # get various info from the structures about the center atom species and indexing
            (strides_by_sp, global_counter, map_by_manager,
             indices_by_sp) = get_index_mappings_sample_per_species(managers, sps)
            selected_ids_by_sp = {key:val[:self.Nselect[key]] for key,val in self.selected_ids_by_sp.items()}
            self.selected_ids = convert_selected_global_index2rascal_sample_per_species(
                managers, selected_ids_by_sp, strides_by_sp, map_by_manager, sps)
            return self.selected_ids
            # # build the pseudo points
            # pseudo_points = PseudoPoints(self._representation)
            # pseudo_points.extend(managers, self.selected_ids)
            # return pseudo_points
        elif self.act_on in ['sample']:
            selected_ids_global = self.selected_ids_global[:self.Nselect]
            strides, _, map_by_manager = get_index_mappings_sample(managers)
            self.selected_ids = convert_selected_global_index2rascal_sample(managers,
                                                        selected_ids_global, strides, map_by_manager)
            return self.selected_ids
            # # build the pseudo points
            # pseudo_points = PseudoPoints(self._representation)
            # pseudo_points.extend(managers, self.selected_ids)
            # return pseudo_points

        elif self.act_on in ['feature']:
            feat_idx2coeff_idx = self._representation.get_feature_index_mapping(managers)
            selected_features = {key:[] for key in feat_idx2coeff_idx[0].keys()}
            for idx in self.selected_ids[:self.Nselect]:
                coef_idx = feat_idx2coeff_idx[idx]
                for key in selected_features.keys():
                    selected_features[key].append(int(coef_idx[key]))
            # keep the global indices for ease of use
            selected_features['selected_features_global_ids'] = self.selected_ids[:self.Nselect].tolist()
            return dict(coefficient_subselection=selected_features)

    def plot(self):
        import matplotlib.pyplot as plt
        if self.fps_minmax_d2_by_sp is None:
            plt.semilogy(self.fps_minmax_d2,label=self.act_on)

        else:
            for sp in self.fps_minmax_d2_by_sp:
                plt.semilogy(self.fps_minmax_d2_by_sp[sp],
                            label='{} species {}'.format(self.act_on, sp))
            plt.legend()
        plt.title('FPSFilter')
        plt.ylabel('fps minmax d^2')

    def fit_transform(self, managers):
        return self.fit(managers).transform(managers)
