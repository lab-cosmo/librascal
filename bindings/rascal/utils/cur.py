from ..utils import BaseIO
from ..utils.filter_utils import (get_index_mappings_sample_per_species,
                                  convert_selected_global_index2perstructure_index_per_species,
                                  get_index_mappings_sample,
                                  convert_selected_global_index2perstructure_index)

from ..models.sparse_points import SparsePoints
import numpy as np
from scipy.sparse.linalg import svds


def do_CUR(X, Nsel, act_on='sample', is_deterministic=False, seed=10, verbose=True):
    """ Apply CUR selection [1] of Nsel rows or columns of the
    given feature matrix X[n_samples, n_features].

    .. [1] Mahoney, M. W., & Drineas, P. (2009). CUR matrix decompositions for
    improved data analysis. Proceedings of the National Academy of Sciences,106(3),
    697–702. https://doi.org/10.1073/pnas.0803205106
    """
    U, _, VT = svds(X, Nsel)
    if 'sample' in act_on:
        weights = np.mean(np.square(U), axis=1)
    elif 'feature' in act_on:
        weights = np.mean(np.square(VT), axis=0)
    if is_deterministic:
        # sorting is smallest to largest hence the minus
        sel = np.argsort(-weights)[:Nsel]
    elif is_deterministic is False:
        np.random.seed(seed)
        # sorting is smallest to largest hence the minus
        sel = np.argsort(np.random.rand(*weights.shape) - weights)[:Nsel]

    if verbose:
        if 'sample' in act_on:
            C = X[sel, :]
            # equivalent to
            # Cp = np.linalg.pinv(C)
            # err = np.sqrt(np.sum((X - np.dot(np.dot(X, Cp), C))**2))
            err = np.sqrt(np.sum((
                X - np.dot(np.linalg.lstsq(C.T, X.T, rcond=None)[0].T, C))**2))

        elif 'feature' in act_on:
            C = X[:, sel]
            # equivalent to
            # Cp = np.linalg.pinv(C)
            # err = np.sqrt(np.sum((X - np.dot(C, np.dot(Cp, X)))**2))
            err = np.sqrt(np.sum((
                X - np.dot(C, np.linalg.lstsq(C, X, rcond=None)[0]))**2))

        print('Reconstruction RMSE={:.3e}'.format(err))

    return sel


class CURFilter(BaseIO):
    """CUR decomposition to select samples or features in a given feature matrix.
    Wrapper around the do_CUR function for convenience.

    Parameters
    ----------
    representation : Calculator
        Representation calculator associated with the kernel

    Nselect: int
        number of points to select. if act_on='sample per species' then it should
        be a dictionary mapping atom type to the number of samples, e.g.
        Nselect = {1:200,6:100,8:50}.

    act_on: string
        Select how to apply the selection. Can be either of 'sample',
        'sample per species','feature'.
        For the moment only 'sample per species' is implemented.

    is_deterministic: bool
        flag to switch between selction criteria

    seed: int
        if is_deterministic==False, seed for the random selection

    """

    def __init__(self, representation, Nselect, act_on='sample per species',
                 is_deterministic=True, seed=10):
        self._representation = representation
        self.Nselect = Nselect
        modes = ['sample', 'sample per species', 'feature']
        if act_on in modes:
            self.act_on = act_on
        else:
            raise ValueError(
                '"act_on" should be either of: "{}", "{}", "{}"'.format(*modes))
        self.is_deterministic = is_deterministic
        self.seed = seed
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

    def select(self, managers):
        """Perform CUR selection of samples/features.

        Parameters
        ----------
        managers : AtomsList
            list of structures containing features computed with representation

        Returns
        -------
        SparsePoints
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
                self.selected_sample_ids_by_sp[sp] = do_CUR(X_by_sp[sp], self.Nselect[sp], self.act_on,
                                                            self.is_deterministic, self.seed)

        elif self.act_on == 'sample':
            self.selected_sample_ids = do_CUR(X, self.Nselect, self.act_on,
                                              self.is_deterministic, self.seed)
        elif self.act_on == 'feature':
            self.selected_feature_ids_global = do_CUR(X, self.Nselect, self.act_on,
                                                      self.is_deterministic, self.seed)
        else:
            raise ValueError("method: {}".format(self.act_on))

        return self

    def filter(self, managers, n_select=None):
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

    def _get_data(self):
        data = super()._get_data()
        data.update(selected_ids=self.selected_ids,
                    selected_sample_ids=self.selected_sample_ids,
                    selected_sample_ids_by_sp=self.selected_sample_ids_by_sp,
                    selected_feature_ids_global=self.selected_feature_ids_global)
        return data

    def _set_data(self, data):
        super()._set_data(data)
        self.selected_ids = data['selected_ids']
        self.selected_sample_ids = data['selected_sample_ids']
        self.selected_sample_ids_by_sp = data['selected_sample_ids_by_sp']
        self.selected_feature_ids_global = data['selected_feature_ids_global']

    def _get_init_params(self):
        return dict(representation=self._representation,
                    Nselect=self.Nselect,
                    act_on=self.act_on,
                    is_deterministic=self.is_deterministic,
                    seed=self.seed,)
