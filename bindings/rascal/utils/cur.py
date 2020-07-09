from ..utils import BaseIO
from ..utils.filter_utils import (get_index_mappings_sample_per_species,
                                  convert_selected_global_index2perstructure_index_per_species,
                                  get_index_mappings_sample,
                                  convert_selected_global_index2perstructure_index)

from ..models.sparse_points import SparsePoints
import numpy as np
from numpy import linalg
import scipy.sparse.linalg as sparselinalg


def get_rank_svals(Svals, tol_ratio=None, max_size=None):
    """Get the rank of a matrix from its singular values

    Parameters:
        Svals   The vector of singular values
    Optional parameters:
        tol_ratio   The cutoff for a singular value to be considered
                    "zero" for the purposes of determining effective
                    rank, as a ratio to the largest singular value in
                    Svals
                    Default is to use the NumPy criterion (eps*size)
        max_size    Largest matrix dimension for determining singular
                    value cutoff
    """
    if max_size is None:
        max_size = max(Svals.shape)
    if tol_ratio is None:
        tol_ratio = np.finfo(Svals.dtype).eps * max_size
    return np.sum(Svals > (np.max(Svals) * tol_ratio))

def do_CUR(X, Nsel, act_on='sample', is_deterministic=False, seed=10,
           use_sparse_svd=False, verbose=True):
    """ Apply CUR selection [1] to the input matrix.

    Return a list of indices to select.

    Parameters:
        X       The matrix to sparsify
        Nsel    The number of features to keep (must be between 1 and
                the size of the matrix along the selected dimension)
    Optional parameters:
        act_on  Whether to act on rows or columns of the matrix.  If
                this string contains 'sample' (the default), selection
                acts on rows.  If it contains 'feature', selection acts
                on columns.  If it contains both or neither, raise a
                ValueError.
        is_deterministic
                Whether to use the "deterministic" selection algorithm
                described in [2], rather than the original probabilistic
                criterion
        seed    Random seed for probabilistic selection (default 10. You
                may want to change this.)
        use_sparse_svd
                Whether to use sparse SVD (i.e., Lanczos iteration) to
                decompose the matrix product, rather than the full
                version.  May be more efficient when selecting a small
                number of entries from a very large number of rows or
                columns.  Default False.
        verbose Print the reconstruction error? (default True)

    .. [1] Mahoney, M. W., & Drineas, P. (2009). CUR matrix decompositions for
    improved data analysis. Proceedings of the National Academy of Sciences,
    106(3), 697–702. https://doi.org/10.1073/pnas.0803205106
    .. [2] Imbalzano, G.; Anelli, A.; Giofré, D.; Klees, S.; Behler, J.;
    Ceriotti, M. Automatic Selection of Atomic Fingerprints and Reference
    Configurations for Machine-Learning Potentials. J. Chem. Phys. 2018, 148
    (24), 241730. https://doi.org/10.1063/1.5024611.
    """
    # TODO sparse SVD does not work for the case M < k < N (or vice versa).
    #      (the statistical "leverage scores" are all the same because there
    #       are at most min(M, N) - 1 singular vectors.  But this doesn't
    #       mean a perfect reconstruction is possible with only min(M, N)
    #       rows / columns of X, does it...?
    if use_sparse_svd:
        max_singular_vals = min(X.shape) - 1
        num_singular_vals = min(max_singular_vals, Nsel)
        U, Svals, VT = sparselinalg.svds(X, num_singular_vals)
    else:
        U, Svals, VT = linalg.svd(X)
        # U = U[:, :Nsel]
        # VT = VT[:Nsel, :]
    X_rank = get_rank_svals(Svals, max_size=max(X.shape))
    num_span_vectors = min(X_rank, Nsel)
    U = U[:, :num_span_vectors]
    VT = VT[:num_span_vectors, :]
    if ('sample' in act_on) and ('feature' in act_on):
        raise ValueError("Must supply either 'sample' or 'feature' in 'act_on'"
                         " string, not both")
    elif 'sample' in act_on:
        weights = np.mean(np.square(U), axis=1)
    elif 'feature' in act_on:
        weights = np.mean(np.square(VT), axis=0)
    else:
        raise ValueError("Must supply either 'sample' or 'feature' in 'act_on'"
                         " string.")
    if is_deterministic:
        # sorting is smallest to largest hence the minus
        sel = np.argsort(-weights)[:Nsel]
    else:
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
            TODO wrong! This only returns self!

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
