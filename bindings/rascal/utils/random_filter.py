from .io import BaseIO
from ..models import SparsePoints
from .filter_utils import (get_index_mappings_sample_per_species,
                           convert_selected_global_index2perstructure_index_per_species,
                           get_index_mappings_sample,
                           convert_selected_global_index2perstructure_index)
import numpy as np


class RandomFilter(BaseIO):
    """Randomly samples or features in a given feature matrix.
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
                 seed=10):
        self._representation = representation
        self.Nselect = Nselect
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
        self.seed = seed

    def select(self, managers):
        return self

    def filter(self, managers, n_select=None):
        if n_select is None:
            n_select = self.Nselect
        np.random.seed(self.seed)

        if self.act_on == 'sample per species':
            sps = list(n_select.keys())
            # get various info from the structures about the center atom species and indexing
            (strides_by_sp, global_counter, map_by_manager,
             indices_by_sp) = get_index_mappings_sample_per_species(managers, sps)

            self.selected_sample_ids_by_sp = {}
            for sp in sps:
                ids = np.arange(global_counter[sp])
                np.random.shuffle(ids)
                self.selected_sample_ids_by_sp[sp] = ids

            selected_ids_by_sp = {key: np.sort(val[:n_select[key]])
                                  for key, val in self.selected_sample_ids_by_sp.items()}
            self.selected_ids = convert_selected_global_index2perstructure_index_per_species(
                managers, selected_ids_by_sp, strides_by_sp, map_by_manager, sps)

            pseudo_points = SparsePoints(self._representation)
            pseudo_points.extend(managers, self.selected_ids)
            return pseudo_points

        elif self.act_on == 'sample':
            strides, _, map_by_manager = get_index_mappings_sample(managers)
            Natoms = strides[-1]
            ids = np.arange(Natoms)
            np.random.shuffle(ids)
            self.selected_sample_ids = ids
            selected_ids_global = np.sort(ids[:n_select])

            self.selected_ids = convert_selected_global_index2perstructure_index(managers,
                                                                            selected_ids_global, strides, map_by_manager)
            return self.selected_ids

        elif self.act_on == 'feature':
            feat_idx2coeff_idx = self._representation.get_feature_index_mapping(
                managers)
            Nfeat = np.max(list(feat_idx2coeff_idx.keys()))
            ids = np.arange(Nfeat)
            np.random.shuffle(ids)
            self.selected_feature_ids_global = ids
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
                    seed=self.seed,)