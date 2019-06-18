from ..lib import RepresentationManager, FeatureManager
from ..utils.pool_worker import FactoryPool
from ..neighbourlist.base import NeighbourListFactory
from ..neighbourlist.structure_manager import convert_to_structure

import numpy as np
import queue

_representations_list = ["sortedcoulomb", "sphericalexpansion", "soap"]
_representations = {}
for k, v in RepresentationManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _representations_list:
            if name in kl:
                _representations[kl] = v

_features_list = ["dense", "blocksparse"]
_features = {}
for k, v in FeatureManager.__dict__.items():
    if "pybind11_builtins.pybind11_type" in str(type(v)):
        kl = k.lower()
        for name in _features_list:
            if name in kl:
                _features[kl] = v


def RepresentationFactory(manager, rep_options):
    name = rep_options['name']
    if name not in _representations:
        raise NameError('The representations factory {} has not been registered. The available combinations are: {}'.format(
            name, list(_representations.keys())))
    return _representations[name](manager, *rep_options['args'])


def FeatureFactory(feature_options):
    name = feature_options['name']
    if name not in _features:
        raise NameError('The features factory {} has not been registered. The available combinations are: {}'.format(
            name, list(_features.keys())))
    return _features[name](*feature_options['args'])


class RepresentationRunner(object):
    def __init__(self, nl_options, rep_options, feature_options, method='thread', n_workers=1, disable_pbar=False):
        self.n_workers = n_workers
        self.managers = [NeighbourListFactory(
            nl_options) for _ in range(self.n_workers)]
        self.representations = [RepresentationFactory(
            manager, rep_options) for manager in self.managers]
        self.feature_options = feature_options
        self.method = method

    def run(self, frames):
        structures = [convert_to_structure(frame) for frame in frames]

        n_atoms = [0]+[len(structure['atom_types'])
                       for structure in structures]
        structure_ids = np.cumsum(n_atoms)[:-1]
        n_centers = np.sum(n_atoms)

        features = FeatureFactory(self.feature_options)
        features.resize(n_centers)

        pool = FactoryPool(self.method, self.n_workers)
        index_q = queue.Queue()
        for ii in range(self.n_workers):
            index_q.put(ii)
        feature_q = queue.Queue()
        feature_q.put(features)

        inputs = [(self.managers, self.representations, it, structure, feature_q, index_q)
                  for it, structure in zip(structure_ids, structures)]

        pool.starmap(compute_wrapper, inputs)

        pool.close()
        pool.join()

        return features


def compute_wrapper(managers, representations, structure_id, structure, feature_q, index_q):
    idx = index_q.get()
    managers[idx].update(**structure)
    representations[idx].compute()
    features = feature_q.get()
    features.assign(structure_id, representations[idx])
    feature_q.put(features)
    index_q.put(idx)
