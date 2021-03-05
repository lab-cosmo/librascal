from rascal.representations import (
    SortedCoulombMatrix,
    SphericalExpansion,
    SphericalInvariants,
)
from rascal.utils import from_dict, to_dict, FPSFilter
from rascal.models import Kernel, SparsePoints
from test_utils import load_json_frame, BoxList, Box, dot
import unittest
import numpy as np
import sys
import os
import json
from copy import copy, deepcopy

rascal_reference_path = "reference_data"
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join(rascal_reference_path, "tests_only")


class TestSortedCoulombRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = os.path.join(inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json")
        self.frame = load_json_frame(fn)

        self.hypers = dict(
            cutoff=3.0,
            sorting_algorithm="row_norm",
            size=50,
            central_decay=0.5,
            interaction_cutoff=3,
            interaction_decay=-1,
        )

    def test_representation_transform(self):

        rep = SortedCoulombMatrix(**self.hypers)

        features = rep.transform([self.frame])

        test = features.get_features(rep)

    def test_serialization(self):
        rep = SortedCoulombMatrix(**self.hypers)

        rep_dict = to_dict(rep)

        rep_copy = from_dict(rep_dict)

        rep_copy_dict = to_dict(rep_copy)

        self.assertTrue(rep_dict == rep_copy_dict)


class TestSphericalExpansionRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fns = [
            os.path.join(inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
            os.path.join(inputs_path, "SiC_moissanite_supercell.json"),
            os.path.join(inputs_path, "methane.json"),
        ]
        self.frames = [load_json_frame(fn) for fn in fns]

        self.hypers = {
            "interaction_cutoff": 6.0,
            "cutoff_smooth_width": 1.0,
            "max_radial": 10,
            "max_angular": 8,
            "gaussian_sigma_type": "Constant",
            "gaussian_sigma_constant": 0.5,
        }

    def test_representation_transform(self):

        rep = SphericalExpansion(**self.hypers)

        features = rep.transform(self.frames)

        test = features.get_features(rep)

    def test_serialization(self):
        rep = SphericalExpansion(**self.hypers)

        rep_dict = to_dict(rep)

        rep_copy = from_dict(rep_dict)

        rep_copy_dict = to_dict(rep_copy)

        self.assertTrue(rep_dict == rep_copy_dict)


class TestSphericalInvariantsRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fns = [
            os.path.join(inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
            os.path.join(inputs_path, "SiC_moissanite_supercell.json"),
            os.path.join(inputs_path, "methane.json"),
        ]
        self.frames = [load_json_frame(fn) for fn in fns]

        global_species = []
        for frame in self.frames:
            global_species.extend(frame["atom_types"])
        self.global_species = list(np.unique(global_species))

        self.hypers = dict(
            soap_type="PowerSpectrum",
            interaction_cutoff=3.5,
            max_radial=2,
            max_angular=1,
            gaussian_sigma_constant=0.4,
            gaussian_sigma_type="Constant",
            cutoff_smooth_width=0.5,
        )

    def test_representation_transform(self):

        rep = SphericalInvariants(**self.hypers)

        features = rep.transform(self.frames)

        test = features.get_features(rep)
        kk_ref = np.dot(test, test.T)

        # test that the feature matrix exported to python in various ways
        # are equivalent
        X_t = features.get_features(rep, self.global_species)
        kk = np.dot(X_t, X_t.T)
        self.assertTrue(np.allclose(kk, kk_ref))

        X_t = features.get_features(rep, self.global_species + [70])
        kk = np.dot(X_t, X_t.T)
        self.assertTrue(np.allclose(kk, kk_ref))

        species = copy(self.global_species)
        species.pop()
        X_t = features.get_features(rep, species)
        kk = np.dot(X_t, X_t.T)
        self.assertFalse(np.allclose(kk, kk_ref))

        X_t = features.get_features_by_species(rep)
        kk = dot(X_t, X_t)
        self.assertTrue(np.allclose(kk, kk_ref))

    def test_representation_gradient(self):
        """
        Test the get_features and get_features_gradient functions by computing
        the linear sparse kernel matrix and check that the exported features
        lead to the same kernel matrix as the reference method.
        """
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True
        rep = SphericalInvariants(**hypers)

        features = rep.transform(self.frames)

        n_sparses = {1: 1, 6: 1, 8: 1, 14: 1, 15: 1, 20: 1, 24: 1}

        compressor = FPSFilter(rep, n_sparses, act_on="sample per species")
        X_pseudo = compressor.select_and_filter(features)

        xs = X_pseudo.get_features()
        n_sparse, n_feat = xs.shape
        masks = {sp: np.zeros(n_sparse, dtype=bool) for sp in n_sparses}
        ii = 0
        for sp, mask in masks.items():
            mask[ii : ii + n_sparses[sp]] = 1
            ii = ii + n_sparses[sp]

        zeta = 1
        kernel = Kernel(
            rep, name="GAP", zeta=zeta, target_type="Structure", kernel_type="Sparse"
        )

        ij = features.get_neighbors_for_gradient()
        n_atoms = len(np.unique(ij[:, 0]))
        n_neigh = ij.shape[0]

        KNM_ref = kernel(features, X_pseudo, (False, False))
        X = features.get_features(rep).reshape((n_atoms, n_feat))
        KNM = np.zeros((len(self.frames), n_sparse))
        ii = 0
        for iff, frame in enumerate(features):
            for at in frame:
                sp = at.atom_type
                KNM[iff, masks[sp]] += np.dot(X[ii], xs[masks[sp]].T)
                ii += 1
        self.assertTrue(np.allclose(KNM_ref, KNM))

        KNM_ref = kernel(features, X_pseudo, (True, False))

        X_der = features.get_features_gradient(rep).reshape((n_neigh, 3, n_feat))

        KNM = np.zeros((n_atoms, 3, n_sparse))
        for ii, (i, j, i_sp, j_sp) in enumerate(ij):
            sp = i_sp
            KNM[j, 0, masks[sp]] += np.dot(X_der[ii, 0], xs[masks[sp]].T)
            KNM[j, 1, masks[sp]] += np.dot(X_der[ii, 1], xs[masks[sp]].T)
            KNM[j, 2, masks[sp]] += np.dot(X_der[ii, 2], xs[masks[sp]].T)

        KNM = KNM.reshape((-1, n_sparse))

        self.assertTrue(np.allclose(KNM_ref, KNM))

    def test_serialization(self):
        rep = SphericalInvariants(**self.hypers)

        rep_dict = to_dict(rep)

        rep_copy = from_dict(rep_dict)

        rep_copy_dict = to_dict(rep_copy)

        self.assertTrue(rep_dict == rep_copy_dict)
