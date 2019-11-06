from rascal.representations import (SortedCoulombMatrix, SphericalExpansion,
                                    SphericalInvariants)
from test_utils import load_json_frame, BoxList, Box, dot
import unittest
import numpy as np
import sys
import json
from copy import copy


class TestSortedCoulombRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = 'reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.hypers = dict(cutoff=3., sorting_algorithm='row_norm',
                           size=50, central_decay=0.5,
                           interaction_cutoff=3, interaction_decay=-1)

    def test_representation_transform(self):

        rep = SortedCoulombMatrix(**self.hypers)

        features = rep.transform([self.frame])

        test = features.get_dense_feature_matrix(rep)


class TestSphericalExpansionRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = 'reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frames = [load_json_frame(fn)]
        fn = 'reference_data/SiC_moissanite_supercell.json'
        self.frames += [load_json_frame(fn)]
        fn = 'reference_data/methane.json'
        self.frames += [load_json_frame(fn)]

        self.hypers = {"interaction_cutoff": 6.0,
                       "cutoff_smooth_width": 1.0,
                       "max_radial": 10,
                       "max_angular": 8,
                       "gaussian_sigma_type": "Constant",
                       "gaussian_sigma_constant": 0.5}

    def test_representation_transform(self):

        rep = SphericalExpansion(**self.hypers)

        features = rep.transform(self.frames)

        test = features.get_dense_feature_matrix(rep)


class TestSphericalInvariantsRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = 'reference_data/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frames = [load_json_frame(fn)]
        fn = 'reference_data/SiC_moissanite.json'
        self.frames += [load_json_frame(fn)]
        fn = 'reference_data/methane.json'
        self.frames += [load_json_frame(fn)]

        global_species = []
        for frame in self.frames:
            global_species.extend(frame['atom_types'])
        self.global_species = list(np.unique(global_species))

        self.hypers = dict(soap_type="PowerSpectrum",
                           interaction_cutoff=3.5,
                           max_radial=6,
                           max_angular=6,
                           gaussian_sigma_constant=0.4,
                           gaussian_sigma_type="Constant",
                           cutoff_smooth_width=0.5,
                           )

    def test_representation_transform(self):

        rep = SphericalInvariants(**self.hypers)

        features = rep.transform(self.frames)

        test = features.get_dense_feature_matrix(rep)
        kk_ref = np.dot(test, test.T)

        # test that the feature matrix exported to python in various ways
        # are equivalent
        X_t = features.get_dense_feature_matrix(rep, self.global_species)
        kk = np.dot(X_t, X_t.T)
        self.assertTrue(np.allclose(kk, kk_ref))

        X_t = features.get_dense_feature_matrix(rep, self.global_species+[70])
        kk = np.dot(X_t, X_t.T)
        self.assertTrue(np.allclose(kk, kk_ref))

        species = copy(self.global_species)
        species.pop()
        X_t = features.get_dense_feature_matrix(rep, species)
        kk = np.dot(X_t, X_t.T)
        self.assertFalse(np.allclose(kk, kk_ref))

        X_t = features.get_sparse_feature_matrix(rep)
        kk = dot(X_t, X_t)
        self.assertTrue(np.allclose(kk, kk_ref))
