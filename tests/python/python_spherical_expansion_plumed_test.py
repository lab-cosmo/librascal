import sys
sys.path.insert(0, "/home/goscinsk/code/librascal-centralmask-noneighbours/build")
import rascal
print(rascal.__file__)
from rascal.representations import (
    SortedCoulombMatrix,
    SphericalExpansion,
    SphericalInvariants,
)
from test_utils import load_json_frame, BoxList, Box, dot
import unittest
import numpy as np
import sys
import os
import json
from copy import copy, deepcopy
from scipy.stats import ortho_group
import pickle

rascal_reference_path = "reference_data"
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join(rascal_reference_path, "tests_only")



class TestSphericalExpansionRepresentation(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """
        self.seed = 18012021

        fns = [
            os.path.join(inputs_path, "CaCrP2O7_mvc-11955_symmetrized.json"),
            os.path.join(inputs_path, "SiC_moissanite_supercell.json"),
            os.path.join(inputs_path, "methane.json"),
        ]
        self.species = [1, 6, 8, 14, 15, 20, 24]
        self.frames = [load_json_frame(fn) for fn in fns]

        self.max_radial = 6
        self.max_angular = 4

        self.hypers = {
            "interaction_cutoff": 6.0,
            "cutoff_smooth_width": 1.0,
            "max_radial": self.max_radial,
            "max_angular": self.max_angular,
            "gaussian_sigma_type": "Constant",
            "gaussian_sigma_constant": 0.5,
        }

    def test_neighbour_mask_all_species_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features(rep)

        hypers_neighbour_mask = deepcopy(self.hypers)
        hypers_neighbour_mask["neighbour_species"] = self.species
        rep = SphericalExpansion(**hypers_neighbour_mask)
        totest_features = rep.transform(self.frames).get_features(rep)

        self.assertTrue(np.allclose(ref_features, totest_features))

    def test_neighbour_mask_oxygen_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features_by_species(rep)

        hypers_neighbour_mask = deepcopy(self.hypers)
        hypers_neighbour_mask["neighbour_species"] = [8]
        rep = SphericalExpansion(**hypers_neighbour_mask)
        totest_features = rep.transform(self.frames).get_features_by_species(rep)

        self.assertTrue(np.allclose(ref_features[(8,)], totest_features[(8,)]))
        for species_key in totest_features.keys():
            if species_key != (8,):
                self.assertTrue(np.allclose(np.zeros(ref_features[species_key].shape), totest_features[species_key]))

    def test_neighbour_mask_all_species_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True

        rep = SphericalExpansion(**hypers)
        ref_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        hypers_neighbour_mask = deepcopy(hypers)
        hypers_neighbour_mask["neighbour_species"] = self.species
        rep = SphericalExpansion(**hypers_neighbour_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        self.assertTrue(np.allclose(ref_features_gradient, totest_features_gradient))

    def test_neighbour_mask_oxygen_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True

        rep = SphericalExpansion(**hypers)
        ref_features = rep.transform(self.frames).get_features_by_species(rep)
        feature_size_per_species = int(hypers["max_radial"] * (hypers["max_angular"]+1)**2)
        #first_key = next(iter(ref_features.keys()))[0]
        #feature_size_per_species = ref_features[first_key].shape[]
        ref_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        hypers_neighbour_mask = deepcopy(hypers)
        hypers_neighbour_mask["neighbour_species"] = [8]
        rep = SphericalExpansion(**hypers_neighbour_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        for i, species_key in enumerate(ref_features.keys()):
            #feature_species_slice = np.index_exp[:, 1:3]
            feature_species_slice= slice(i*feature_size_per_species, (i+1)*feature_size_per_species)
            if species_key == (8,):
                self.assertTrue(np.allclose(totest_features_gradient[:, feature_species_slice], totest_features_gradient[:, feature_species_slice]))
            else:
                self.assertTrue(np.allclose(np.zeros(totest_features_gradient[:, feature_species_slice].shape), totest_features_gradient[:, feature_species_slice]))

    #def test_radial_dimension_reduction_test(self):
    #    rep = SphericalExpansion(**self.hypers)
    #    features_ref = rep.transform(self.frames).get_features(rep)
    #    K_ref = features_ref.dot(features_ref.T)

    #    # identity test
    #    hypers = deepcopy(self.hypers)
    #    projection_matrices = [
    #        np.eye(self.max_radial).tolist() for _ in range(self.max_angular + 1)
    #    ]

    #    hypers["optimization"] = {
    #        "Spline": {"accuracy": 1e-8},
    #        "RadialDimReduction": {
    #            "projection_matrices": {sp: projection_matrices for sp in self.species}
    #        },
    #    }
    #    rep = SphericalExpansion(**hypers)
    #    features_test = rep.transform(self.frames).get_features(rep)
    #    self.assertTrue(np.allclose(features_ref, features_test))

    #    # random orthogonal matrix test,
    #    # for angular_l=0 we can do the projection on the python side
    #    hypers = deepcopy(self.hypers)
    #    np.random.seed(self.seed)
    #    projection_matrices = {
    #        sp: [
    #            ortho_group.rvs(self.max_radial)[: self.max_radial - 1, :].tolist()
    #            for _ in range(self.max_angular + 1)
    #        ]
    #        for sp in self.species
    #    }

    #    hypers["max_radial"] = self.max_radial - 1
    #    hypers["optimization"] = {
    #        "Spline": {"accuracy": 1e-8},
    #        "RadialDimReduction": {"projection_matrices": projection_matrices},
    #    }
    #    rep = SphericalExpansion(**hypers)
    #    features_test = rep.transform(self.frames).get_features(rep)
    #    features_test = features_test.reshape(
    #        len(features_test),
    #        len(self.species),
    #        self.max_radial - 1,
    #        (self.max_angular + 1) ** 2,
    #    )

    #    hypers["max_radial"] = self.max_radial
    #    hypers["optimization"] = {
    #        "Spline": {"accuracy": 1e-8},
    #    }
    #    rep = SphericalExpansion(**hypers)
    #    features_ref = rep.transform(self.frames).get_features(rep)
    #    features_ref = features_ref.reshape(
    #        len(features_ref),
    #        len(self.species),
    #        self.max_radial,
    #        (self.max_angular + 1) ** 2,
    #    )

    #    for a, species in enumerate(self.species):
    #        self.assertTrue(
    #            np.allclose(
    #                features_ref[:, a, :, 0]
    #                @ np.array(projection_matrices[species][0]).T,
    #                features_test[:, a, :, 0],
    #            )
    #        )

    #    # checks if error is raised if wrong number of species is given
    #    with self.assertRaises(RuntimeError):
    #        species = [1]
    #        hypers["optimization"] = {
    #            "Spline": {"accuracy": 1e-8},
    #            "RadialDimReduction": {
    #                "projection_matrices": {sp: projection_matrices for sp in species}
    #            },
    #        }
    #        rep = SphericalExpansion(**hypers)
    #        features_test = rep.transform(self.frames).get_features(rep)

