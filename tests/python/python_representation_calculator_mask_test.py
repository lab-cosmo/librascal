from rascal.representations import (
    SortedCoulombMatrix,
    SphericalExpansion,
    SphericalInvariants,
)
from test_utils import load_json_frame, BoxList, Box, dot
import unittest
import numpy as np
from scipy.spatial.distance import cdist
import sys
import os
import json
from copy import copy, deepcopy
from scipy.stats import ortho_group
import pickle

rascal_reference_path = "reference_data"
inputs_path = os.path.join(rascal_reference_path, "inputs")
dump_path = os.path.join(rascal_reference_path, "tests_only")



class TestSphericalExpansionMask(unittest.TestCase):
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
        feature_size_per_species = int(hypers["max_radial"] * (hypers["max_angular"]+1)**2)
        ref_features = rep.transform(self.frames).get_features_by_species(rep)
        ref_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        hypers_neighbour_mask = deepcopy(hypers)
        hypers_neighbour_mask["neighbour_species"] = [8]
        rep = SphericalExpansion(**hypers_neighbour_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        for i, species_key in enumerate(ref_features.keys()):
            feature_species_slice = slice(i*feature_size_per_species, (i+1)*feature_size_per_species)
            if species_key == (8,):
                self.assertTrue(np.allclose(totest_features_gradient[:, feature_species_slice], totest_features_gradient[:, feature_species_slice]))
            else:
                self.assertTrue(np.allclose(np.zeros(totest_features_gradient[:, feature_species_slice].shape), totest_features_gradient[:, feature_species_slice]))

    def test_central_mask_all_species_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features(rep)

        hypers_central_mask = deepcopy(self.hypers)
        hypers_central_mask["central_species"] = self.species
        rep = SphericalExpansion(**hypers_central_mask)
        totest_features = rep.transform(self.frames).get_features(rep)

        self.assertTrue(np.allclose(ref_features, totest_features))

    def test_central_mask_oxygen_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features(rep)

        hypers_central_mask = deepcopy(self.hypers)
        hypers_central_mask["central_species"] = [8]
        rep = SphericalExpansion(**hypers_central_mask)
        totest_features = rep.transform(self.frames).get_features(rep)
        
        central_species = np.concatenate([frame["atom_types"] for frame in self.frames]).flatten()
        central_oxygen_mask = central_species == 8

        self.assertTrue(np.allclose(ref_features[central_oxygen_mask], totest_features[central_oxygen_mask]))
        self.assertTrue(np.allclose(np.zeros(ref_features[~central_oxygen_mask].shape), totest_features[~central_oxygen_mask]))


    def test_central_mask_all_species_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True

        rep = SphericalExpansion(**hypers)
        ref_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        hypers_central_mask = deepcopy(hypers)
        hypers_central_mask["central_species"] = self.species
        rep = SphericalExpansion(**hypers_central_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        self.assertTrue(np.allclose(ref_features_gradient, totest_features_gradient))
        
    def test_central_mask_oxygen_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True
    
        rep = SphericalExpansion(**hypers)
        manager = rep.transform(self.frames)
        ref_features_gradient = manager.get_features_gradient(rep)
        gradients_info = manager.get_gradients_info()
        # chooses all pair gradients that have oxygen as center or neighbour
        computed_gradients = np.logical_or(gradients_info[:, 3] == 8, gradients_info[:, 4] == 8)
        computed_gradients = np.concatenate([[pair_mask]*3 for pair_mask in computed_gradients])
    
        hypers_central_mask = deepcopy(hypers)
        hypers_central_mask["central_species"] = [8]
        rep = SphericalExpansion(**hypers_central_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)
        
        self.assertTrue(np.allclose(ref_features_gradient[computed_gradients], totest_features_gradient[computed_gradients]))
        self.assertTrue(np.allclose(np.zeros(ref_features_gradient[~computed_gradients].shape), totest_features_gradient[~computed_gradients]))


    def test_neighbour_and_central_mask_all_species_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features(rep)

        hypers_mask = deepcopy(self.hypers)
        hypers_mask["neighbour_species"] = self.species
        hypers_mask["central_species"] = self.species
        rep = SphericalExpansion(**hypers_mask)
        totest_features = rep.transform(self.frames).get_features(rep)

        self.assertTrue(np.allclose(ref_features, totest_features))

    def test_neighbour_and_central_mask_oxygen_features(self):
        rep = SphericalExpansion(**self.hypers)
        ref_features = rep.transform(self.frames).get_features_by_species(rep)

        hypers_mask = deepcopy(self.hypers)
        hypers_mask["neighbour_species"] = [8]
        hypers_mask["central_species"] = [8]
        rep = SphericalExpansion(**hypers_mask)
        totest_features = rep.transform(self.frames).get_features_by_species(rep)

        central_species = np.concatenate([frame["atom_types"] for frame in self.frames]).flatten()
        central_oxygen_mask = central_species == 8

        self.assertTrue(np.allclose(ref_features[(8,)][central_oxygen_mask], totest_features[(8,)][central_oxygen_mask]))
        for species_key in totest_features.keys():
            if species_key != (8,):
                self.assertTrue(np.allclose(np.zeros(ref_features[species_key].shape), totest_features[species_key]))

            self.assertTrue(np.allclose(np.zeros(ref_features[species_key][~central_oxygen_mask].shape), totest_features[species_key][~central_oxygen_mask]))

    def test_neighbour_and_central_mask_all_species_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True

        rep = SphericalExpansion(**hypers)
        ref_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        hypers_mask = deepcopy(hypers)
        hypers_mask["neighbour_species"] = self.species
        hypers_mask["central_species"] = self.species
        rep = SphericalExpansion(**hypers_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        self.assertTrue(np.allclose(ref_features_gradient, totest_features_gradient))

    def test_neighbour_central_mask_oxygen_features_gradient(self):
        hypers = deepcopy(self.hypers)
        hypers["compute_gradients"] = True

        rep = SphericalExpansion(**hypers)
        feature_size_per_species = int(hypers["max_radial"] * (hypers["max_angular"]+1)**2)
        manager = rep.transform(self.frames)
        ref_features = manager.get_features_by_species(rep)
        ref_features_gradient = manager.get_features_gradient(rep)

        gradients_info = manager.get_gradients_info()
        # chooses all pair gradients that have oxygen as center or neighbour
        computed_gradients = np.logical_or(gradients_info[:, 3] == 8, gradients_info[:, 4] == 8)
        computed_gradients = np.concatenate([[pair_mask]*3 for pair_mask in computed_gradients])

        hypers_mask = deepcopy(hypers)
        hypers_mask["neighbour_species"] = [8]
        hypers_mask["central_species"] = [8]
        rep = SphericalExpansion(**hypers_mask)
        totest_features_gradient = rep.transform(self.frames).get_features_gradient(rep)

        self.assertTrue(np.allclose(np.zeros(ref_features_gradient[~computed_gradients].shape), totest_features_gradient[~computed_gradients]))

        for i, species_key in enumerate(ref_features.keys()):
            feature_species_slice = slice(i*feature_size_per_species, (i+1)*feature_size_per_species)
            if species_key == (8,):
                self.assertTrue(np.allclose(totest_features_gradient[computed_gradients][:, feature_species_slice], totest_features_gradient[computed_gradients][:, feature_species_slice]))
            else:
                self.assertTrue(np.allclose(np.zeros(totest_features_gradient[computed_gradients][:, feature_species_slice].shape), totest_features_gradient[computed_gradients][:, feature_species_slice]))
