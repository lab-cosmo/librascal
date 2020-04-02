from rascal.representations import SphericalInvariants
from rascal.models import Kernel
from rascal.utils import from_dict, to_dict
from test_utils import load_json_frame, BoxList, Box
import unittest
import numpy as np
import sys
import json


class TestCosineKernel(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = 'reference_data/inputs/CaCrP2O7_mvc-11955_symmetrized.json'
        self.frame = load_json_frame(fn)

        self.hypers = dict(soap_type="PowerSpectrum",
                           interaction_cutoff=3.5,
                           max_radial=6,
                           max_angular=6,
                           gaussian_sigma_constant=0.4,
                           gaussian_sigma_type="Constant",
                           cutoff_smooth_width=0.5,
                           )

    def test_model_call(self):

        rep = SphericalInvariants(**self.hypers)

        features = rep.transform([self.frame])

        for target_type in ["Atom", "Structure"]:
            cosine_kernel = Kernel(
                rep, name="Cosine", target_type=target_type, zeta=2)
            cosine_kernel(features)

        # wrong name
        with self.assertRaises(RuntimeError):
            Kernel(rep, name="WrongName", target_type="Structure", zeta=2)
        with self.assertRaises(RuntimeError):
            Kernel(rep, name="cosine", target_type="Structure", zeta=2)
        # wrong target_type
        with self.assertRaises(RuntimeError):
            Kernel(rep, name="Cosine", target_type="WrongType", zeta=2)
        with self.assertRaises(RuntimeError):
            Kernel(rep, name="Cosine", target_type="structure", zeta=2)
        with self.assertRaises(RuntimeError):
            Kernel(rep, name="Cosine", target_type="atom", zeta=2)
        # wrong zeta
        with self.assertRaises(ValueError):
            Kernel(rep, name="Cosine", target_type="Structure", zeta=2.5)
        with self.assertRaises(ValueError):
            Kernel(rep, name="Cosine", target_type="Structure", zeta=-2)

    def test_serialization(self):
        rep = SphericalInvariants(**self.hypers)

        for target_type in ["Atom", "Structure"]:
            cosine_kernel = Kernel(
                rep, name="Cosine", target_type=target_type, zeta=2)

            cosine_kernel_dict = to_dict(cosine_kernel)
            cosine_kernel_copy = from_dict(cosine_kernel_dict)
            cosine_kernel_copy_dict = to_dict(cosine_kernel_copy)

            self.assertTrue(cosine_kernel_dict == cosine_kernel_copy_dict)
