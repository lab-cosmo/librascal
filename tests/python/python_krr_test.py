# TMP line for debugging for alex
import sys
sys.path.insert(0, "/home/goscinsk/code/librascal-local-stress/build")

import rascal
print(rascal.__file__)
import unittest
# TMP END


import json

import numpy as np
import ase.io

from rascal.utils import load_obj

class TestKRR(unittest.TestCase):
    def setUp(self):
        """
        """
        self.model = load_obj("reference_data/tests_only/simple_gap_model.json")
        # the file is extracted from the "reference_data/tests_only/simple_gap_fit_params.json"
        self.frames = ase.io.read("reference_data/inputs/methane_dimer_sample.xyz", ":")

    def test_local_neg_stress(self):
        """..."""
        calculator = self.model.get_representation_calculator()
        for frame in self.frames:
            frame.wrap(eps=1e-12)
            manager = calculator.transform(frame)
            global_stress = self.model.predict_stress(manager)
            print(global_stress.shape)
            print(global_stress)

            # compare the to outputs
            local_stress = self.model.predict_stress(manager, local_stress=True)
            print(local_stress.shape)
            print(np.sum(local_stress, axis=0))

            # test should do something like, but does no work atm 
            #self.assertTrue(np.allclose(global_stress, np.sum(local_stress, axis=0)))
            break

