# TMP line for debugging for alex
import sys
sys.path.insert(0, "/scratch/tisi/Programs/librascal/")
sys.path.insert(0, "/scratch/tisi/Programs/librascal/build")

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
        self.model = load_obj("../../reference_data/tests_only/simple_gap_model.json")
        # the file is extracted from the "reference_data/tests_only/simple_gap_fit_params.json"
        self.frames = ase.io.read("../../reference_data/inputs/methane_dimer_sample.xyz", ":")
        for frame in self.frames:
            frame.wrap(eps=1e-12)
        # librascal stacks all central local stress computation for each
        # structure into one dimension, so we need to keep track of the
        # staring index for each structure
        self.structure_starting_indices = np.cumsum([0] + [len(frame) for frame in self.frames])

    def test_summed_local_stress_is_equal_to_global_stress(self):
        """ Tests if the local stress contribution summed up match the global one"""
        calculator = self.model.get_representation_calculator()
        manager = calculator.transform(self.frames)
        global_stress = self.model.predict_stress(manager)
        print(type(self.model))

        print(self.model.kernel.kernel_type)
        print("type( self.model.X_tra._sparse_points",type(self.model.X_train._sparse_points))
        #breakpoint()
        #self.model.kernel._kernel.compute(manager.managers,self.model.X_train._sparse_points)
        self.model.kernel._kernel.compute_local(
            self.model.kernel._representation,
            manager.managers,
            self.model.X_train._sparse_points)

        # compare to local stress 
        local_stress = self.model.predict_stress(manager, local_stress=True)
        local_stress_matrix = local_stress.reshape(-1, 3, 3)
        # voight indicies of stress matrix as stored in the global stress in librascal
        #                 xx      yy     zz     yz     xz    xy
        voigt_indices = [(0,0), (1,1), (2,2), (1,2), (0,2), (0,1)]
        local_stress_voigt = np.concatenate([local_stress_matrix[:, i, j].reshape(-1,1) for (i,j) in voigt_indices], axis=1)
        for i in range(len(self.frames)):
            print(global_stress[i],np.sum(local_stress_voigt[self.structure_starting_indices[i]:self.structure_starting_indices[i+1]], axis=0))
            self.assertTrue(
                np.allclose(global_stress[i],
                np.sum(local_stress_voigt[self.structure_starting_indices[i]:self.structure_starting_indices[i+1]], axis=0)))

