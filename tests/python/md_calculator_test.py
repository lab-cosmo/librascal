"""Test for the generic MD interface, chiefly used to work with i-PI"""

import unittest

import ase.io
import numpy as np

from rascal.models.genericmd import GenericMDCalculator


class TestGenericMD(unittest.TestCase):

    """Test the initialization and use of a GenericMDCalculator"""

    def setUp(self):
        self.calculator = GenericMDCalculator(
            "reference_data/tests_only/simple_gap_model.json",
            "reference_data/inputs/methane_dimer_sample.xyz",
            assume_pbc=False,
        )
        self.test_structure = ase.io.read(
            "reference_data/inputs/methane_dimer_sample.xyz"
        )
        self.test_positions = self.test_structure.get_positions()
        # Shift the positions to match the ASE/rascal cell convention
        self.test_positions += self.test_structure.get_cell()[np.diag_indices(3)] / 2.0
        self.test_cell = self.test_structure.get_cell()

    def test_calculate_multiple(self):
        """Test two subsequent calls -- once to initialize the manager,
        once to make sure it still works"""
        energy, forces, stress = self.calculator.calculate(
            self.test_positions, self.test_cell
        )
        self.assertFalse(np.isnan(energy))
        self.assertFalse(np.any(np.isnan(forces)))
        self.assertFalse(np.any(np.isnan(stress)))
        self.assertEqual(forces.shape, self.test_positions.shape)
        self.assertEqual(stress.shape, (3, 3))
        new_energy, new_forces, new_stress = self.calculator.calculate(
            self.test_positions, self.test_cell
        )
        self.assertTrue(np.allclose(energy, new_energy))
        self.assertTrue(np.allclose(forces, new_forces))
        self.assertTrue(np.allclose(stress, new_stress))

    def test_positions_update(self):
        energy, forces, stress = self.calculator.calculate(
            self.test_positions, self.test_cell
        )
        shift = np.array([5.0, 5.0, 5.0])
        new_energy, new_forces, new_stress = self.calculator.calculate(
            self.test_positions + shift, self.test_cell
        )
        # The descriptors (and model) should be invariant to uniform translation
        self.assertTrue(np.allclose(energy, new_energy))
        self.assertTrue(np.allclose(forces, new_forces))
        self.assertTrue(np.allclose(stress, new_stress))

    def test_cell_update(self):
        new_cell = np.eye(3) * 10.0  # Works for the first methane dimer structure, YMMV
        energy, forces, stress = self.calculator.calculate(
            self.test_positions, self.test_cell
        )
        new_energy, new_forces, new_stress = self.calculator.calculate(
            self.test_positions, new_cell
        )
        self.assertNotEqual(energy, new_energy)
        self.assertFalse(np.allclose(forces, new_forces))
        self.assertFalse(np.allclose(stress, new_stress))

    #TODO test non-periodic
    #TODO test wrong # atoms, cell shape
