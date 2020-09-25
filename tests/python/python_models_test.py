from rascal.representations import SphericalInvariants
from rascal.models import Kernel, SparsePoints
from rascal.models.kernels import compute_numerical_kernel_gradients
from rascal.utils import from_dict, to_dict
from test_utils import load_json_frame, BoxList, Box
from ase.calculators.lj import LennardJones
import ase.io
import unittest
import numpy as np
import sys
import copy
import json


def displace_strain_tensor(frame, alpha_index, beta_index, h_disp):
    shift = np.eye(3)
    shift[alpha_index, beta_index] += h_disp
    original_cell = copy.deepcopy(frame.cell.array)
    displaced_cell = np.dot(original_cell, shift)
    frame.set_cell(displaced_cell)
    # adapt scaled positions of atoms
    M = np.linalg.inv(original_cell).dot(displaced_cell)
    frame.positions = np.dot(frame.positions, M)
    return frame

class TestNumericalKernel(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """
        self.verbose = True
        self.calc = LennardJones()
        self.treshold = 1e-8
        self.kernel_input_filename = 'reference_data/tests_only/numerical_gradient_stress_kernel_inputs.json'

        self.voigt_ids = [ (0,0), (1,1), (2,2), (1,2), (0,2), (0,1)]

    def test_displace_strain_tensor(self):
        """ Tests if displace_strain_tensor functions works properly, by
        checking if the stress tensor can be computed properly for Lennard Jones
        energy calculations with ase"""

        with open(self.kernel_input_filename, 'r') as f:
            kernel_inputs = json.load(f)

        for kernel_input in kernel_inputs:
            structures_filename = kernel_input["filename"]
            frames = ase.io.read(structures_filename,
                    ":"+str(kernel_input["n_structures"]))
            h_disp = kernel_input["h"]

            for frame in frames:
                frame.calc = self.calc
                numerical_stress = np.zeros(6)
                for i in range(6):
                    frame_displaced_plus  = displace_strain_tensor(
                            copy.deepcopy(frame), self.voigt_ids[i][0], 
                            self.voigt_ids[i][1], h_disp)
                    frame_displaced_minus = displace_strain_tensor(
                            copy.deepcopy(frame), self.voigt_ids[i][0],
                            self.voigt_ids[i][1], -h_disp)
                    e_plus  = frame_displaced_plus.get_total_energy()
                    e_minus = frame_displaced_minus.get_total_energy()
                    numerical_stress[i] = (e_plus - e_minus) / (2*h_disp)
                numerical_stress /= frame.get_volume()
                analytical_stress = frame.get_stress()

                absolute_error = np.abs(numerical_stress-analytical_stress)
                relative_error = np.abs((numerical_stress-analytical_stress)/analytical_stress)
                passes_test = ( np.all(absolute_error < self.treshold) or
                        np.all(relative_error< self.treshold) )
                if (self.verbose and not(passes_test)):
                    print("structures_filename:", structures_filename)
                    print("structure index:", frame.index)
                    print()
                    print("numerical_stress:",numerical_stress)
                    print("analytical_stress:",analytical_stress)
                    print()
                    print("absolute_error:", absolute_error)
                    print("relative_error:", relative_error)

                self.assertTrue(passes_test)

    def test_numerical_stress_gradient(self):
        """ Tests if numerical stress gradient on the c++ site agrees with the
        one on the python site. The python site numerical stress gradient is
        verified with Lennard Jones energy calculations by the
        `test_displace_strain_tensor` test. """

        with open(self.kernel_input_filename, 'r') as f:
            kernel_inputs = json.load(f)

        for kernel_input in kernel_inputs:
            structures_filename = kernel_input["filename"]
            frames = ase.io.read(structures_filename,
                    ":"+str(kernel_input["n_structures"]))
            h_disp = kernel_input["h"]

            selected_ids = kernel_input["selected_ids"]
            hypers = kernel_input["calculator"]
            # TODO(alex) the cutoff function is a bit hard coded
            #            a general function transformation c++ parameters to
            #            python would be more suitable here
            calculator = SphericalInvariants(
                    soap_type=hypers["soap_type"],
                    radial_basis=hypers["radial_contribution"]["type"],
                    max_radial=hypers["max_radial"],
                    max_angular=hypers["max_angular"],
                    cutoff_function_type=hypers["cutoff_function"]["type"],
                    interaction_cutoff=hypers["cutoff_function"]["cutoff"]["value"],
                    cutoff_smooth_width=hypers["cutoff_function"]["smooth_width"]["value"],
                    gaussian_sigma_type=hypers["gaussian_density"]["type"],
                    gaussian_sigma_constant=hypers["gaussian_density"]["gaussian_sigma"]["value"],
                    compute_gradients=hypers["compute_gradients"],
                    normalize=hypers["normalize"])
            kernel = Kernel(calculator, kernel_type='Sparse',
                    **kernel_input["kernel"])
            for j in range(len(frames)):
                # we do this frame by frame to be able to use the function
                # `displace_strain_tensor` as in the
                # `test_displace_strain_tensor` test
                frame = frames[j]
                selected_id = selected_ids[j]
                managers = calculator.transform([frame])
                sparse_points = SparsePoints(calculator)
                sparse_points.extend(managers, [selected_id])

                cpp_site_stress = compute_numerical_kernel_gradients(kernel,
                        calculator, managers, sparse_points, h_disp, True)[-6:]
                python_site_stress = np.zeros( (6, len(selected_id)) )

                for i in range(6):
                    frame_displaced_plus = displace_strain_tensor(
                            copy.deepcopy(frame), self.voigt_ids[i][0],
                            self.voigt_ids[i][1], h_disp)
                    managers = calculator.transform([frame_displaced_plus])
                    kernel_plus = kernel(managers, sparse_points)

                    frame_displaced_minus = displace_strain_tensor(
                            copy.deepcopy(frame), self.voigt_ids[i][0],
                            self.voigt_ids[i][1], -h_disp)
                    managers = calculator.transform([frame_displaced_minus])
                    kernel_minus = kernel(managers, sparse_points)

                    python_site_stress[i] = np.sum( (kernel_plus - kernel_minus)
                            / (2*h_disp), axis=0 )

                absolute_error = np.abs(python_site_stress-cpp_site_stress)
                relative_error = (python_site_stress-cpp_site_stress)
                relative_error[cpp_site_stress != 0] /= cpp_site_stress[cpp_site_stress != 0]
                relative_error = np.abs(relative_error)
                # change to logic or elementwise
                passes_test = ( np.all(absolute_error < self.treshold) or
                        np.all(relative_error< self.treshold) )
                if (self.verbose and not(passes_test)):
                    np.set_printoptions(suppress=True)
                    print("structures_filename:", structures_filename)
                    print("structure index:", j)
                    print()
                    print("python_site_stress:\n",python_site_stress)
                    print("cpp_site_stress:\n",cpp_site_stress)
                    print()
                    print("absolute_error:\n", absolute_error)
                    print("relative_error:\n", relative_error)

                self.assertTrue(passes_test)

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
