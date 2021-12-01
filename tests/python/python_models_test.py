from rascal.representations import SphericalInvariants
from rascal.models import Kernel
from rascal.models.sparse_points import SparsePoints
from rascal.models.kernels import compute_numerical_kernel_gradients
from rascal.utils import from_dict, to_dict
from test_utils import load_json_frame, BoxList, Box, compute_relative_error
from ase.calculators.lj import LennardJones
import ase.io
import unittest
import numpy as np
import sys
import copy
import json
import pickle


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


class TestNumericalKernelGradient(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """
        self.calc = LennardJones()
        self.error_threshold = 1e-8
        # test file contains a reduced selection of
        self.kernel_input_filename = (
            "reference_data/tests_only/sparse_kernel_inputs.json"
        )
        # only some inputs are selected from the test inputs to reduce test time
        self.selected_test_indices = [0, 1]
        self.matrix_indices_in_voigt_notation = [
            (0, 0),
            (1, 1),
            (2, 2),
            (1, 2),
            (0, 2),
            (0, 1),
        ]

    def test_numerical_stress(self):
        """Tests if the numerical stress tensor on python site matches the one
        from ase library."""

        with open(self.kernel_input_filename, "r") as f:
            kernel_inputs = json.load(f)
        kernel_inputs = [kernel_inputs[i] for i in self.selected_test_indices]
        for kernel_input in kernel_inputs:
            structures_filename = kernel_input["filename"]
            frames = ase.io.read(
                structures_filename, ":" + str(kernel_input["n_structures"])
            )
            h_disp = kernel_input["h"]

            for frame in frames:
                frame.calc = self.calc

                ase_numerical_stress = self.calc.calculate_numerical_stress(
                    frame, d=h_disp, voigt=True
                )

                def compute_numerical_stress():
                    numerical_stress = np.zeros(6)
                    for i in range(6):
                        frame_displaced_plus = displace_strain_tensor(
                            copy.deepcopy(frame),
                            self.matrix_indices_in_voigt_notation[i][0],
                            self.matrix_indices_in_voigt_notation[i][1],
                            h_disp,
                        )
                        frame_displaced_minus = displace_strain_tensor(
                            copy.deepcopy(frame),
                            self.matrix_indices_in_voigt_notation[i][0],
                            self.matrix_indices_in_voigt_notation[i][1],
                            -h_disp,
                        )
                        e_plus = frame_displaced_plus.get_total_energy()
                        e_minus = frame_displaced_minus.get_total_energy()
                        numerical_stress[i] = (e_plus - e_minus) / (2 * h_disp)
                    return numerical_stress / frame.get_volume()

                numerical_stress = compute_numerical_stress()

                relative_error = compute_relative_error(
                    ase_numerical_stress, numerical_stress
                )
                absolute_error = np.abs(ase_numerical_stress - numerical_stress)
                passes_test = np.all(
                    np.logical_or(
                        relative_error < self.error_threshold,
                        absolute_error < self.error_threshold,
                    )
                )
                if not (passes_test):
                    print("structures_filename:", structures_filename)
                    print()
                    print("relative_error:", relative_error)
                    print()
                    print("ase_numerical_stress:", ase_numerical_stress)
                    print("numerical_stress:", numerical_stress)
                self.assertTrue(passes_test)

    def test_numerical_kernel_stress(self):
        """Tests if the numerical kernel stress on the python site matches the one
        on the cpp site."""

        with open(self.kernel_input_filename, "r") as f:
            kernel_inputs = json.load(f)

        kernel_inputs = [kernel_inputs[i] for i in self.selected_test_indices]
        for kernel_input in kernel_inputs:
            structures_filename = kernel_input["filename"]
            frames = ase.io.read(
                structures_filename, ":" + str(kernel_input["n_structures"])
            )
            h_disp = kernel_input["h"]

            selected_ids = kernel_input["selected_ids"]
            hypers = kernel_input["calculator"]
            # TODO(alex) the cutoff function is kind of hard coded
            #            a general function transformation c++ parameters to
            #            python would be more suitable here
            #            future work
            calculator = SphericalInvariants(**hypers)
            kernel = Kernel(calculator, kernel_type="Sparse", **kernel_input["kernel"])
            for j in range(len(frames)):
                # we do this frame by frame to be able to use the function
                # `displace_strain_tensor` as in the
                # `test_displace_strain_tensor` test
                frame = frames[j]
                selected_id = selected_ids[j]
                managers = calculator.transform([frame])
                sparse_points = SparsePoints(calculator)
                sparse_points.extend(managers, [selected_id])

                # the binded cpp function; the minus is because the function
                # returns the negative stress
                cpp_site_stress = -compute_numerical_kernel_gradients(
                    kernel, calculator, managers, sparse_points, h_disp, True
                )[-6:]

                def compute_numerical_kernel_gradient_on_python_site():
                    python_site_stress = np.zeros((6, len(selected_id)))
                    for i in range(6):
                        frame_displaced_plus = displace_strain_tensor(
                            copy.deepcopy(frame),
                            self.matrix_indices_in_voigt_notation[i][0],
                            self.matrix_indices_in_voigt_notation[i][1],
                            h_disp,
                        )
                        managers = calculator.transform([frame_displaced_plus])
                        kernel_plus = kernel(managers, sparse_points)

                        frame_displaced_minus = displace_strain_tensor(
                            copy.deepcopy(frame),
                            self.matrix_indices_in_voigt_notation[i][0],
                            self.matrix_indices_in_voigt_notation[i][1],
                            -h_disp,
                        )
                        managers = calculator.transform([frame_displaced_minus])
                        kernel_minus = kernel(managers, sparse_points)

                        python_site_stress[i] = np.sum(
                            (kernel_plus - kernel_minus) / (2 * h_disp), axis=0
                        )
                    return python_site_stress / frame.get_volume()

                python_site_stress = compute_numerical_kernel_gradient_on_python_site()

                relative_error = compute_relative_error(
                    python_site_stress, cpp_site_stress
                )
                absolute_error = np.abs(python_site_stress - cpp_site_stress)
                passes_test = np.all(
                    np.logical_or(
                        relative_error < self.error_threshold,
                        absolute_error < self.error_threshold,
                    )
                )
                if not (passes_test):
                    np.set_printoptions(suppress=True)
                    print("structures_filename:", structures_filename)
                    print("structure index:", j)
                    print()
                    print("relative_error:\n", relative_error)
                    print()
                    print("python_site_stress:\n", python_site_stress)
                    print("cpp_site_stress:\n", cpp_site_stress)

                self.assertTrue(passes_test)


class TestCosineKernel(unittest.TestCase):
    def setUp(self):
        """
        builds the test case. Test the order=1 structure manager implementation
        against a triclinic crystal.
        """

        fn = "reference_data/inputs/CaCrP2O7_mvc-11955_symmetrized.json"
        self.frame = load_json_frame(fn)

        self.hypers = dict(
            soap_type="PowerSpectrum",
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
            cosine_kernel = Kernel(rep, name="Cosine", target_type=target_type, zeta=2)
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

    def test_pickle(self):
        rep = SphericalInvariants(**self.hypers)
        cosine_kernel = Kernel(rep, name="Cosine", target_type="Structure", zeta=2)
        serialized = pickle.dumps(cosine_kernel)
        cosine_kernel_ = pickle.loads(serialized)
        self.assertTrue(to_dict(cosine_kernel) == to_dict(cosine_kernel_))

    def test_serialization(self):
        rep = SphericalInvariants(**self.hypers)

        for target_type in ["Atom", "Structure"]:
            cosine_kernel = Kernel(rep, name="Cosine", target_type=target_type, zeta=2)

            cosine_kernel_dict = to_dict(cosine_kernel)
            cosine_kernel_copy = from_dict(cosine_kernel_dict)
            cosine_kernel_copy_dict = to_dict(cosine_kernel_copy)

            self.assertTrue(cosine_kernel_dict == cosine_kernel_copy_dict)
