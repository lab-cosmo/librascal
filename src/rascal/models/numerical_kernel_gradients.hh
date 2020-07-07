/**
 * @file   rascal/models/numerical_kernel_gradients.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   9 April 2020
 *
 * @brief Implementation of numerical evaluation of kernel gradients
 *
 * Copyright 2019 Felix Musil COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_MODELS_NUMERICAL_KERNEL_GRADIENTS_HH_
#define SRC_RASCAL_MODELS_NUMERICAL_KERNEL_GRADIENTS_HH_

#include "rascal/math/utils.hh"
#include "rascal/models/kernels.hh"
#include "rascal/models/sparse_kernels.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/structure_managers/make_structure_manager.hh"

namespace rascal {

  /**
   * Compute the sparse kernel associated to an atomic structure
   * displacing atom `i_atom` by `disp`.
   *
   * @see compute_numerical_kernel_gradients
   * @param i_atom index of the atom to displace
   * @param disp 3D vector used to displace the position of i_atom
   *
   */
  template <class KernelImpl, class Calculator, class Manager,
            class SparsePoints, typename Derived>
  math::Matrix_t compute_displaced_kernel(
      KernelImpl & kernel, Calculator & calculator, Manager & manager,
      const SparsePoints & sparse_points, const size_t & i_atom,
      const Eigen::MatrixBase<Derived> & disp) {
    // get a copy of the atomic_structure object
    auto manager_root = extract_underlying_manager<0>(manager);
    json structure_copy = manager_root->get_atomic_structure();
    auto atomic_structure =
        structure_copy.template get<AtomicStructure<ThreeD>>();
    atomic_structure.displace_position(i_atom, disp);
    // make sure all atoms are in the unit cell
    atomic_structure.wrap();
    manager->update(atomic_structure);
    calculator.compute(manager);
    std::vector<std::remove_const_t<Manager>> managers{};
    managers.emplace_back(manager);
    math::Matrix_t KNM = kernel.compute(calculator, managers, sparse_points);
    // reset neighborlist to the original structure
    manager->update(structure_copy.template get<AtomicStructure<ThreeD>>());
    return KNM;
  }

  /**
   * Compute finite-difference gradient of the kernel of a sparse GPR model
   * w.r.t. atomic positions for an atomic structure.
   *
   * @see compute_numerical_kernel_gradients
   */
  template <class KernelImpl, class Calculator, class Manager,
            class SparsePoints>
  math::Matrix_t compute_numerical_kernel_gradient(
      KernelImpl & kernel, Calculator & calculator, Manager & manager,
      const SparsePoints & sparse_points, const double & h_disp) {
    Eigen::Matrix3d disps = h_disp * Eigen::Matrix3d::Identity();
    size_t n_sparse_points{sparse_points.size()};
    math::Matrix_t KNM{manager->size() * ThreeD, n_sparse_points};
    KNM.setZero();
    for (size_t i_atom{0}; i_atom < manager->size(); ++i_atom) {
      for (int i_der{0}; i_der < ThreeD; ++i_der) {
        // use centered finite difference to estimate gradient
        math::Matrix_t KNM_p =
            compute_displaced_kernel(kernel, calculator, manager, sparse_points,
                                     i_atom, disps.row(i_der));
        math::Matrix_t KNM_m =
            compute_displaced_kernel(kernel, calculator, manager, sparse_points,
                                     i_atom, -disps.row(i_der));
        KNM.row(i_atom * ThreeD + i_der) =
            (KNM_p - KNM_m).colwise().sum() / (2 * h_disp);
      }
    }
    return KNM;
  }

  /**
   * Compute finite-difference gradient of the kernel of a sparse GPR model
   * w.r.t. atomic positions for a collection of atomic structures
   * using centered finite differences.
   *
   * @param kernel a sparse kernel
   * @param calculator a representation of the atomic neighborhood
   * @param managers a collection of structure managers
   * @param sparse_points basis points used in the sparse GPR model
   * @param h_disp displacement used for the centered finite difference
   *
   */
  template <class KernelImpl, class Calculator, class Managers,
            class SparsePoints>
  math::Matrix_t compute_numerical_kernel_gradients(
      KernelImpl & kernel, Calculator & calculator, Managers & managers,
      const SparsePoints & sparse_points, double h_disp = 1e-5) {
    size_t n_centers{0};
    for (const auto & manager : managers) {
      n_centers += manager->size() * ThreeD;
    }
    size_t n_sparse_points{sparse_points.size()};
    math::Matrix_t KNM{n_centers, n_sparse_points};
    KNM.setZero();
    size_t i_centers{0};
    for (const auto & manager : managers) {
      KNM.block(i_centers, 0, manager->size() * ThreeD, n_sparse_points) =
          compute_numerical_kernel_gradient(kernel, calculator, manager,
                                            sparse_points, h_disp);
      i_centers += manager->size() * ThreeD;
    }
    return KNM;
  }
}  // namespace rascal

#endif  // SRC_RASCAL_MODELS_NUMERICAL_KERNEL_GRADIENTS_HH_
