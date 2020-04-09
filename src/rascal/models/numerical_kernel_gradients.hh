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
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/models/kernels.hh"
#include "rascal/models/sparse_kernels.hh"

namespace rascal {



  template<class KernelImpl, class Calculator, class Manager, class SparsePoints, typename Derived>
  math::Matrix_t compute_displaced_kernel(const KernelImpl& kernel, const Calculator&  calculator, const Manager& manager, const SparsePoints& sparse_points, const size_t& i_atom, const Eigen::MatrixBase<Derived>& disp) {
    // get a copy of the atomic_structure object
    auto manager_root = extract_underlying_manager<0>(manager);
    json structure_copy = manager_root->get_atomic_structure();
    auto atomic_structure = structure_copy.template get<AtomicStructure<ThreeD>>();
    // displace the atom
    atomic_structure.displace_position(i_atom, disp);
    // make sure all atoms are in the unit cell
    atomic_structure.wrap();
    // update the neighborlist
    manager->update(atomic_structure);
    // compute the representation of the modified structure
    calculator.compute(manager);
    // compute the kernel of the modified structure
    math::Matrix_t KNM = kernel.compute(calculator, manager, sparse_points);
    // reset neighborlist to the original structure
    manager->update(structure_copy.template get<AtomicStructure<ThreeD>>());
    return KNM;
  }

  template<class KernelImpl, class Calculator, class Manager, class SparsePoints>
  math::Matrix_t compute_numerical_kernel_gradient(const KernelImpl& kernel,
                                                    const Calculator& calculator,
                                                    const Manager& manager,
                                                    const SparsePoints& sparse_points,
                                                    const double& h_disp) {
    auto disps = h_disp*Eigen::Matrix3d::Identity();
    size_t n_sparse_points{sparse_points.size()};
    math::Matrix_t KNM{manager.size() * ThreeD, n_sparse_points};
    for (size_t i_atom{0}; i_atom < manager.size(); ++i_atom) {
      for (int i_der{0}; i_der < ThreeD; ++i_der) {
        // use centered finite difference to estimate gradient
        math::Matrix_t KNM_p = compute_displaced_kernel(kernel, calculator, manager, sparse_points, i_atom, disps.row(i_der));
        math::Matrix_t KNM_m = compute_displaced_kernel(kernel, calculator, manager, sparse_points, i_atom, -disps.row(i_der));
        KNM.row(i_atom * ThreeD + i_der) =
            (KNM_p - KNM_m).rowwise().sum() / (2*h_disp);
      }
    }
  }

  template<class KernelImpl, class Calculator, class Managers, class SparsePoints>
  math::Matrix_t compute_numerical_kernel_gradients(const KernelImpl& kernel,
                                                    const Calculator& calculator,
                                                    const Managers& managers,
                                                    const SparsePoints& sparse_points,
                                                    double h_disp = 1e-5) {
    size_t n_centers{0};
    for (const auto& manager : managers) {
      n_centers += manager.size() * ThreeD;
    }
    size_t n_sparse_points{sparse_points.size()};
    math::Matrix_t KNM{n_centers, n_sparse_points};
    size_t i_centers{0};
    for (const auto& manager : managers) {
      KNM.block(i_centers, 0, manager.size() * ThreeD, n_sparse_points) =
        compute_numerical_kernel_gradient(kernel, calculator, manager, sparse_points, h_disp);
      i_centers += manager.size() * ThreeD;
    }
  }
}  // namespace rascal

#endif  // SRC_RASCAL_MODELS_NUMERICAL_KERNEL_GRADIENTS_HH_