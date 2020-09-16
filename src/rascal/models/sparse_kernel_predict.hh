/**
 * @file   rascal/models/sparse_kernel_predict.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   15 Sept 2020
 *
 * @brief Implementation of similarity kernels for sparse methods
 *
 * Copyright 2020 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_MODELS_SPARSE_KERNEL_PREDICT_HH_
#define SRC_RASCAL_MODELS_SPARSE_KERNEL_PREDICT_HH_

#include "rascal/math/utils.hh"
#include "rascal/models/kernels.hh"
#include "rascal/models/sparse_kernels.hh"
#include "rascal/structure_managers/property_block_sparse.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils/json_io.hh"

namespace rascal {

  /**
   * Compute the partial gradients of a structure w.r.t atomic positions
   * using the GAP model (see
   * @ref SparseKernelImpl<internal::SparseKernelType::GAP>::compute_derivative)
   *
   * The point of this function is to provide a faster prediction routine
   * compared to computing the kernel element and then multiplying them with
   * the weights of the model.
   *
   * The partial forces are attached to the input manager in a Property of
   * Order 2.
   * partial forces [N_{neighbor}, 3]
   *
   * @tparam StructureManagers should be an iterable over shared pointer
   *          of structure managers like ManagerCollection
   * @param sparse_points a SparsePoints* class
   * @param managers a ManagerCollection or similar collection of
   * structure managers
   * @param zeta exponent of the GAP kernel
   * @param weights of the GAP model
   * @param representation_name name used to get the representation in the
   * managers
   * @param representation_grad_name name used to get representation gradient
   * the in the managers
   * @param fij_name name used to get/register the partial gradients in the
   * managers
   * @return
   */
  template <class Property_t, class PropertyGradient_t, class StructureManager,
            class SparsePoints>
  void
  compute_partial_gradients_gap(StructureManager & manager,
                                SparsePoints & sparse_points,
                                math::Vector_t & weights, const size_t zeta,
                                const std::string & representation_name,
                                const std::string & representation_grad_name,
                                const std::string & fij_name) {
    using Manager_t = typename StructureManager::element_type;
    using Keys_t = typename SparsePoints::Keys_t;
    using Key_t = typename SparsePoints::Key_t;

    auto && prop{
        *manager->template get_property<Property_t>(representation_name, true)};
    auto && prop_grad{*manager->template get_property<PropertyGradient_t>(
        representation_grad_name, true)};

    const int inner_size{prop.get_nb_comp()};

    bool do_block_by_key_dot{false};
    if (prop_grad.are_keys_uniform()) {
      do_block_by_key_dot = true;
    }

    std::set<int> unique_species{};
    for (auto center : manager) {
      unique_species.insert(center.get_atom_type());
    }

    // find shared central atom species
    std::set<int> species_intersect{
        internal::set_intersection(unique_species, sparse_points.species())};

    // find offsets alongs the sparse points direction
    std::map<int, int> offsets{sparse_points.get_offsets()};
    Keys_t rep_keys{prop_grad.get_keys()};
    std::map<int, Keys_t> keys_intersect{};
    for (const int & sp : species_intersect) {
      keys_intersect[sp] =
          internal::set_intersection(rep_keys, sparse_points.keys_sp.at(sp));
    }

    // attach partial forces array to manager
    auto && f_ij{
        *manager
             ->template get_property<Property<double, 2, Manager_t, 1, ThreeD>>(
                 fij_name, true, true)};
    // don't recompute the partial forces if already up to date
    if (f_ij.is_updated()) {
      return;
    }

    f_ij.resize();
    f_ij.setZero();

    if (species_intersect.empty()) {
      return;
    }

    math::Vector_t weights_scaled(weights.size());
    size_t i_row{0};
    for (auto center : manager) {
      const int a_sp{center.get_atom_type()};
      // \alpha_n^{scaled} = \alpha_n * [z* (X_j \dot T_n)^{z-1}]
      weights_scaled =
          (weights.array() *
           (zeta *
            internal::pow_zeta(sparse_points.dot(a_sp, prop[center]), zeta - 1))
               .transpose()
               .array())
              .matrix();
      // \sum_n \alpha_n^{scaled} T_n
      SparsePoints sparse_point_scaled{sparse_points.dot(a_sp, weights_scaled)};
      const size_t n_neigh{center.pairs_with_self_pair().size()};
      // sparse_point_scaled \dot dX_i/dr_j
      if (do_block_by_key_dot) {
        auto rep_grads = prop_grad.get_raw_data_view();
        const auto & values_by_sp = sparse_point_scaled.values.at(a_sp);

        for (const Key_t & key : keys_intersect.at(a_sp)) {
          const auto & values_by_sp_key = values_by_sp.at(key);
          auto spts = Eigen::Map<const math::Vector_t>(
              values_by_sp_key.data(), static_cast<Eigen::Index>(inner_size));

          math::Vector_t fij_block(n_neigh);

          int col_st{prop_grad.get_gradient_col_by_key(key)};
          for (int i_der{0}; i_der < ThreeD; i_der++) {
            fij_block = (rep_grads.block(i_row, col_st + i_der * inner_size,
                                         n_neigh, inner_size) *
                         spts.transpose())
                            .transpose();

            int i_row_{0};
            for (auto neigh : center.pairs_with_self_pair()) {
              f_ij[neigh](i_der) += fij_block(i_row_);
              i_row_++;
            }  // neigh
          }    // i_der
        }      // key
        i_row += n_neigh;
      } else {
        for (auto neigh : center.pairs_with_self_pair()) {
          f_ij[neigh] =
              sparse_point_scaled.dot_derivative(a_sp, prop_grad[neigh]);
        }
      }
    }  // center
  }

  /**
   * Compute the gradients of a structure w.r.t atomic positions
   * using a sparse GPR model. Only GAP model is implemented at the moment (see
   * @ref SparseKernelImpl<internal::SparseKernelType::GAP>::compute_derivative)
   *
   * The point of this function is to provide a faster prediction routine
   * compared to computing the kernel elements and then multiplying them with
   * the weights of the model.
   *
   * The forces are attached to the input manager in a Property of
   * Order 1 with the name
   * 'representation_grad_name+" forces; weight_hash:"+weight_hash'
   * forces [N_{atoms}, 3]
   *
   * @tparam StructureManagers should be an iterable over shared pointer
   *          of structure managers like ManagerCollection
   * @param sparse_points a SparsePoints* class
   * @param managers a ManagerCollection or similar collection of
   * structure managers
   * @param weights regression weights of the sparse GPR model
   * @return name used to register the forces in the managers
   */
  template <class Calculator, class StructureManagers, class SparsePoints>
  std::string compute_sparse_kernel_forces(const Calculator & calculator,
                                           SparseKernel & kernel,
                                           StructureManagers & managers,
                                           SparsePoints & sparse_points,
                                           math::Vector_t & weights) {
    using Manager_t = typename StructureManagers::Manager_t;
    using Property_t = typename Calculator::template Property_t<Manager_t>;
    using PropertyGradient_t =
        typename Calculator::template PropertyGradient_t<Manager_t>;
    auto && representation_name{calculator.get_name()};
    const auto representation_grad_name{calculator.get_gradient_name()};
    size_t n_centers{0};
    // find the total number of forces
    for (const auto & manager : managers) {
      n_centers += manager->size();
    }

    internal::Hash<math::Vector_t, double> hasher{};
    auto kernel_type_str = kernel.parameters.at("name").get<std::string>();
    std::string weight_hash = std::to_string(hasher(weights));
    std::string fij_name = representation_grad_name +
                           std::string(" partial gradients; weight_hash:") +
                           weight_hash;
    std::string force_name = representation_grad_name +
                             std::string(" forces; weight_hash:") + weight_hash;

    for (const auto & manager : managers) {
      if (kernel_type_str == "GAP") {
        const auto zeta = kernel.parameters.at("zeta").get<size_t>();

        compute_partial_gradients_gap<Property_t, PropertyGradient_t>(
            manager, sparse_points, weights, zeta, representation_name,
            representation_grad_name, fij_name);
      }

      auto && forces{*manager->template get_property<
          Property<double, 1, Manager_t, 1, ThreeD>>(force_name, true, true,
                                                     true)};
      forces.resize();
      forces.setZero();

      auto && f_ij{*manager->template get_property<
          Property<double, 2, Manager_t, 1, ThreeD>>(fij_name, true)};
      for (auto center : manager) {
        // accumulate partial gradients onto forces
        for (auto neigh : center.pairs_with_self_pair()) {
          forces[neigh.get_atom_j()] -= f_ij[neigh];
        }
      }
    }  // center
    return force_name;
  }
}  // namespace rascal
#endif  // SRC_RASCAL_MODELS_SPARSE_KERNEL_PREDICT_HH_
