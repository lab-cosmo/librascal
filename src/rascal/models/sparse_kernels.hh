/**
 * @file   rascal/models/sparse_kernels.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   8 Jan 2020
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

#ifndef SRC_RASCAL_MODELS_SPARSE_KERNELS_HH_
#define SRC_RASCAL_MODELS_SPARSE_KERNELS_HH_

#include "rascal/math/utils.hh"
#include "rascal/models/kernels.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils/json_io.hh"

namespace rascal {

  namespace internal {
    enum class SparseKernelType { GAP };

    template <internal::SparseKernelType Type>
    struct SparseKernelImpl {};

    /**
     * Implementation of the sparse kernel used in GAP. The kernel between
     * local environments \f$X_i^a\f$ and \f$X_j^b\f$ with central atom types
     * \f$a\f$ and \f$b\f$ is:
     *
     * @f[
     *    k(X_i^a, X_j^b) = \delta_{ab} k(X_i, X_j),
     * @f]
     *
     * where \f$k(X_i, X_j)\f$ is the cosine kernel. When building a model for
     * properties associated with the structure we assume that the training
     * will be done on the property itself (not divided by the number of
     * atoms).
     *
     * It is particularly designed to build the \f$K_{MM}\f$ and \f$K_{NM}\f$
     * matrices needed by the sparse kernel formulation where the kernel matrix
     * is given by:
     *
     * @f[
     *   K = K_{MM} + K_{MN} \Lambda^{-2} K_{NM}
     * @f]
     *
     * with \f$\Lambda\f$ the regularization matrix.
     */
    template <>
    struct SparseKernelImpl<internal::SparseKernelType::GAP> : KernelImplBase {
      using Hypers_t = typename KernelImplBase::Hypers_t;

      //! exponent of the cosine kernel
      size_t zeta{1};
      constexpr static const int SpatialDims{3};

      SparseKernelImpl() = default;

      explicit SparseKernelImpl(const Hypers_t & hypers) : KernelImplBase{} {
        this->set_hyperparmeters(hypers);
      }

      void set_hyperparmeters(const Hypers_t & hypers) {
        if (hypers.count("zeta") == 1) {
          zeta = hypers["zeta"].get<size_t>();
        } else {
          throw std::runtime_error(
              R"(zeta should be specified for the GAP kernel)");
        }
      }

      /**
       * Compute the kernel between a set of structure(s) and a set of pseudo
       * points, per structure.
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *          of structure managers like ManagerCollection
       * @param sparse_points a SparsePoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation data
       * has been registered in the elements of managers and managers_b
       * @return kernel matrix
       */
      template <
          class Property_t, internal::TargetType Type,
          std::enable_if_t<Type == internal::TargetType::Structure, int> = 0,
          class StructureManagers, class SparsePoints>
      math::Matrix_t compute(StructureManagers & managers,
                             SparsePoints & sparse_points,
                             const std::string & representation_name) {
        math::Matrix_t KNM(managers.size(), sparse_points.size());
        KNM.setZero();
        size_t ii_A{0};
        for (auto & manager : managers) {
          auto && propA{*manager->template get_property<Property_t>(
              representation_name, true)};
          for (auto center : manager) {
            int sp = center.get_atom_type();
            // only the pseudo points of species sp contribute
            KNM.row(ii_A) +=
                pow_zeta(sparse_points.dot(sp, propA[center]), this->zeta)
                    .transpose();
          }
          ++ii_A;
        }
        return KNM;
      }

      /**
       * Compute the kernel between a set of pseudo points.
       *
       * @tparam SparsePoints should be a set a set of SparsePoints
       * @param sparse_points a SparsePoints* class
       * @return kernel matrix MxM
       */
      template <class SparsePoints>
      math::Matrix_t compute(SparsePoints & sparse_points) {
        math::Matrix_t KMM(sparse_points.size(), sparse_points.size());
        KMM.setZero();
        int start{0};
        // loop over the species
        for (const int & sp : sparse_points.species()) {
          // only the pseudo points of the same species contribute
          auto block_size = sparse_points.size_by_species(sp);
          KMM.block(start, start, block_size, block_size) =
              pow_zeta(sparse_points.self_dot(sp), this->zeta);
          start += block_size;
        }
        return KMM;
      }

      /**
       * Compute the kernel between a set of structure(s) and a set of pseudo
       * points.
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *          of structure managers like ManagerCollection
       * @param sparse_points a SparsePoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation data
       * has been registered in the elements of managers and managers_b
       * @return kernel matrix
       */
      template <class Property_t, internal::TargetType Type,
                std::enable_if_t<Type == internal::TargetType::Atom, int> = 0,
                class StructureManagers, class SparsePoints>
      math::Matrix_t compute(const StructureManagers & managers,
                             const SparsePoints & sparse_points,
                             const std::string & representation_name) {
        size_t n_centersA{0};
        for (const auto & manager : managers) {
          n_centersA += manager->size();
        }
        size_t nb_sparse_points{sparse_points.size()};
        math::Matrix_t KNM(n_centersA, nb_sparse_points);
        size_t ii_A{0};
        for (auto & manager : managers) {
          auto && propA{*manager->template get_property<Property_t>(
              representation_name, true)};
          for (auto center : manager) {
            int sp = center.get_atom_type();
            KNM.row(ii_A) =
                pow_zeta(sparse_points.dot(sp, propA[center]), this->zeta)
                    .transpose();
            ii_A++;
          }
        }
        return KNM;
      }

      /**
       * This documentation contains specific information for the GAP
       * implementation See
       * [SparseKernel::compute_derivative](./cpp.html#_CPPv4I000EN6rascal12SparseKernel18compute_derivativeEN4math8Matrix_tERK10CalculatorRK17StructureManagersRK12SparsePointsKb)
       * for general kernel gradient derivation.
       *
       * The storage of the kernel forces would scale with indices \f$m\f$ x
       * \f$i\f$ x \f$j\f$, but we don't need to store individual neighbour
       * contributions and can therefore contract over all neighbors of a
       * particular center. This routine returns the already contracted kernel
       * derivative for center \f$X_i\f$
       *
       * @f[
       *      \sum_{j \in A_i} \vec{\nabla}_i k(X_j, T_m)
       * @f]
       *
       * so the resulting matrix has the shape of the number of centers (index
       * \f$i\f$) times 3 (three Cartesian dimensions) times the number of
       * sparse points (index \f$m\f$). The number of centers and Cartesian
       * dimensions are merged to one dimension, such that the x, y, z spatial
       * dimensions appear in successive order in the row dimension.
       *
       * The derivative of the kernel of environment \f$X_j\f$ with respect to
       * the position of \f$\mathbf{r}_i\f$ atom \f$i\f$ and the sparse point
       * \f$T_m\f$ is obtained by the chain rule
       *
       * @f[
       *      \vec{\nabla}_i k(X_j, T_m) = \frac{ \partial k(X_j, T_m) }{
       * \partial \mathbf{r}_i } = \frac{ \partial k(X_j, T_m) }{ \partial X_{j}
       * } \cdot \frac{ \partial X_j }{ \partial \mathbf{r}_i },
       * @f]
       *
       * where \f$\cdot\f$ corresponds to the contraction of the chain rule.
       * For the derivative of the GAP kernel we have:
       *
       * @f[
       *      \frac{ \partial k(X_{j}, T_m) }{ \partial X_{j}} =
       *        \zeta (\sum_n X_{jn} T_{mn})^{\zeta-1}  T_m
       * @f]
       *
       * Besides that, the negative stress tensor in the kernel formulation is
       * computed by also using the chain rule and contracting over the number
       * of atoms dimension as
       *
       * @f[
       *      \frac{ \partial k(X_j, T_m) }{ \partial\eta_{\alpha\beta} } =
       *          \sum_{i\in A} \sum_{j\in i} [\mathbf{r}_{ji}
       *                  \otimes \vec{\nabla}_j k(X_i, T_m)]_{\alpha\beta}.
       * @f]
       *
       * See https://aip.scitation.org/doi/full/10.1063/1.2214719 for more
       * details with a short range potential.
       * The elements are stored at the end of the returned matrix using the
       * Voigt format (xx, yy, zz, zy, zx, yx) in the order of managers.
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *         of structure managers like ManagerCollection
       * @param sparse_points a SparsePoints* class
       * @param managers a ManagerCollection or similar collection of
       *        structure managers
       * @param representation_name name under which the representation
       *        gradient data has been registered in the elements of managers
       * @param compute_neg_stress the computed negative stresses are appended
       *        at the end of the kernel matrix
       * @return kernel matrix
       */
      template <class Property_t, class PropertyGradient_t,
                internal::TargetType Type,
                std::enable_if_t<Type == internal::TargetType::Atom, int> = 0,
                class StructureManagers, class SparsePoints>
      math::Matrix_t
      compute_derivative(StructureManagers & managers,
                         SparsePoints & sparse_points,
                         const std::string & representation_name,
                         const std::string & representation_grad_name,
                         const bool compute_neg_stress) {
        using Manager_t = typename StructureManagers::Manager_t;
        using Keys_t = typename SparsePoints::Keys_t;
        using Key_t = typename SparsePoints::Key_t;
        // the nb of rows of the kernel matrix consist of:
        // - 3*nb_centers rows for each center for each spatial_dim
        // - 6 rows at the end for the stress tensor in voigt notation if
        //   `compute_neg_stress` is true
        size_t nb_kernel_gradient_rows{0};
        // find the total number of rows the matrix block should have
        for (const auto & manager : managers) {
          nb_kernel_gradient_rows += manager->size() * SpatialDims;
        }
        // the stress terms are stored at the bottom of the KNM in a block
        size_t row_idx_stress{nb_kernel_gradient_rows};
        if (compute_neg_stress) {
          nb_kernel_gradient_rows += 2 * SpatialDims * managers.size();
        }
        // Voigt order is xx, yy, zz, yz, xz, xy. To compute xx, yy, zz
        // and yz, xz, xy in one loop over the three spatial dimensions
        // dK/dr_{x,y,z}, we fill the off-diagonals yz, xz, xy by computing
        // xz, yx, zy exploiting the symmetry of the stress tensor
        // thus yx=xy and zx=xz
        // array accessed by voigt_idx and returns spatial_dim_idx
        const std::array<std::array<int, 2>, SpatialDims>
            voigt_id_to_spatial_dim = {{// voigt_idx,  spatial_dim_idx
                                        {{4, 2}},    //    xz,            z
                                        {{5, 0}},    //    xy,            x
                                        {{3, 1}}}};  //    yz,            y

        const size_t nb_sparse_points{sparse_points.size()};
        math::Matrix_t KNM_der(nb_kernel_gradient_rows, nb_sparse_points);
        KNM_der.setZero();
        const size_t zeta{this->zeta};
        size_t idx_center{0};
        // loop over the structures
        for (auto & manager : managers) {
          auto && prop_repr{*manager->template get_property<Property_t>(
              representation_name, true)};
          auto && prop_repr_grad{
              *manager->template get_property<PropertyGradient_t>(
                  representation_grad_name, true)};
          // dk/dr_i this is col major hence dim order
          Property<double, 1, Manager_t, Eigen::Dynamic, SpatialDims> dkdr{
              *manager, "dkdr", true};
          dkdr.set_nb_row(nb_sparse_points);
          dkdr.resize();
          dkdr.setZero();

          // dk/dX without sparse point factor T
          Property<double, 1, Manager_t, Eigen::Dynamic, 1> dkdX_missing_T{
              *manager, "dkdX without sparse point factor T", true};
          dkdX_missing_T.set_nb_row(nb_sparse_points);
          dkdX_missing_T.resize();
          if (zeta > 1) {
            dkdX_missing_T.setZero();
            // zeta * (X * T)**(zeta-1)
            for (auto center : manager) {
              int a_species{center.get_atom_type()};
              dkdX_missing_T[center] =
                  zeta *
                  pow_zeta(sparse_points.dot(a_species, prop_repr[center]),
                           zeta - 1);
            }
          }
          //
          const int block_size{prop_repr.get_nb_comp()};

          bool do_block_by_key_dot{false};
          if (prop_repr_grad.are_keys_uniform()) {
            do_block_by_key_dot = true;
          }

          std::set<int> unique_species{};
          for (auto center : manager) {
            unique_species.insert(center.get_atom_type());
          }

          // find shared central atom species
          std::set<int> species_intersect{internal::set_intersection(
              unique_species, sparse_points.species())};

          if (species_intersect.size() == 0) {
            return KNM_der;
          }

          // find offsets along the sparse points spatial_dim
          std::map<int, int> offsets{sparse_points.get_offsets()};
          Keys_t repr_keys{prop_repr_grad.get_keys()};
          std::map<int, Keys_t> keys_intersect{};
          for (const int & species : species_intersect) {
            keys_intersect[species] = internal::set_intersection(
                repr_keys, sparse_points.keys_sp.at(species));
          }
          // compute dX/dr * T * k_{zeta-1} * zeta
          if (do_block_by_key_dot) {
            size_t idx_row{0};
            auto repr_grads = prop_repr_grad.get_raw_data_view();
            for (auto center : manager) {
              int a_species{center.get_atom_type()};
              Eigen::Vector3d r_i = center.get_position();
              if (species_intersect.count(a_species) == 0) {
                continue;
              }
              auto dkdX_i_missing_T = dkdX_missing_T[center];
              const int offset = offsets.at(a_species);
              const auto & values_by_sp = sparse_points.values.at(a_species);
              const auto & indices_by_sp = sparse_points.indices.at(a_species);
              const size_t nb_rows{center.pairs_with_self_pair().size()};
              for (const Key_t & key : keys_intersect.at(a_species)) {
                const auto & indices_by_sp_key = indices_by_sp.at(key);
                const auto & values_by_sp_key = values_by_sp.at(key);
                auto sparse_points_block = Eigen::Map<const math::Matrix_t>(
                    values_by_sp_key.data(),
                    static_cast<Eigen::Index>(indices_by_sp_key.size()),
                    static_cast<Eigen::Index>(block_size));
                assert(indices_by_sp_key.size() * block_size ==
                       values_by_sp_key.size());
                math::Matrix_t KNM_der_block(nb_rows, indices_by_sp_key.size());
                // copy subset of k_{zeta-1} * zeta that matches a_species and
                // key
                math::Vector_t dkdX_i_missing_T_a_species_block{
                    indices_by_sp_key.size()};
                if (zeta > 1) {
                  for (size_t idx_col{0}; idx_col < indices_by_sp_key.size();
                       idx_col++) {
                    dkdX_i_missing_T_a_species_block(idx_col) =
                        dkdX_i_missing_T(offset + indices_by_sp_key[idx_col]);
                  }  // sparse points
                }
                int block_start_col_idx{
                    prop_repr_grad.get_gradient_col_by_key(key)};
                for (int idx_spatial_dim{0}; idx_spatial_dim < SpatialDims;
                     idx_spatial_dim++) {
                  // dX/dr * T
                  KNM_der_block =
                      repr_grads.block(idx_row,
                                       block_start_col_idx +
                                           idx_spatial_dim * block_size,
                                       nb_rows, block_size) *
                      sparse_points_block.transpose();
                  // * zeta * (X * T)**(zeta-1)
                  if (zeta > 1) {
                    KNM_der_block *=
                        dkdX_i_missing_T_a_species_block.asDiagonal();
                  }
                  int idx_neigh{0};
                  for (auto neigh : center.pairs_with_self_pair()) {
                    auto dkdr_ji{dkdr[neigh.get_atom_j()]};
                    for (int idx_col{0}; idx_col < KNM_der_block.cols();
                         idx_col++) {
                      dkdr_ji(offset + indices_by_sp_key[idx_col],
                              idx_spatial_dim) +=
                          KNM_der_block(idx_neigh, idx_col);
                    }  // sparse points
                    idx_neigh++;
                  }  // neigh
                  if (compute_neg_stress) {
                    const auto & voigt =
                        voigt_id_to_spatial_dim[idx_spatial_dim];
                    idx_neigh = 0;
                    for (auto neigh : center.pairs_with_self_pair()) {
                      Eigen::Vector3d r_ji = r_i - neigh.get_position();
                      for (int idx_col{0}; idx_col < KNM_der_block.cols();
                           idx_col++) {
                        // computes in order xx, yy, zz
                        KNM_der(row_idx_stress + idx_spatial_dim,
                                offset + indices_by_sp_key[idx_col]) +=
                            r_ji(idx_spatial_dim) *
                            KNM_der_block(idx_neigh, idx_col);
                        // computes in order xz, xy, yz
                        KNM_der(row_idx_stress + voigt[0],
                                offset + indices_by_sp_key[idx_col]) +=
                            r_ji(voigt[1]) * KNM_der_block(idx_neigh, idx_col);
                      }  // sparse points
                      idx_neigh++;
                    }  // neigh
                  }    // if compute_neg_stress
                }      // idx_spatial_dim
              }        // key
              idx_row += nb_rows;
            }  // center
          } else {
            for (auto center : manager) {
              int a_species{center.get_atom_type()};
              Eigen::Vector3d r_i = center.get_position();
              auto dkdX_i_missing_T = dkdX_missing_T[center];
              for (auto neigh : center.pairs_with_self_pair()) {
                // T * dX/dr
                auto T_times_dXdr = sparse_points.dot_derivative(
                    a_species, prop_repr_grad[neigh]);
                // * zeta * (X * T)**(zeta-1)
                if (zeta > 1) {
                  T_times_dXdr.transpose() *= dkdX_i_missing_T.asDiagonal();
                }
                dkdr[neigh.get_atom_j()] += T_times_dXdr;
                if (compute_stress) {
                  Eigen::Vector3d r_ji = r_i - neigh.get_position();
                  for (int i_der{0}; i_der < SpatialDims; i_der++) {
                    const auto & voigt = voigt_id_to_spatial_dim[i_der];
                    // computes in order xx, yy, zz
                    KNM_der.row(row_idx_stress + i_der) +=
                        r_ji(i_der) * T_times_dXdr.transpose().row(i_der);
                    // computes in order xz, xy, yz
                    KNM_der.row(row_idx_stress + voigt[0]) +=
                        r_ji(voigt[1]) * T_times_dXdr.transpose().row(i_der);
                  }
                }
              }  // neigh
            }    // center
          }      // if do_block_by_key_dot

          // copy the data to the kernel matrix
          for (auto center : manager) {
            KNM_der.block(idx_center, 0, SpatialDims, nb_sparse_points) =
                dkdr[center].transpose();
            idx_center += SpatialDims;
          }
          if (compute_neg_stress) {
            // TODO(alex) when we established how we deal with
            // `get_atomic_structure` method for other root managers
            // replace this part
            auto manager_root = extract_underlying_manager<0>(manager);
            json structure_copy = manager_root->get_atomic_structure();
            auto atomic_structure =
                structure_copy.template get<AtomicStructure<SpatialDims>>();
            KNM_der.block(row_idx_stress, 0, 6, KNM_der.cols()) /=
                atomic_structure.get_volume();
            row_idx_stress += SpatialDims * 2;
          }
        }  // managers
        return KNM_der;
      }
    };
  }  // namespace internal

  template <internal::SparseKernelType Type, class Hypers>
  std::shared_ptr<internal::KernelImplBase>
  make_sparse_kernel_impl(const Hypers & hypers) {
    return std::static_pointer_cast<internal::KernelImplBase>(
        std::make_shared<internal::SparseKernelImpl<Type>>(hypers));
  }

  template <internal::SparseKernelType Type>
  std::shared_ptr<internal::SparseKernelImpl<Type>> downcast_sparse_kernel_impl(
      std::shared_ptr<internal::KernelImplBase> & kernel_impl) {
    return std::static_pointer_cast<internal::SparseKernelImpl<Type>>(
        kernel_impl);
  }

  class SparseKernel {
   public:
    using Hypers_t = typename internal::KernelImplBase::Hypers_t;

    explicit SparseKernel(const Hypers_t & hypers) {
      using internal::SparseKernelType;
      using internal::TargetType;

      this->parameters = hypers;

      if (hypers.count("target_type") == 1) {
        if (hypers["target_type"] == "Structure") {
          this->target_type = TargetType::Structure;
        } else if (hypers["target_type"] == "Atom") {
          this->target_type = TargetType::Atom;
        } else {
          throw std::runtime_error("Given target_type " +
                                   hypers["target_type"].get<std::string>() +
                                   " is not known."
                                   " It is either 'Structure' or 'Atom')");
        }
      } else {
        throw std::runtime_error(
            R"(No target_type given. It is either 'Structure' or Atom')");
      }

      auto kernel_type_str = hypers.at("name").get<std::string>();
      if (kernel_type_str == "GAP") {
        this->kernel_type = SparseKernelType::GAP;
        this->kernel_impl =
            make_sparse_kernel_impl<SparseKernelType::GAP>(hypers);
      } else {
        throw std::logic_error("Requested SparseKernel \'" + kernel_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'GAP\'.");
      }
    }

    //! Move constructor
    SparseKernel(SparseKernel && other) noexcept
        : identifiers{std::move(other.identifiers)}, parameters{},
          target_type{std::move(other.target_type)},
          kernel_type{std::move(other.kernel_type)}, kernel_impl{std::move(
                                                         other.kernel_impl)} {
      this->parameters = std::move(other.parameters);
    }

    /**
     * The root compute kernel function. It computes the kernel between 2 set of
     * structures for a given representation specified by the calculator.
     *
     * @param calculator the calculator which has been used to calculate
     * the representation on the two managers
     * has been registered in the elements of managers and managers_b
     * @param managers a ManagerCollection or similar collection of
     * structure managers
     * @param managers_b a ManagerCollection or similar collection of
     * structure managers or a set of pseudo points
     */
    template <class Calculator, class StructureManagers, class SparsePoints>
    math::Matrix_t compute(const Calculator & calculator,
                           const StructureManagers & managers,
                           const SparsePoints & sparse_points) {
      using ManagerPtr_t = typename StructureManagers::value_type;
      using Manager_t = typename ManagerPtr_t::element_type;
      using Property_t = typename Calculator::template Property_t<Manager_t>;
      auto && representation_name{calculator.get_name()};
      using internal::TargetType;

      switch (this->target_type) {
      case TargetType::Structure:
        return this->compute_helper<Property_t, TargetType::Structure>(
            representation_name, managers, sparse_points);
      case TargetType::Atom:
        return this->compute_helper<Property_t, TargetType::Atom>(
            representation_name, managers, sparse_points);
      default:
        throw std::logic_error(
            "Given target_type " +
            this->parameters["target_type"].get<std::string>() +
            " is not known."
            " It is either 'Structure' or 'Atom')");
      }
    }

    template <class Property_t, internal::TargetType Type,
              class StructureManagers, class SparsePoints>
    math::Matrix_t compute_helper(const std::string & representation_name,
                                  const StructureManagers & managers,
                                  const SparsePoints & sparse_points) {
      using internal::SparseKernelType;

      if (this->kernel_type == SparseKernelType::GAP) {
        auto kernel =
            downcast_sparse_kernel_impl<SparseKernelType::GAP>(kernel_impl);
        return kernel->template compute<Property_t, Type>(
            managers, sparse_points, representation_name);
      } else {
        throw std::logic_error(
            "Given kernel_type " +
            this->parameters["kernel_type"].get<std::string>() +
            " is not known."
            " It is 'GAP'");
      }
    }

    template <class SparsePoints>
    math::Matrix_t compute(const SparsePoints & sparse_points) {
      using internal::SparseKernelType;

      if (this->kernel_type == SparseKernelType::GAP) {
        auto kernel =
            downcast_sparse_kernel_impl<SparseKernelType::GAP>(kernel_impl);
        return kernel->compute(sparse_points);
      } else {
        throw std::logic_error(
            "Given kernel_type " +
            this->parameters["kernel_type"].get<std::string>() +
            " is not known."
            " It is 'GAP'");
      }
    }

    /**
     * The root compute kernel function. It computes the kernel between the
     * representation gradients of a set of structures with the set of pseudo
     * points.
     * The energy is defined as a sum of atomic contributions
     * @f[
     *      E(A) = \sum_{i \in A} \epsilon_{a_i}(X_i),
     * @f]
     *
     * where \f$\epsilon_{a_i}\f$ is an atomic energy function for atoms of
     * type \f$a_i\f$ modeled with a GPR model with sparse points \f$T_m\f$
     *
     * @f[
     *      \epsilon_{a_i} (X_i) =
     *        \sum_m \alpha_m \delta_{b_m a_i} k(X_i,T_m),
     * @f]
     *
     * where \f$b_m\f$ is the atom type associated with the sparse point
     * \f$T_m\f$. To summarize our model, the energy is given by
     *
     * @f[
     *      E(A) = \sum_m \alpha_m \sum_{i \in A} \delta_{b_m a_i} k(X_i,T_m),
     * @f]
     *
     *
     * thus the forces with respect to the position of atom \f$i\f$ is in the
     kernel formulation
     *
     * @f[
     *      \vec{\nabla}_i E = \sum_m \alpha_m \sum_{i \in A} \delta_{b_m a_i}
     \sum_{j \in i} \vec{\nabla}_i k(X_j,T_m),
     * @f]
     *
     * and the virial negative stress
     *
     * @f[
     *      \frac{ \partial E }{ \partial\eta_{\alpha\beta} } =
     *            \sum_m  \alpha_m \sum_{i \in A}  \delta_{b_m a_i}
     *            \frac{ \partial k(X_i, T_m) }{ \partial\eta_{\alpha\beta} },
     * @f]
     *
     * where \f$\eta_{\alpha\beta}\f$ correspond to the 3x3 deformation (strain)
     tensor
     * in the \f$\alpha\beta \in \{xx, yy, zz, yz, xz, xy\}\f$ spatial
     dimensions (following Voigt convention).

     *
     * @param calculator the calculator which has been used to calculate
     *        the representation on the two managers
     *        has been registered in the elements of managers and managers_b
     * @param sparse_points class of pseudo points
     * @param managers_b a ManagerCollection or similar collection of
     *        structure managers
     */
    template <class Calculator, class StructureManagers, class SparsePoints>
    math::Matrix_t compute_derivative(const Calculator & calculator,
                                      const StructureManagers & managers,
                                      const SparsePoints & sparse_points,
                                      const bool compute_neg_stress) {
      using ManagerPtr_t = typename StructureManagers::value_type;
      using Manager_t = typename ManagerPtr_t::element_type;
      using Property_t = typename Calculator::template Property_t<Manager_t>;
      using PropertyGradient_t =
          typename Calculator::template PropertyGradient_t<Manager_t>;
      if (not calculator.does_gradients()) {
        throw std::runtime_error(
            "This representation does not compute gradients.");
      }
      const auto & representation_name{calculator.get_name()};
      const auto representation_grad_name{calculator.get_gradient_name()};
      using internal::SparseKernelType;
      using internal::TargetType;

      if (this->kernel_type == SparseKernelType::GAP) {
        auto kernel =
            downcast_sparse_kernel_impl<SparseKernelType::GAP>(kernel_impl);
        return kernel->template compute_derivative<
            Property_t, PropertyGradient_t, TargetType::Atom>(
            managers, sparse_points, representation_name,
            representation_grad_name, compute_neg_stress);
      } else {
        throw std::logic_error(
            "Given kernel_type " +
            this->parameters["kernel_type"].get<std::string>() +
            " is not known."
            " It is 'GAP'");
      }
    }

    //! list of names identifying the properties that should be used
    //! to compute the kernels
    std::vector<std::string> identifiers{};
    //! parameters of the kernel
    Hypers_t parameters{};
    /**
     * Defines if the similarity is defined on the level of structures or per
     * atom
     */
    internal::TargetType target_type{};

    internal::SparseKernelType kernel_type{};
    std::shared_ptr<internal::KernelImplBase> kernel_impl{};
  };

}  // namespace rascal

namespace nlohmann {
  /**
   * Special specialization of the json serialization for non default
   * constructible type.
   */
  template <>
  struct adl_serializer<rascal::SparseKernel> {
    static rascal::SparseKernel from_json(const json & j) {
      return rascal::SparseKernel{j};
    }

    static void to_json(json & j, const rascal::SparseKernel & t) {
      j = t.parameters;
    }
  };
}  // namespace nlohmann

#endif  // SRC_RASCAL_MODELS_SPARSE_KERNELS_HH_
