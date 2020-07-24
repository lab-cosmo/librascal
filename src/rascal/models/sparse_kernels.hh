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
     * local environments $X_i^a$ and $X_j^b$ with central atom types $a$ and
     * $b$ is:
     *
     * $$k(X_i^a, X_j^b) = \delta_{ab} k(X_i, X_j),$$
     *
     * where $k(X_i, X_j)$ is the cosine kernel. When building a model for
     * properties associated with the structure we assume that the training
     * will be done on the property itself (not divided by the number of
     * atoms).
     *
     * It is particularly designed to build the $K_{MM}$ and $K_{NM}$ matrices
     * needed by the sparse kernel formulation where the kernel matrix is given
     * by:
     *
     * $$K = K_{MM} + K_{MN} \Lambda^{-2} K_{NM}$$
     *
     * with $\Lambda$ the regularization matrix.
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
        size_t n_sparse_points{sparse_points.size()};
        math::Matrix_t KNM(n_centersA, n_sparse_points);
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
       * Compute the kernel between the representation gradients of structure(s)
       * and a set of pseudo points.
       * The derivative of the kernel of environment $X_j^c$ with respect to
       * the position of $r_i$ atom $i$ of type $c$ and the pseudo point T_m^b
       * is:
       *
       * $$dk(X_j^c, T_m^b)/dr_i =
       *              \sum_{p} dk(X_j^c, T_m^b)/dX_{jp}^c dX_j^c/dr_i$$,
       *
       * where $p$ is an index running over all the representation elements.
       * This object has the dimension of the number of neighbors (so quite
       * large) but to compute forces we only need the contraction over all
       * neighbors of a particular center so this routine return the already
       * contracted kernel derivative for center $X_i^a$:
       *
       * $$dk(X_i^a, T_m^b)/dr_i =
       *            \sum_{(c,j) \in X_i^a} dk(X_j^c, T_m^b)/dr_i$$,
       *
       * so the resulting matrix has the dimension of the number of centers
       * times 3 (the three cartesian dimensions) and the number of sparse
       * points.
       *
       * For the GAP kernel we have:
       *
       * $$dk(X_j^c, T_m^b)/dX_{jp}^c =
       *   \delta_{cb} \zeta (\sum_n X_{jn}^c T_{mn}^b)^{\zeta-1}  T_{mp}^b$$
       *
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *          of structure managers like ManagerCollection
       * @param sparse_points a SparsePoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation
       * gradient data has been registered in the elements of managers
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
                         const std::string & representation_grad_name) {
        using Manager_t = typename StructureManagers::Manager_t;
        size_t n_centers{0};
        // find the total number of rows the matrix block should have
        for (const auto & manager : managers) {
          n_centers += manager->size() * SpatialDims;
        }
        size_t n_sparse_points{sparse_points.size()};
        math::Matrix_t KNM(n_centers, n_sparse_points);
        math::Matrix_t spts{sparse_points.get_features()};
        const int n_features{static_cast<int>(spts.cols())};
        std::map<int, std::array<int, 2>> spts_slices{
            sparse_points.get_blocks_info()};

        size_t i_center{0};
        // loop over the structures
        for (auto & manager : managers) {
          Property<double, 1, Manager_t, SpatialDims, Eigen::Dynamic> dKdr{
              *manager, "no metadata", true};
          auto && prop{*manager->template get_property<Property_t>(
              representation_name, true)};
          dKdr.set_nb_col(n_sparse_points);
          dKdr.resize();
          dKdr.setZero();
          auto && prop_grad{*manager->template get_property<PropertyGradient_t>(
              representation_grad_name, true)};
          // compute dX/dr * T
          // sparse_points.dot_derivative(prop_grad, manager, dKdr);

          math::Matrix_t dXdr =
              prop_grad.get_features_gradient(sparse_points.keys);
          int i_row{0};
          for (auto center : manager) {
            auto a_sp{center.get_atom_type()};
            auto n_rows{center.pairs_with_self_pair().size() * ThreeD};
            const auto & spts_slice{spts_slices[a_sp]};
            math::Matrix_t KNM_block =
                dXdr.block(i_row, 0, n_rows, n_features) *
                spts.block(spts_slice[0], 0, spts_slice[1], n_features)
                    .transpose();

            if (this->zeta > 1) {
              // dk/dX without the pseudo point factor
              Eigen::Matrix<double, 1, Eigen::Dynamic> rep =
                  this->zeta * pow_zeta(sparse_points.dot(a_sp, prop[center]),
                                        this->zeta - 1)
                                   .transpose();
              KNM_block *=
                  rep.segment(spts_slice[0], spts_slice[1]).asDiagonal();
            }

            int i_row_{0};
            for (auto neigh : center.pairs_with_self_pair()) {
              dKdr[neigh.get_atom_j()].block(0, spts_slice[0], ThreeD,
                                             spts_slice[1]) +=
                  KNM_block.block(i_row_, 0, ThreeD, spts_slice[1]);
              i_row_ += ThreeD;
            }
            i_row += n_rows;
          }

          // copy the data to the kernel matrix
          for (auto center : manager) {
            KNM.block(i_center, 0, SpatialDims, n_sparse_points) = dKdr[center];
            i_center += SpatialDims;
          }
        }
        return KNM;
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

    /*
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

    /*
     * The root compute kernel function. It computes the kernel between the
     * representation gradients of a set of structures with the set of pseudo
     * points.
     *
     * @param calculator the calculator which has been used to calculate
     * the representation on the two managers
     * has been registered in the elements of managers and managers_b
     * @param sparse_points class of pseudo points
     * @param managers_b a ManagerCollection or similar collection of
     * structure managers
     */
    template <class Calculator, class StructureManagers, class SparsePoints>
    math::Matrix_t compute_derivative(const Calculator & calculator,
                                      const StructureManagers & managers,
                                      const SparsePoints & sparse_points) {
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
            representation_grad_name);
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
