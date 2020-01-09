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
#include "rascal/structure_managers/structure_manager_collection.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/models/kernels.hh"


namespace rascal {

  namespace internal {
    enum class SparseKernelType {SparseGAP};

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
    struct SparseKernelImpl<internal::SparseKernelType::SparseGAP> : KernelImplBase {
      using Hypers_t = typename KernelImplBase::Hypers_t;

      //! exponent of the cosine kernel
      size_t zeta{1};
      constexpr static int n_spatial_dimensions{3};

      SparseKernelImpl() = default;

      explicit SparseKernelImpl(const Hypers_t & hypers) : KernelImplBase{} {
        this->set_hyperparmeters(hypers);
      }

      void set_hyperparmeters(const Hypers_t & hypers) {
        if (hypers.count("zeta") == 1) {
          zeta = hypers["zeta"].get<size_t>();
        } else {
          throw std::runtime_error(
              R"(zeta should be specified for the SparseGAP kernel)");
        }
      }

      /**
       * Compute the kernel between a set of structure(s) and a set of pseudo
       * points, structure wise.
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *          of structure managers like ManagerCollection
       * @param pseudo_points a PseudoPoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation data
       * has been registered in the elements of managers and managers_b
       * @return kernel matrix
       */
      template <
          class Property_t, internal::TargetType Type,
          std::enable_if_t<Type == internal::TargetType::Structure, int> = 0,
          class StructureManagers, class PseudoPoints>
      math::Matrix_t compute(StructureManagers & managers,
                             PseudoPoints & pseudo_points,
                             const std::string & representation_name) {
        math::Matrix_t KNM(managers.size(), pseudo_points.size());
        KNM.setZero();
        size_t ii_A{0};
        for (auto & manager : managers) {
          auto && propA{*manager->template get_property<Property_t>(
              representation_name, true)};
          for (auto center : manager) {
            int sp = center.get_atom_type();
            // only the pseudo points of species sp contribute
            KNM.row(ii_A) += pow_zeta(pseudo_points.dot(sp, propA[center]), this->zeta);
          }
          ++ii_A;
        }
        return KNM;
      }

      /**
       * Compute the kernel between a set of pseudo points.
       *
       * @tparam PseudoPoints should be a set a set of PseudoPoints
       * @param pseudo_points a PseudoPoints* class
       * @return kernel matrix MxM
       */
      template <class PseudoPoints>
      math::Matrix_t compute(PseudoPoints & pseudo_points) {
        math::Matrix_t KMM(pseudo_points.size(), pseudo_points.size());
        KMM.setZero();
        int start{0};
        // loop over the species
        for (const int& sp : pseudo_points.species()) {
          // only the pseudo points of the same species contribute
          auto KMM_by_sp = pseudo_points.self_dot(sp);
          auto block_size = KMM_by_sp.rows();
          KMM.block(start, start, block_size, block_size) = pow_zeta(KMM_by_sp, this->zeta);
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
       * @param pseudo_points a PseudoPoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation data
       * has been registered in the elements of managers and managers_b
       * @return kernel matrix
       */
      template <class Property_t, internal::TargetType Type,
                std::enable_if_t<Type == internal::TargetType::Atom, int> = 0,
                class StructureManagers, class PseudoPoints>
      math::Matrix_t compute(const StructureManagers & managers,
                             const PseudoPoints & pseudo_points,
                             const std::string & representation_name) {
        size_t n_centersA{0};
        for (const auto & manager : managers) {
          n_centersA += manager->size();
        }
        size_t n_pseudo_points{pseudo_points.size()};
        math::Matrix_t KNM(n_centersA, n_pseudo_points);
        size_t ii_A{0};
        for (auto & manager : managers) {
          auto && propA{*manager->template get_property<Property_t>(
              representation_name, true)};
          for (auto center : manager) {
            int sp = center.get_atom_type();
            KNM.row(ii_A) = pow_zeta(pseudo_points.dot(sp, propA[center]), this->zeta);
            ii_A++;
          }
        }
        return KNM;
      }


      /**
       * Compute the kernel between the representation gradients of structure(s)
       * and a set of pseudo points.
       *
       * @tparam StructureManagers should be an iterable over shared pointer
       *          of structure managers like ManagerCollection
       * @param pseudo_points a PseudoPoints* class
       * @param managers a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation
       * gradient data has been registered in the elements of managers
       * @return kernel matrix
       */
      template <
          class Property_t, class PropertyGradient_t, internal::TargetType Type,
          std::enable_if_t<Type == internal::TargetType::Atom, int> = 0,
          class StructureManagers, class PseudoPoints>
      math::Matrix_t compute_derivative(StructureManagers & managers,
                             PseudoPoints & pseudo_points,
                             const std::string & representation_name,
                             const std::string & representation_grad_name) {
        size_t n_centers{0};
        for (const auto & manager : managers) {
          n_centers += manager->size() * n_spatial_dimensions;
        }
        size_t n_pseudo_points{pseudo_points.size()};
        math::Matrix_t KNM(n_centers, n_pseudo_points);
        KNM.setZero();
        size_t ii_A{0};
        for (auto & manager : managers) {
          auto && prop{*manager->template get_property<Property_t>(
              representation_name, true)};
          auto && prop_grad{*manager->template get_property<PropertyGradient_t>(
              representation_grad_name, true)};
          for (auto center : manager) {
            for (auto neigh : center.pairs_with_self_pair()) {
              int sp = neigh.get_atom_type();
              auto atom_j = neigh.get_atom_j();
              auto row = this->zeta * pow_zeta(pseudo_points.dot(sp, prop[atom_j]), this->zeta - 1).transpose();

              KNM.block(ii_A, 0, n_spatial_dimensions, n_pseudo_points) += row.array() * pseudo_points.dot_grad(sp, prop_grad[neigh]).transpose().array().rowwise();
            }
            ii_A += n_spatial_dimensions;
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
    return std::static_pointer_cast<internal::SparseKernelImpl<Type>>(kernel_impl);
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
      if (kernel_type_str == "SparseGAP") {
        this->kernel_type = SparseKernelType::SparseGAP;
        this->kernel_impl = make_sparse_kernel_impl<SparseKernelType::SparseGAP>(hypers);
      } else {
        throw std::logic_error("Requested SparseKernel \'" + kernel_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'SparseGAP\'.");
      }
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
        throw std::logic_error("The combination of parameter is not handled.");
      }
    }

    template <class Property_t, internal::TargetType Type,
              class StructureManagers, class SparsePoints>
    math::Matrix_t compute_helper(const std::string & representation_name,
                                  const StructureManagers & managers,
                                  const SparsePoints & sparse_points) {
      using internal::SparseKernelType;

      if (this->kernel_type == SparseKernelType::SparseGAP) {
        auto kernel = downcast_sparse_kernel_impl<SparseKernelType::SparseGAP>(kernel_impl);
        return kernel->template compute<Property_t, Type>(
            managers, sparse_points, representation_name);
      } else {
        throw std::logic_error("The combination of parameter is not handled.");
      }
    }

    template <class Calculator, class SparsePoints>
    math::Matrix_t compute(const Calculator & calculator,
                           const SparsePoints & sparse_points) {

      using internal::TargetType;

      switch (this->target_type) {
      case TargetType::Atom:
        return this->compute_helper<TargetType::Atom>(sparse_points);
      default:
        throw std::logic_error("The combination of parameter is not handled.");
      }
    }

    template <internal::TargetType Type,
              class SparsePoints>
    math::Matrix_t compute_helper(const SparsePoints & sparse_points) {
      using internal::SparseKernelType;

      if (this->kernel_type == SparseKernelType::SparseGAP) {
        auto kernel = downcast_sparse_kernel_impl<SparseKernelType::SparseGAP>(kernel_impl);
        return kernel->template compute<Type>(sparse_points);
      } else {
        throw std::logic_error("The combination of parameter is not handled.");
      }
    }

    /*
     * The root compute kernel function. It computes the kernel the
     * representation gradients of a set of structures and the set of pseudo
     * points.
     *
     * @param calculator the calculator which has been used to calculate
     * the representation on the two managers
     * has been registered in the elements of managers and managers_b
     * @param pseudo_points class of pseudo points
     * @param managers_b a ManagerCollection or similar collection of
     * structure managers
     */
    template <class Calculator, class PseudoPoints, class StructureManagers>
    math::Matrix_t compute_derivative(const Calculator & calculator,
                           const StructureManagers & managers,
                           const PseudoPoints & pseudo_points) {
      using ManagerPtr_t = typename StructureManagers::value_type;
      using Manager_t = typename ManagerPtr_t::element_type;
      using Property_t = typename Calculator::template Property_t<Manager_t>;
      using PropertyGradient_t = typename Calculator::template PropertyGradient_t<Manager_t>;
      auto && representation_name{calculator.get_name()};
      auto && representation_grad_name{calculator.get_gradient_name()};
      using internal::TargetType;
      using internal::SparseKernelType;

      if (this->kernel_type == SparseKernelType::SparseGAP) {
        auto kernel = downcast_sparse_kernel_impl<SparseKernelType::SparseGAP>(kernel_impl);
        return kernel->template compute<Property_t, PropertyGradient_t, TargetType::Atom>(
            managers, pseudo_points, representation_name, representation_grad_name);
      } else {
        throw std::logic_error("The combination of parameter is not handled.");
      }
    }


    //! list of names identifying the properties that should be used
    //! to compute the kernels
    std::vector<std::string> identifiers{};
    //! parameters of the kernel
    Hypers_t parameters{};
    /**
     * Defines if the similarity if defined structure or atom wise
     */
    internal::TargetType target_type{};

    internal::SparseKernelType kernel_type{};
    std::shared_ptr<internal::KernelImplBase> kernel_impl{};
  };

}  // namespace rascal

#endif  // SRC_RASCAL_MODELS_SPARSE_KERNELS_HH_
