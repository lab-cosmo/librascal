/**
 * file   kernels.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   16 Jun 2019
 *
 * @brief Implementation of similarity kernels
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

#ifndef SRC_MODELS_KERNELS_HH_
#define SRC_MODELS_KERNELS_HH_

#include "math/math_utils.hh"
#include "structure_managers/structure_manager_collection.hh"
#include "json_io.hh"

namespace rascal {

  namespace internal {
    enum class KernelType { Cosine };

    enum class TargetType { Structure, Atom };

    struct KernelImplBase {
      using Hypers_t = json;
    };

    template <internal::KernelType Type>
    struct KernelImpl {};

    template <>
    struct KernelImpl<internal::KernelType::Cosine> : KernelImplBase {
      using Hypers_t = typename KernelImplBase::Hypers_t;

      //! exponent of the cosine kernel
      size_t zeta{1};

      KernelImpl() = default;

      explicit KernelImpl(const Hypers_t & hypers) : KernelImplBase{} {
        this->set_hyperparmeters(hypers);
      }

      void set_hyperparmeters(const Hypers_t & hypers) {
        if (hypers.count("zeta") == 1) {
          zeta = hypers["zeta"].get<size_t>();
        } else {
          throw std::runtime_error(
              R"(zeta should be specified for the cosine kernel)");
        }
      }

      template <
          class Property_t, internal::TargetType Type,
          std::enable_if_t<Type == internal::TargetType::Structure, int> = 0,
          class StructureManagers>
      inline math::Matrix_t compute(StructureManagers & managers_a,
                                    StructureManagers & managers_b,
                                    const std::string & representation_name) {
        math::Matrix_t kernel(managers_a.size(), managers_b.size());
        auto integer_power{math::MakePositiveIntegerPower<double>(this->zeta)};
        size_t ii_A{0};
        for (auto & manager_a : managers_a) {
          size_t ii_B{0};
          auto && propA{
              manager_a->template get_validated_property_ref<Property_t>(
                  representation_name)};
          for (auto & manager_b : managers_b) {
            auto && propB{
                manager_b->template get_validated_property_ref<Property_t>(
                    representation_name)};

            kernel(ii_A, ii_B) =
                propA.dot(propB).unaryExpr(integer_power).mean();
            ++ii_B;
          }
          ++ii_A;
        }
        return kernel;
      }

      /**
       * Compute the kernel between 2 set of structures for a given representation specified by the name.
       *
       * @param managers_a a ManagerCollection or similar collection of
       * structure managers
       * @param managers_b a ManagerCollection or similar collection of
       * structure managers
       * @param representation_name name under which the representation data
       * has been registered in the elements of managers_a and managers_b
       *
       * @return kernel matrix
       */
      template <class Property_t, internal::TargetType Type,
                std::enable_if_t<Type == internal::TargetType::Atom, int> = 0,
                class StructureManagers>
      inline math::Matrix_t compute(const StructureManagers & managers_a,
                                    const StructureManagers & managers_b,
                                    const std::string & representation_name) {
        size_t n_centersA{0};
        for (const auto & manager_a : managers_a) {
          n_centersA += manager_a->size();
        }
        size_t n_centersB{0};
        for (const auto & manager_b : managers_b) {
          n_centersB += manager_b->size();
        }
        math::Matrix_t kernel(n_centersA, n_centersB);
        auto integer_power{math::MakePositiveIntegerPower<double>(this->zeta)};
        size_t ii_A{0};
        for (auto & manager_a : managers_a) {
          size_t ii_B{0};
          auto && propA{
              manager_a->template get_validated_property_ref<Property_t>(
                  representation_name)};
          for (auto & manager_b : managers_b) {
            auto && propB{
                manager_b->template get_validated_property_ref<Property_t>(
                    representation_name)};

            kernel.block(ii_A, ii_B, manager_a->size(), manager_b->size()) =
                propA.dot(propB).unaryExpr(integer_power);
            ii_B += manager_b->size();
          }
          ii_A += manager_a->size();
        }
        return kernel;
      }
    };
  }  // namespace internal

  template <internal::KernelType Type, class Hypers>
  decltype(auto) make_kernel_impl(const Hypers & hypers) {
    return std::static_pointer_cast<internal::KernelImplBase>(
        std::make_shared<internal::KernelImpl<Type>>(hypers));
  }

  template <internal::KernelType Type>
  decltype(auto) downcast_kernel_impl(
      std::shared_ptr<internal::KernelImplBase> & kernel_impl) {
    return std::static_pointer_cast<internal::KernelImpl<Type>>(kernel_impl);
  }

  class Kernel {
   public:
    using Hypers_t = typename internal::KernelImplBase::Hypers_t;

    explicit Kernel(const Hypers_t & hypers) {
      using internal::KernelType;
      using internal::TargetType;

      this->parameters = hypers;

      if (hypers.count("target_type") == 1) {
        if (hypers["target_type"] == "Structure") {
          this->target_type = TargetType::Structure;
        } else if (hypers["target_type"] == "Atom") {
          this->target_type = TargetType::Atom;
        }
      } else {
        throw std::runtime_error(R"(target_type is either structure or atom)");
      }

      auto kernel_type_str = hypers.at("name").get<std::string>();
      if (kernel_type_str.compare("Cosine") == 0) {
        this->kernel_type = KernelType::Cosine;
        this->kernel_impl = make_kernel_impl<KernelType::Cosine>(hypers);
      } else {
        throw std::logic_error("Requested Kernel \'" + kernel_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'Cosine\'.");
      }
    }


    /*
     * The root compute kernel function. It computes the kernel between 2 set of structures for a given representation specified by the calculator.
     *
     * @param calculator the calculator which has been used to calculate
     * the representation on the two managers
     * has been registered in the elements of managers_a and managers_b
     * @param managers_a a ManagerCollection or similar collection of
     * structure managers
     * @param managers_b a ManagerCollection or similar collection of
     * structure managers
     */
    template <class Calculator, class StructureManagers>
    math::Matrix_t compute(const Calculator & calculator,
                           const StructureManagers & managers_a,
                           const StructureManagers & managers_b) {
      using ManagerPtr_t = typename StructureManagers::value_type;
      using Manager_t = typename ManagerPtr_t::element_type;
      using Property_t = typename Calculator::template Property_t<Manager_t>;
      auto && representation_name{calculator.get_name()};
      using internal::TargetType;

      switch (this->target_type) {
      case TargetType::Structure: {
        return this->compute_helper<Property_t, TargetType::Structure>(
            representation_name, managers_a, managers_b);
        break;
      }
      case TargetType::Atom: {
        return this->compute_helper<Property_t, TargetType::Atom>(
            representation_name, managers_a, managers_b);
        break;
      }
      default:
        throw std::logic_error("The combination of parameter is not handdled.");
        break;
      }
    }

    template <class Property_t, internal::TargetType Type,
              class StructureManagers>
    inline math::Matrix_t
    compute_helper(const std::string & representation_name,
                   const StructureManagers & managers_a,
                   const StructureManagers & managers_b) {
      using internal::KernelType;

      switch (this->kernel_type) {
      case KernelType::Cosine: {
        auto kernel = downcast_kernel_impl<KernelType::Cosine>(kernel_impl);
        return kernel->template compute<Property_t, Type>(
            managers_a, managers_b, representation_name);
        break;
      }
      default:
        throw std::logic_error("The combination of parameter is not handdled.");
        break;
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

    internal::KernelType kernel_type{};
    std::shared_ptr<internal::KernelImplBase> kernel_impl{};
  };

}  // namespace rascal

#endif  // SRC_MODELS_KERNELS_HH_
