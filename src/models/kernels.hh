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
    enum class KernelType {Cosine};

    enum class TargetType {Structure, Atom};
  }


  struct KernelBase {
    using Hypers_t = json;

    KernelBase(const Hypers_t& hypers) {
      this->parameters = hypers;
      if (hypers.count("representation_identifiers") == 1) {
        for (auto& id : hypers["representation_identifiers"]) {
          identifiers.emplace_back(id.get<std::string>());
        }
      }

      if (hypers.count("target_type") == 1) {
        if (hypers["target_type"] == "structure") {
          this->target_type = internal::TargetType::Structure;
        } else if (hypers["target_type"] == "atom") {
          this->target_type = internal::TargetType::Atom;
        }
      } else {
        throw std::runtime_error(R"(target_type is either structure or atom)");
      }
    }

    /**
     * Computes the kernel associated with a given identifier
     * This should be implemented but can't be virtual since
     * ManagerCollection is templated.
     */
    // virtual math::Matrix_t compute(std::string property_identifier,
    //                      ManagerCollection<>, ManagerCollection<>) = 0;

    //! list of names identifying the properties that should be used
    //! to compute the kernels
    std::vector<std::string> identifiers{};
    //! parameters of the kernel
    Hypers_t parameters{};
    /**
     * Defines if the similarity if defined structure or atom wise
     */
    internal::TargetType target_type{};
  };

  template<internal::KernelType Type>
  struct Kernel {};

  template<>
  struct Kernel<internal::KernelType::Cosine> : KernelBase {
    using Parent = KernelBase;
    using Hypers_t = typename Parent::Hypers_t;

    //! exponent of the cosine kernel
    size_t zeta{1};

    Kernel(const Hypers_t& hypers) : Parent{hypers} {
      if (hypers.count("zeta") == 1) {
        zeta = hypers["zeta"].get<size_t>();
      } else {
        throw std::runtime_error(R"(zeta should be specified for the global cosine kernel)");
      }
    }

    template<class Calculator, class StructureManagers>
    math::Matrix_t compute(const Calculator& calculator,
                        const StructureManagers& managersA,
                        const StructureManagers& managersB) {
      using ManagerPtr_t = typename StructureManagers::value_type;
      using Manager_t = typename ManagerPtr_t::element_type;
      using Property_t = typename Calculator::template Property_t<Manager_t>;
      auto&& representation_name{calculator.get_name()};

      switch (this->target_type) {
      case internal::TargetType::Structure: {
        return this->compute<Property_t, internal::TargetType::Structure>(managersA, managersB, representation_name);
        break;
      }
      case internal::TargetType::Atom: {
        return this->compute<Property_t, internal::TargetType::Atom>(managersA, managersB, representation_name);
        break;
      }
      default:
        throw std::logic_error("The combination of parameter is not handdled.");
        break;
      }
    }

    /**
     *
     */
    template<class Property_t, internal::TargetType Type, std::enable_if_t<Type == internal::TargetType::Structure, int> = 0, class StructureManagers>
    math::Matrix_t compute( StructureManagers& managersA,
                            StructureManagers& managersB,
                           const std::string& representation_name) {

      math::Matrix_t kernel(managersA.size(), managersB.size());
      auto integer_power{math::MakePositiveIntegerPower<double>(this->zeta)};
      size_t ii_A{0};
      for ( auto& managerA : managersA) {
        size_t ii_B{0};
        auto&& propA{managerA->template get_validated_property_ref<Property_t>(representation_name)};
        for ( auto& managerB : managersB) {
          auto&& propB{managerB->template get_validated_property_ref<Property_t>(representation_name)};

          kernel(ii_A, ii_B) = propA.dot(propB).unaryExpr(integer_power).mean();
          ++ii_B;
        }
        ++ii_A;
      }
      return kernel;
    }

    template<class Property_t, internal::TargetType Type, std::enable_if_t<Type == internal::TargetType::Atom, int> = 0, class StructureManagers>
    math::Matrix_t compute(const StructureManagers& managersA,
                           const StructureManagers& managersB,
                           const std::string& representation_name) {
      size_t n_centersA{0};
      for (const auto& managerA : managersA) {
        n_centersA += managerA->size();
      }
      size_t n_centersB{0};
      for (const auto& managerB : managersB) {
        n_centersB += managerB->size();
      }
      math::Matrix_t kernel(n_centersA, n_centersB);
      auto integer_power{math::MakePositiveIntegerPower<double>(this->zeta)};
      size_t ii_A{0};
      for ( auto& managerA : managersA) {
        size_t ii_B{0};
        auto&& propA{managerA->template get_validated_property_ref<Property_t>(representation_name)};
        for ( auto& managerB : managersB) {
          auto&& propB{managerB->template get_validated_property_ref<Property_t>(representation_name)};

          kernel.block(ii_A, ii_B, managerA->size(), managerB->size()) = propA.dot(propB).unaryExpr(integer_power);
          ii_B += managerB->size();
        }
        ii_A += managerA->size();
      }
      return kernel;
    }

  };


}




#endif  // SRC_MODELS_KERNELS_HH_
