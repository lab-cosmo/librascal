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

  enum class KernelType {Cosine};

  enum class TargetType {Structure, Atom};

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

        } else if (hypers["target_type"] == "atom") {

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
    TargetType target_type{};
  };

  template<KernelType Type>
  struct Kernel : {};

  template<>
  struct Kernel<KernelType::Cosine> : KernelBase {
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

    template<class StructureManager>
    math::Matrix_t compute(const std::string& representation_name, const ManagerCollection<StructureManager>& managersA, const ManagerCollection<StructureManager>& managersB) {
      math::Matrix_t kernel(managersA.size(), managersB.size());

      for (const auto& managerA : managersA) {
        for (const auto& managerB : managersB) {

        }
      }
    }

  };


}




#endif  // SRC_MATH_KERNELS_HH_
