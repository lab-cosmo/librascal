/**
 * @file   test_utils.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date    2019
 *
 * @brief Test some utility functions
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TESTS_TEST_UTILS_HH_
#define TESTS_TEST_UTILS_HH_

#include "tests.hh"
#include "rascal_utility.hh"
#include "structure_managers/cluster_ref_key.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/structure_manager_centers.hh"

namespace rascal {
  namespace internal {
    template <typename T>
    struct PrettyFunction {
      /**
       * Get the type name
       * @returns std::string name of the type
       */
      static const std::string get() {
#define FUNCTION_MACRO __PRETTY_FUNCTION__
        const size_t funcNameLength{sizeof(FUNCTION_MACRO) - 1u};
        std::string typeName{FUNCTION_MACRO, funcNameLength};
        return typeName;
#undef FUNCTION_MACRO
      }
    };
  }  // namespace internal
}  // namespace rascal

#endif  // TESTS_TEST_UTILS_HH_
