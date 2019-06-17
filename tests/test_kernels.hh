/**
 * file test_kernels.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 June 2019
 *
 * @brief test the implementation of similarity kernel classes
 *
 * @section LICENSE
 *
 * Copyright  2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TESTS_TEST_KERNELS_HH_
#define TESTS_TEST_KERNELS_HH_

#include "tests.hh"
#include "test_adaptor.hh"
#include "test_manager_collection.hh"
#include "test_calculator.hh"
#include "models/kernels.hh"

namespace rascal {

  template<internal::KernelType Type, class CalculatorFixture, class CollectionFixture>
  struct KernelFixture : CollectionFixture, CalculatorFixture {
    using ManagerCollection_t = typename CollectionFixture::ManagerCollection_t;
    using Manager_t = typename ManagerCollection_t::Manager_t;
    using Calculator_t = typename CalculatorFixture::Representation_t;
    using Property_t = typename Calculator_t::template Property_t<Manager_t>;

    KernelFixture() :CollectionFixture{}, CalculatorFixture{} {
      this->collection.set_adaptor_inputs(this->adaptors);
      this->collection.add_structures(this->filename, 0, 20);
      for (auto& hyper : this->hypers) {
        this->calculators.push_back(hyper);
        this->calculators.back().compute(this->collection);
      }

      for (auto& hyper : this->kernel_hypers) {
        this->kernels.push_back(hyper);
      }

    }

    std::vector<json> kernel_hypers{{{"zeta", 2},
                                      {"target_type", "structure"},
                                      {"zeta", 2},
                                      {"target_type", "atom"}}};

    std::vector<Kernel<Type>> kernels{};

    ~KernelFixture() = default;

    std::vector<Calculator_t> calculators{};

    bool verbose{false};

  };

}  // namespace rascal

#endif  // TESTS_TEST_KERNELS_HH_