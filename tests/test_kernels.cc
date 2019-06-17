/**
 * file test_kernels.cc
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


#include "test_kernels.hh"

namespace rascal {
  BOOST_AUTO_TEST_SUITE(kernels_test);

  using multiple_fixtures = boost::mpl::list<
      KernelFixture<internal::KernelType::Cosine,MultipleStructureSphericalInvariant, ManagerNLCollectionFixture>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_kernel_compute_test, Fix,
                                   multiple_fixtures, Fix) {
    auto& kernels = Fix::kernels;
    auto& calculators = Fix::calculators;

    for (auto & calculator : calculators) {
      for (auto & kernel : kernels) {
        kernel.compute(calculator, this->collection, this->collection);
      }
    }

  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
