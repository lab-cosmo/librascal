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

  struct StrictNLKernelFixture : StrictNLCollectionFixture {
    using Parent = StrictNLCollectionFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalInvariants;

    StrictNLKernelFixture() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~StrictNLKernelFixture() = default;

    std::vector<json> representation_hypers{};

    std::vector<json> fc_hypers{
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "A"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "A"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "A"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{{{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}}};
  };

  /**
   * BaseFixture is expected to be similar to
   * StrictNLKernelFixture
   */
  template<internal::KernelType Type, class BaseFixture>
  struct KernelFixture : CollectionFixture<BaseFixture>, CalculatorFixture<BaseFixture> {
    using ParentA = CollectionFixture<BaseFixture>;
    using ParentB = CalculatorFixture<BaseFixture>;

    using ManagerCollection_t = typename ParentA::ManagerCollection_t;
    using Manager_t = typename ManagerCollection_t::Manager_t;
    using Calculator_t = typename ParentB::Representation_t;
    using Property_t = typename Calculator_t::template Property_t<Manager_t>;

    KernelFixture() :ParentA{}, ParentB{} {
      for (auto& collection : this->collections) {
        collection.add_structures(this->ParentA::filename, this->ParentA::start, this->ParentA::lenght);
        for (auto& hyper : this->ParentB::representation_hypers) {
          this->representations.push_back(hyper);
          this->representations.back().compute(collection);
        }
      }

      for (auto& hyper : this->kernel_hypers) {
        this->kernels.push_back(hyper);
      }
    }

    std::vector<json> kernel_hypers{{ {"zeta", 2},
                                      {"target_type", "structure"},
                                      {"name", "Cosine"}},
                                    { {"zeta", 2},
                                      {"target_type", "atom"},
                                      {"name", "Cosine"}}
                                    };

    std::vector<Kernel> kernels{};

    ~KernelFixture() = default;

    bool verbose{false};

  };

}  // namespace rascal

#endif  // TESTS_TEST_KERNELS_HH_