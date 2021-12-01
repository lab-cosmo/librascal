/**
 * @file test_sparse_kernels.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 Jan 2020
 *
 * @brief test the implementation of sparse similarity kernel classes
 *
 * @section LICENSE
 *
 * Copyright  2020 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TESTS_TEST_SPARSE_KERNELS_HH_
#define TESTS_TEST_SPARSE_KERNELS_HH_

#include "test_adaptor.hh"
#include "test_calculator.hh"
#include "test_manager_collection.hh"

#include "rascal/models/numerical_kernel_gradients.hh"
#include "rascal/models/sparse_kernel_predict.hh"
#include "rascal/models/sparse_kernels.hh"
#include "rascal/models/sparse_points.hh"

namespace rascal {

  struct StrictNLSparseKernelFixture : StrictNLCCCollectionFixture {
    using Parent = StrictNLCCCollectionFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    using Representation_t = CalculatorSphericalInvariants;

    StrictNLSparseKernelFixture() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function_type"] = fc_hyp["cutoff_function_type"];
              rep_hyp["interaction_cutoff"] = fc_hyp["interaction_cutoff"];
              rep_hyp["cutoff_smooth_width"] = fc_hyp["cutoff_smooth_width"];
              if (fc_hyp.count("cutoff_function_parameters")) {
                rep_hyp["cutoff_function_parameters"] =
                    fc_hyp["cutoff_function_parameters"];
              }

              rep_hyp["gaussian_sigma_type"] = sig_hyp["gaussian_sigma_type"];
              rep_hyp["gaussian_sigma_constant"] =
                  sig_hyp["gaussian_sigma_constant"];

              rep_hyp["radial_basis"] = ri_hyp["radial_basis"];
              this->representation_hypers.push_back(rep_hyp);
            }
          }
        }
      }
    }
    ~StrictNLSparseKernelFixture() = default;

    std::vector<json> representation_hypers{};

    std::vector<json> fc_hypers{{{"cutoff_function_type", "ShiftedCosine"},
                                 {"interaction_cutoff", 2.0},
                                 {"cutoff_smooth_width", 0.5}}};

    std::vector<json> density_hypers{{{"gaussian_sigma_type", "Constant"},
                                      {"gaussian_sigma_constant", 0.4}}};
    std::vector<json> radial_contribution_hypers{{{"radial_basis", "GTO"}}};
    std::vector<json> rep_hypers{
        {{"max_radial", 3},
         {"max_angular", 0},
         {"soap_type", "RadialSpectrum"},
         {"expansion_by_species_method", "environment wise"},
         {"compute_gradients", false},
         {"normalize", true}},
        {{"max_radial", 2},
         {"max_angular", 2},
         {"soap_type", "PowerSpectrum"},
         {"expansion_by_species_method", "environment wise"},
         {"compute_gradients", false},
         {"normalize", true}}};

    std::vector<json> kernel_hypers{
        {{"zeta", 2}, {"target_type", "Structure"}, {"name", "GAP"}},
        {{"zeta", 2}, {"target_type", "Atom"}, {"name", "GAP"}}};
  };

  /**
   * BaseFixture is expected to be similar to
   * StrictNLSparseKernelFixture
   */
  template <class BaseFixture>
  struct SparseKernelFixture : CollectionFixture<BaseFixture>,
                               CalculatorFixture<BaseFixture> {
    using ParentA = CollectionFixture<BaseFixture>;
    using ParentB = CalculatorFixture<BaseFixture>;

    using ManagerCollection_t = typename ParentA::ManagerCollection_t;
    using Manager_t = typename ManagerCollection_t::Manager_t;
    using Calculator_t = typename ParentB::Representation_t;
    using Property_t = typename Calculator_t::template Property_t<Manager_t>;

    SparseKernelFixture() : ParentA{}, ParentB{} {
      for (auto & collection : this->collections) {
        collection.add_structures(this->ParentA::filename, this->ParentA::start,
                                  this->ParentA::length);
        for (auto & hypers : this->ParentB::representation_hypers) {
          this->representations.emplace_back(hypers);
          this->representations.back().compute(collection);
        }
      }

      for (auto & hypers : this->ParentA::kernel_hypers) {
        this->kernels.emplace_back(hypers);
      }
    }

    std::vector<SparseKernel> kernels{};

    ~SparseKernelFixture() = default;

    bool verbose{false};
  };

}  // namespace rascal

#endif  // TESTS_TEST_SPARSE_KERNELS_HH_
