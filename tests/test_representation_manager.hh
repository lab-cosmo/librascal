/**
 * file   test_representation_manager_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef TESTS_TEST_REPRESENTATION_MANAGER_HH_
#define TESTS_TEST_REPRESENTATION_MANAGER_HH_

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"
#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap.hh"
#include "representations/feature_manager_block_sparse.hh"

#include "json_io.hh"
#include "rascal_utility.hh"

#include <tuple>

namespace rascal {

  struct TestData {
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorStrict>;
    TestData() = default;

    void get_ref(const std::string & ref_filename) {
      std::vector<std::uint8_t> ref_data_ubjson;
      internal::read_binary_file(ref_filename, ref_data_ubjson);
      this->ref_data = json::from_ubjson(ref_data_ubjson);
      auto filenames =
          this->ref_data.at("filenames").get<std::vector<std::string>>();
      auto cutoffs = this->ref_data.at("cutoffs").get<std::vector<double>>();

      for (auto && filename : filenames) {
        for (auto && cutoff : cutoffs) {
          // std::cout << filename << " " << cutoff << std::endl;
          json parameters;
          json structure{{"filename", filename}};
          json adaptors;
          json ad1{
              {"name", "AdaptorNeighbourList"},
              {"initialization_arguments",
               {{"cutoff", cutoff},
                {"consider_ghost_neighbours", consider_ghost_neighbours}}}};
          json ad2{{"name", "AdaptorStrict"},
                   {"initialization_arguments", {{"cutoff", cutoff}}}};
          adaptors.emplace_back(ad1);
          adaptors.emplace_back(ad2);

          parameters["structure"] = structure;
          parameters["adaptors"] = adaptors;

          this->factory_args.emplace_back(parameters);
        }
      }
    }

    ~TestData() = default;

    const bool consider_ghost_neighbours{false};
    json ref_data{};
    json factory_args{};
  };

  struct MultipleStructureSOAP : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSOAP() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~MultipleStructureSOAP() = default;

    std::vector<json> hypers{};

    std::vector<json> fc_hypers{
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "A"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "A"}}}},
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "A"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "A"}}}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.2}, {"unit", "A"}}}},
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.4}, {"unit", "A"}}}}};
    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};
    std::vector<json> rep_hypers{{{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 0},
                                  {"soap_type", "RadialSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}},
                                 {{"max_radial", 6},
                                  {"max_angular", 6},
                                  {"soap_type", "PowerSpectrum"},
                                  {"normalize", true}}};
  };

  struct SOAPTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;
    SOAPTestData() : Parent{} { this->get_ref(this->ref_filename); }
    ~SOAPTestData() = default;
    std::string ref_filename{"reference_data/soap_reference.ubjson"};
  };

  struct MultipleStructureSphericalExpansion
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSphericalExpansion() : Parent{} {
      for (auto & ri_hyp : this->radial_contribution_hypers) {
        for (auto & fc_hyp : this->fc_hypers) {
          for (auto & sig_hyp : this->density_hypers) {
            for (auto & rep_hyp : this->rep_hypers) {
              rep_hyp["cutoff_function"] = fc_hyp;
              rep_hyp["gaussian_density"] = sig_hyp;
              rep_hyp["radial_contribution"] = ri_hyp;
              this->hypers.push_back(rep_hyp);
            }
          }
        }
      }
    };
    ~MultipleStructureSphericalExpansion() = default;

    std::vector<json> hypers{};

    std::vector<json> fc_hypers{
        {{"type", "Cosine"},
         {"cutoff", {{"value", 3.0}, {"unit", "A"}}},
         {"smooth_width", {{"value", 0.5}, {"unit", "A"}}}},
        {{"type", "Cosine"},
         {"cutoff", {{"value", 2.0}, {"unit", "A"}}},
         {"smooth_width", {{"value", 1.0}, {"unit", "A"}}}}};

    std::vector<json> radial_contribution_hypers{{{"type", "GTO"}}};

    std::vector<json> density_hypers{
        {{"type", "Constant"},
         {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}}};

    std::vector<json> rep_hypers{{{"max_radial", 10}, {"max_angular", 8}}};
  };

  struct SphericalExpansionTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    SphericalExpansionTestData() : Parent{} {
      this->get_ref(this->ref_filename);
    }
    ~SphericalExpansionTestData() = default;
    std::string ref_filename{
        "reference_data/spherical_expansion_reference.ubjson"};
  };

  struct MultipleStructureSortedCoulomb
      : MultipleStructureManagerNLStrictFixture {
    using Parent = MultipleStructureManagerNLStrictFixture;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    MultipleStructureSortedCoulomb() : Parent{} {};
    ~MultipleStructureSortedCoulomb() = default;

    std::vector<json> hypers{{{"central_decay", 0.5},
                              {"interaction_cutoff", 10.},
                              {"interaction_decay", 0.5},
                              {"size", 120},
                              {"sorting_algorithm", "distance"}},
                             {{"central_decay", 0.5},
                              {"interaction_cutoff", 10.},
                              {"interaction_decay", 0.5},
                              {"size", 120},
                              {"sorting_algorithm", "row_norm"}}};
  };

  struct SortedCoulombTestData : TestData {
    using Parent = TestData;
    using ManagerTypeHolder_t = typename Parent::ManagerTypeHolder_t;

    SortedCoulombTestData() : Parent{} { this->get_ref(this->ref_filename); }
    ~SortedCoulombTestData() = default;

    // name of the file containing the reference data. it has been generated
    // with the following python code:
    // script/generate_sorted_coulomb_ref_data.py

    const bool consider_ghost_neighbours{false};
    std::string ref_filename{"reference_data/sorted_coulomb_reference.ubjson"};
  };

  template <class BaseFixture, template <class> class RepresentationManager>
  struct RepresentationFixture : MultipleStructureFixture<BaseFixture> {
    using Parent = MultipleStructureFixture<BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = RepresentationManager<Manager_t>;

    RepresentationFixture() : Parent{} {}
    ~RepresentationFixture() = default;

    std::list<Representation_t> representations{};
  };

  /* ---------------------------------------------------------------------- */

}  // namespace rascal

#endif  // TESTS_TEST_REPRESENTATION_MANAGER_HH_
