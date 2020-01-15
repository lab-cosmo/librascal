/**
 * @file   test_behler_feature.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  testing Behler-Parinello G-functions
 *
 * Copyright © 2019 Till Junge
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
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

// #include "rascal/representations/behler_feature.hh"
#include "behler_fixtures.hh"
#include "test_structure.hh"

#include "rascal/representations/behler_feature.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/permutation.hh"

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>

namespace rascal {

  template <SymmetryFunctionType MySymFunType,
            SymmetryFunctionType... SymFunTypes>
  struct BehlerFeatureFixture {
    BehlerFeatureFixture() {}
    using CutFun_t = CutoffFunction<InlCutoffFunctionType::Cosine>;
    const double r_cut{1.1};
    const UnitStyle unit_style{units::metal};
    std::shared_ptr<CutFun_t> cut_fun{std::make_shared<CutFun_t>(
        unit_style,
        json{{"params", {}}, {"r_cut", {{"value", r_cut}, {"unit", "Å"}}}})};
    json raw_params{{"type", "Gaussian"},
                    {"index", 0},
                    {"unit", "eV"},
                    {"params",
                     {{"eta", {{"value", 0.1}, {"unit", "(Å)^(-2)"}}},
                      {"r_s", {{"value", 0.6}, {"unit", "Å"}}}}},
                    {"species", {"Mg", "Si"}},
                    {"r_cut", {{"value", 1.1}, {"unit", "Å"}}}};
    BehlerFeature<MySymFunType, SymFunTypes...> bf{cut_fun, unit_style,
                                                   raw_params};
  };

  using Features =
      boost::mpl::list<BehlerFeatureFixture<SymmetryFunctionType::Gaussian,
                                            SymmetryFunctionType::Gaussian>>;

  BOOST_AUTO_TEST_SUITE(behler_parinello_feature_tests);
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Features, Fix) {}

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test, Fix, Features, Fix) {
    ManagerFixture<StructureManagerLammps> manager_fix{};
    auto & manager{
        *make_adapted_manager<AdaptorStrict>(manager_fix.manager, Fix::r_cut)};
    manager.update();

    using GVals_t =
        Property<double, AtomOrder, AdaptorStrict<StructureManagerLammps>>;
    auto G_vals{std::make_shared<GVals_t>(manager)};

    // Yes, the pairs in this manager do not have the correct species, but this
    // doesn't interfere with testing the compute algo
    Fix::bf.template compute<RepeatedSpecies::All, Permutation<2, 0, 1>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::Not, Permutation<2, 0, 1>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::FirstTwo, Permutation<2, 0, 1>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::SecondTwo, Permutation<2, 0, 1>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::OuterTwo, Permutation<2, 0, 1>>(
        manager, G_vals);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
