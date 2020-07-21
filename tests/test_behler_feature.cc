/**
 * @file   test_behler_feature.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  testing Behler-Parrinello G-functions
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
    json make_params() {
      json retval{};
      switch (MySymFunType) {
      case SymmetryFunctionType::Gaussian: {
        retval = json{{"type", "Gaussian"},
                      {"index", 0},
                      {"unit", "eV"},
                      {"params",
                       {{"eta", {{"value", 0.1}, {"unit", "(Å)^(-2)"}}},
                        {"r_s", {{"value", 0.6}, {"unit", "Å"}}}}},
                      {"species", {"Mg", "Si"}},
                      {"r_cut", {{"value", this->r_cut}, {"unit", "Å"}}}};
        break;
      }
      case SymmetryFunctionType::AngularNarrow: {
        retval = json{{"type", "AngularNarrow"},
                      {"index", 1},
                      {"unit", "eV"},
                      {"params",
                       {{"eta", {{"value", 0.1}, {"unit", "(Å)^(-2)"}}},
                        {"zeta", {{"value", 0.6}, {"unit", "-"}}},
                        {"lambda", {{"value", 0.6}, {"unit", "-"}}}}},
                      {"species", {"Mg", "Si", "Si"}},
                      {"r_cut", {{"value", this->r_cut}, {"unit", "Å"}}}};
        break;
      }
      default:
        throw std::runtime_error(
            "This SymmetryFunctionType is not yet implemented");
      }
      return retval;
    }

    BehlerFeatureFixture() : raw_params(make_params()) {}

    using CutFun_t = CutoffFunction<InlCutoffFunctionType::Cosine>;
    const double r_cut{1.42};
    const UnitStyle unit_style{units::metal};
    static constexpr auto SymFunType() { return MySymFunType; }

    std::shared_ptr<CutFun_t> cut_fun{std::make_shared<CutFun_t>(
        unit_style,
        json{{"params", {}}, {"r_cut", {{"value", r_cut}, {"unit", "Å"}}}})};
    json raw_params;
    static constexpr auto Order{SymmetryFunction<MySymFunType>::Order};
    using BehlerFeature_t =
        BehlerFeatureOrderSelector_t<Order, MySymFunType, SymFunTypes...>;
    BehlerFeature_t bf{this->cut_fun, this->unit_style, this->raw_params};
  };

  // list of all tested BehlerFeatures
  using Features =
      boost::mpl::list<BehlerFeatureFixture<SymmetryFunctionType::AngularNarrow,
                                            SymmetryFunctionType::AngularNarrow,
                                            SymmetryFunctionType::Gaussian>,
                       BehlerFeatureFixture<SymmetryFunctionType::Gaussian,
                                            SymmetryFunctionType::AngularNarrow,
                                            SymmetryFunctionType::Gaussian>>;

  BOOST_AUTO_TEST_SUITE(behler_parrinello_feature_tests);
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(constructor_test, Fix, Features, Fix) {}

  // list of all tested BehlerFeatures defined on pairs
  using PairFeatures =
      boost::mpl::list<BehlerFeatureFixture<SymmetryFunctionType::Gaussian,
                                            SymmetryFunctionType::Gaussian>>;
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test, Fix, PairFeatures, Fix) {
    ManagerFixture<StructureManagerLammpsMinimal> manager_fix{};
    auto strict_manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, Fix::r_cut)};
    auto & manager{*strict_manager_ptr};
    manager.update();

    constexpr auto Order{Fix::BehlerFeature_t::SymmetryFunction_t::Order};
    using GVals_t = Property<double, AtomOrder,
                             AdaptorStrict<StructureManagerLammpsMinimal>,
                             nb_distances(Order)>;
    auto G_vals{std::make_shared<GVals_t>(manager)};

    // Yes, the pairs in this manager do not have the correct species, but this
    // doesn't interfere with testing the compute algo
    Fix::bf.template compute<RepeatedSpecies::All, Permutation<Order, 0, 1>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::Not, Permutation<Order, 0, 1>>(
        manager, G_vals);

    auto throw_unknown_species_rep{[&manager, &G_vals, this]() {
      this->bf.template compute<RepeatedSpecies::FirstTwo,
                                Permutation<Order, 0, 1>>(manager, G_vals);
    }};

    BOOST_CHECK_THROW(throw_unknown_species_rep(), std::runtime_error);
  }

  // list of all tested BehlerFeatures defined on triplets
  using TripletFeatures = boost::mpl::list<
      BehlerFeatureFixture<SymmetryFunctionType::AngularNarrow,
                           SymmetryFunctionType::AngularNarrow>>;

  BOOST_AUTO_TEST_CASE(TripFeaturetest) {
    const double r_cut{1.42};
    static constexpr auto Order{TripletOrder};
    json params{{"type", "AngularNarrow"},
                {"index", 1},
                {"unit", "eV"},
                {"params",
                 {{"eta", {{"value", 0.1}, {"unit", "(Å)^(-2)"}}},
                  {"zeta", {{"value", 0.6}, {"unit", "-"}}},
                  {"lambda", {{"value", 0.6}, {"unit", "-"}}}}},
                {"species", {"Mg", "Si", "Si"}},
                {"r_cut", {{"value", r_cut}, {"unit", "Å"}}}};
    const UnitStyle unit_style{units::metal};
    std::shared_ptr<CutoffFunction<InlCutoffFunctionType::Cosine>> cut_fun{
        std::make_shared<CutoffFunction<InlCutoffFunctionType::Cosine>>(
            unit_style, json{{"params", {}},
                             {"r_cut", {{"value", r_cut}, {"unit", "Å"}}}})};
    using BehlerFeature_t =
        BehlerFeatureOrderSelector_t<Order, SymmetryFunctionType::AngularNarrow,
                                     SymmetryFunctionType::AngularNarrow>;
    BehlerFeature_t bf{cut_fun, unit_style, params};
    BehlerFeatureFixture<SymmetryFunctionType::AngularNarrow,
                         SymmetryFunctionType::AngularNarrow>
        bff{};
  }
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test_triplet_const, Fix,
                                   TripletFeatures, Fix) {}
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test_triplet, Fix, TripletFeatures,
                                   Fix) {
    ManagerFixture<StructureManagerLammpsMinimal> manager_fix{};
    auto strict_manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, Fix::r_cut)};
    auto triplet_manager_ptr{
        make_adapted_manager<AdaptorMaxOrder>(strict_manager_ptr)};
    auto & manager{*triplet_manager_ptr};
    manager.update();
    constexpr auto Order{Fix::BehlerFeature_t::SymmetryFunction_t::Order};
    using GVals_t =
        Property<double, AtomOrder,
                 AdaptorMaxOrder<AdaptorStrict<StructureManagerLammpsMinimal>>>;
    auto G_vals{std::make_shared<GVals_t>(manager)};

    int pair_counter{0};
    for (auto && atom : manager) {
      for (auto && pair : atom.pairs()) {
        auto pair_offset{pair.get_global_index()};
        std::cout << "pair (" << atom.get_atom_tag() << ", "
                  << pair.get_atom_tag() << "), pair_counter = " << pair_counter
                  << ", pair_offset = " << pair_offset << std::endl;
        ++pair_counter;
      }
    }
    // Yes, the triplets in this manager do not have the correct species, but
    // this doesn't interfere with testing the compute algo
    Fix::bf.template compute<RepeatedSpecies::All, Permutation<Order, 0, 1, 2>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::Not, Permutation<Order, 0, 1, 2>>(
        manager, G_vals);

    // auto throw_unknown_species_rep{[&manager, &G_vals, this]() {
    //   this->bf.template compute<RepeatedSpecies::FirstTwo,
    //                             Permutation<Order, 2, 0, 1>>(manager,
    //                             G_vals);
    // }};

    // BOOST_CHECK_NO_THROW(throw_unknown_species_rep());
  }

  /*--------------------------------------------------------------------------*/
  /**
   * Tests evaluation of a pair behler feature including permutation, and
   * compares results to independently computed reference values on a half
   * neighbour list
   */
  using GaussianSymFun = BehlerFeatureFixture<SymmetryFunctionType::Gaussian,
                                              SymmetryFunctionType::Gaussian>;

  BOOST_FIXTURE_TEST_CASE(pairmutation_test, GaussianSymFun) {
    ManagerFixture<StructureManagerLammpsMinimal> manager_fix{};
    auto manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, this->r_cut)};
    auto & manager{*manager_ptr};
    manager.update();
    using GVals_t =
        Property<double, AtomOrder, AdaptorStrict<StructureManagerLammpsMinimal>>;

    using dGSelfVals_t =
        Property<double, AtomOrder, AdaptorStrict<StructureManagerLammpsMinimal>,
                 ThreeD>;
    using dGOtherVals_t =
        Property<double, PairOrder, AdaptorStrict<StructureManagerLammpsMinimal>,
                 ThreeD, 2>;

    // results without permutation
    auto G01_vals{std::make_shared<GVals_t>(manager)};
    // results with permutation
    auto G10_vals{std::make_shared<GVals_t>(manager)};
    // results with equal species
    auto G11_vals{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto G01_vals2{std::make_shared<GVals_t>(manager)};
    // results with permutation
    auto G10_vals2{std::make_shared<GVals_t>(manager)};
    // results with equal species
    auto G11_vals2{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dGSelf01_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    // results with permutation
    auto dGSelf10_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    // results with equal species
    auto dGSelf11_derivatives{std::make_shared<dGSelfVals_t>(manager)};

    // results with derivative without permutation
    auto dGOther01_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with permutation
    auto dGOther10_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with equal species
    auto dGOther11_derivatives{std::make_shared<dGOtherVals_t>(manager)};

    // manual without permutation
    auto G01_ref{std::make_shared<GVals_t>(manager)};
    // manual with permutation
    auto G10_ref{std::make_shared<GVals_t>(manager)};
    // manual with equal species
    auto G11_ref{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dGSelf01_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    // results with permutation
    auto dGSelf10_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    // results with equal species
    auto dGSelf11_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};

    // results with derivative without permutation
    auto dGOther01_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with permutation
    auto dGOther10_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with equal species
    auto dGOther11_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};

    // calculate all behler feature values: symmetry and cutoff function
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 0, 1>>(
        manager, G01_vals);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 1, 0>>(
        manager, G10_vals);
    this->bf.template compute<RepeatedSpecies::All, Permutation<2, 1, 0>>(
        manager, G11_vals);

    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 0, 1>>(
        manager, G01_vals2, dGSelf01_derivatives, dGOther01_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 1, 0>>(
        manager, G10_vals2, dGSelf10_derivatives, dGOther10_derivatives);
    this->bf.template compute<RepeatedSpecies::All, Permutation<2, 1, 0>>(
        manager, G11_vals2, dGSelf11_derivatives, dGOther11_derivatives);

    const double eta{this->raw_params.at("params")
                         .at("eta")
                         .at("value")
                         .template get<double>()};
    const double r_s{this->raw_params.at("params")
                         .at("r_s")
                         .at("value")
                         .template get<double>()};

    // calculate reference values, directly here
    G01_ref->resize();
    G10_ref->resize();
    G11_ref->resize();

    dGSelf01_ref_derivatives->resize();
    dGSelf10_ref_derivatives->resize();
    dGSelf11_ref_derivatives->resize();

    dGOther01_ref_derivatives->resize();
    dGOther10_ref_derivatives->resize();
    dGOther11_ref_derivatives->resize();

    for (auto && atom : manager) {
      for (auto && pair : atom.pairs()) {
        double r_ij{manager.get_distance(pair)};
        auto & dir_ij{manager.get_direction_vector(pair)};
        // '_c' for cutoff, '_s' for symmetry function
        double f_c{.5 * (std::cos(math::PI * r_ij / this->r_cut) + 1)};
        double f_s{std::exp(-eta * (r_ij - r_s) * (r_ij - r_s))};
        double G_incr{f_s * f_c};
        G01_ref->operator[](atom) += G_incr;
        G10_ref->operator[](pair) += G_incr;
        G11_ref->operator[](atom) += G_incr;
        G11_ref->operator[](pair) += G_incr;

        double df_c{-.5 * (math::PI * r_ij / this->r_cut) *
                    std::sin(math::PI * r_ij / this->r_cut)};

        double df_s{-2 * eta * (r_ij - r_s) * f_s};

        auto && dG_incr{(df_s * f_c + f_s * df_c) * dir_ij};

        dGSelf01_ref_derivatives->operator[](atom) += dG_incr;
        dGSelf10_ref_derivatives->operator[](pair) -= dG_incr;

        dGSelf11_ref_derivatives->operator[](atom) += dG_incr;
        dGSelf11_ref_derivatives->operator[](pair) -= dG_incr;

        dGOther01_ref_derivatives->operator[](pair).col(0) -= dG_incr;

        dGOther10_ref_derivatives->operator[](pair).col(1) -= -dG_incr;

        dGOther11_ref_derivatives->operator[](pair).col(0) -= dG_incr;
        dGOther11_ref_derivatives->operator[](pair).col(1) -= -dG_incr;
      }
    }

    // Test eval of just the values
    double rel_error{(G01_vals->eigen() - G01_ref->eigen()).norm() /
                     G01_ref->eigen().norm()};
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (G10_vals->eigen() - G10_ref->eigen()).norm() / G10_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (G11_vals->eigen() - G11_ref->eigen()).norm() / G11_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    // Test eval of values when both values and derivatives are computed
    rel_error = (G01_vals2->eigen() - G01_ref->eigen()).norm() /
                G01_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G10_vals2->eigen() - G10_ref->eigen()).norm() /
                G10_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G11_vals2->eigen() - G11_ref->eigen()).norm() /
                G11_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    // Test eval of derivatives
    rel_error =
        (dGSelf01_derivatives->eigen() - dGSelf01_ref_derivatives->eigen())
            .norm() /
        dGSelf01_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGSelf10_derivatives->eigen() - dGSelf10_ref_derivatives->eigen())
            .norm() /
        dGSelf10_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGSelf11_derivatives->eigen() - dGSelf11_ref_derivatives->eigen())
            .norm() /
        dGSelf11_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGOther01_derivatives->eigen() - dGOther01_ref_derivatives->eigen())
            .norm() /
        dGOther01_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGOther10_derivatives->eigen() - dGOther10_ref_derivatives->eigen())
            .norm() /
        dGOther10_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGOther11_derivatives->eigen() - dGOther11_ref_derivatives->eigen())
            .norm() /
        dGOther11_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
  }

  /*--------------------------------------------------------------------------*/
  /**
   * Tests evaluation of a triplet behler feature including permutation, and
   * compares results to independently computed reference values on a half
   * neighbour list
   */
  using TripletSymFunNarrow =
      BehlerFeatureFixture<SymmetryFunctionType::AngularNarrow,
                           SymmetryFunctionType::AngularNarrow>;

  BOOST_FIXTURE_TEST_CASE(tripletmutation_test, TripletSymFunNarrow) {
    ManagerFixture<StructureManagerLammpsMinimal> manager_fix{};
    auto strict_manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, this->r_cut)};
    auto triplet_manager_ptr{
        make_adapted_manager<AdaptorMaxOrder>(strict_manager_ptr)};
    auto & manager{*triplet_manager_ptr};
    manager.update();

    using TripletManager_t = AdaptorMaxOrder<
        AdaptorStrict<StructureManagerLammpsMinimal>>;

    using GVals_t = Property<double, AtomOrder, TripletManager_t>;

    using dGSelfVals_t = Property<double, AtomOrder, TripletManager_t, ThreeD>;
    using dGOtherVals_t =
        Property<double, PairOrder, TripletManager_t, ThreeD, 2>;


    // for (const auto RepSpecies : AllRepSpecies) {
    //   auto orderings{
    //       Permutation<3, 0, 1, 2>::get_triplet_orderings<RepSpecies>()};
    // }

    // auto orderings_not{
    //     Permutation<3, 0, 1,
    //                 2>::template
    //                 get_triplet_orderings<RepeatedSpecies::Not>()};
    // auto orderings_first_two{
    //     Permutation<3, 0, 1,
    //                 2>::template
    //                 get_triplet_orderings<RepeatedSpecies::FirstTwo>()};
    // auto orderings_outer_two{
    //     Permutation<3, 0, 1,
    //                 2>::template
    //                 get_triplet_orderings<RepeatedSpecies::OuterTwo>()};
    // auto orderings_all{
    //     Permutation<3, 0, 1,
    //                 2>::template
    //                 get_triplet_orderings<RepeatedSpecies::All>()};

    /**
     * The used manager has not the correct species combination for the given
     * symmetry function, but we permute atoms anyways to test numerics.
     */

    // results without permutation
    auto G012_vals{std::make_shared<GVals_t>(manager)};
    // results with permutations
    auto G120_vals{std::make_shared<GVals_t>(manager)};
    auto G201_vals{std::make_shared<GVals_t>(manager)};

    // results with equal species (??) ::RepeatedSpecies::All
    auto GAAA_vals{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto G012_vals2{std::make_shared<GVals_t>(manager)};
    // results with permutation
    auto G120_vals2{std::make_shared<GVals_t>(manager)};
    auto G201_vals2{std::make_shared<GVals_t>(manager)};

    // results with equal species (??)
    auto GAAA_vals2{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dGSelf012_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto dGOther012_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with permutation
    auto dGSelf120_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto dGOther120_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    auto dGSelf201_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto dGOther201_derivatives{std::make_shared<dGOtherVals_t>(manager)};

    // results with equal species (??)
    auto dGSelfAAA_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto dGOtherAAA_derivatives{std::make_shared<dGOtherVals_t>(manager)};

    // manual without permutation
    auto G012_ref{std::make_shared<GVals_t>(manager)};
    // manual with permutation
    auto G120_ref{std::make_shared<GVals_t>(manager)};
    auto G201_ref{std::make_shared<GVals_t>(manager)};

    // manual with equal species (??)
    auto GAAA_ref{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dGSelf012_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto
    dGOther012_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with permutation
    auto dGSelf120_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto
    dGOther120_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    auto dGSelf201_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto
    dGOther201_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};
    // results with equal species
    auto dGSelfAAA_ref_derivatives{std::make_shared<dGSelfVals_t>(manager)};
    auto
    dGOtherAAA_ref_derivatives{std::make_shared<dGOtherVals_t>(manager)};

    // calculate all behler feature values: symmetry and cutoff function
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 0, 1, 2>>(
        manager, G012_vals);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 1, 2, 0>>(
        manager, G120_vals);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 2, 0, 1>>(
        manager, G201_vals);
    this->bf.template compute<RepeatedSpecies::All, Permutation<3, 0, 1, 2>>(
        manager, GAAA_vals);

    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 0, 1, 2>>(
        manager, G012_vals2, dGSelf012_derivatives, dGOther012_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 1, 2, 0>>(
        manager, G120_vals2, dGSelf120_derivatives, dGOther120_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 2, 0, 1>>(
        manager, G201_vals2, dGSelf201_derivatives, dGOther201_derivatives);
    this->bf.template compute<RepeatedSpecies::All, Permutation<3, 0, 1, 2>>(
        manager, GAAA_vals2, dGSelfAAA_derivatives, dGOtherAAA_derivatives);

    const double eta{this->raw_params.at("params")
                         .at("eta")
                         .at("value")
                         .template get<double>()};
    const double zeta{this->raw_params.at("params")
                          .at("zeta")
                          .at("value")
                          .template get<double>()};
    const double lambda{this->raw_params.at("params")
                            .at("lambda")
                            .at("value")
                            .template get<double>()};

    // calculate reference values, directly here
    G012_ref->resize();
    G120_ref->resize();
    G201_ref->resize();
    GAAA_ref->resize();

    dGSelf012_ref_derivatives->resize();
    dGSelf120_ref_derivatives->resize();
    dGSelf201_ref_derivatives->resize();
    dGSelfAAA_ref_derivatives->resize();

    dGOther012_ref_derivatives->resize();
    dGOther120_ref_derivatives->resize();
    dGOther201_ref_derivatives->resize();
    dGOtherAAA_ref_derivatives->resize();

    auto & triplet_distances{manager.get_triplet_distance()};
    auto & triplet_directions{manager.get_triplet_direction_vectors()};
    auto & cos_angles{get_cos_angles(manager)};
    auto & neigh_to_i_atom{
        manager.template get_neighbours_to_i_atoms<TripletOrder>()};
    auto && pairs_container{
        manager.template get_sub_clusters<PairOrder, TripletOrder>()};

    for (auto && atom : manager) {
      for (auto && triplet : atom.triplets()) {
        auto && atom_cluster_indices{neigh_to_i_atom[triplet]};

        auto && trip_dist{triplet_distances[triplet]};
        auto && trip_cos{cos_angles[triplet]};

        auto && r_ij{trip_dist[0]};
        auto && r_jk{trip_dist[1]};
        auto && r_ki{trip_dist[2]};

        double prefactor{math::pow(2., 1 - zeta)};

        double f_cut_ij{.5 * (std::cos(math::PI * r_ij / this->r_cut) + 1)};
        double f_cut_jk{.5 * (std::cos(math::PI * r_jk / this->r_cut) + 1)};
        double f_cut_ki{.5 * (std::cos(math::PI * r_ki / this->r_cut) + 1)};

        auto && cos_theta_ijk{trip_cos[0]};
        auto && cos_theta_jki{trip_cos[1]};
        auto && cos_theta_kij{trip_cos[2]};

        double f_cut_val{f_cut_ij * f_cut_jk * f_cut_ki};

        double lam_cos_theta_1_ijk{1. + lambda * cos_theta_ijk};
        double lam_cos_theta_1_jki{1. + lambda * cos_theta_jki};
        double lam_cos_theta_1_kij{1. + lambda * cos_theta_kij};

        double angular_contrib_ijk{math::pow(lam_cos_theta_1_ijk, zeta)};
        double angular_contrib_jki{math::pow(lam_cos_theta_1_jki, zeta)};
        double angular_contrib_kij{math::pow(lam_cos_theta_1_kij, zeta)};

        double exp_contrib{exp(-eta * trip_dist.squaredNorm())};

        double f_sym_i{prefactor * angular_contrib_ijk * exp_contrib};
        double f_sym_j{prefactor * angular_contrib_jki * exp_contrib};
        double f_sym_k{prefactor * angular_contrib_kij * exp_contrib};

        double G_incr_i{f_sym_i * f_cut_val};
        double G_incr_j{f_sym_j * f_cut_val};
        double G_incr_k{f_sym_k * f_cut_val};

        auto && atom_i{manager[atom_cluster_indices[0]]};
        auto && atom_j{manager[atom_cluster_indices[1]]};
        auto && atom_k{manager[atom_cluster_indices[2]]};

        auto && triplet_pairs{pairs_container[triplet]};
        auto && pair_ij{triplet_pairs[0]};
        auto && pair_jk{triplet_pairs[1]};
        auto && pair_ki{triplet_pairs[2]};

        G012_ref->operator[](atom_i) += G_incr_i;
        G120_ref->operator[](atom_j) += G_incr_j;
        G201_ref->operator[](atom_k) += G_incr_k;
        GAAA_ref->operator[](atom_i) += G_incr_i + G_incr_j + G_incr_k;
        GAAA_ref->operator[](atom_j) += G_incr_i + G_incr_j + G_incr_k;
        GAAA_ref->operator[](atom_k) += G_incr_i + G_incr_j + G_incr_k;

        double cutoff_der_ij{-.5 * (math::PI * r_ij / this->r_cut) *
                             std::sin(math::PI * r_ij / this->r_cut)};
        double cutoff_der_jk{-.5 * (math::PI * r_jk / this->r_cut) *
                             std::sin(math::PI * r_jk / this->r_cut)};
        double cutoff_der_ki{-.5 * (math::PI * r_ki / this->r_cut) *
                             std::sin(math::PI * r_ki / this->r_cut)};

        double d_f_cut_ij{cutoff_der_ij * f_cut_jk * f_cut_ki};
        double d_f_cut_jk{cutoff_der_jk * f_cut_ki * f_cut_ij};
        double d_f_cut_ki{cutoff_der_ki * f_cut_ij * f_cut_jk};

        // components of derivatives for symmetry function centered on
        // atom 1
        double d_f_sym_i_ij{
            -2 * eta * r_ij * f_sym_i +
            zeta * (lambda / r_ki - lambda / r_ij * cos_theta_ijk) * f_sym_i
            /
                lam_cos_theta_1_ijk};
        double d_f_sym_i_ki{
            -2 * eta * r_ki * f_sym_i +
            zeta * (lambda / r_ij - lambda / r_ki * cos_theta_ijk) * f_sym_i
            /
                lam_cos_theta_1_ijk};
        double d_f_sym_i_jk{-lambda * r_jk * zeta * f_sym_i /
                                (r_ij * r_ki * lam_cos_theta_1_ijk) -
                            2 * eta * r_jk * f_sym_i};

        // components of derivatives for symmetry function centered on
        // atom 2
        double d_f_sym_j_jk{
            -2 * eta * r_jk * f_sym_j +
            zeta * (lambda / r_ij - lambda / r_jk * cos_theta_jki) * f_sym_j
            /
                lam_cos_theta_1_jki};
        double d_f_sym_j_ij{
            -2 * eta * r_ij * f_sym_j +
            zeta * (lambda / r_jk - lambda / r_ij * cos_theta_jki) * f_sym_j
            /
                lam_cos_theta_1_jki};
        double d_f_sym_j_ki{-lambda * r_ki * zeta * f_sym_j /
                                (r_jk * r_ij * lam_cos_theta_1_jki) -
                            2 * eta * r_ki * f_sym_j};

        // components of derivatives for symmetry function centered on
        // atom 3
        double d_f_sym_k_ki{
            -2 * eta * r_ki * f_sym_k +
            zeta * (lambda / r_jk - lambda / r_ki * cos_theta_kij) * f_sym_k
            /
                lam_cos_theta_1_kij};
        double d_f_sym_k_jk{
            -2 * eta * r_jk * f_sym_k +
            zeta * (lambda / r_ki - lambda / r_jk * cos_theta_kij) * f_sym_k
            /
                lam_cos_theta_1_kij};
        double d_f_sym_k_ij{-lambda * r_ij * zeta * f_sym_k /
                                (r_ki * r_jk * lam_cos_theta_1_kij) -
                            2 * eta * r_ij * f_sym_k};

        auto && dir_ij{triplet_directions[triplet].col(0)};
        auto && dir_jk{triplet_directions[triplet].col(1)};
        auto && dir_ki{triplet_directions[triplet].col(2)};

        // dG_i/d_direction
        auto && dG_incr_i_ij{
            (d_f_sym_i_ij * f_cut_val + f_sym_i * d_f_cut_ij) * dir_ij};
        auto && dG_incr_i_jk{
            (d_f_sym_i_jk * f_cut_val + f_sym_i * d_f_cut_jk) * dir_jk};
        auto && dG_incr_i_kj{-dG_incr_i_jk};
        auto && dG_incr_i_ki{
            (d_f_sym_i_ki * f_cut_val + f_sym_i * d_f_cut_ki) * dir_ki};
        auto && dG_incr_i_ik{-dG_incr_i_ki};
        // dG_j/d_direction
        auto && dG_incr_j_ij{
            (d_f_sym_j_ij * f_cut_val + f_sym_j * d_f_cut_ij) * dir_ij};
        auto && dG_incr_j_ji{-dG_incr_j_ij};
        auto && dG_incr_j_jk{
            (d_f_sym_j_jk * f_cut_val + f_sym_j * d_f_cut_jk) * dir_jk};

        auto && dG_incr_j_ki{
            (d_f_sym_j_ki * f_cut_val + f_sym_j * d_f_cut_ki) * dir_ki};
        auto && dG_incr_j_ik{-dG_incr_j_ki};
        // dG_k/d_direction
        auto && dG_incr_k_ij{
            (d_f_sym_k_ij * f_cut_val + f_sym_k * d_f_cut_ij) * dir_ij};
        auto && dG_incr_k_jk{
            (d_f_sym_k_jk * f_cut_val + f_sym_k * d_f_cut_jk) * dir_jk};
        auto && dG_incr_k_kj{-dG_incr_k_jk};
        auto && dG_incr_k_ki{
            (d_f_sym_k_ki * f_cut_val + f_sym_k * d_f_cut_ki) * dir_ki};
        auto && dG_incr_k_ik{-dG_incr_k_ki};

        dGSelf012_ref_derivatives->operator[](atom_i) +=
            dG_incr_i_ij + dG_incr_i_ik;
        dGOther012_ref_derivatives->operator[](pair_ij).col(0) -=
            dG_incr_i_ij + dG_incr_i_kj;
        dGOther012_ref_derivatives->operator[](pair_ki).col(1) -=
            dG_incr_i_ik + dG_incr_i_jk;

        dGSelf120_ref_derivatives->operator[](atom_j) +=
            dG_incr_j_ji + dG_incr_j_jk;
        dGOther120_ref_derivatives->operator[](pair_jk).col(0) -=
            dG_incr_j_jk + dG_incr_j_ik;
        dGOther120_ref_derivatives->operator[](pair_ij).col(1) -=
            dG_incr_j_ji + dG_incr_j_ki;

        dGSelf201_ref_derivatives->operator[](atom_k) +=
            -dG_incr_k_ki + dG_incr_k_kj;
        dGOther201_ref_derivatives->operator[](pair_jk).col(1) -=
            dG_incr_k_kj + dG_incr_k_ij;
        dGOther201_ref_derivatives->operator[](pair_ki).col(0) -=
            dG_incr_k_ik + dG_incr_k_jk;

        /* ------------------------------------------------------------------
        */ dGSelfAAA_ref_derivatives->operator[](atom_i) +=
            dG_incr_i_ij + dG_incr_i_ik;
        dGOtherAAA_ref_derivatives->operator[](pair_ij).col(0) -=
            dG_incr_i_ij + dG_incr_i_kj;
        dGOtherAAA_ref_derivatives->operator[](pair_ki).col(1) -=
            dG_incr_i_ik + dG_incr_i_jk;

        dGSelfAAA_ref_derivatives->operator[](atom_j) +=
            dG_incr_j_ji + dG_incr_j_jk;
        dGOtherAAA_ref_derivatives->operator[](pair_jk).col(0) -=
            dG_incr_j_jk + dG_incr_j_ik;
        dGOtherAAA_ref_derivatives->operator[](pair_ij).col(1) -=
            dG_incr_j_ji + dG_incr_j_ki;

        dGSelfAAA_ref_derivatives->operator[](atom_k) +=
            -dG_incr_k_ki + dG_incr_k_kj;
        dGOtherAAA_ref_derivatives->operator[](pair_jk).col(1) -=
            dG_incr_k_kj + dG_incr_k_ij;
        dGOtherAAA_ref_derivatives->operator[](pair_ki).col(0) -=
            dG_incr_k_ik + dG_incr_k_jk;
      }
    }

    // Test eval of just the values
    std::cout << "G012_vals->eigen(): " << std::endl
              << G012_vals->eigen() << std::endl;
    std::cout << "G012_ref->eigen():  " << std::endl
              << G012_ref->eigen() << std::endl;

    double rel_error{(G012_vals->eigen() - G012_ref->eigen()).norm() /
                     G012_ref->eigen().norm()};
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G120_vals->eigen() - G120_ref->eigen()).norm() /
                G120_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G201_vals->eigen() - G201_ref->eigen()).norm() /
                G201_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (GAAA_vals->eigen() - GAAA_ref->eigen()).norm() /
                GAAA_ref->eigen().norm();
    std::cout << "GAAA_vals->eigen(): " << std::endl
              << GAAA_vals->eigen() << std::endl;
    std::cout << "GAAA_ref->eigen():  " << std::endl
              << GAAA_ref->eigen() << std::endl;
    BOOST_CHECK_EQUAL(rel_error, 0);

    // Test eval of values when both values and derivatives are
    // computed
    rel_error = (G012_vals2->eigen() - G012_ref->eigen()).norm() /
                G012_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G120_vals2->eigen() - G120_ref->eigen()).norm() /
                G120_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (G201_vals2->eigen() - G201_ref->eigen()).norm() /
                G201_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error = (GAAA_vals2->eigen() - GAAA_ref->eigen()).norm() /
                GAAA_ref->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    // Test eval of derivatives
    rel_error =
        (dGSelf012_derivatives->eigen() -
        dGSelf012_ref_derivatives->eigen()).norm() /
        dGSelf012_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
    rel_error =
        (dGOther012_derivatives->eigen() -
        dGOther012_ref_derivatives->eigen()).norm() /
        dGOther012_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGSelf120_derivatives->eigen() -
        dGSelf120_ref_derivatives->eigen()).norm() /
        dGSelf120_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
    rel_error =
        (dGOther120_derivatives->eigen() -
        dGOther120_ref_derivatives->eigen()).norm() /
        dGOther120_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGSelf201_derivatives->eigen() -
        dGSelf201_ref_derivatives->eigen()).norm() /
        dGSelf201_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
    rel_error =
        (dGOther201_derivatives->eigen() -
        dGOther201_ref_derivatives->eigen()).norm() /
        dGOther201_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGSelfAAA_derivatives->eigen() - dGSelfAAA_ref_derivatives->eigen())
            .norm() /
        dGSelfAAA_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
    rel_error =
        (dGOtherAAA_derivatives->eigen() -
        dGOtherAAA_ref_derivatives->eigen())
            .norm() /
        dGOtherAAA_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

   } 

  /**
   * todo(markus) same evaluation and permutation test as above but
   * for angular wide
   */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
