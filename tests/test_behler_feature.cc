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
    ManagerFixture<StructureManagerLammps> manager_fix{};
    auto manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, Fix::r_cut)};
    auto & manager{*manager_ptr};
    manager.update();

    constexpr auto Order{Fix::BehlerFeature_t::SymmetryFunction_t::Order};
    using GVals_t =
        Property<double, AtomOrder, AdaptorStrict<StructureManagerLammps>,
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

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(eval_test_triplet, Fix, TripletFeatures,
                                   Fix) {
    ManagerFixture<StructureManagerLammps> manager_fix{};
    auto strict_manager_ptr{
        make_adapted_manager<AdaptorStrict>(manager_fix.manager, Fix::r_cut)};
    auto triplet_manager_ptr{
        make_adapted_manager<AdaptorMaxOrder>(strict_manager_ptr)};
    auto & manager{*triplet_manager_ptr};
    manager.update();
    constexpr auto Order{Fix::BehlerFeature_t::SymmetryFunction_t::Order};
    using GVals_t =
        Property<double, AtomOrder,
                 AdaptorMaxOrder<AdaptorStrict<StructureManagerLammps>>>;
    auto G_vals{std::make_shared<GVals_t>(manager)};

    // Yes, the triplets in this manager do not have the correct species, but
    // this doesn't interfere with testing the compute algo
    Fix::bf.template compute<RepeatedSpecies::All, Permutation<Order, 0, 1, 2>>(
        manager, G_vals);
    Fix::bf.template compute<RepeatedSpecies::Not, Permutation<Order, 0, 1, 2>>(
        manager, G_vals);

    auto throw_unknown_species_rep{[&manager, &G_vals, this]() {
      this->bf.template compute<RepeatedSpecies::FirstTwo,
                                Permutation<Order, 2, 0, 1>>(manager, G_vals);
    }};

    BOOST_CHECK_NO_THROW(throw_unknown_species_rep());
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
    ManagerFixture<StructureManagerLammps> manager_fix{};
    auto half_list_ptr{
        make_adapted_manager<AdaptorHalfList>(manager_fix.manager)};
    auto manager_ptr{
        make_adapted_manager<AdaptorStrict>(half_list_ptr, this->r_cut)};
    auto & manager{*manager_ptr};
    manager.update();
    using GVals_t =
        Property<double, AtomOrder,
                 AdaptorStrict<AdaptorHalfList<StructureManagerLammps>>>;

    constexpr auto Order{SymmetryFunction<SymFunType()>::Order};
    using dGVals_t =
        Property<double, Order,
                 AdaptorStrict<AdaptorHalfList<StructureManagerLammps>>,
                 nb_distances(Order)>;

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
    auto dG01_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with permutation
    auto dG10_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with equal species
    auto dG11_derivatives{std::make_shared<dGVals_t>(manager)};

    // manual without permutation
    auto G01_ref{std::make_shared<GVals_t>(manager)};
    // manual with permutation
    auto G10_ref{std::make_shared<GVals_t>(manager)};
    // manual with equal species
    auto G11_ref{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dG01_ref_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with permutation
    auto dG10_ref_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with equal species
    auto dG11_ref_derivatives{std::make_shared<dGVals_t>(manager)};

    // calculate all behler feature values: symmetry and cutoff function
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 0, 1>>(
        manager, G01_vals);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 1, 0>>(
        manager, G10_vals);
    this->bf.template compute<RepeatedSpecies::All, Permutation<2, 1, 0>>(
        manager, G11_vals);

    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 0, 1>>(
        manager, G01_vals2, dG01_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<2, 1, 0>>(
        manager, G10_vals2, dG10_derivatives);
    this->bf.template compute<RepeatedSpecies::All, Permutation<2, 1, 0>>(
        manager, G11_vals2, dG11_derivatives);

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

    dG01_ref_derivatives->resize();
    dG10_ref_derivatives->resize();
    dG11_ref_derivatives->resize();
    for (auto && atom : manager) {
      for (auto && pair : atom.pairs()) {
        double r_ij{manager.get_distance(pair)};
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

        auto && dG_incr{df_s * f_c + f_s * df_c};

        dG01_ref_derivatives->operator[](pair) = dG_incr;
        dG10_ref_derivatives->operator[](pair) = dG_incr;
        dG11_ref_derivatives->operator[](pair) = dG_incr;
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
        (dG01_derivatives->eigen() - dG01_ref_derivatives->eigen()).norm() /
        dG01_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dG10_derivatives->eigen() - dG10_ref_derivatives->eigen()).norm() /
        dG10_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dG11_derivatives->eigen() - dG11_ref_derivatives->eigen()).norm() /
        dG11_ref_derivatives->eigen().norm();
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
    ManagerFixture<StructureManagerLammps> manager_fix{};
    auto half_list_ptr{
        make_adapted_manager<AdaptorHalfList>(manager_fix.manager)};
    auto strict_manager_ptr{
        make_adapted_manager<AdaptorStrict>(half_list_ptr, this->r_cut)};
    auto triplet_manager_ptr{
        make_adapted_manager<AdaptorMaxOrder>(strict_manager_ptr)};
    auto & manager{*triplet_manager_ptr};
    manager.update();

    using TripletManager_t =
        AdaptorMaxOrder<AdaptorStrict<AdaptorHalfList<StructureManagerLammps>>>;

    using GVals_t = Property<double, AtomOrder, TripletManager_t>;

    // todo(jungestricker): this makes no sense for triplets - property should
    // be Atom, since it is in general not given that all pairs exist constexpr
    const auto Order{SymmetryFunction<SymFunType()>::Order};
    using dGVals_t =
        Property<double, Order, TripletManager_t, nb_distances(Order)>;
    using DerivativeCast_t =
        Eigen::Map<Eigen::Matrix<double, nb_distances(Order), 1>>;

    /**
     * temporary comment (markus): triplets need to test for permutation: 012,
     * 120, 201 as well as RepeatedSpecies:Not RepeatedSpecies::FirstTwo,
     * RepeatedSpecies::SecondTwo, RepeatedSpecies::All
     */

    // get all triplet orderings Not, FirstTwo, OuterTwo, All

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
     * todo(markus) for now assuming no repetition, but permuting anyways
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
    auto dG012_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with permutation
    auto dG120_derivatives{std::make_shared<dGVals_t>(manager)};
    auto dG201_derivatives{std::make_shared<dGVals_t>(manager)};

    // results with equal species (??)
    auto dGAAA_derivatives{std::make_shared<dGVals_t>(manager)};

    // manual without permutation
    auto G012_ref{std::make_shared<GVals_t>(manager)};
    // manual with permutation
    auto G120_ref{std::make_shared<GVals_t>(manager)};
    auto G201_ref{std::make_shared<GVals_t>(manager)};

    // manual with equal species (??)
    auto GAAA_ref{std::make_shared<GVals_t>(manager)};

    // results with derivative without permutation
    auto dG012_ref_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with permutation
    auto dG120_ref_derivatives{std::make_shared<dGVals_t>(manager)};
    auto dG201_ref_derivatives{std::make_shared<dGVals_t>(manager)};
    // results with equal species
    auto dGAAA_ref_derivatives{std::make_shared<dGVals_t>(manager)};

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
        manager, G012_vals2, dG012_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 1, 2, 0>>(
        manager, G120_vals2, dG120_derivatives);
    this->bf.template compute<RepeatedSpecies::Not, Permutation<3, 2, 0, 1>>(
        manager, G201_vals2, dG201_derivatives);
    this->bf.template compute<RepeatedSpecies::All, Permutation<3, 0, 1, 2>>(
        manager, GAAA_vals2, dGAAA_derivatives);

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

    dG012_ref_derivatives->resize();
    dG120_ref_derivatives->resize();
    dG201_ref_derivatives->resize();
    dGAAA_ref_derivatives->resize();

    auto & triplet_distances{manager.get_triplet_distance()};
    auto & cos_angles{get_cos_angles(manager)};
    auto & neigh_to_i_atom{
        manager.template get_neighbours_to_i_atoms<TripletOrder>()};

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

        double f_sym_ij{prefactor * angular_contrib_ijk * exp_contrib};
        double f_sym_jk{prefactor * angular_contrib_jki * exp_contrib};
        double f_sym_ki{prefactor * angular_contrib_kij * exp_contrib};

        double G_incr_ij{f_sym_ij * f_cut_val};
        double G_incr_jk{f_sym_jk * f_cut_val};
        double G_incr_ki{f_sym_ki * f_cut_val};

        auto && atom_i{manager[atom_cluster_indices[0]]};
        auto && atom_j{manager[atom_cluster_indices[1]]};
        auto && atom_k{manager[atom_cluster_indices[2]]};

        G012_ref->operator[](atom_i) += G_incr_ij;
        G120_ref->operator[](atom_j) += G_incr_jk;
        G201_ref->operator[](atom_k) += G_incr_ki;
        GAAA_ref->operator[](atom_i) += G_incr_ij + G_incr_jk + G_incr_ki;
        GAAA_ref->operator[](atom_j) += G_incr_ij + G_incr_jk + G_incr_ki;
        GAAA_ref->operator[](atom_k) += G_incr_ij + G_incr_jk + G_incr_ki;

        double cutoff_der_ij{-.5 * (math::PI * r_ij / this->r_cut) *
                             std::sin(math::PI * r_ij / this->r_cut)};
        double cutoff_der_jk{-.5 * (math::PI * r_jk / this->r_cut) *
                             std::sin(math::PI * r_jk / this->r_cut)};
        double cutoff_der_ki{-.5 * (math::PI * r_ki / this->r_cut) *
                             std::sin(math::PI * r_ki / this->r_cut)};

        double d_f_cut_ij{cutoff_der_ij * f_cut_jk * f_cut_ki};
        double d_f_cut_jk{cutoff_der_jk * f_cut_ki * f_cut_ij};
        double d_f_cut_ki{cutoff_der_ki * f_cut_ij * f_cut_jk};

        // components of derivatives for symmetry function centered on atom 1
        double d_f_sym_i_ij{
            -2 * eta * r_ij * f_sym_ij +
            zeta * (lambda / r_ki - lambda / r_ij * cos_theta_ijk) * f_sym_ij /
                lam_cos_theta_1_ijk};
        double d_f_sym_i_ki{
            -2 * eta * r_ki * f_sym_ij +
            zeta * (lambda / r_ij - lambda / r_ki * cos_theta_ijk) * f_sym_ij /
                lam_cos_theta_1_ijk};
        double d_f_sym_i_jk{-lambda * r_jk * zeta * f_sym_ij /
                                (r_ij * r_ki * lam_cos_theta_1_ijk) -
                            2 * eta * r_jk * f_sym_ij};

        // components of derivatives for symmetry function centered on atom 2
        double d_f_sym_j_jk{
            -2 * eta * r_jk * f_sym_jk +
            zeta * (lambda / r_ij - lambda / r_jk * cos_theta_jki) * f_sym_jk /
                lam_cos_theta_1_jki};
        double d_f_sym_j_ij{
            -2 * eta * r_ij * f_sym_jk +
            zeta * (lambda / r_jk - lambda / r_ij * cos_theta_jki) * f_sym_jk /
                lam_cos_theta_1_jki};
        double d_f_sym_j_ki{-lambda * r_ki * zeta * f_sym_jk /
                                (r_jk * r_ij * lam_cos_theta_1_jki) -
                            2 * eta * r_ki * f_sym_jk};

        // components of derivatives for symmetry function centered on atom 3
        double d_f_sym_k_ki{
            -2 * eta * r_ki * f_sym_ki +
            zeta * (lambda / r_jk - lambda / r_ki * cos_theta_kij) * f_sym_ki /
                lam_cos_theta_1_kij};
        double d_f_sym_k_jk{
            -2 * eta * r_jk * f_sym_ki +
            zeta * (lambda / r_ki - lambda / r_jk * cos_theta_kij) * f_sym_ki /
                lam_cos_theta_1_kij};
        double d_f_sym_k_ij{-lambda * r_ij * zeta * f_sym_ki /
                                (r_ki * r_jk * lam_cos_theta_1_kij) -
                            2 * eta * r_ij * f_sym_ki};

        // dG_i/d_direction
        auto && dG_incr_i_ij{d_f_sym_i_ij * f_cut_val + f_sym_ij * d_f_cut_ij};
        auto && dG_incr_i_jk{d_f_sym_i_jk * f_cut_val + f_sym_ij * d_f_cut_jk};
        auto && dG_incr_i_ki{d_f_sym_i_ki * f_cut_val + f_sym_ij * d_f_cut_ki};
        // dG_j/d_direction
        auto && dG_incr_j_ij{d_f_sym_j_ij * f_cut_val + f_sym_jk * d_f_cut_ij};
        auto && dG_incr_j_jk{d_f_sym_j_jk * f_cut_val + f_sym_jk * d_f_cut_jk};
        auto && dG_incr_j_ki{d_f_sym_j_ki * f_cut_val + f_sym_jk * d_f_cut_ki};
        // dG_k/d_direction
        auto && dG_incr_k_ij{d_f_sym_k_ij * f_cut_val + f_sym_ki * d_f_cut_ij};
        auto && dG_incr_k_jk{d_f_sym_k_jk * f_cut_val + f_sym_ki * d_f_cut_jk};
        auto && dG_incr_k_ki{d_f_sym_k_ki * f_cut_val + f_sym_ki * d_f_cut_ki};

        std::array<double, 3> tmp{dG_incr_i_ij, dG_incr_i_jk, dG_incr_i_ki};
        dG012_ref_derivatives->operator[](triplet) =
            DerivativeCast_t{tmp.data()};

        tmp = {dG_incr_j_ij, dG_incr_j_jk, dG_incr_j_ki};
        dG120_ref_derivatives->operator[](triplet) =
            DerivativeCast_t{tmp.data()};
        tmp = {dG_incr_k_ij, dG_incr_k_jk, dG_incr_k_ki};
        dG201_ref_derivatives->operator[](triplet) =
            DerivativeCast_t{tmp.data()};

        tmp = {dG_incr_i_ij + dG_incr_j_ij + dG_incr_k_ij,
               dG_incr_i_jk + dG_incr_j_jk + dG_incr_k_jk,
               dG_incr_i_ki + dG_incr_j_ki + dG_incr_k_ki};
        dGAAA_ref_derivatives->operator[](triplet) =
            DerivativeCast_t{tmp.data()};
      }
    }

    // Test eval of just the values
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
    BOOST_CHECK_EQUAL(rel_error, 0);

    // Test eval of values when both values and derivatives are computed
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
        (dG012_derivatives->eigen() - dG012_ref_derivatives->eigen()).norm() /
        dG012_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dG120_derivatives->eigen() - dG120_ref_derivatives->eigen()).norm() /
        dG120_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dG201_derivatives->eigen() - dG201_ref_derivatives->eigen()).norm() /
        dG201_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);

    rel_error =
        (dGAAA_derivatives->eigen() - dGAAA_ref_derivatives->eigen()).norm() /
        dGAAA_ref_derivatives->eigen().norm();
    BOOST_CHECK_EQUAL(rel_error, 0);
  }  // namespace rascal

  /**
   * todo(markus) same evaluation and permutation test as above but for angular
   * wide
   */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
