/**
 * file   test_properties.hh
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  tests for cluster-related properties
 *
 * @section LICENSE
 *
 * Copyright  2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TESTS_TEST_PROPERTIES_HH_
#define TESTS_TEST_PROPERTIES_HH_

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"

#include "structure_managers/property.hh"
#include "structure_managers/property_block_sparse.hh"

#include <random>
#include <set>

namespace rascal {

  /* These Fixtures allow flexible stacking of different adaptors for templated
   * tests. The manager is initialized in the initialization list of the
   * constructor. This is required for the usage with the PropertyFixtures,
   * since references have to be initialized in the intialization list of the
   * constructor and every kind of Property object has a reference to its
   * corresponding ManagerImplementation. For an usage example, see in this
   * file the PropertyFixtures (e.g. AtomPropertyFixture).
   */
  struct StructureManagerCentersStackFixture {
    using Manager_t = StructureManagerCenters;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using ManagerTypeHolder_t =
        StructureManagerTypeHolder<StructureManagerCenters>;
    const std::string filename{
        "reference_data/CaCrP2O7_mvc-11955_symmetrized.json"};
    StructureManagerCentersStackFixture()
        : manager{make_structure_manager<StructureManagerCenters>()} {
      manager->update(filename);
    }
    ManagerPtr_t manager;
  };

  template <class StackFixture, bool consider_ghost_neighbours_>
  struct AdaptorNeighbourListStackFixture : StackFixture {
    using Parent = StackFixture;
    using Manager_t = AdaptorNeighbourList<typename Parent::Manager_t>;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    const double consider_ghost_neighbours{consider_ghost_neighbours_};
    const double cutoff{1.};
    AdaptorNeighbourListStackFixture()
        : manager{make_adapted_manager<AdaptorNeighbourList>(
              StackFixture::manager, cutoff, consider_ghost_neighbours)} {
      manager->update();
    }
    ManagerPtr_t manager;
  };

  template <class StackFixture>
  struct AdaptorHalfListStackFixture : StackFixture {
    using Parent = StackFixture;
    using Manager_t = AdaptorHalfList<typename Parent::Manager_t>;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    AdaptorHalfListStackFixture()
        : manager{
              make_adapted_manager<AdaptorHalfList>(StackFixture::manager)} {
      manager->update();
    }

    ManagerPtr_t manager;
  };

  template <class StackFixture>
  struct AdaptorStrictStackFixture : StackFixture {
    using Parent = StackFixture;
    using Manager_t = AdaptorStrict<typename Parent::Manager_t>;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    const double cutoff{2.};
    AdaptorStrictStackFixture()
        : manager{make_adapted_manager<AdaptorStrict>(StackFixture::manager,
                                                      cutoff)} {
      manager->update();
    }

    ManagerPtr_t manager;
  };

  template <class StackFixture>
  struct AdaptorMaxOrderStackFixture : StackFixture {
    using Parent = StackFixture;
    using Manager_t = AdaptorMaxOrder<typename Parent::Manager_t>;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    AdaptorMaxOrderStackFixture()
        : manager{
              make_adapted_manager<AdaptorMaxOrder>(StackFixture::manager)} {
      manager->update();
    }

    ManagerPtr_t manager;
  };

  // #BUG8486@(markus), #BUG8486@(felix) I made new Fixtures, Multiple*Fixture
  // does not match this use case because I need the manager already initialized
  // in the initialization list of the constructor, an the old PairFixtures
  // where not versitile enough to allow different stackings. I shadow the
  // manager member variable in each stack such that it is usable in the new
  // stack and replaces it for higher stacks. For the property tests what test
  // data is load into the manager is not important, since we make our own
  // properties.
  /* ---------------------------------------------------------------------- */
  /*
   * A fixture for testing properties. The StackFixture, which it inherits from,
   * provides the ManagerImplementation for the property.
   */

  template <class StackFixture>
  struct AtomPropertyFixture : StackFixture {
    using Parent = StackFixture;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using AtomScalarProperty_t =
        typename Manager_t::template Property_t<double, 1>;
    using AtomVectorProperty_t =
        typename Manager_t::template Property_t<double, 1, 3, 1>;
    using AtomDynamicUnitProperty_t =
        typename Manager_t::template TypedProperty_t<size_t, 1>;
    using AtomDynamicProperty_t =
        typename Manager_t::template TypedProperty_t<double, 1>;
    using SparseAtomScalarProperty_t =
        typename Manager_t::template BlockSparseProperty_t<double, 1>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string atom_vector_property_metadata{"positions"};
    std::string atom_dynamic_vector_unit_property_metadata{
        "arbitrary counters"};
    std::string atom_dynamic_vector_property_metadata{"distances"};

    AtomPropertyFixture()
        : StackFixture{}, atom_scalar_property{*this->manager},
          atom_vector_property{*this->manager, atom_vector_property_metadata},
          atom_dynamic_scalar_property{
              *this->manager, 1, 1, atom_dynamic_vector_unit_property_metadata},
          atom_dynamic_vector_unit_property{
              *this->manager, DynSize(), 1,
              atom_dynamic_vector_unit_property_metadata},
          atom_dynamic_vector_property{*this->manager, DynSize(), 1,
                                       atom_dynamic_vector_property_metadata},
          sparse_atom_scalar_property{*this->manager} {}

    AtomScalarProperty_t atom_scalar_property;
    AtomVectorProperty_t atom_vector_property;
    AtomDynamicProperty_t atom_dynamic_scalar_property;
    AtomDynamicUnitProperty_t atom_dynamic_vector_unit_property;
    AtomDynamicProperty_t atom_dynamic_vector_property;
    SparseAtomScalarProperty_t sparse_atom_scalar_property;
  };

  template <class StackFixture>
  struct PairPropertyFixture : AtomPropertyFixture<StackFixture> {
    using Parent = AtomPropertyFixture<StackFixture>;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using PairScalarProperty_t =
        typename Manager_t::template Property_t<double, 2>;

    PairPropertyFixture() : Parent{}, pair_property{*this->manager} {}
    PairScalarProperty_t pair_property;
  };

  template <class StackFixture>
  struct TriplePropertyFixture : PairPropertyFixture<StackFixture> {
    using Parent = PairPropertyFixture<StackFixture>;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using TripleScalarProperty_t =
        typename Manager_t::template Property_t<double, 3>;
    TriplePropertyFixture() : Parent{}, triple_property{*this->manager} {}

    TripleScalarProperty_t triple_property;
  };

  template <bool consider_ghost_neighbours>
  struct ANL_SMC_StackFixture_Helper {
    using type =
        AdaptorNeighbourListStackFixture<StructureManagerCentersStackFixture,
                                         consider_ghost_neighbours>;
  };

  /* Helper struct to concat tuple type. The same functionality using directly
   * a boost list does not work. Therefore we do the workaround with tuples.
   */
  template <typename ta, typename tb>
  struct tuple_cat;

  template <typename... a, typename... b>
  struct tuple_cat<std::tuple<a...>, std::tuple<b...>> {
    typedef std::tuple<a..., b...> type;
  };

  /* Helper struct to unpack the tuple arguments and pack them into a
   * boost::mpl::list.
   */
  template <typename ta>
  struct pack_into_list;

  template <typename... a>
  struct pack_into_list<std::tuple<a...>> {
    using type = boost::mpl::list<a...>;
  };

  /* Helper usings because the direct usage of type with template arguments
   * returns compiler errors.
   */
  using ANLWithGhosts_SMC_StackFixture =
      AdaptorNeighbourListStackFixture<StructureManagerCentersStackFixture,
                                       true>;
  using ANLWithoutGhosts_SMC_StackFixture =
      AdaptorNeighbourListStackFixture<StructureManagerCentersStackFixture,
                                       false>;

  /* The general procedure here is, that the tuple_some_description contains all
   * the types in a tuple, and type_some_description contains the boost list so
   * it can be used for the templated tests of boost. See test_properties.cc for
   * usage of these structs (e.g. atom_vector_property_fixtures_with_ghosts).
   */
  struct OrderOnePropertyBoostList {
    template <bool consider_ghost_neighbours>
    struct OrderTwoFixtureStacksTuple {
      using ANL_SMC_TemplatedStackFixture =
          AdaptorNeighbourListStackFixture<StructureManagerCentersStackFixture,
                                           consider_ghost_neighbours>;
      using type = std::tuple<
          AtomPropertyFixture<ANL_SMC_TemplatedStackFixture>,
          AtomPropertyFixture<
              AdaptorHalfListStackFixture<ANL_SMC_TemplatedStackFixture>>,
          AtomPropertyFixture<
              AdaptorStrictStackFixture<ANL_SMC_TemplatedStackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANL_SMC_TemplatedStackFixture>>>>;
    };
    using tuple_order_3 = std::tuple<
        AtomPropertyFixture<
            AdaptorMaxOrderStackFixture<ANLWithGhosts_SMC_StackFixture>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
        AtomPropertyFixture<
            AdaptorMaxOrderStackFixture<AdaptorStrictStackFixture<
                AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>>>;

    using tuple_order_2_with_ghosts = OrderTwoFixtureStacksTuple<true>::type;
    using tuple_with_ghosts =
        tuple_cat<tuple_order_2_with_ghosts, tuple_order_3>::type;
    using type_with_ghosts = pack_into_list<tuple_with_ghosts>::type;
    using tuple_without_ghosts = OrderTwoFixtureStacksTuple<false>::type;
    using type_without_ghosts = pack_into_list<tuple_without_ghosts>::type;
    using tuple = tuple_cat<tuple_with_ghosts, tuple_without_ghosts>::type;
    using type = pack_into_list<tuple>::type;
  };

  struct OrderTwoPropertyBoostList {
    template <bool consider_ghost_neighbours>
    struct OrderTwoFixtureStacksTuple {
      using ANL_SMC_TemplatedStackFixture =
          AdaptorNeighbourListStackFixture<StructureManagerCentersStackFixture,
                                           consider_ghost_neighbours>;
      using type = std::tuple<
          PairPropertyFixture<ANL_SMC_TemplatedStackFixture>,
          PairPropertyFixture<
              AdaptorHalfListStackFixture<ANL_SMC_TemplatedStackFixture>>,
          PairPropertyFixture<
              AdaptorStrictStackFixture<ANL_SMC_TemplatedStackFixture>>,
          PairPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANL_SMC_TemplatedStackFixture>>>>;
    };

    using tuple_order_3 = std::tuple<
        PairPropertyFixture<
            AdaptorMaxOrderStackFixture<ANLWithGhosts_SMC_StackFixture>>,
        PairPropertyFixture<
            AdaptorMaxOrderStackFixture<AdaptorStrictStackFixture<
                AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>>>;

    using tuple_order_2_with_ghosts = OrderTwoFixtureStacksTuple<true>::type;
    using tuple_with_ghosts =
        tuple_cat<tuple_order_2_with_ghosts, tuple_order_3>::type;
    using type_with_ghosts = pack_into_list<tuple_with_ghosts>::type;

    using tuple_without_ghosts = OrderTwoFixtureStacksTuple<false>::type;
    using type_without_ghosts = pack_into_list<tuple_without_ghosts>::type;
    using tuple = tuple_cat<tuple_without_ghosts, tuple_with_ghosts>::type;
    using type = pack_into_list<tuple>::type;
  };

  struct OrderThreePropertyBoostList {
    using tuple = std::tuple<
        TriplePropertyFixture<
            AdaptorMaxOrderStackFixture<ANLWithGhosts_SMC_StackFixture>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
        TriplePropertyFixture<
            AdaptorMaxOrderStackFixture<AdaptorStrictStackFixture<
                AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>>>;
    using type = pack_into_list<tuple>::type;
  };
}  // namespace rascal

#endif  // TESTS_TEST_PROPERTIES_HH_
