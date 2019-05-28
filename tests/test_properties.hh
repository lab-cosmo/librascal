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
#include "structure_managers/property_block_sparse.hh"

#include <random>
#include <set>


namespace rascal {

  /* ---------------------------------------------------------------------- */
  /*
   * A fixture for testing properties. It is based on the PairFixture, which is
   * basically a fixture which provides a pair neighbour list based on positions
   * which are initialized in the tests.
   */ 

  // #BUG8486@(markus), #BUG8486@(felix) I made new Fixtures, Multiple*Fixture
  // does not match this use case because I need the manager already initialized
  // in the initialization list of the constructor, an the old PairFixtures
  // where not versitile enough to allow different stackings. I shadow the
  // manager member variable in each stack such that it is usable in the new
  // stack and replaces it for higher stacks. For the property tests what test
  // data is load into the manager is not important, since we make our own
  // properties.
  
  template <class StackFixture>
  struct AtomPropertyFixture: StackFixture {
    using Parent = StackFixture;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using AtomScalarProperty_t =
        typename Manager_t::template Property_t<double, 1>;
    using AtomVectorProperty_t =
        typename Manager_t::template Property_t<double, 1, 3, 1>;
    using AtomDynamicProperty_t =
        typename Manager_t::template TypedProperty_t<size_t, 1>;
    using AtomDynamicProperty2_t =
        typename Manager_t::template TypedProperty_t<double, 1>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string atom_property_metadata{"positions"};
    std::string dynamic_property_metadata{"arbitrary counters"};
    std::string dynamic_property2_metadata{"distances"};

    AtomPropertyFixture()
          : StackFixture{}, scalar_atom_property{*this->manager},
          atom_property{*this->manager},
          atom_dynamic_property{*this->manager, 1, 1,
                           dynamic_property_metadata},
          dynamic_property{*this->manager, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->manager, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t atom_dynamic_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  template <class StackFixture>
  struct PairPropertyFixture: StackFixture {
    using Parent = StackFixture;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using AtomScalarProperty_t =
        typename Manager_t::template Property_t<double, 1>;
    using PairScalarProperty_t =
        typename Manager_t::template Property_t<double, 2>;
    using AtomVectorProperty_t =
        typename Manager_t::template Property_t<double, 1, 3, 1>;
    using AtomDynamicProperty_t =
        typename Manager_t::template TypedProperty_t<size_t, 1>;
    using AtomDynamicProperty2_t =
        typename Manager_t::template TypedProperty_t<double, 1>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string atom_property_metadata{"positions"};
    std::string dynamic_property_metadata{"arbitrary counters"};
    std::string dynamic_property2_metadata{"distances"};

    PairPropertyFixture()
          : StackFixture{}, scalar_atom_property{*this->manager}, 
          atom_property{*this->manager, atom_property_metadata},
          pair_property{*this->manager},
          atom_dynamic_property{*this->manager, 1, 1,
                           dynamic_property_metadata},
          dynamic_property{*this->manager, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->manager, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    AtomVectorProperty_t atom_property;
    PairScalarProperty_t pair_property;
    AtomDynamicProperty_t atom_dynamic_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  template <class StackFixture>
  struct TriplePropertyFixture: StackFixture {
    using Parent = StackFixture;
    using Manager_t = typename Parent::Manager_t;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;

    using AtomScalarProperty_t =
        typename Manager_t::template Property_t<double, 1>;
    using PairScalarProperty_t =
        typename Manager_t::template Property_t<double, 2>;
    using TripleScalarProperty_t =
        typename Manager_t::template Property_t<double, 3>;
    using AtomVectorProperty_t =
        typename Manager_t::template Property_t<double, 1, 3, 1>;
    using AtomDynamicProperty_t =
        typename Manager_t::template TypedProperty_t<size_t, 1>;
    using AtomDynamicProperty2_t =
        typename Manager_t::template TypedProperty_t<double, 1>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string atom_property_metadata{"positions"};
    std::string dynamic_property_metadata{"arbitrary counters"};
    std::string dynamic_property2_metadata{"distances"};

    TriplePropertyFixture()
          : StackFixture{}, scalar_atom_property{*this->manager}, 
          atom_property{*this->manager, atom_property_metadata},
          pair_property{*this->manager},
          triple_property{*this->manager},
          dynamic_property{*this->manager, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->manager, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    AtomVectorProperty_t atom_property;
    PairScalarProperty_t pair_property;
    TripleScalarProperty_t triple_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };
  
  template<bool consider_ghost_neighbours>
  struct ANL_SMC_StackFixture_Helper {
    using type = AdaptorNeighbourListStackFixture<
      StructureManagerCentersStackFixture, consider_ghost_neighbours>;
  };

  // Helper types to concat template lists
  template< typename ta, typename tb >
  struct type_cat;

  template< typename ... a, typename ... b >
  struct type_cat< std::tuple< a ... >, std::tuple< b ... > >
    { typedef std::tuple< a ..., b ... > type; };

  // Helper types to uppack the tuple template arguments and pack them into a
  //  boost::mpl::list
  template< typename ta>
  struct pack_into_list;

  template<typename ... a>
  struct pack_into_list< std::tuple < a ... > >
    { using list = boost::mpl::list<a ...>;};

  using ANLWithGhosts_SMC_StackFixture = AdaptorNeighbourListStackFixture<
    StructureManagerCentersStackFixture, true>;
  using ANLWithoutGhosts_SMC_StackFixture = AdaptorNeighbourListStackFixture<
    StructureManagerCentersStackFixture, false>;

  struct CommonStacksBoostList {
    template<bool consider_ghost_neighbours>
    struct  CommonStacksBoostListTemplated {
      using ANL_SMC_TemplatedStackFixture = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, consider_ghost_neighbours>;      
      using type = std::tuple<
        AtomPropertyFixture<ANL_SMC_TemplatedStackFixture>,
        AtomPropertyFixture<AdaptorHalfListStackFixture<
            ANL_SMC_TemplatedStackFixture>>,
        AtomPropertyFixture<AdaptorStrictStackFixture<
            ANL_SMC_TemplatedStackFixture>>,
        AtomPropertyFixture<AdaptorStrictStackFixture<
            AdaptorHalfListStackFixture<ANL_SMC_TemplatedStackFixture>>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            ANLWithGhosts_SMC_StackFixture>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>,
        AtomPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>>
        >;
    };
  
    using type_with_ghosts = CommonStacksBoostListTemplated<true>::type;
    using list_with_ghosts = pack_into_list<type_with_ghosts>::list;
    using type_without_ghosts = CommonStacksBoostListTemplated<false>::type;
    using list_without_ghosts = pack_into_list<type_without_ghosts>::list;
    using type = type_cat<type_with_ghosts, type_without_ghosts>::type;
    using list = pack_into_list<type>::list;
  };

  struct CommonOrderTwoStacksBoostList {
    template<bool consider_ghost_neighbours>
    struct  CommonOrderTwoStacksBoostListTemplated {
      using ANL_SMC_TemplatedStackFixture = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, consider_ghost_neighbours>;      
      using type = std::tuple<
        PairPropertyFixture<ANLWithGhosts_SMC_StackFixture>,
        PairPropertyFixture<AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>,
        PairPropertyFixture<AdaptorStrictStackFixture<
            ANLWithGhosts_SMC_StackFixture>>,
        PairPropertyFixture<AdaptorStrictStackFixture<
            AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
        PairPropertyFixture<AdaptorMaxOrderStackFixture<
            ANLWithGhosts_SMC_StackFixture>>,
        PairPropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>>
        >;
    };
    using type_with_ghosts = CommonOrderTwoStacksBoostListTemplated<true>::type;
    using list_with_ghosts = pack_into_list<type_with_ghosts>::list;

    using type_without_ghosts = CommonOrderTwoStacksBoostListTemplated<false>::type;
    using list_without_ghosts = pack_into_list<type_without_ghosts>::list;
    using type = type_cat<type_with_ghosts, type_without_ghosts>::type;
    using list = pack_into_list<type>::list;
  };

  struct CommonOrderThreeStacksBoostList {
    using type = std::tuple<
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<ANLWithGhosts_SMC_StackFixture>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<AdaptorStrictStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<AdaptorStrictStackFixture<
            AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>>,
        TriplePropertyFixture<AdaptorMaxOrderStackFixture<
            AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
            ANLWithGhosts_SMC_StackFixture>>>>
        >;
    using list = pack_into_list<type>::list;
  };
}  // namespace rascal

#endif  // TESTS_TEST_PROPERTIES_HH_
