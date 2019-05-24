/**
 * file   test_properties.cc
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

#include "tests.hh"
#include "test_structure.hh"
#include "test_adaptor.hh"
#include "structure_managers/property_block_sparse.hh"

#include <random>
#include <set>

namespace rascal {
  // TODO(felix) test dynamic sized Property
  BOOST_AUTO_TEST_SUITE(Property_tests);

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

    constexpr static Dim_t DynSize() { return 3; }

    std::string atom_property_metadata{"positions"};
    std::string dynamic_property_metadata{"arbitrary counters"};
    std::string dynamic_property2_metadata{"distances"};

    AtomPropertyFixture()
          : StackFixture{}, scalar_atom_property{*this->manager} {}
     AtomScalarProperty_t scalar_atom_property;
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
          pair_property{*this->manager},
          atom_property{*this->manager, atom_property_metadata},
          dynamic_property{*this->manager, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->manager, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  // the old one is still in usage
  template <class ManagerImplementation>
  struct PropertyFixture: public PairFixture<ManagerImplementation> {
    using Manager_t = AdaptorNeighbourList<ManagerImplementation>;

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

    PropertyFixture(bool consider_ghost_neighbours = false)
        : PairFixture<ManagerImplementation>{consider_ghost_neighbours},
          scalar_atom_property{*this->pair_manager},
          pair_property{*this->pair_manager},
          atom_property{*this->pair_manager, atom_property_metadata},
          dynamic_property{*this->pair_manager, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->pair_manager, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };



  template <class ManagerImplementation>
  struct PropertyFixtureWithGhosts : public PropertyFixture<ManagerImplementation> {
    PropertyFixtureWithGhosts()
        : PropertyFixture<ManagerImplementation>{true} {}
  };
  
  template <class ManagerImplementation>
  struct PropertyFixtureStrict: public PairFixtureStrict<ManagerImplementation> {
    using Manager_t = AdaptorStrict<AdaptorNeighbourList<ManagerImplementation>>;

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

    PropertyFixtureStrict(bool consider_ghost_neighbours = false)
        : PairFixtureStrict<ManagerImplementation>{consider_ghost_neighbours},
          scalar_atom_property{*this->adaptor_strict},
          pair_property{*this->adaptor_strict},
          atom_property{*this->adaptor_strict, atom_property_metadata},
          dynamic_property{*this->adaptor_strict, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->adaptor_strict, DynSize(), 1,
                            dynamic_property2_metadata} {}

    AtomScalarProperty_t scalar_atom_property;
    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  template <class ManagerImplementation>
  struct PropertyFixtureStrictWithGhosts : public PairFixtureStrictWithGhosts<ManagerImplementation> {
    using Manager_t = AdaptorStrict<AdaptorNeighbourList<ManagerImplementation>>;

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

    PropertyFixtureStrictWithGhosts()
        : PairFixtureStrictWithGhosts<ManagerImplementation>{},
          pair_property{*this->adaptor_strict},
          atom_property{*this->adaptor_strict, atom_property_metadata},
          dynamic_property{*this->adaptor_strict, DynSize(), 1,
                           dynamic_property_metadata},
          dynamic_property2{*this->adaptor_strict, DynSize(), 1,
                            dynamic_property2_metadata} {}

    PairScalarProperty_t pair_property;
    AtomVectorProperty_t atom_property;
    AtomDynamicProperty_t dynamic_property;
    AtomDynamicProperty2_t dynamic_property2;
  };

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          PropertyFixture<StructureManagerCenters>) {}

  /* ---------------------------------------------------------------------- */
  /*
   * checks if the properties associated with atoms and pairs can be filled
   */
  BOOST_FIXTURE_TEST_CASE(fill_test_simple,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();
    atom_property.resize();
    int pair_property_counter{};
    for (auto atom : pair_manager) {
      atom_property[atom] = atom.get_position();
      for (auto pair : atom) {
        pair_property[pair] = ++pair_property_counter;
      }
    }

    pair_property_counter = 0;
    for (auto atom : pair_manager) {
      auto error = (atom_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(pair_property[pair], ++pair_property_counter);
      }
    }
  }

  //TODO(alex) BEGIN
  /* ---------------------------------------------------------------------- */
  // TODO(alex) tests: AHF<NL<SMC; AS<AHF<ANL
  // TODO(alex) tests:NL<SMC no ghost
  /* ---------------------------------------------------------------------- */
  /*
   * This test assumes that atom index is equal to the cluster index of Order=1
   * Layer=0 using the StructureManagerCenters.
   * 
   */
    template<bool consider_ghost_neighbours>
    struct ANL_SMC_StackFixture_Helper {
      using type = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, consider_ghost_neighbours>;
    };
    template<bool consider_ghost_neighbours>
    struct CommonStacksBoostList {
      using ANL_SMC_StackFixture = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, consider_ghost_neighbours>;
      using type = boost::mpl::list<
          AtomPropertyFixture<StructureManagerCentersStackFixture>,
          AtomPropertyFixture<ANL_SMC_StackFixture>,
          AtomPropertyFixture<ANL_SMC_StackFixture>,
          AtomPropertyFixture<AdaptorHalfListStackFixture<
              ANL_SMC_StackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANL_SMC_StackFixture>>>,
          AtomPropertyFixture<AdaptorMaxOrderStackFixture<
              AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
              ANL_SMC_StackFixture>>>>
          >;
    };

    struct CommonOrderTwoStacksBoostList {
      using ANLWithGhosts_SMC_StackFixture = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, true>;
      using ANLWithoutGhosts_SMC_StackFixture = AdaptorNeighbourListStackFixture<
        StructureManagerCentersStackFixture, false>;
      using type_with_ghosts = boost::mpl::list<
          AtomPropertyFixture<ANLWithGhosts_SMC_StackFixture>,
          AtomPropertyFixture<ANLWithGhosts_SMC_StackFixture>,
          AtomPropertyFixture<AdaptorHalfListStackFixture<
              ANLWithGhosts_SMC_StackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
          AtomPropertyFixture<AdaptorMaxOrderStackFixture<
              AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
              ANLWithGhosts_SMC_StackFixture>>>>
          >;
      using type_without_ghosts = boost::mpl::list<
          AtomPropertyFixture<ANLWithoutGhosts_SMC_StackFixture>,
          AtomPropertyFixture<ANLWithoutGhosts_SMC_StackFixture>,
          AtomPropertyFixture<AdaptorHalfListStackFixture<
              ANLWithoutGhosts_SMC_StackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANLWithoutGhosts_SMC_StackFixture>>>,
          AtomPropertyFixture<AdaptorMaxOrderStackFixture<
              AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
              ANLWithoutGhosts_SMC_StackFixture>>>>
          >;
      // TODO(alex) is there a way to concatenate two lists? boost::mpl::joint_view
      // seems to not work
      // using type = <type_with_ghosts, type_without_ghosts>;
      using type= boost::mpl::list<
          AtomPropertyFixture<ANLWithGhosts_SMC_StackFixture>,
          AtomPropertyFixture<ANLWithGhosts_SMC_StackFixture>,
          AtomPropertyFixture<AdaptorHalfListStackFixture<
              ANLWithGhosts_SMC_StackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANLWithGhosts_SMC_StackFixture>>>,
          AtomPropertyFixture<AdaptorMaxOrderStackFixture<
              AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
              ANLWithGhosts_SMC_StackFixture>>>>,
          AtomPropertyFixture<ANLWithoutGhosts_SMC_StackFixture>,
          AtomPropertyFixture<ANLWithoutGhosts_SMC_StackFixture>,
          AtomPropertyFixture<AdaptorHalfListStackFixture<
              ANLWithoutGhosts_SMC_StackFixture>>,
          AtomPropertyFixture<AdaptorStrictStackFixture<
              AdaptorHalfListStackFixture<ANLWithoutGhosts_SMC_StackFixture>>>,
          AtomPropertyFixture<AdaptorMaxOrderStackFixture<
              AdaptorStrictStackFixture<AdaptorHalfListStackFixture<
              ANLWithoutGhosts_SMC_StackFixture>>>>
          >;
    };

    //TODO(alex)  
    
    using atom_property_fixtures_with_ghosts = CommonStacksBoostList<true>::type;
    using atom_property_fixtures_without_ghosts = CommonStacksBoostList<true>::type;

    //using atom_property_fixtures_order_one = boost::mpl::list<
    //  AtomPropertyFixture<StructureManagerCentersStackFixture>>;

  // If consider_ghost_neighbours is trie the atoms index should always
  // correspond to the cluster index of order 1 using SMC as root implementation
  // and no filtering of order 1
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_property_fixtures_tests, Fix, atom_property_fixtures_with_ghosts, Fix) {
    bool verbose{false};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    size_t cluster_index{0};
    for (auto atom : Fix::manager) {
      BOOST_CHECK_EQUAL(Fix::manager->get_cluster_index(atom.get_atom_index()), cluster_index);
      cluster_index++;
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  //TODO(alex) replace PropertyFixture with PropertyFixtureNew
  
  BOOST_FIXTURE_TEST_CASE_TEMPLATE(atom_propery_access_with_pair_tests, Fix, CommonOrderTwoStacksBoostList::type, Fix) {
    bool verbose{true};
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " starts now." << std::endl;
    }
    Fix::scalar_atom_property.resize();
    // initalize the positions
    std::vector<int> atom_indices{}; // TODO(alex) change to get_manager_atom_indices()
    atom_indices.reserve(Fix::manager->get_size());
    for (auto atom : Fix::manager) {
        Fix::scalar_atom_property[atom] = 0;
        atom_indices.push_back(atom.get_atom_index());
    }
    std::vector<size_t> counter{};
    counter.reserve(atom_indices.size());
    for (size_t i{0}; i<atom_indices.size(); i++) {
      counter.push_back(0);
    }
    // add the position to the atom and count how often this happens
    for (auto atom : Fix::manager) {
      for (auto pair : atom) {
        Fix::scalar_atom_property[pair]++;
        counter.at(Fix::manager->get_cluster_index(pair.get_internal_neighbour_atom_index()))++;
      }
    }
    for (auto atom : Fix::manager) {
      BOOST_CHECK_EQUAL(Fix::scalar_atom_property[atom],
          counter.at(Fix::manager->get_cluster_index(atom.get_atom_index())));
    }
    if (verbose) {
      std::cout << ">> Test for manager ";
      std::cout << Fix::manager->get_name();
      std::cout << " finished." << std::endl;
    }
  }
  // TODO(alex)make this test for new functionality and then delete
  // AHF<ANL<SMC TODO
  //BOOST_FIXTURE_TEST_CASE(atom_index_equal_order_one_cluster_for_smc_strict_test4,
  //                        PropertyFixtureStrictWithGhosts<StructureManagerCenters>) {
  //  atom_property.resize();
  //  // initalize the positions
  //  for (auto atom : adaptor_strict) {
  //      atom_property[atom] = atom.get_position();
  //  }
  //  std::vector<int> atom_indices = adaptor_strict->get_manager_atom_indices();
  //  std::vector<size_t> counter{};
  //  counter.reserve(atom_indices.size());
  //  for (size_t i{0}; i<atom_indices.size(); i++) {
  //    counter.push_back(1);
  //  }
  //  // add the position to the atom and count how often this happens
  //  for (auto atom : adaptor_strict) {
  //    for (auto pair : atom) {
  //      counter.at(pair.get_atom_index())++;
  //      atom_property[pair] += adaptor_strict->get_position(pair.get_atom_index());
  //    }
  //  }
  //  for (auto atom : adaptor_strict) {
  //    auto error = (atom_property[atom] - counter.at(atom.get_atom_index())*atom.get_position()).norm();
  //    BOOST_CHECK_LE(error, tol * 100);
  //  }
  //}

  //TODO(alex) END
  /* ---------------------------------------------------------------------- */

  /* ---------------------------------------------------------------------- */
  /*
   * The access of an order one property with the atom itself
   * and the pair with the atom as neighbour should be the same.
   */ 
  BOOST_FIXTURE_TEST_CASE(fill_test_simple_order_one_property,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();
    atom_property.resize();
    for (auto atom : pair_manager) {
      atom_property[atom] = atom.get_position();
    }

    for (auto atom : pair_manager) {
      for (auto atom2 : pair_manager) {
        for (auto pair : atom) {
          if (atom.back() == pair.back()) {
            auto error = (atom_property[atom] - atom_property[pair]).norm();
            BOOST_CHECK_LE(error, tol * 100);
          }
        }
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE(fill_test_simple_order_one_property_strict,
                          PropertyFixture<StructureManagerCenters>) {
    auto adaptor_strict{
        make_adapted_manager<AdaptorStrict>(pair_manager, 2)};
    adaptor_strict->update();
    pair_property.resize();
    atom_property.resize();
    for (auto atom : adaptor_strict) {
      atom_property[atom] = atom.get_position();
    }

    for (auto atom : adaptor_strict) {
      for (auto atom2 : adaptor_strict) {
        for (auto pair : atom) {
          if (atom.back() == pair.back()) {
            auto error = (atom_property[atom] - atom_property[pair]).norm();
            BOOST_CHECK_LE(error, tol * 100);
          }
        }
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test, if metadata can be assigned to properties
   */
  BOOST_FIXTURE_TEST_CASE(meta_data_test,
                          PropertyFixture<StructureManagerCenters>) {
    auto atom_metadata = atom_property.get_metadata();
    auto dynamic_metadata = dynamic_property.get_metadata();
    auto dynamic_metadata2 = dynamic_property2.get_metadata();

    BOOST_CHECK_EQUAL(atom_metadata, atom_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata, dynamic_property_metadata);
    BOOST_CHECK_EQUAL(dynamic_metadata2, dynamic_property2_metadata);
  }

  /* ---------------------------------------------------------------------- */
  /*
   * test filling statically and dynamically sized properties with actual data
   * and comparing if retrieval of those is consistent with the data that was
   * put in
   */
  BOOST_FIXTURE_TEST_CASE(fill_test_complex,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();
    atom_property.resize();
    dynamic_property.resize();
    dynamic_property2.resize();

    BOOST_CHECK_THROW(
        AtomVectorProperty_t ::check_compatibility(dynamic_property),
        std::runtime_error);

    BOOST_CHECK_NO_THROW(
        AtomVectorProperty_t ::check_compatibility(atom_property));

    int pair_property_counter{};
    size_t counter{};
    for (auto atom : pair_manager) {
      atom_property[atom] = atom.get_position();
      dynamic_property2[atom] = atom.get_position();

      dynamic_property[atom] << counter++, counter, counter;
      for (auto pair : atom) {
        pair_property[pair] = ++pair_property_counter;
      }
    }

    auto & FakeSizedProperty{
        AtomVectorProperty_t::check_compatibility(dynamic_property2)};

    pair_property_counter = 0;
    counter = 0;
    for (auto atom : pair_manager) {
      auto error = (atom_property[atom] - atom.get_position()).norm();
      BOOST_CHECK_LE(error, tol * 100);
      Eigen::Matrix<size_t, DynSize(), Eigen::Dynamic> tmp(DynSize(), 1);
      tmp << counter++, counter, counter;

      auto ierror{(tmp - dynamic_property[atom]).norm()};
      BOOST_CHECK_EQUAL(ierror, 0);

      error = (atom_property[atom] - dynamic_property2[atom]).norm();
      BOOST_CHECK_LE(error, tol * 100);
      error = (atom_property[atom] - FakeSizedProperty[atom]).norm();
      BOOST_CHECK_LE(error, tol * 100);
      for (auto pair : atom) {
        BOOST_CHECK_EQUAL(pair_property[pair], ++pair_property_counter);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  /*
   * test for retrieval of information from property: is it the same that was
   * put in?
   */
  BOOST_FIXTURE_TEST_CASE(compute_distances,
                          PropertyFixture<StructureManagerCenters>) {
    pair_property.resize();

    for (auto atom : pair_manager) {
      for (auto pair : atom) {
        pair_property[pair] =
            (atom.get_position() - pair.get_position()).norm();
      }
    }

    for (auto atom : pair_manager) {
      for (auto pair : atom) {
        auto dist{(atom.get_position() - pair.get_position()).norm()};
        auto error{pair_property[pair] - dist};
        BOOST_CHECK_LE(error, tol / 100);
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

  /* ---------------------------------------------------------------------- */

  /*
   * A fixture for testing partially sparse properties.
   * TODO(felix) use MultipleStructureManagerCentersFixture instead of NL
   */
  struct BlockSparsePropertyFixture
      : public MultipleStructureFixture<MultipleStructureManagerNLFixture> {
    using Parent = MultipleStructureFixture<MultipleStructureManagerNLFixture>;
    using ManagerTypeList_t = typename Parent::ManagerTypeHolder_t::type_list;

    using Key_t = std::vector<int>;
    using BlockSparseProperty_t = BlockSparseProperty<double, 1, 0>;
    using Dense_t = typename BlockSparseProperty_t::Dense_t;
    using InputData_t = typename BlockSparseProperty_t::InputData_t;
    using test_data_t = std::vector<InputData_t>;

    constexpr static Dim_t DynSize() { return 3; }

    std::string sparse_features_desc{"some atom centered sparse features"};

    BlockSparsePropertyFixture() : Parent{} {
      std::random_device rd;
      std::mt19937 gen{rd()};
      auto size_dist{std::uniform_int_distribution<size_t>(1, 10)};
      auto key_dist{std::uniform_int_distribution<int>(1, 100)};
      // size_t i_manager{0};
      for (auto & manager : managers) {
        sparse_features.emplace_back(*manager, sparse_features_desc);
        this->keys_list.emplace_back();
        test_data_t test_data{};
        for (auto atom : manager) {
          // set up random unique keys
          auto set_max_size{size_dist(gen)};
          std::set<Key_t> keys{};
          for (size_t ii{0}; ii < set_max_size; ii++) {
            keys.emplace(key_dist(gen));
          }

          // set up the data to fill the property later
          InputData_t datas{};
          for (auto & key : keys) {
            auto data = Dense_t::Random(21, 8);
            datas.emplace(std::make_pair(key, data));
          }
          this->keys_list.back().push_back(keys);
          test_data.push_back(datas);
        }
        this->test_datas.push_back(test_data);
      }
    }

    std::vector<std::vector<std::set<Key_t>>> keys_list{};
    std::vector<test_data_t> test_datas{};
    std::vector<BlockSparseProperty_t> sparse_features{};
  };

  BOOST_AUTO_TEST_SUITE(Property_partially_sparse_tests);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test, BlockSparsePropertyFixture) {}

  /* ---------------------------------------------------------------------- */
  /*
   * checks if the partially sparse properties associated with centers can be
   * filled and that the data can be accessed consistently.
   */
  BOOST_FIXTURE_TEST_CASE(fill_test_simple, BlockSparsePropertyFixture) {
    bool verbose{false};
    // fill the property structures
    auto i_manager{0};
    for (auto & manager : managers) {
      auto i_center{0};
      sparse_features[i_manager].set_shape(21, 8);
      for (auto center : manager) {
        sparse_features[i_manager].push_back(test_datas[i_manager][i_center]);
        i_center++;
      }
      i_manager++;
    }

    i_manager = 0;
    for (auto & manager : managers) {
      if (verbose)
        std::cout << "manager: " << i_manager << std::endl;
      auto i_center{0};
      for (auto center : manager) {
        if (verbose)
          std::cout << "center: " << i_center << std::endl;

        auto data = sparse_features[i_manager].get_dense_row(center);
        size_t key_id{0};
        double error1{0};
        for (auto & key : sparse_features[i_manager].get_keys(center)) {
          auto & value{test_datas[i_manager][i_center][key]};
          Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> test_data(
              value.data(), value.size());
          error1 += (data.col(key_id) - test_data).squaredNorm();
          key_id += 1;
        }
        if (verbose)
          std::cout << "error1: " << error1 << std::endl;
        BOOST_CHECK_LE(std::sqrt(error1), tol * 100);
        for (auto & key : keys_list[i_manager][i_center]) {
          auto error2 = (sparse_features[i_manager](center, key) -
                         test_datas[i_manager][i_center][key])
                            .norm();
          if (verbose)
            std::cout << "error2: " << error2 << std::endl;
          BOOST_CHECK_LE(error2, tol * 100);
        }
        i_center++;
      }
      i_manager++;
    }
  }

  /* ---------------------------------------------------------------------- */
  /**
   * test, if metadata can be assigned to properties
   */
  BOOST_FIXTURE_TEST_CASE(meta_data_test, BlockSparsePropertyFixture) {
    for (auto & sparse_feature : sparse_features) {
      auto sparse_feature_metadata = sparse_feature.get_metadata();
      BOOST_CHECK_EQUAL(sparse_feature_metadata, sparse_features_desc);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
