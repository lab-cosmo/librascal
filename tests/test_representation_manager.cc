/**
 * file   test_representation_manager.cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  test representation managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_representation_manager.hh"


namespace rascal {
  BOOST_AUTO_TEST_SUITE(representation_sorted_coulomb_test);
  /* ---------------------------------------------------------------------- */
  BOOST_TEST_CASE_TEMPLATE_FUNCTION(test_internals,T)
  {
    bool verbose{false};
    typedef std::vector<double>::const_iterator myiter;
    typedef Eigen::Matrix<double,10,1> vec;

    vec  numbers = vec::Random();
    vec  numbers2 = vec::Random();
    std::vector<double> index{};
    std::vector<double> index2{};

    if (verbose) std::cout << "Test simple sorting " << std::endl;
    for (int ii{0}; ii < numbers.size(); ii++){
        index.push_back(numbers(ii));
        index2.push_back(numbers2(ii));
        if (verbose) std::cout << numbers(ii) << ", " ;
    }
    if (verbose) std::cout << std::endl;

    std::vector<std::pair<size_t, myiter> > order(index.size());

    size_t n{0};
    for (myiter it = index.begin(); it != index.end(); ++it, ++n)
        {order[n] = make_pair(n, it);}

    std::sort(order.begin(), order.end(), internal::ordering());
    auto sorted_index = internal::sort_from_ref(index,order);
    auto sorted_index2 = internal::sort_from_ref(index2,order);
    for (size_t ii{0}; ii < sorted_index.size(); ii++){
        if (verbose) std::cout << sorted_index[ii] << ", " ;
    }
    if (verbose) std::cout << std::endl;
    for (size_t ii{0}; ii < sorted_index2.size(); ii++){
        if (verbose) std::cout << index2[ii] << ", " ;
    }
    if (verbose) std::cout << std::endl;

    for (size_t ii{0}; ii < sorted_index2.size(); ii++){
        if (verbose) std::cout << sorted_index2[ii] << ", " ;
    }
    if (verbose) std::cout << std::endl;

    if (verbose) std::cout << "Test upper diag sorting " << std::endl;
    typedef Eigen::Matrix<double,5,5> matrix;
    typedef Eigen::Matrix<double,5*(5+1)/2,1> lin_mat;

    matrix  mat0 = matrix::Random();
    lin_mat  mat1 = lin_mat::Ones();
    std::vector<double> dists{{2,4,0,1,3}};

    internal::sort_coulomb_matrix(mat0,mat1,dists);
    for (int ii{0}; ii < mat0.rows(); ii++){
      for (int jj{0}; jj < mat0.cols(); jj++){
         if (verbose) std::cout << mat0(ii,jj) << ",\t" ;
      }
      if (verbose) std::cout << std::endl;
    }

    for (int jj{0}; jj < mat1.rows(); jj++){
        if (verbose) std::cout << mat1(jj)<< ", " ;
    }
    if (verbose) std::cout << std::endl;

  }
  /* ---------------------------------------------------------------------- */


  // TODO define more test that could be streamlined
  // gets a list of fixtures for all the different possible structure managers
  using multiple_fixtures = boost::mpl::list<
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSortedCoulomb,
                          MultipleStructureSortedCoulomb,
                          Option::CMSortDistance>,
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSortedCoulomb,
                          MultipleStructureSortedCoulomb,
                          Option::CMSortRowNorm>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_constructor_test,
            Fix, multiple_fixtures, Fix) {
    auto & managers = Fix::managers_strict;
    auto& representations = Fix::representations;
    auto& hypers = Fix::hypers;

    for (auto& manager : managers) {
      for (auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(multiple_compute_test,
            Fix, multiple_fixtures, Fix) {
    auto & managers = Fix::managers_strict;
    auto& representations = Fix::representations;
    auto& hypers = Fix::hypers;
    for (auto& manager : managers) {
      for (auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
        representations.back().compute();
      }
    }
  }


  BOOST_AUTO_TEST_SUITE_END();

  /* ---------------------------------------------------------------------- */

  BOOST_AUTO_TEST_SUITE(representation_spherical_expansion_test);

  using multiple_fixtures = boost::mpl::list<
    RepresentationFixture<StructureManagerCenters,
                          RepresentationManagerSphericalExpansion,
                          MultipleStructureSphericalExpansion>>;

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(
      multiple_constructor_test, Fix, multiple_fixtures, Fix) {

    auto& managers = Fix::managers_strict;
    auto& representations = Fix::representations;
    auto& hypers = Fix::hypers;

    for (auto& manager : managers) {
      for (auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(
      multiple_precompute_test, Fix, multiple_fixtures, Fix) {
    auto& managers = Fix::managers_strict;
    auto& representations = Fix::representations;
    const auto& hypers = Fix::hypers;

    for (auto& manager : managers) {
      for (const auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
        BOOST_TEST(representations.back().get_is_precomputed() == false);
        representations.back().precompute();
        BOOST_TEST(representations.back().get_is_precomputed() == true);
        // And now test automatic precomputation
        representations.emplace_back(manager, hyper);
        BOOST_TEST(representations.back().get_is_precomputed() == false);
        representations.back().compute();
        BOOST_TEST(representations.back().get_is_precomputed() == true);
      }
    }
  }

  BOOST_FIXTURE_TEST_CASE_TEMPLATE(
      multiple_compute_test, Fix, multiple_fixtures, Fix) {

    auto& managers = Fix::managers_strict;
    auto& cutoffs = Fix::cutoffs;
    auto& representations = Fix::representations;
    const auto& filenames = Fix::filenames;
    const auto& hypers = Fix::hypers;
    bool verbose{true};

    auto filename_it{filenames.begin()};
    auto cutoff_it{cutoffs.begin()};

    for (auto& manager : managers) {
      if (verbose) {
        if (filename_it != filenames.end()) {
          std::cout << "Structure: " << *filename_it;
          std::cout << " with cutoff " << *cutoff_it << std::endl;
          ++cutoff_it;
          if (cutoff_it == cutoffs.end()) {
            cutoff_it = cutoffs.begin();
            ++filename_it;
          }
        } else {
          std::cout << "Structure: Filename unknown" << std::endl;
        }
      }
      for (const auto& hyper : hypers) {
        representations.emplace_back(manager, hyper);
        // Should be done automatically in compute()
        // TODO(max-veit) make that its own test case
        //representations.back().precompute();
        representations.back().compute();
        // Check dimensions of the storage array
        size_t max_radial = hyper.at("max_radial");
        size_t max_angular =  hyper.at("max_angular");
        BOOST_TEST(representations.back().get_feature_size() ==
          max_radial * (max_angular + 1) * (max_angular + 1));
        if (verbose) {
          size_t center_idx{0};
          for (auto center : manager) {
            std::cout << "Soap vector for center: " << center_idx++;
            std::cout << std::endl;
            representations.back().print_soap_vector(center, std::cout);
            //std::cout << std::endl;
          }
        }
      }
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
} // RASCAL
