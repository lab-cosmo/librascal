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
  BOOST_FIXTURE_TEST_CASE(internal_test,
  RepresentationFixture<StructureManagerCenters>)
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

    std::vector<std::pair<size_t, myiter> > order_mat(dists.size());

    size_t n_{0};
    for (myiter it = dists.begin(); it != dists.end(); ++it, ++n_)
        {order_mat[n_] = make_pair(n_, it);}

    std::sort(order_mat.begin(), order_mat.end(), internal::ordering());
    internal::sort_coulomb_matrix(mat0,mat1,order_mat);
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
  BOOST_FIXTURE_TEST_CASE(constructor_test,
  RepresentationFixture<StructureManagerCenters>)
  {

    AdaptorNeighbourList<StructureManagerCenters> nl{manager,cutoff_max};
    nl.update();
    AdaptorStrict<AdaptorNeighbourList<
                              StructureManagerCenters>> strict_nl{nl,cutoff_max};
    strict_nl.update();

    using Representation_t = RepresentationManagerSortedCoulomb<
                   AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

    Representation_t representation{strict_nl,central_decay,
                                    interaction_cutoff,interaction_decay,size};

  }
    /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(compute_test,
  RepresentationFixture<StructureManagerCenters>)
  {
    bool verbose{false};
    AdaptorNeighbourList<StructureManagerCenters> nl{manager,cutoff_max};
    nl.update();
    AdaptorStrict<AdaptorNeighbourList<
                              StructureManagerCenters>> strict_nl{nl,cutoff_max};
    strict_nl.update();

    using Representation_t = RepresentationManagerSortedCoulomb<
                   AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

    Representation_t representation{strict_nl,central_decay,
                                    interaction_cutoff,interaction_decay,size};
    representation.compute();

    auto rep = representation.get_representation_full();
    if (verbose){
        std::cout << rep.size() <<", "<< rep.cols() <<", "<< rep.rows()<< std::endl;
        for (auto ii{0}; ii < rep.cols(); ++ii){
            for (auto jj{0}; jj < rep.rows(); ++jj){
                std::cout << rep(jj,ii) << ", ";
            }
            std::cout << std::endl;
        }
    }

  }

  BOOST_AUTO_TEST_SUITE_END();
} // RASCAL