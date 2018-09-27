/**
 * file   test_adaptor_strict.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   04 Jun 2018
 *
 * @brief  tests the implementation of the strict structure adaptor
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
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
#include "test_structure.hh"
#include <vector>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(strict_adaptor_test);

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(constructor_test,
                          ManagerFixture<StructureManagerCenters>) {
    double cut_off{3.5};
    std::vector<std::vector<int>> neigh_ids;
    std::vector<std::vector<double>> neigh_dist;

    std::vector<std::vector<int>> neigh_ids_strict;
    std::vector<std::vector<double>> neigh_dist_strict;
    
    AdaptorMaxOrder<StructureManagerCenters> pair_manager{manager, cutoff};
    pair_manager.update();
    std::cout << "Setting up strict manager" << std::endl;
    AdaptorStrict<AdaptorMaxOrder<StructureManagerCenters>> 
                                      adaptor_strict{pair_manager, cut_off};
    adaptor_strict.update();

    // std::cout << "Testing adaptor_strict" << std::endl;
    for (auto center : pair_manager) {
      std::vector<int> indices;
      std::vector<double> distances;
      // std::cout << "strict atom out " << atom.get_index() << std::endl;
      for (auto neigh : center) {
        double distance{(center.get_position()
                         - neigh.get_position()).norm()};
        if (distance <= cut_off) {              
          indices.push_back(neigh.get_atom_index());
          distances.push_back(distance);
        }
        // std::cout << "  strict pair out " << pair.get_index() << std::endl;
        // auto dist = adaptor.get_distance(pair);
        // std::cout << "distance " << dist << std::endl;
      }
      neigh_ids.push_back(indices);
      neigh_dist.push_back(distances);
    }

    for (auto center : adaptor_strict) {
      std::vector<int> indices;
      std::vector<double> distances;
      // std::cout << "strict atom out " << atom.get_index() << std::endl;
      for (auto neigh : center) {
        double distance{(center.get_position()
                         - neigh.get_position()).norm()};
        if (distance <= cut_off) {              
          indices.push_back(neigh.get_atom_index());
          distances.push_back(distance);
        }
      }
      neigh_ids_strict.push_back(indices);
      neigh_dist_strict.push_back(distances);
    }

    
    BOOST_CHECK_EQUAL(neigh_ids.size(),neigh_ids_strict.size());
    for (int ii{0};ii<neigh_ids.size();++ii){
      BOOST_CHECK_EQUAL(neigh_ids[ii].size(),neigh_ids_strict[ii].size());
      for (int jj{0};jj<neigh_ids[ii].size();++jj){
        int a0{neigh_ids[ii][jj]};
        int a1{neigh_ids_strict[ii][jj]};
        double d0{neigh_dist[ii][jj]};
        double d1{neigh_dist_strict[ii][jj]};
        BOOST_CHECK_EQUAL(a0,a1);
        BOOST_CHECK_EQUAL(d0,d1);
      }
    }


  }


  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
