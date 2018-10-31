/**
 * file   test_nl.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   2018
 *
 * @brief Example for Neighbour list
 *
 * Copyright © 2018 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "basic_types.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>

using namespace std;

using Manager_t = rascal::StructureManagerCenters;
constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

using Representation_t = rascal::RepresentationManagerSortedCoulomb<
                   rascal::AdaptorStrict<rascal::AdaptorNeighbourList<rascal::StructureManagerCenters>>>;

int main() {
  Manager_t manager{};
  Eigen::MatrixXd positions(3,8);
  Eigen::Matrix<int, Eigen::Dynamic, 1> numbers(8);
  Eigen::MatrixXd cell(3, 3);
  std::array<int, 3> pbc{{true,true,true}};
  bool verbose{false};
  bool verbose_rep{false};

  cell << 2., 0., 0.,
        0., 2., 0.,
        0., 0., 2.;

  positions << 0.4, 1.4, 0.4, 1.4, 0.4, 1.4, 0.4, 1.4,
    0.4, 0.4, 1.4, 1.4, 0.4, 0.4, 1.4, 1.4,
    0.4, 0.4, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4;

  numbers << 1,2,3,4,5,6,7,8;

  // setting up the manager
  manager.update(positions, numbers, cell,
                 Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});

  double cut_off{2};

  int mult{4};
  double rc_max{mult*0.5 + cut_off};
  // build a neighbourlist
  rascal::AdaptorNeighbourList<rascal::StructureManagerCenters>
    pair_manager{manager, rc_max };
  // update to execute the build
  pair_manager.update();

  for (auto i{0}; i < mult; ++i) {
    auto cutoff_tmp = i*0.5 + cut_off;
    std::vector<std::vector<int>> neigh_ids;
    std::vector<std::vector<double>> neigh_dist;

    std::vector<std::vector<int>> neigh_ids_strict;
    std::vector<std::vector<double>> neigh_dist_strict;

    std::cout << "Setting up strict manager with rc="<<cutoff_tmp << std::endl;
    // make strict neighbour list
    rascal::AdaptorStrict<
      rascal::AdaptorNeighbourList<rascal::StructureManagerCenters>>
      adaptor_strict{pair_manager, cutoff_tmp};
    // execute
    adaptor_strict.update();


    if (verbose) std::cout << "Setting get adaptor_strict info" << std::endl;
    for (auto center : adaptor_strict) {
      auto icenter{center.get_index()};
      std::vector<int> indices_{};
      std::vector<double> distances_{};

      if (verbose) {
        // get_index returns iteration index
        std::cout << "strict atom out " << center.get_index();
        // get_atom_index returns index from
        std::cout << " " << center.get_atom_index() << " " ;

        for (int ii{0};ii<3;++ii){
          std::cout << center.get_position()[ii] << " ";
        }
        std::cout << " " << center.get_atom_type() << std::endl;
      }

      int Nneigh{0};
      for (auto neigh : center) {
        double distance{(center.get_position()
                         - neigh.get_position()).norm()};
        Nneigh += 1;
        indices_.push_back(neigh.get_atom_index());
        distances_.push_back(adaptor_strict.get_distance(neigh));

        if (verbose) {
          std::cout << "strict neigh out " << neigh.get_index();
          std::cout << " " << neigh.get_atom_index() << "\t " ;

          for (int ii{0};ii<3;++ii){
            std::cout << neigh.get_position()[ii] << ", ";
          }
          std::cout << "\t dist=" << distance;
          std::cout << "\t " << neigh.get_atom_type() << std::endl;
        }

      }

      std::cout << "Number of Neighbourg: " << Nneigh << std::endl;
      neigh_ids_strict.push_back(indices_);
      neigh_dist_strict.push_back(distances_);
      if (icenter > 1) break;
    }
  }

  rascal::AdaptorStrict<
    rascal::AdaptorNeighbourList<rascal::StructureManagerCenters>>
    adaptor_strict{pair_manager, rc_max};
  // execute
  adaptor_strict.update();
  double central_decay{10};
  double interaction_cutoff{10};
  double interaction_decay{10};
  size_t size{50};
  Representation_t representation{adaptor_strict,central_decay,
                                  interaction_cutoff,interaction_decay,size};
  representation.compute();

  auto rep = representation.get_representation_full();
  if (verbose_rep){
      std::cout << rep.size() <<", "<< rep.cols() <<", "<< rep.rows()<< std::endl;
      for (auto ii{0}; ii < rep.cols(); ++ii){
          for (auto jj{0}; jj < rep.rows(); ++jj){
              std::cout << rep(jj,ii) << ", ";
          }
          std::cout << std::endl;
      }
  }



  return(0);
}
