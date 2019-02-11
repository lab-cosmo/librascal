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

#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/make_structure_manager.hh"
#include "rascal_utility.hh"
#include "representations/representation_manager_sorted_coulomb.hh"
#include "basic_types.hh"

#include <Eigen/StdVector>

#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>

//using namespace std;
using namespace rascal; // NOLINT

//using Manager_t = StructureManagerCenters;
constexpr static int dim{3};
using Vector_t = Eigen::Matrix<double, dim, 1>;

using Representation_t = RepresentationManagerSortedCoulomb<
              AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;





int main() {

  Eigen::MatrixXd positions(22, 3);
  Eigen::VectorXi atom_types(22);
  Eigen::MatrixXd cell(3, 3);
  std::array<int, 3> pbc{{true, true, true}};
  cell << 6.19, 2.41, 0.21, 0.00, 6.15, 1.02, 0.00, 0.00, 7.31;

  // clang-format off
  positions << 3.689540159937393, 5.123016813620886, 1.994119731169116,
      6.818437242389163, 2.630056617829216, 6.182500355729062,
      2.114977334498767, 6.697579639059512, 1.392155450018263,
      7.420401523540017, 2.432242071439904, 6.380314902118375,
      1.112656394115962, 7.699900579442317, 3.569715877854675,
      5.242841095703604, 3.122826344932127, 5.689730628626151,
      3.248684682453303, 5.563872291104976, 2.608353462112637,
      6.204203511445642, 5.035681855581504, 2.134827911489532,
      0.946910011088814, 6.223599755982222, 4.168634519120968,
      3.001875247950068, 1.980327734683430, 5.190182032387606,
      2.943861424421339, 4.226648342649697, 5.457161501166098,
      1.713348265904937, 1.501663178733906, 5.668846588337130,
      5.208365510425203, 1.962144256645833, 2.728127406527150,
      4.442382360543885, 2.839975217222644, 4.330534549848392,
      0.744216089807768, 6.426293677263268, 4.643695520786083,
      2.662204050783991, 1.250682335857938, 6.055217235712136,
      0.860905287815103, 6.444994283754972, 4.536108843695142,
      2.769790727874932, 5.609177455068640, 1.696722116501434,
      6.703053268421970, 0.602846303148105, 3.487609972580834,
      3.818289598989240, 1.436734374347541, 5.869165197222533,
      1.054504320562138, 6.251395251007936, 3.998423858825871,
      3.307475712744203, 5.323662899811682, 1.982236671758393;
  // clang-format on
  positions.transposeInPlace();
  atom_types << 20, 20, 24, 24, 15, 15, 15, 15, 8, 8, 8, 8, 8, 8, 8, 8, 8,
      8, 8, 8, 8, 8;

  using PBC_t = Eigen::Map<Eigen::Matrix<int, 3, 1>>;
  // manager->update(positions, atom_types, cell, PBC_t{pbc.data()});

  double cutoff{2.};
  bool verbose{false};
  int mult = 6;
  double rc_max{mult * 0.5 + cutoff};
  auto manager{make_structure_manager<StructureManagerCenters>()};
  auto pair_manager{make_adapted_manager<AdaptorNeighbourList>(manager, rc_max)};

  for (auto i{0}; i < mult; ++i) {
    auto cutoff_tmp = i * 0.5 + cutoff;
    std::vector<std::vector<int>> neigh_ids{};
    std::vector<std::vector<double>> neigh_dist{};
    std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec{};
    std::vector<std::vector<int>> neigh_ids_strict{};
    std::vector<std::vector<double>> neigh_dist_strict{};
    std::vector<std::vector<std::array<double, 3>>> neigh_dir_vec_strict{};

    if (verbose) {
      std::cout << "Setting up strict manager with rc = " << cutoff_tmp
                << std::endl;
    }

    auto adaptor_strict{make_adapted_manager<AdaptorStrict>(pair_manager, cutoff_tmp)};
    adaptor_strict->update(positions, atom_types, cell, PBC_t{pbc.data()});

    if (verbose) {
      std::cout << "Setting up comparison list with rc = " << cutoff_tmp
                << std::endl;
    }
    for (auto center : pair_manager) {
      std::vector<int> indices{};
      std::vector<double> distances{};
      std::vector<std::array<double, 3>> dir_vecs{};

      if (verbose) {
        // get_index returns iteration index
        std::cout << "cell atom out " << center.get_index();
        // get_atom_index returns index from
        std::cout << " " << center.get_atom_index() << " ";

        for (int ii{0}; ii < 3; ++ii) {
          std::cout << center.get_position()[ii] << " ";
        }
        std::cout << " " << center.get_atom_type() << std::endl;
      }

      for (auto neigh : center) {
        double distance{
            (center.get_position() - neigh.get_position()).norm()};
        if (distance <= cutoff_tmp) {
          indices.push_back(neigh.get_atom_index());
          distances.push_back(distance);
          auto dir_vec{
              (neigh.get_position() - center.get_position()).array() /
              distance};
          std::array<double, 3> aa{{dir_vec(0), dir_vec(1), dir_vec(2)}};
          dir_vecs.push_back(aa);
          if (verbose) {
            std::cout << "cell neigh out " << neigh.get_index();
            std::cout << " " << neigh.get_atom_index() << " ";

            for (int ii{0}; ii < 3; ++ii) {
              std::cout << neigh.get_position()[ii] << " ";
            }
            std::cout << " " << neigh.get_atom_type() << std::endl;
          }
        }
      }
      neigh_ids.push_back(indices);
      neigh_dist.push_back(distances);
      neigh_dir_vec.push_back(dir_vecs);
      // break;
    }

    if (verbose) {
      std::cout << "Setting get adaptor_strict info" << std::endl;
    }
    for (auto center : adaptor_strict) {
      std::vector<int> indices_{};
      std::vector<double> distances_{};
      std::vector<std::array<double, 3>> dir_vecs_{};

      if (verbose) {
        // get_index returns iteration index
        std::cout << "strict atom out " << center.get_index();
        // get_atom_index returns index from
        std::cout << " " << center.get_atom_index() << " ";

        for (int ii{0}; ii < 3; ++ii) {
          std::cout << center.get_position()[ii] << " ";
        }
        std::cout << " " << center.get_atom_type() << std::endl;
      }

      for (auto neigh : center) {
        double distance{adaptor_strict->get_distance(neigh)};

        indices_.push_back(neigh.get_atom_index());
        distances_.push_back(distance);
        auto dir_vec{adaptor_strict->get_direction_vector(neigh)};
        std::array<double, 3> bb{{dir_vec(0), dir_vec(1), dir_vec(2)}};
        dir_vecs_.push_back(bb);

        if (verbose) {
          std::cout << "strict neigh out " << neigh.get_index();
          std::cout << " " << neigh.get_atom_index() << "\t ";

          for (int ii{0}; ii < 3; ++ii) {
            std::cout << neigh.get_position()[ii] << ", ";
          }
          std::cout << "\t dist=" << distance;
          std::cout << "\t " << neigh.get_atom_type() << std::endl;
        }
      }

      if (verbose) {
        std::cout << "Number of Neighbourg: " << indices_.size() << std::endl;
      }

      neigh_ids_strict.push_back(indices_);
      neigh_dist_strict.push_back(distances_);
      neigh_dir_vec_strict.push_back(dir_vecs_);
    }


    for (size_t ii{0}; ii < neigh_ids.size(); ++ii) {

      for (size_t jj{0}; jj < neigh_ids[ii].size(); ++jj) {
        int a0{neigh_ids[ii][jj]};
        int a1{neigh_ids_strict[ii][jj]};
        double d0{neigh_dist[ii][jj]};
        double d1{neigh_dist_strict[ii][jj]};
        for (size_t kk{0}; kk < neigh_dir_vec[ii][jj].size(); ++kk) {
          double dv0{neigh_dir_vec[ii][jj][kk]};
          double dv1{neigh_dir_vec_strict[ii][jj][kk]};
        }
      }
    }
  }

  return(0);
}
