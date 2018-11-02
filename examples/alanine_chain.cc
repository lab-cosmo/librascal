/**
 * file   alanine_chain.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   18 Jun 2018
 *
 * @brief Example for NeighbourhoodManagerChain
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include <structure_managers/structure_manager_chain.hh>
#include <structure_managers/property.hh>

#include "json.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using Manager_t = rascal::StructureManagerChain;

int main() {

  // integer lists containing the definition of the quadruplets for the
  // calculation of the dihedral angles in the alanine-unit.
  std::vector<int> phi{4, 3, 19, 1};
  std::vector<int> psi{19, 1, 0, 18};
  // put them in a container for iteration
  std::vector<std::vector<int>> quadruplets{};
  quadruplets.push_back(phi);
  quadruplets.push_back(psi);

  // solution vector for the dihedral angles
  std::vector<double> dihedral_angles{};

  // initialize the manager
  Manager_t manager;
  double cutoff{1.0};

  // read atomic structure from the JSON file
  manager.read_structure_from_json("alanine-X.json");
  manager.update(cutoff);

  // Get the positions to work with. The return type is an
  // <code>Eigen::Map</code> to the underlying array. This means, that
  // all the <code>Eigen</code> magic works on the structure, e.g. dot
  // and cross products or norms.
  auto positions = manager.get_positions();

  // Loop over the defined quadruplets and calculate the respective
  // angles with atan2 and cosine definition.
  for (auto q : quadruplets) {
    auto b1 = positions.col(q[1]) - positions.col(q[0]);
    auto b2 = positions.col(q[1]) - positions.col(q[2]);
    auto b3 = positions.col(q[3]) - positions.col(q[2]);

    auto na = b1.cross(b2);
    auto nb = b2.cross(b3);

    auto arg1 = ((b1.cross(b2)).cross(b2.cross(b3))).dot(b2 / b2.norm());
    auto arg2 = (b1.cross(b2)).dot(b2.cross(b3));

    auto angle = std::atan2(arg1, arg2);

    dihedral_angles.push_back(angle);

    std::cout << "atan2 " << angle << std::endl;
    std::cout << "cos " << std::acos(na.dot(nb) / na.norm() / nb.norm())
              << std::endl;
  }

  std::cout << "Dihedral angles in alanine unit \n";
  for (auto a : dihedral_angles) {
    std::cout << a << " ";
  }
  std::cout << std::endl;
} // end main
