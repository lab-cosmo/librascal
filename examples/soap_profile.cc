/**
 * file   soap_profile.cc
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   7 May 2019
 *
 * @brief  Example for profiling the spherical expansion and SOAP
 *
 * Copyright © 2018 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
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
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap.hh"
#include "representations/feature_manager_dense.hh"
#include "basic_types.hh"


#include <iostream>
#include <basic_types.hh>
#include <cmath>
#include <list>
#include <functional>
#include <string>
#include <initializer_list>

// using namespace std;
using namespace rascal;  // NOLINT

const int N_ITERATIONS = 10000; // it's over 9000

using Representation_t = RepresentationManagerSOAP<
    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>>;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide atomic structure json filename as argument";
    std::cerr << std::endl;
  }

  //TODO(max) put these in a file so they can be varied systematically
  //maybe together with the filename and iteration count
  std::string filename{argv[1]};
  json hypers{{"interaction_cutoff", 2.0},
              {"cutoff_smooth_width", 0.0},
              {"max_radial", 6},
              {"max_angular", 6},
              {"gaussian_sigma_type", "Constant"},
              {"gaussian_sigma_constant", 0.2},
              {"soap_type", "PowerSpectrum"}};
  json structure{{"filename", filename}};

  double cutoff{hypers["interaction_cutoff"]};
  json adaptors;
        json ad1{{"name", "AdaptorNeighbourList"},
                {"initialization_arguments",
                  {{"cutoff", cutoff},
                  {"consider_ghost_neighbours", false}}}};
        json ad2{{"name", "AdaptorStrict"},
                {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);
  auto manager =
    make_structure_manager_stack<StructureManagerCenters, AdaptorNeighbourList,
                                 AdaptorStrict>(structure, adaptors);

  Representation_t representation{manager, hypers};

  // This is the part that should get profiled
  for (size_t looper{0}; looper < N_ITERATIONS; looper++) {
    representation.compute();
  }
}
