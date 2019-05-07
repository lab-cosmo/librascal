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

#include <string>
#include <memory>

#include "structure_managers/make_structure_manager.hh"
#include "structure_managers/structure_manager_centers.hh"
#include "structure_managers/adaptor_neighbour_list.hh"
#include "structure_managers/adaptor_strict.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap.hh"

#include "json_io.hh"


//using namespace rascal;
const int N_ITERATIONS = 10000; // it's over 9000


int main(int argc, char* argv[]) {
  using rascal::StructureManagerTypeHolder;
  using rascal::StructureManagerCenters;
  using rascal::AdaptorNeighbourList;
  using rascal::AdaptorStrict;
  using rascal::RepresentationManagerSOAP;
  using rascal::make_structure_manager_stack_with_hypers_and_typeholder;
  //using json::json;
  using ManagerTypeHolder_t =
     StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
                                AdaptorStrict>;
  using ManagerTypeList_t = typename ManagerTypeHolder_t::type_list;
  using Manager_t = typename ManagerTypeHolder_t::type;
  using ManagerPtr_t = std::shared_ptr<Manager_t>;

  //double cutoff{3.};
  if (argc < 1) {
    std::cerr << "Must provide atomic structure json filename as argument";
    std::cerr << std::endl;
  }
  std::string filename{argv[0]};
  //std::string filename{"../tests/reference_data/methane.json"};
  //std::string filename{"../tests/reference_data/small_molecule.json"};
  //TODO(max) put these in a file so they can be varied systematically
  //maybe together with the filename and iteration count
  json hypers{{"interaction_cutoff", 3.0},
              {"cutoff_smooth_width", 0.5},
              {"max_radial", 6},
              {"max_angular", 6},
              {"gaussian_sigma_type", "Constant"},
              {"gaussian_sigma_constant", 0.3},
              {"soap_type", "PowerSpectrum"}};
  double cutoff{hypers["interaction_cutoff"]};
  json structure{{"filename", filename}};
  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments",
            {{"cutoff", cutoff},
             {"consider_ghost_neighbours", false}}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad2);

  auto manager = make_structure_manager_stack_with_hypers_and_typeholder<
    ManagerTypeList_t>::apply(structure, adaptors);
  RepresentationManagerSOAP<ManagerPtr_t> representation{manager, hypers};

  // This is the part that should get profiled
  for (size_t looper{0}; looper < N_ITERATIONS; looper++) {
    representation.compute();
  }
}
