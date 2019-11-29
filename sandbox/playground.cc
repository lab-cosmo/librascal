/**
 * @file   examples/spherical_invariants_example.cc
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   26 June 2019
 *
 * @brief  Example for computing the spherical invariants (SOAP)
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

#include "rascal/atomic_structure.hh"
#include "rascal/basic_types.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_half_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/utils.hh"

#include <chrono>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <string>

using namespace rascal;  // NOLINT

using Representation_t = CalculatorSphericalExpansion;
using Manager_t = AdaptorStrict<
    AdaptorCenterContribution<
                AdaptorNeighbourList<StructureManagerCenters>>>;
using Prop_t = typename CalculatorSphericalExpansion::Property_t<Manager_t>;
using ManagerHalf_t = AdaptorStrict<
    AdaptorCenterContribution<AdaptorHalfList<
                AdaptorNeighbourList<StructureManagerCenters>>>>;
using PropHalf_t = typename CalculatorSphericalExpansion::Property_t<ManagerHalf_t>;
using PropGradHalf_t =
    typename CalculatorSphericalExpansion::PropertyGradient_t<ManagerHalf_t>;

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide atomic structure json filename as argument";
    std::cerr << std::endl;
    return -1;
  }

  std::string filename{argv[1]};

  double cutoff{4.};
  json hypers{{"max_radial", 1},
              {"max_angular", 1},
              {"compute_gradients", true},
              {"soap_type", "PowerSpectrum"},
              {"normalize", true}};

  json fc_hypers{{"type", "ShiftedCosine"},
                 {"cutoff", {{"value", cutoff}, {"unit", "AA"}}},
                 {"smooth_width", {{"value", 0.5}, {"unit", "AA"}}}};
  json sigma_hypers{{"type", "Constant"},
                    {"gaussian_sigma", {{"value", 0.4}, {"unit", "AA"}}}};

  hypers["cutoff_function"] = fc_hypers;
  hypers["gaussian_density"] = sigma_hypers;
  hypers["radial_contribution"] = {{"type", "GTO"}};

  json structure{{"filename", filename}};
  json adaptors;
  json adaptors_half;
  json ad1a{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  json ad1b{{"name", "AdaptorHalfList"},
            {"initialization_arguments", {}}};
  json ad1c{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1a);
  adaptors.emplace_back(ad1c);
  adaptors.emplace_back(ad2);

  adaptors_half.emplace_back(ad1a);
  adaptors_half.emplace_back(ad1b);
  adaptors_half.emplace_back(ad1c);
  adaptors_half.emplace_back(ad2);
  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structure, adaptors);
  auto manager_half =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList, AdaptorHalfList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structure, adaptors_half);

  Representation_t representation{hypers};

  representation.compute(manager);

  representation.compute(manager_half);

  constexpr size_t n_centers_print{5};
  constexpr size_t n_neigh_print{1000};
  size_t center_count{0};

  // Print the first few elements and gradients, so we know we're getting
  // something
  std::cout << "Expansion of first " << n_centers_print << " centers:";
  std::cout << std::endl;
  std::cout << "Note that the coefficients are printed with species pairs along"
               " the columns and n-n'-l along the rows."
            << std::endl;
  std::cout << "Gradients are printed with: First Cartesian component, "
               "then species pairs, along the columns; n-n'-l along the rows.";
  std::cout << std::endl;


  auto && soap_vectors{
      *manager->template get_property_ptr<Prop_t>(representation.get_name())};
  auto && soap_vectors_half{
      *manager_half->template get_property_ptr<PropHalf_t>(representation.get_name())};
  // auto && soap_vector_gradients{*manager->template get_property_ptr<PropGrad_t>(
  //     representation.get_gradient_name())};

  for (auto center : manager) {
    // if (center_count >= n_centers_print) {
    //   break;
    // }
    size_t n_species_center{soap_vectors.get_keys(center).size()};
    std::cout << "============================" << std::endl;
    std::cout << "Center " << center.get_index();
    std::cout << " of type " << center.get_atom_type() << std::endl;
    if (std::abs((soap_vectors.get_dense_row(center)-soap_vectors_half.get_dense_row(center)).mean()) > 1e-12) {
      std::cout << "Ref: " << std::endl<< soap_vectors.get_dense_row(center);
      std::cout << std::endl;
      std::cout << "Test: " << std::endl<< soap_vectors_half.get_dense_row(center);
      std::cout << std::endl;
    }
    auto keys_center = soap_vectors[center].get_keys();
    std::cout << "Center data keys: ";
    for (auto key : keys_center) {
      std::cout << "(";
      for (auto key_sp : key) {
        std::cout << key_sp << ", ";
      }
      std::cout << "\b\b) ";
    }
    std::cout << std::endl;
    auto ii_pair = center.get_atom_ii();

    auto half_it = manager_half->get_iterator_at(center_count, 0);
    auto half_center = *(half_it);
    std::cout << "Tags:  (";
    for (auto neigh : half_center) {
      auto neigh_type = neigh.get_atom_tag();
      std::cout << neigh_type << ", ";
    }
    std::cout << "\b\b) ";
    std::cout << std::endl;

    std::cout << "Types: (";
    for (auto neigh : half_center) {
      auto neigh_type = neigh.get_atom_type();
      std::cout << neigh_type << ", ";
    }
    std::cout << "\b\b) ";
    std::cout << std::endl;
    // auto keys_grad_center = soap_vector_gradients[ii_pair].get_keys();
    // std::cout << "Center gradient keys: ";
    // for (auto key : keys_grad_center) {
    //   std::cout << "(";
    //   for (auto key_sp : key) {
    //     std::cout << key_sp << ", ";
    //   }
    //   std::cout << "\b\b) ";
    // }
    // std::cout << std::endl;

    // size_t neigh_count{0};
    // for (auto neigh : center) {
    //   if (neigh_count >= n_neigh_print) {
    //     break;
    //   }
      // auto neigh_type = neigh.get_atom_type();
      // auto tags = neigh.get_atom_tag_list();
      // std::cout << "Neighbour "<< neigh_type<<" tags: ";

      // std::cout << "(";
      // for (auto tag : tags) {
      //   std::cout << tag << ", ";
      // }
      // std::cout << "\b\b) ";

      // std::cout << std::endl;
      // auto keys_neigh = soap_vector_gradients[neigh].get_keys();
      // std::cout << "Neighbour "<< neigh_type<<" keys: ";
      // for (auto key : keys_neigh) {
      //   std::cout << "(";
      //   for (auto key_sp : key) {
      //     std::cout << key_sp << ", ";
      //   }
      //   std::cout << "\b\b) ";
      // }
      // std::cout << std::endl;
    //   ++neigh_count;
    // }
    ++center_count;
  }
}
