/**
 * @file   sandbox/playground.cc
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   26 June 2019
 *
 * @brief an executable to test ideas
 *
 * Copyright © 2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "rascal/basic_types.hh"
#include "rascal/models/kernels.hh"
#include "rascal/utils.hh"
#include "rascal/representations/calculator_sorted_coulomb.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/calculator_spherical_invariants.hh"
#include "rascal/structure_managers/adaptor_increase_maxorder.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_half_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <algorithm>
#include <iterator>

using namespace rascal;  // NOLINT

using Representation_t = CalculatorSphericalExpansion;
using Manager_t = AdaptorStrict<
    AdaptorCenterContribution<
                AdaptorNeighbourList<StructureManagerCenters>>>;
using Prop_t = typename CalculatorSphericalExpansion::Property_t<Manager_t>;
using PropGrad_t =
    typename CalculatorSphericalExpansion::PropertyGradient_t<Manager_t>;
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

  double cutoff{3.};

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

   auto manager =
       make_structure_manager_stack<StructureManagerCenters,
                                    AdaptorNeighbourList,
                                    AdaptorCenterContribution, AdaptorStrict>(
           structure, adaptors);
//  auto manager =
//      make_structure_manager_stack<StructureManagerCenters,
//                                   AdaptorNeighbourList, AdaptorStrict>(
//          structure, adaptors);



   std::cout << "n_centers: " << manager->size() << std::endl;
   for (auto center : manager) {
     auto ctag = center.get_atom_tag();
     std::cout << "Center: " << ctag << " n. neighbors " << center.pairs().size()
               << std::endl;

     for (auto neigh : center.pairs()) {
       auto tag_list = neigh.get_atom_tag_list();

       auto atom_j = neigh.get_atom_j();
       auto atom_j_tag = atom_j.get_atom_tag_list();
       auto atom_j_ids = atom_j.get_cluster_indices();
       std::cout << "neigh: " << tag_list[0] << ", " << tag_list[1] << ", "
                 << " tag_j: " << atom_j_tag[0] << ", " << atom_j_ids[0]
                 << " -- global index " << neigh.get_global_index()
                 <<std::endl;
     }
     // for (auto triplet : center.triplets()) {
     //   std::cout << "triplet: " << std::endl;
     // }
   }


  auto triplet_manager{make_adapted_manager<AdaptorMaxOrder>(manager)};
  triplet_manager->update();

//  for (auto center : triplet_manager) {
//    auto proxy = center.pairs();
//    auto it = proxy.begin();
//    auto neigh = *it;
//    auto pos = neigh.get_position();
//  }

  for (auto center : triplet_manager) {
    auto ctag = center.get_atom_tag();
    auto it = center.triplets();
    auto size{it.size()};
    std::cout << "Center: " << ctag << " n. neighbors " << size
              << std::endl;
//     std::cout << "Center: " << ctag << " n. neighbors "
//               << std::endl;

    for (auto triplet : center.triplets()) {
      auto tags = triplet.get_atom_tag_list();
      std::cout << center.get_atom_tag() << " triplet ("
                << tags[0] << ", " << tags[1] << ", " << tags[2]
                << ") global index " << triplet.get_global_index()
              << ") index " << triplet.get_index()
                << std::endl;
    }
  }

  representation.compute(manager_half);

  constexpr size_t n_centers_print{5};
  constexpr size_t n_neigh_print{1000};
  size_t center_count{0};

  std::vector<int> new_tag_list_s{{1,6,7,8}};
  std::sort(new_tag_list_s.begin(), new_tag_list_s.end());
  auto any_equal = std::adjacent_find(new_tag_list_s.begin(), new_tag_list_s.end());
  std::cout << std::boolalpha
            << (any_equal == new_tag_list_s.end())
            << std::endl;
  std::cout << "Gradients are printed with: First Cartesian component, "
               "then species pairs, along the columns; n-n'-l along the rows.";
  std::cout << std::endl;


  auto && soap_vectors{
      *manager->template get_property_ptr<Prop_t>(representation.get_name())};
  auto && soap_vectors_half{
      *manager_half->template get_property_ptr<PropHalf_t>(representation.get_name())};
  auto && soap_vector_gradients{*manager->template get_property_ptr<PropGrad_t>(
      representation.get_gradient_name())};
  auto && soap_vector_gradients_half{*manager_half->template get_property_ptr<PropGradHalf_t>(
      representation.get_gradient_name())};

  for (auto center : manager) {
    // if (center_count >= n_centers_print) {
    //   break;
    // }
    size_t n_species_center{soap_vectors.get_keys(center).size()};
    std::cout << "============================" << std::endl;
    std::cout << "Center " << center.get_index();
    std::cout << " of type " << center.get_atom_type() << std::endl;
    int maxRow, maxCol;
    auto diff_rep_m{math::relative_error(soap_vectors.get_dense_row(center), soap_vectors_half.get_dense_row(center), 1e-15)};
    double diff_rep = diff_rep_m.maxCoeff(&maxRow, &maxCol);
    std::cout << "max error: " << diff_rep << " ref val: " << soap_vectors.get_dense_row(center)(maxRow, maxCol) << " Nb_neigh: "<< center.size() << std::endl;
    if (diff_rep > 1e-13) {
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

    auto diff_ii_m{math::relative_error(soap_vector_gradients.get_dense_row(ii_pair), soap_vector_gradients_half.get_dense_row(half_center.get_atom_ii()), 1e-15)};
    double diff_ii = diff_ii_m.maxCoeff(&maxRow, &maxCol);
    std::cout << "max error: " << diff_ii << " ref val: " << soap_vector_gradients.get_dense_row(ii_pair)(maxRow, maxCol) << " Nb_neigh: "<< center.size() << std::endl;
    if (diff_ii > 1e-12) {
      std::cout << "Ref: " << std::endl<< soap_vector_gradients.get_dense_row(ii_pair).transpose();
      std::cout << std::endl;
      std::cout << "Test: " << std::endl<< soap_vector_gradients_half.get_dense_row(half_center.get_atom_ii()).transpose();
      std::cout << std::endl;
    }
    auto keys_grad_center = soap_vector_gradients[ii_pair].get_keys();
    std::cout << "Center gradient keys: ";
    for (auto key : keys_grad_center) {
      std::cout << "(";
      for (auto key_sp : key) {
        std::cout << key_sp << ", ";
      }
      std::cout << "\b\b) ";
    }
    std::cout << std::endl;

    size_t neigh_count{0};
    for (auto neigh : center) {
      // if (neigh_count >= n_neigh_print) {
      //   break;
      // }
      // auto neigh_type = neigh.get_atom_type();
      // auto tags = neigh.get_atom_tag_list();
      // if (tags[1] <= tags[0]) {continue;}
      // auto half_neigh_it = half_center.begin();
      // for (size_t ii{0}; ii < neigh_count; ii++) {
      //   ++half_neigh_it;
      // }
      // auto half_neigh = *(half_neigh_it);
      // std::cout << "Neighbour "<< neigh_type <<" tags: ";

      // std::cout << "(";
      // for (auto tag : tags) {
      //   std::cout << tag << ", ";
      // }
      // std::cout << "\b\b) ";
      // std::cout << std::endl;

      // auto keys_neigh = soap_vector_gradients[neigh].get_keys();
      // std::cout << "Neighbour "<< neigh_type <<" keys: ";
      // for (auto key : keys_neigh) {
      //   std::cout << "(";
      //   for (auto key_sp : key) {
      //     std::cout << key_sp << ", ";
      //   }
      //   std::cout << "\b\b) ";
      // }
      // std::cout << std::endl;
      // auto keys_neigh_half = soap_vector_gradients_half[half_neigh].get_keys();
      // std::cout << " /// ";
      // for (auto key : keys_neigh_half) {
      //   std::cout << "(";
      //   for (auto key_sp : key) {
      //     std::cout << key_sp << ", ";
      //   }
      //   std::cout << "\b\b) ";
      // }
      // std::cout << std::endl;
      // double diff_ij{math::relative_error(soap_vector_gradients.get_dense_row(neigh), soap_vector_gradients_half.get_dense_row(half_neigh))};
      // if (diff_ij > 1e-12) {
      //   std::cout << "Ref: " << std::endl<< soap_vector_gradients.get_dense_row(neigh).transpose();
      //   std::cout << std::endl;
      //   std::cout << "Test: " << std::endl<< soap_vector_gradients_half.get_dense_row(half_neigh).transpose();
      //   std::cout << std::endl;
      // }

      ++neigh_count;
    }
    ++center_count;
  }
}
