/**
 * file   representation_manager_soap.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Michael Willatt <michael.willatt@epfl.ch>
 * @author Andrea Grisafi <andrea.grisafi@epfl.ch>
 *
 * @date   12 March 2019
 *
 * @brief  Compute the spherical harmonics expansion of the local atom density
 *
 * Copyright © 2019 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_

#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "structure_managers/property_block_sparse.hh"
#include "rascal_utility.hh"
#include "math/math_utils.hh"

#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unordered_set>

namespace rascal {

  namespace internal {
    enum class SOAPType { RadialSpectrum, PowerSpectrum };
  }  // namespace internal

  template <class StructureManager>
  class RepresentationManagerSOAP : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using Key_t = std::vector<int>;
    using SparseProperty_t = BlockSparseProperty<double, 1, 0>;
    using Data_t = typename SparseProperty_t::Data_t;

    RepresentationManagerSOAP(ManagerPtr_t sm, const Hypers_t & hyper)
        : soap_vectors{*sm}, structure_manager{sm}, rep_expansion{std::move(sm),
                                                                  hyper} {
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    RepresentationManagerSOAP(const RepresentationManagerSOAP & other) = delete;

    //! Move constructor
    RepresentationManagerSOAP(RepresentationManagerSOAP && other) = default;

    //! Destructor
    virtual ~RepresentationManagerSOAP() = default;

    //! Copy assignment operator
    RepresentationManagerSOAP &
    operator=(const RepresentationManagerSOAP & other) = delete;

    //! Move assignment operator
    RepresentationManagerSOAP &
    operator=(RepresentationManagerSOAP && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      this->normalize = hypers.at("normalize").get<bool>();
      this->soap_type_str = hypers.at("soap_type").get<std::string>();

      if (this->soap_type_str.compare("PowerSpectrum") == 0) {
        this->soap_type = internal::SOAPType::PowerSpectrum;
      } else if (this->soap_type_str.compare("RadialSpectrum") == 0) {
        this->soap_type = internal::SOAPType::RadialSpectrum;
        if (this->max_angular > 0) {
          throw std::logic_error("max_angular should be 0 with RadialSpectrum");
        }
      } else {
        throw std::logic_error("Requested SOAP type \'" + this->soap_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'PowerSpectrum or RadialSpectrum\'.");
      }
    }

    std::vector<Precision_t> & get_representation_raw_data() {
      return this->dummy;
    }

    Data_t & get_representation_sparse_raw_data() {
      return this->soap_vectors.get_raw_data();
    }

    size_t get_feature_size() { return this->soap_vectors.get_nb_comp(); }

    size_t get_center_size() { return this->soap_vectors.get_nb_item(); }

    auto get_representation_full() {
      return this->soap_vectors.get_dense_rep();
    }

    //! compute representation
    void compute();

    //! compute representation \nu == 1
    void compute_radialspectrum();

    //! compute representation \nu == 2
    void compute_powerspectrum();

    SparseProperty_t soap_vectors;

   protected:
    size_t max_radial{};
    size_t max_angular{};
    bool normalize{};
    ManagerPtr_t structure_manager;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;
    internal::SOAPType soap_type{};
    std::string soap_type_str{};
    std::vector<Precision_t> dummy{};
  };

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute() {
    using internal::SOAPType;
    switch (this->soap_type) {
    case SOAPType::RadialSpectrum:
      this->compute_radialspectrum();
      break;
    case SOAPType::PowerSpectrum:
      this->compute_powerspectrum();
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  /*//////////////////////////////////////////////////////////////////////////////

    template <class Mngr>
    void RepresentationManagerSOAP<Mngr>::compute_bispectrum() {
      rep_expansion.compute();
      auto& expansions_coefficients{rep_expansion.expansions_coefficients};

      size_t n_row{pow(this->max_radial, 3)}
      size_t n_col{*pow(this->max_angular, 3)*pow((2*this->max_angular + 1),
    3)};

      this->soap_vectors.clear();
      this->soap_vectors.set_shape(n_row, n_col);
      this->soap_vectors.resize();

      for (auto center : this->structure_manager) {
        auto& coefficients{expansions_coefficients[center]};
        auto& soap_vector{this->soap_vectors[center]};
        Key_t triplet_type{0, 0, 0};
        for (const auto& el1: coefficients) {
          triplet_type[0] = el1.first[0];
          auto& coef1{el1.second};
          for (const auto& el2: coefficients) {
            triplet_type[1] = el2.first[0];
            auto& coef2{el2.second};
            for (const auto& el2: coefficients) {
              triplet_type[2] = el3.first[0];
              auto& coef3{el3.second};

              if (soap_vector.count(triplet_type) == 0) {
                soap_vector[triplet_type] = Dense_t::Zero(n_row, n_col);
              }
            }
          }
        }
      }

    }

  */  //////////////////////////////////////////////////////////////////////////////

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_powerspectrum() {
    using math::pow;
    // Compute the spherical expansions of the current structure
    rep_expansion.compute();
    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    size_t n_row{static_cast<size_t>(pow(this->max_radial, 2))};
    size_t n_col{this->max_angular + 1};

    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      Key_t pair_type{0, 0};

      std::unordered_set<Key_t, internal::Hash<Key_t>> pair_list{};
      auto & center_type{center.get_atom_type()};
      pair_list.insert({center_type, center_type});
      for (auto neigh1 : center) {
        auto && neigh1_type{neigh1.get_atom_type()};
        pair_list.insert({center_type, neigh1_type});
        for (auto neigh2 : center) {
          auto && neigh2_type{neigh2.get_atom_type()};
          if (neigh1_type <= neigh2_type) {
            pair_list.insert({neigh1_type, neigh2_type});
          }
        }
      }
      // initialize the power spectrum to 0 with the proper dimension
      soap_vector.resize(pair_list, n_row, n_col);

      for (const auto & el1 : coefficients) {
        pair_type[0] = el1.first[0];
        auto & coef1{el1.second};
        for (const auto & el2 : coefficients) {
          pair_type[1] = el2.first[0];
          auto & coef2{el2.second};
          /* avoid computing p^{ab} and p^{ba} since p^{ab} = p^{ba}^T
           */
          if (pair_type[0] > pair_type[1]) {
            continue;
          }

          auto && soap_vector_by_pair{soap_vector[pair_type]};

          size_t lm{0};
          for (size_t l = 0; l < this->max_angular + 1; l++) {
            double l_factor{pow(std::sqrt(2 * l + 1), -1)};
            for (size_t m = 0; m < 2 * l + 1; m++) {
              size_t nn{0};
              for (size_t n1 = 0; n1 < this->max_radial; n1++) {
                for (size_t n2 = 0; n2 < this->max_radial; n2++) {
                  double val{l_factor * coef1(n1, lm) * coef2(n2, lm)};
                  soap_vector_by_pair(nn, l) += val;
                  nn++;
                }
              }
              lm++;
            }
          }
        }
      }

      // the SQRT_TWO factor comes from the fact that
      // the upper diagonal of the species is not considered
      for (const auto & el : soap_vector) {
        auto && pair_type{el.first};
        if (pair_type[0] != pair_type[1]) {
          soap_vector[pair_type] *= math::SQRT_TWO;
        }
      }

      // normalize the soap vector
      if (this->normalize) {
        double norm{0.};
        for (const auto & el : soap_vector) {
          auto && pair_type{el.first};
          norm += soap_vector[pair_type].squaredNorm();
        }
        norm = std::sqrt(norm);
        for (const auto & el : soap_vector) {
          auto && pair_type{el.first};
          soap_vector[pair_type] /= norm;
        }
      }
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_radialspectrum() {
    rep_expansion.compute();
    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    size_t n_row{this->max_radial};
    size_t n_col{1};

    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      Key_t element_type{0};

      std::unordered_set<Key_t, internal::Hash<Key_t>> keys{};
      for (auto neigh : center) {
        keys.insert({neigh.get_atom_type()});
      }
      keys.insert({center.get_atom_type()});
      // initialize the radial spectrum to 0 and the proper size
      soap_vector.resize(keys, n_row, n_col);

      std::vector<Key_t> element_list{};
      for (const auto & el : coefficients) {
        element_type[0] = el.first[0];
        auto & coef{el.second};
        soap_vector[element_type] += coef;
        element_list.push_back(element_type);
      }

      // normalize the soap vector
      if (this->normalize) {
        double norm{0.};
        for (const auto & element_type : element_list) {
          norm += soap_vector[element_type].squaredNorm();
        }
        norm = std::sqrt(norm);
        for (const auto & element_type : element_list) {
          soap_vector[element_type] /= norm;
        }
      }
    }
  }
}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
