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


namespace rascal {

  namespace internal {
    enum class SOAPType { RadialSpectrum, PowerSpectrum, BiSpectrum };
  }  // namespace internal

  template <class StructureManager>
  class RepresentationManagerSOAP : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using key_t = std::vector<int>;
    using SparseProperty_t = BlockSparseProperty<double, 1, 0>;
    using data_t = typename SparseProperty_t::data_t;

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
      this->soap_type_str = hypers.at("soap_type").get<std::string>();

      if (this->soap_type_str.compare("PowerSpectrum") == 0) {
        this->soap_type = internal::SOAPType::PowerSpectrum;
      } else if (this->soap_type_str.compare("RadialSpectrum") == 0) {
        this->soap_type = internal::SOAPType::RadialSpectrum;
      } else if (this->soap_type_str.compare("BiSpectrum") == 0) {
        this->soap_type = internal::SOAPType::BiSpectrum;
        if (hypers.find("inversion_symmetry") != hypers.end()) {
          this->inversion_symmetry = hypers.at("inversion_symmetry");
        }
      } else {
        throw std::logic_error("Requested SOAP type \'" + this->soap_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'PowerSpectrum\'.");
      }
    }

    std::vector<Precision_t> & get_representation_raw_data() {
      return this->dummy;
    }

    data_t & get_representation_sparse_raw_data() {
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

    //! compute representation \nu == 3
    void compute_bispectrum();

    //! precompute the Wigner 3j symbols
    void precompute_w3js();

    SparseProperty_t soap_vectors;

   protected:
    size_t max_radial{};
    size_t max_angular{};
    ManagerPtr_t structure_manager;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;
    internal::SOAPType soap_type{};
    std::string soap_type_str{};
    std::vector<Precision_t> dummy{};
    bool is_precomputed{false};
    bool inversion_symmetry{true};
    std::vector<double> w3js{};
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
    case SOAPType::BiSpectrum:
      this->compute_bispectrum();
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  /** Compute Wigner 3j symbols */
  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::precompute_w3js() {
    //2*lmax and Wigner symbol type (3)
    wig_table_init(2*(this->max_angular + 1), 3);
    wig_temp_init(2*(this->max_angular + 1));
    for (size_t l1 = 0; l1 < this->max_angular+1; l1++) {
      for (size_t l2 = 0; l2 < this->max_angular+1; l2++) {
        for (size_t l3 = 0; l3 < this->max_angular+1; l3++) {
          if (l1 < (size_t)std::abs<int>(l2 - l3) || l1 > l2 + l3) { continue; }
          if (this->inversion_symmetry == true) {
            if ((l1 + l2 + l3) % 2 == 1) { continue; }
          }
          for (size_t m1 = 0; m1 < 2*l1 + 1; m1++) {
          int m1s = m1 - l1;
          for (size_t m2 = 0; m2 < 2*l2 + 1; m2++) {
          int m2s = m2 - l2;
          for (size_t m3 = 0; m3 < 2*l3 + 1; m3++) {
          int m3s = m3 - l3;
          if (m1s + m2s + m3s != 0) { continue; }
          this->w3js.push_back(wig3jj(2*l1, 2*l2, 2*l3, 2*m1s, 2*m2s, 2*m3s));
          }
          }
          }
        }
      }
    }
    wig_temp_free();
    wig_table_free();
    this->is_precomputed = true;
  }

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_bispectrum() {
    rep_expansion.compute();
    auto& expansions_coefficients{rep_expansion.expansions_coefficients};
    using complex = std::complex<double>;
    size_t n_row{(size_t)pow(this->max_radial, 3)};
    size_t n_col{0};
    //size_t n_col{pow((this->max_angular + 1), 3)};
    if (this->inversion_symmetry == false) {
      /* sum of the next lmax + 1 natural numbers */
      n_col = (size_t)(1.0 + 2.0*(double)this->max_angular + \
              1.5*(double)pow(this->max_angular, 2) + \
              (double)pow(this->max_angular, 3)/2.0);
    }
    else {
      /* number of 3x3 symmetric matrices with non-negative integer entries,
       * such that the sum of every row and column equal lmax */
      n_col = (size_t)std::floor((double)((pow(this->max_angular + 1, 2) + 1)* \
              (2*(this->max_angular + 1) + 3))/8.0);
    }

    double mult{1.0};

    if (this->is_precomputed == false) { this->precompute_w3js(); }
    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    for (auto center : this->structure_manager) {
      auto& coefficients{expansions_coefficients[center]};
      auto& soap_vector{this->soap_vectors[center]};
      key_t triplet_type{0, 0, 0};
      for (const auto& el1: coefficients) {
        triplet_type[0] = el1.first[0];
        auto& coef1{el1.second};
        for (const auto& el2: coefficients) {
          triplet_type[1] = el2.first[0];
          auto& coef2{el2.second};
          for (const auto& el3: coefficients) {
            triplet_type[2] = el3.first[0];
            auto& coef3{el3.second};

            //triplet multiplicity
            //all the same
            if (triplet_type[0] == triplet_type[1] && \
                triplet_type[1] == triplet_type[2]) {
              mult = 1.0;
            }
            //two the same
            else if (triplet_type[0] == triplet_type[1] || \
                     triplet_type[0] == triplet_type[2] || \
                     triplet_type[1] == triplet_type[2]) {
              mult = std::sqrt(3.0);
            }
            //all different
            else {
              mult = std::sqrt(6.0);
            }

            if (soap_vector.count(triplet_type) == 0) {
              soap_vector[triplet_type] = dense_t::Zero(n_row, n_col);

              size_t nn{0};
              for (size_t n1 = 0; n1 < this->max_radial; n1++) {
                for (size_t n2 = 0; n2 < this->max_radial; n2++) {
                  for (size_t n3 = 0; n3 < this->max_radial; n3++) {
                    size_t l0{0};
                    int count{0};
                    for (size_t l1 = 0; l1 < this->max_angular+1; l1++) {
                      for (size_t l2 = 0; l2 < this->max_angular+1; l2++) {
                        for (size_t l3 = 0; l3 < this->max_angular+1; l3++) {
                          if (l1 < (size_t)std::abs<int>(l2 - l3) || l1 > l2 + l3) {
                            continue;
                          }
                          if ((l1 + l2 + l3) % 2 == 1) { continue; }
                          for (size_t m1 = 0; m1 < 2*l1 + 1; m1++) {
                          int m1s = m1 - l1;
                          int lm1 = std::pow(l1, 2) + m1;
                          for (size_t m2 = 0; m2 < 2*l2 + 1; m2++) {
                          int m2s = m2 - l2;
                          int lm2 = std::pow(l2, 2) + m2;
                          for (size_t m3 = 0; m3 < 2*l3 + 1; m3++) {
                          int m3s = m3 - l3;
                          if (m1s + m2s + m3s != 0) { continue; }
                          int lm3 = std::pow(l3, 2) + m3;
                          double w3j = w3js[count];
                          complex coef1c, coef2c, coef3c;

                          //usual formulae for converting from real to complex
                          if (m1s > 0) {
                            coef1c = std::pow(-1.0, m1s)* \
                                     complex(coef1(n1, lm1), \
                                     coef1(n1, lm1 - 2*m1s));
                          }
                          else if (m1s == 0) {
                            coef1c = complex(coef1(n1, lm1), 0.0)* \
                                     std::sqrt(2.0);
                          }
                          else if (m1s < 0) {
                            coef1c = complex(coef1(n1, lm1 - 2*m1s), \
                                     -coef1(n1, lm1));
                          }
                          if (m2s > 0) {
                            coef2c = std::pow(-1.0, m2s)* \
                                     complex(coef2(n2, lm2), \
                                     coef2(n2, lm2 - 2*m2s));
                          }
                          else if (m2s == 0) {
                            coef2c = complex(coef2(n2, lm2), 0.0)* \
                                     std::sqrt(2.0);
                          }
                          else if (m2s < 0) {
                            coef2c = complex(coef2(n2, lm2 - 2*m2s), \
                                     -coef2(n2, lm2));
                          }
                          if (m3s > 0) {
                            coef3c = std::pow(-1.0, m3s)* \
                                     complex(coef3(n3, lm3), \
                                     coef3(n3, lm3 - 2*m3s));
                          }
                          else if (m3s == 0) {
                            coef3c = complex(coef3(n3, lm3), 0.0)* \
                                     std::sqrt(2.0);
                          }
                          else if (m3s < 0) {
                            coef3c = complex(coef3(n3, lm3 - 2*m3s), \
                                     -coef3(n3, lm3));
                          }
                          coef1c /= std::sqrt(2.0);
                          coef2c /= std::sqrt(2.0);
                          coef3c /= std::sqrt(2.0);

                          soap_vector[triplet_type](nn, l0) += \
                            w3j*mult*(coef1c*coef2c*coef3c).real();
                          count++;
                          }
                          }
                          }
                          l0++;
                        }
                      }
                    }
                    nn++;
                  }
                }
              }

            }

          }
        }
      }
    }

  }

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_powerspectrum() {
    rep_expansion.compute();
    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    size_t n_row{(size_t)pow(this->max_radial, 2)};
    size_t n_col{this->max_angular + 1};

    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      key_t pair_type{0, 0};

      for (const auto & el1 : coefficients) {
        pair_type[0] = el1.first[0];
        auto & coef1{el1.second};
        for (const auto & el2 : coefficients) {
          pair_type[1] = el2.first[0];
          auto & coef2{el2.second};
          /* avoid computing p^{ab} and p^{ba} since p^{ab} = p^{ba}^T
           */
          if (soap_vector.count(pair_type) == 0) {
            soap_vector[pair_type] = dense_t::Zero(n_row, n_col);
            size_t nn{0};
            for (size_t n1 = 0; n1 < this->max_radial; n1++) {
              for (size_t n2 = 0; n2 < this->max_radial; n2++) {
                size_t lm{0};
                for (size_t l = 0; l < this->max_angular + 1; l++) {
                  // TODO(andrea) pre compute l_factor
                  double l_factor{1 / std::sqrt(2 * l + 1)};
                  for (size_t m = 0; m < 2 * l + 1; m++) {
                    soap_vector[pair_type](nn, l) +=
                        l_factor * coef1(n1, lm) * coef2(n2, lm);
                    lm++;
                  }
                }
                nn++;
              }
            }
          }
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
      key_t element_type{0};
      for (const auto & el : coefficients) {
        element_type[0] = el.first[0];
        soap_vector[element_type] = dense_t::Zero(n_row, n_col);
        auto & coef{el.second};
        soap_vector[element_type] += coef;
      }
    }
  }
}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
