 /**
 * file   representation_manager_soap_covariant.hh
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

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_COVARIANT_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_COVARIANT_HH_

#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "representations/representation_manager_soap_invariant.hh"
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
    enum class SOAPCovariantType { LambdaSpectrum, End_};

    template <SOAPCovariantType SpectrumType>
    struct SOAPCovariantPrecomputation {};

    template <>
    struct SOAPCovariantPrecomputation<SOAPCovariantType::LambdaSpectrum> :SOAPPrecomputationBase {
      using Parent = SOAPPrecomputationBase;
      using Hypers_t = typename SOAPPrecomputationBase::Hypers_t;

      SOAPCovariantPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular").get<size_t>();
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
        this->lambda = hypers.at("lam").get<size_t>();

        // get the number of non zero elements in the w3j
        int n_elements{0};
        size_t l3{this->lambda};
        for (size_t l1{0}; l1 < this->max_angular+1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular+1; ++l2) {
            if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) { continue; }
            if (this->inversion_symmetry == true) {
              if ((l1 + l2 + l3) % 2 == 1) { continue; }
            }
            for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
            int m1s{m1 - l1};
            for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
            int m2s{m2 - l2};
            for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
            int m3s{m3 - l3};
            if (m1s + m2s + m3s != 0 && m1s + m2s - m3s != 0) { continue; }
            ++n_elements;
            }
            }
            }
          }
        }

        this->w3js.resize(n_elements);
        n_elements = 0;
        wig_table_init(2*(this->max_angular + 1), 3);
        wig_temp_init(2*(this->max_angular + 1));
        for (size_t l1{0}; l1 < this->max_angular+1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular+1; ++l2) {
            if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) { continue; }
            if (this->inversion_symmetry == true) {
              if ((l1 + l2 + l3) % 2 == 1) { continue; }
            }
            for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
            int m1s{m1 - l1};
            for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
            int m2s{m2 - l2};
            for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
            int m3s{m3 - l3};
            if (m1s + m2s + m3s != 0 && m1s + m2s - m3s != 0) { continue; }
            double w3j1{wig3jj(2*l1, 2*l2, 2*l3, 2*m1s, 2*m2s, 2*m3s)};
            double w3j2{wig3jj(2*l1, 2*l2, 2*l3, 2*m1s, 2*m2s, -2*m3s)};
            if (m3s > 0) {
              this->w3js(n_elements) = (w3j2 + math::pow(-1, m3s)*w3j1)/std::sqrt(2.0);
            } else if (m3s < 0) {
              this->w3js(n_elements) = ((w3j1 - math::pow(-1, m3s)*w3j2))/std::sqrt(2.0));
            } else if (m3s == 0) {
              this->w3js(n_elements) = w3j1;
            }
            //*/
            //change to the following for agreement with SOAPFAST
            //(different definition of the real spherical harmonics)
            /*
            if (m3s > 0) { this->w3js.push_back((w3j1 + pow(-1, m3s)*w3j2)/sqrt(2.0)); }
            else if (m3s == 0) { this->w3js.push_back(w3j1); }
            else if (m3s < 0) { this->w3js.push_back(((w3j2 - pow(-1, m3s)*w3j1))/sqrt(2.0)); }
            */
            ++n_elements;
            }
            }
            }
          }
        }
        wig_temp_free();
        wig_table_free();
      }

      size_t max_angular{0};
      bool inversion_symmetry{false};
      Eigen::ArrayXd w3js{};
      size_t lambda{2};
    };

  }  // namespace internal

  template <internal::SOAPCovariantType Type, class Hypers>
  decltype(auto) make_soap_covariant_precompute(const Hypers & hypers) {
    return std::static_pointer_cast<internal::SOAPPrecomputationBase>(
        std::make_shared<internal::SOAPCovariantPrecomputation<Type>>(hypers));
  }

  template <internal::SOAPCovariantType Type>
  decltype(auto) downcast_soap_covariant_precompute(
      const std::shared_ptr<internal::SOAPPrecomputationBase> &
          soap_precompute) {
    return std::static_pointer_cast<internal::SOAPCovariantPrecomputation<Type>>(
        soap_precompute);
  }

  template <class StructureManager>
  class RepresentationManagerSOAPCovariant : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using Key_t = std::vector<int>;
    using SparseProperty_t =
        BlockSparseProperty<double, 1, 0, Manager_t, Key_t>;
    using Data_t = typename SparseProperty_t::Data_t;

    RepresentationManagerSOAPCovariant(ManagerPtr_t sm, const Hypers_t & hyper)
        : soap_vectors{*sm}, structure_manager{sm}, rep_expansion{std::move(sm),
                                                                  hyper} {
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    RepresentationManagerSOAPCovariant(const RepresentationManagerSOAPCovariant & other) = delete;

    //! Move constructor
    RepresentationManagerSOAPCovariant(RepresentationManagerSOAPCovariant && other) = default;

    //! Destructor
    virtual ~RepresentationManagerSOAPCovariant() = default;

    //! Copy assignment operator
    RepresentationManagerSOAPCovariant &
    operator=(const RepresentationManagerSOAPCovariant & other) = delete;

    //! Move assignment operator
    RepresentationManagerSOAPCovariant &
    operator=(RepresentationManagerSOAPCovariant && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::SOAPCovariantType;
      this->max_radial = hypers.at("max_radial").get<size_t>();
      this->max_angular = hypers.at("max_angular").get<size_t>();
      this->soap_type_str = hypers.at("soap_type").get<std::string>();
      this->lambda = hypers.at("lam").get<size_t>();
      this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();

      if (this->soap_type_str.compare("LambdaSpectrum") == 0) {
        this->soap_type = SOAPCovariantType::LambdaSpectrum;
        this->precompute_soap[enumValue(SOAPCovariantType::LambdaSpectrum)] =
            make_soap_covariant_precompute<SOAPCovariantType::LambdaSpectrum>(hypers);

      } else {
        throw std::logic_error("Requested SOAP type \'" + this->soap_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'LambdaSpectrum\'.");
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

    //! compute tensor representation with \nu == 2
    void compute_lambdaspectrum();

    void initialize_percenter_lambda_soap_vectors();

    SparseProperty_t soap_vectors;

   protected:
    size_t max_radial{};
    size_t max_angular{};
    ManagerPtr_t structure_manager;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;
    internal::SOAPCovariantType soap_type{};
    std::array<std::shared_ptr<internal::SOAPPrecomputationBase>,
               internal::enumSize<internal::SOAPCovariantType>()>
        precompute_soap{};
    std::string soap_type_str{};
    std::vector<Precision_t> dummy{};
    bool inversion_symmetry{false};
    size_t lambda{0};
  };

  template <class Mngr>
  void RepresentationManagerSOAPCovariant<Mngr>::compute() {
    using internal::SOAPCovariantType;
    switch (this->soap_type) {
    case SOAPCovariantType::LambdaSpectrum:
      this->compute_lambdaspectrum();
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAPCovariant<
      Mngr>::initialize_percenter_lambda_soap_vectors() {
    using math::pow;

    size_t n_row{pow(this->max_radial, 2)};
    // number of combinations of l1 and l2 satisfying the triangle constraint
    size_t n_col{0};
    if (this->inversion_symmetry == false) {
      n_col = static_cast<size_t>((2 + this->lambda - 3*pow(this->lambda, 2) + 2*this->max_angular + 4*this->lambda*this->max_angular)/2*(2*this->lambda + 1));
    }
    else {
      n_col = static_cast<size_t>(
              std::ceil(pow(this->max_angular + 1, 2)/2.0) -
              pow(1.0 + std::floor((this->lambda - 1)/2.0), 2) -
              std::floor(pow(this->max_angular + 1 - this->lambda, 2)/2.0)*
              (this->lambda % 2) -
              (std::ceil(pow(this->max_angular + 1 - this->lambda, 2)/2.0) -
              (this->max_angular - this->lambda + 1))*(1.0 - this->lambda % 2))* (2*this->lambda + 1);
      if (this->lambda % 2 == 1) {
        n_col = static_cast<size_t>(
                0.5*(2.0 + this->lambda - 3.0*pow(this->lambda, 2) +
                2*this->max_angular + 4*this->lambda*this->max_angular)*
                (2*this->lambda + 1) - static_cast<int>(n_col));
      }
    }

    // clear the data container and resize it
    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    // identify the species in each environment and initialize soap_vectors
    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      internal::Sorted<false> is_not_sorted{};

      std::vector<internal::SortedKey<Key_t>> pair_list{};
      auto & center_type{center.get_atom_type()};
      Key_t pair_type{center_type, center_type};
      // TODO(felix) optimize this loop
      Key_t pair_type{0, 0};
      for (const auto& el1: coefficients) {
        pair_type[0] = el1.first[0];
        for (const auto& el2: coefficients) {
          pair_type[1] = el2.first[0];
          pair_list.emplace_back(is_not_sorted, pair_type);
          }
        }
      }
      // initialize the power spectrum with the proper dimension
      soap_vector.resize(pair_list, n_row, n_col);
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAPCovariant<Mngr>::compute_lambdaspectrum() {
    using math::pow;
    using internal::SOAPCovariantType;
    using complex = std::complex<double>;

    rep_expansion.compute();
    auto& expansions_coefficients{rep_expansion.expansions_coefficients};

    this->initialize_percenter_lambda_soap_vectors();

    auto precomputation{downcast_soap_covariant_precompute<SOAPCovariantType::LambdaSpectrum>(
        this->precompute_soap[enumValue(SOAPCovariantType::LambdaSpectrum)])};
    auto & w3js{precomputation->w3js};

    for (auto center : this->structure_manager) {
      auto& coefficients{expansions_coefficients[center]};
      auto& soap_vector{this->soap_vectors[center]};
      Key_t pair_type{0, 0};
      for (const auto& el1: coefficients) {
        pair_type[0] = el1.first[0];
        auto& coef1{el1.second};
        for (const auto& el2: coefficients) {
          pair_type[1] = el2.first[0];
          auto& coef2{el2.second};

          size_t& l3{this->lambda};
          if (soap_vector.count(pair_type) == 0) {
            auto&& soap_vector_by_type{soap_vector[pair_type]};

            size_t nn{0};
            for (size_t n1{0}; n1 < this->max_radial; n1++) {
              for (size_t n2{0}; n2 < this->max_radial; n2++) {
                size_t l0{0};
                int count{0};
                for (size_t l1{0}; l1 < this->max_angular+1; l1++) {
                  for (size_t l2{0}; l2 < this->max_angular+1; l2++) {
                    if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) { continue; }
                    if (this->inversion_symmetry == true) {
                      if ((l1 + l2 + l3) % 2 == 1) { continue; }
                    }
                    for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
                      int m3s{m3 - l3};
                      for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
                        int m1s{m1 - l1};
                        int lm1{pow(l1, 2) + m1};
                        for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
                          int m2s{m2 - l2};
                          if (m1s + m2s + m3s != 0 && m1s + m2s - m3s != 0) { continue; }
                          int lm2{pow(l2, 2) + m2};
                          complex coef1c, coef2c;
                          double& w3j = w3js[count];
                          // usual formulae for converting from real to complex
                          if (m1s > 0) {
                            coef1c = pow(-1.0, m1s)*
                                     complex(coef1(n1, lm1),
                                     coef1(n1, lm1 - 2*m1s));
                          } else if (m1s == 0) {
                            coef1c = complex(coef1(n1, lm1), 0.0)*
                                     std::sqrt(2.0);
                          } else if (m1s < 0) {
                            coef1c = complex(coef1(n1, lm1 - 2*m1s),
                                     -coef1(n1, lm1));
                          }
                          if (m2s > 0) {
                            coef2c = pow(-1.0, m2s)*
                                     complex(coef2(n2, lm2),
                                     coef2(n2, lm2 - 2*m2s));
                          } else if (m2s == 0) {
                            coef2c = complex(coef2(n2, lm2), 0.0)*
                                     std::sqrt(2.0);
                          } else if (m2s < 0) {
                            coef2c = complex(coef2(n2, lm2 - 2*m2s),
                                     -coef2(n2, lm2));
                          }
                          coef1c /= std::sqrt(2.0);
                          coef2c /= std::sqrt(2.0);
                          //combine the coefficients with Wigner 3j symbols
                          complex ci = complex(0.0, 1.0);
                          if ((l1 + l2 + l3) % 2 == 0) {
                            if (m3s < 0) {
                              soap_vector_by_type(nn, l0) += w3j*(ci*coef1c*coef2c).real();
                            } else {
                              soap_vector_by_type(nn, l0) += w3j*(coef1c*coef2c).real();
                            }
                          } else if (this->inversion_symmetry == false) {
                            if (m3s < 0) {
                              soap_vector_by_type(nn, l0) += w3j*(ci*coef1c*coef2c).imag();
                            } else {
                              soap_vector_by_type(nn, l0) += w3j*(coef1c*coef2c).imag();
                            }
                          }
                          count++;
                        } // m2
                      } // m1
                      l0++;
                    }  // m3
                  } // l2
                } // l1
                nn++;
              } // n2
            } // n1
          } // if pair_type count
        } // coef1
      } // coef1

      // the SQRT_TWO factor comes from the fact that
      // the upper diagonal of the species is not considered
      soap_vector.multiply_offdiagonal_elements_by(math::SQRT_TWO);

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }

    } // center
  } // compute_lambdaspectrum

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_COVARIANT_HH_
