 /**
 * file   representation_manager_soap_invariant.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Michael Willatt <michael.willatt@epfl.ch>
 * @author Andrea Grisafi <andrea.grisafi@epfl.ch>
 *
 * @date   12 March 2019
 *
 * @brief  compute spherical invariants
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

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_INVARIANT_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_INVARIANT_HH_

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
    enum class SOAPType { RadialSpectrum, PowerSpectrum, BiSpectrum, End_ };

    /**
     * Base class for the specification of the atomic smearing.
     */
    struct SOAPPrecomputationBase {
      //! Constructor
      SOAPPrecomputationBase() = default;
      //! Destructor
      virtual ~SOAPPrecomputationBase() = default;
      //! Copy constructor
      SOAPPrecomputationBase(const SOAPPrecomputationBase & other) = delete;
      //! Move constructor
      SOAPPrecomputationBase(SOAPPrecomputationBase && other) = default;
      //! Copy assignment operator
      SOAPPrecomputationBase &
      operator=(const SOAPPrecomputationBase & other) = delete;
      //! Move assignment operator
      SOAPPrecomputationBase &
      operator=(SOAPPrecomputationBase && other) = default;

      using Hypers_t = RepresentationManagerBase::Hypers_t;
    };

    template <SOAPType SpectrumType>
    struct SOAPPrecomputation {};

    template <>
    struct SOAPPrecomputation<SOAPType::RadialSpectrum>
        : SOAPPrecomputationBase {
      using Hypers_t = typename SOAPPrecomputationBase::Hypers_t;
      explicit SOAPPrecomputation(const Hypers_t &) {}
    };

    template <>
    struct SOAPPrecomputation<SOAPType::PowerSpectrum>
        : SOAPPrecomputationBase {
      using Hypers_t = typename SOAPPrecomputationBase::Hypers_t;
      explicit SOAPPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular");
        this->l_factors.resize(math::pow(this->max_angular + 1, 2_z));

        size_t lm{0};
        for (size_t l{0}; l < this->max_angular + 1; ++l) {
          double l_factor{math::pow(std::sqrt(2 * l + 1), -1)};
          for (size_t m{0}; m < 2 * l + 1; ++m) {
            this->l_factors(lm) = l_factor;
            ++lm;
          }
        }
      }

      size_t max_angular{0};
      //! factor of 1 / sqrt(2*l+1) in front of the powerspectrum
      Eigen::VectorXd l_factors{};
    };

    template <>
    struct SOAPPrecomputation<SOAPType::BiSpectrum>
        : SOAPPrecomputationBase {
      using Hypers_t = typename SOAPPrecomputationBase::Hypers_t;
      explicit SOAPPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular").get<size_t>();
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
        // get the number of non zero elements in the w3j
        int n_elements{0};
        for (size_t l1{0}; l1 < this->max_angular+1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular+1; ++l2) {
            for (size_t l3{0}; l3 < this->max_angular+1; ++l3) {
              if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) { continue; }
              if (this->inversion_symmetry == true) {
                if ((l1 + l2 + l3) % 2 == 1) { continue; }
              }
              for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
              int m1s{static_cast<int>(m1 - l1)};
              for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
              int m2s{static_cast<int>(m2 - l2)};
              for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
              int m3s{static_cast<int>(m3 - l3)};
              if (m1s + m2s + m3s != 0) { continue; }
              ++n_elements;
              }
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
            for (size_t l3{0}; l3 < this->max_angular+1; ++l3) {
              if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) { continue; }
              if (this->inversion_symmetry == true) {
                if ((l1 + l2 + l3) % 2 == 1) { continue; }
              }
              for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
              int m1s{static_cast<int>(m1 - l1)};
              for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
              int m2s{static_cast<int>(m2 - l2)};
              for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
              int m3s{static_cast<int>(m3 - l3)};
              if (m1s + m2s + m3s != 0) { continue; }
              this->w3js(n_elements) = wig3jj(2*l1, 2*l2, 2*l3, 2*m1s, 2*m2s, 2*m3s);
              ++n_elements;
              }
              }
              }
            }
          }
        }
        wig_temp_free();
        wig_table_free();
      }

      size_t max_angular{0};
      bool inversion_symmetry{};
      Eigen::ArrayXd w3js{};
    };

  }  // namespace internal

  template <internal::SOAPType Type, class Hypers>
  decltype(auto) make_soap_precompute(const Hypers & hypers) {
    return std::static_pointer_cast<internal::SOAPPrecomputationBase>(
        std::make_shared<internal::SOAPPrecomputation<Type>>(hypers));
  }

  template <internal::SOAPType Type>
  decltype(auto) downcast_soap_precompute(
      const std::shared_ptr<internal::SOAPPrecomputationBase> &
          soap_precompute) {
    return std::static_pointer_cast<internal::SOAPPrecomputation<Type>>(
        soap_precompute);
  }

  template <class StructureManager>
  class RepresentationManagerSOAPInvariant : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using Key_t = std::vector<int>;
    using SparseProperty_t =
        BlockSparseProperty<double, 1, 0, Manager_t, Key_t>;
    using Data_t = typename SparseProperty_t::Data_t;

    RepresentationManagerSOAPInvariant(ManagerPtr_t sm, const Hypers_t & hyper)
        : soap_vectors{*sm}, structure_manager{sm}, rep_expansion{std::move(sm),
                                                                  hyper} {
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    RepresentationManagerSOAPInvariant(const RepresentationManagerSOAPInvariant & other) = delete;

    //! Move constructor
    RepresentationManagerSOAPInvariant(RepresentationManagerSOAPInvariant && other) = default;

    //! Destructor
    virtual ~RepresentationManagerSOAPInvariant() = default;

    //! Copy assignment operator
    RepresentationManagerSOAPInvariant &
    operator=(const RepresentationManagerSOAPInvariant & other) = delete;

    //! Move assignment operator
    RepresentationManagerSOAPInvariant &
    operator=(RepresentationManagerSOAPInvariant && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::enumValue;
      using internal::SOAPType;

      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      this->normalize = hypers.at("normalize").get<bool>();
      this->soap_type_str = hypers.at("soap_type").get<std::string>();

      if (this->soap_type_str.compare("PowerSpectrum") == 0) {
        this->soap_type = SOAPType::PowerSpectrum;
        this->precompute_soap[enumValue(SOAPType::PowerSpectrum)] =
            make_soap_precompute<SOAPType::PowerSpectrum>(hypers);
      } else if (this->soap_type_str.compare("RadialSpectrum") == 0) {
        this->soap_type = SOAPType::RadialSpectrum;
        this->precompute_soap[enumValue(SOAPType::RadialSpectrum)] =
            make_soap_precompute<SOAPType::RadialSpectrum>(hypers);
        if (this->max_angular > 0) {
          throw std::logic_error("max_angular should be 0 with RadialSpectrum");
        }
      } else if (this->soap_type_str.compare("BiSpectrum") == 0) {
        this->soap_type = internal::SOAPType::BiSpectrum;
        this->inversion_symmetry = hypers.at("inversion_symmetry");

      } else {
        throw std::logic_error("Requested SOAP type \'" + this->soap_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'RadialSpectrum\', \'PowerSpectrum\', " +
                               "\'BiSpectrum\'.");
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

    //! compute representation \nu == 3
    void compute_bispectrum();

    SparseProperty_t soap_vectors;

    //! initialize the soap vectors with only the keys needed for each center
    void initialize_percenter_powerspectrum_soap_vectors();

    void initialize_percenter_radialspectrum_soap_vectors();

    void initialize_percenter_bispectrum_soap_vectors();

   protected:
    size_t max_radial{};
    size_t max_angular{};
    bool normalize{};
    ManagerPtr_t structure_manager;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;
    internal::SOAPType soap_type{};
    //! collection of precomputation for the different body order
    std::array<std::shared_ptr<internal::SOAPPrecomputationBase>,
               internal::enumSize<internal::SOAPType>()>
        precompute_soap{};
    std::string soap_type_str{};
    std::vector<Precision_t> dummy{};
    bool inversion_symmetry{false};
  };

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<Mngr>::compute() {
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

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<Mngr>::compute_powerspectrum() {
    using internal::enumValue;
    using internal::SOAPType;
    using math::pow;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_soap_precompute<SOAPType::PowerSpectrum>(
        this->precompute_soap[enumValue(SOAPType::PowerSpectrum)])};
    auto & l_factors{precomputation->l_factors};

    // Compute the spherical expansions of the current structure
    rep_expansion.compute();
    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    this->initialize_percenter_powerspectrum_soap_vectors();

    Key_t pair_type{0, 0};
    // use special container to tell that there is not need to sort when
    // using operator[] of soap_vector
    internal::SortedKey<Key_t> spair_type{pair_type};

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};

      for (const auto & el1 : coefficients) {
        spair_type[0] = el1.first[0];

        // multiply with the precomputed factors
        auto coef1{el1.second * l_factors.asDiagonal()};

        for (const auto & el2 : coefficients) {
          // avoid computing p^{ab} and p^{ba} since p^{ab} = p^{ba}^T
          if (spair_type[0] > el2.first[0]) {
            continue;
          }
          spair_type[1] = el2.first[0];
          auto & coef2{el2.second};
          auto && soap_vector_by_pair{soap_vector[spair_type]};

          size_t n1n2{0};
          size_t pos{0}, size{0};
          for (size_t n1{0}; n1 < this->max_radial; ++n1) {
            for (size_t n2{0}; n2 < this->max_radial; ++n2) {
              soap_vector_by_pair(n1n2, 0) = coef1(n1, 0) * coef2(n2, 0);
              pos = 1;
              for (size_t l{1}; l < this->max_angular + 1; ++l) {
                size = 2 * l + 1;
                // do the reduction over m (with vectorization)
                soap_vector_by_pair(n1n2, l) =
                    (coef1.block(n1, pos, 1, size).array() *
                     coef2.block(n2, pos, 1, size).array())
                        .sum();
                pos += size;
              }
              ++n1n2;
            }
          }
        }  // for coefficients
      }    // for coefficients

      // the SQRT_TWO factor comes from the fact that
      // the upper diagonal of the species is not considered
      soap_vector.multiply_offdiagonal_elements_by(math::SQRT_TWO);

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<Mngr>::compute_radialspectrum() {
    rep_expansion.compute();
    using math::pow;

    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    this->initialize_percenter_radialspectrum_soap_vectors();
    Key_t element_type{0};

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};

      for (const auto & el : coefficients) {
        element_type[0] = el.first[0];
        auto & coef{el.second};
        soap_vector[element_type] += coef;
      }

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }
    }
  }


  template <class Mngr>
  void RepresentationManagerSOAPInvariant<Mngr>::compute_bispectrum() {
    using internal::SOAPType;
    using internal::enumValue;

    auto precomputation{downcast_soap_precompute<SOAPType::BiSpectrum>(
        this->precompute_soap[enumValue(SOAPType::BiSpectrum)])};
    auto & w3js{precomputation->w3js};

    rep_expansion.compute();

    this->initialize_percenter_bispectrum_soap_vectors();

    auto& expansions_coefficients{rep_expansion.expansions_coefficients};

    using complex = std::complex<double>;

    // factor that takes into acount the missing equivalent off diagonal
    // element with respect to the key (or species) index
    double mult{1.0};

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
          for (const auto& el3: coefficients) {
            triplet_type[2] = el3.first[0];
            auto& coef3{el3.second};
            // triplet multiplicity
            // all the same
            if (triplet_type[0] == triplet_type[1] && \
                triplet_type[1] == triplet_type[2]) {
              mult = 1.0;
            }
            // two the same
            else if (triplet_type[0] == triplet_type[1] || \
                     triplet_type[0] == triplet_type[2] || \
                     triplet_type[1] == triplet_type[2]) {
              mult = std::sqrt(3.0);
            }
            // all different
            else {
              mult = std::sqrt(6.0);
            }

            if (soap_vector.count(triplet_type) == 0) {
              auto && soap_vector_by_type{soap_vector[triplet_type]};

              size_t nn{0};
              for (size_t n1{0}; n1 < this->max_radial; n1++) {
                for (size_t n2{0}; n2 < this->max_radial; n2++) {
                  for (size_t n3{0}; n3 < this->max_radial; n3++) {
                    size_t l0{0};
                    int count{0};
                    for (size_t l1{0}; l1 < this->max_angular+1; l1++) {
                      for (size_t l2{0}; l2 < this->max_angular+1; l2++) {
                        for (size_t l3{0}; l3 < this->max_angular+1; l3++) {
                          if (this->inversion_symmetry == true) {
                            if ((l1 + l2 + l3) % 2 == 1) { continue; }
                          }

                          if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) || l1 > l2 + l3) {
                            continue;
                          }

                          for (size_t m1{0}; m1 < 2*l1 + 1; m1++) {
                          int m1s{static_cast<int>(m1 - l1)};
                          size_t lm1{math::pow(l1, 2_z) + m1};
                          for (size_t m2{0}; m2 < 2*l2 + 1; m2++) {
                          int m2s{static_cast<int>(m2 - l2)};
                          size_t lm2{math::pow(l2, 2_z) + m2};
                          for (size_t m3{0}; m3 < 2*l3 + 1; m3++) {
                          int m3s{static_cast<int>(m3 - l3)};
                          if (m1s + m2s + m3s != 0) { continue; }
                          size_t lm3{math::pow(l3, 2_z) + m3};
                          double w3j = w3js[count];
                          complex coef1c, coef2c, coef3c;
                          // usual formulae for converting from real to complex
                          if (m1s > 0) {
                            coef1c = math::pow(-1, m1s)* \
                                     complex(coef1(n1, lm1), \
                                     coef1(n1, lm1 - 2*m1s));
                          } else if (m1s == 0) {
                            coef1c = complex(coef1(n1, lm1), 0.0)* \
                                     std::sqrt(2.0);
                          } else if (m1s < 0) {
                            coef1c = complex(coef1(n1, lm1 - 2*m1s), \
                                     -coef1(n1, lm1));
                          }
                          if (m2s > 0) {
                            coef2c = math::pow(-1.0, m2s)* \
                                     complex(coef2(n2, lm2), \
                                     coef2(n2, lm2 - 2*m2s));
                          } else if (m2s == 0) {
                            coef2c = complex(coef2(n2, lm2), 0.0)* \
                                     std::sqrt(2.0);
                          } else if (m2s < 0) {
                            coef2c = complex(coef2(n2, lm2 - 2*m2s), \
                                     -coef2(n2, lm2));
                          }
                          if (m3s > 0) {
                            coef3c = math::pow(-1.0, m3s)* \
                                     complex(coef3(n3, lm3), \
                                     coef3(n3, lm3 - 2*m3s));
                          } else if (m3s == 0) {
                            coef3c = complex(coef3(n3, lm3), 0.0)* \
                                     std::sqrt(2.0);
                          } else if (m3s < 0) {
                            coef3c = complex(coef3(n3, lm3 - 2*m3s), \
                                     -coef3(n3, lm3));
                          }
                          coef1c /= std::sqrt(2.0);
                          coef2c /= std::sqrt(2.0);
                          coef3c /= std::sqrt(2.0);
                          // The descriptor components are purely real or
                          // imaginary
                          // depending on the divisibility of l1 + l2 +l3 by 2.
                          if ((l1 + l2 + l3) % 2 == 0) {
                            soap_vector_by_type(nn, l0) += \
                              w3j*mult*(coef1c*coef2c*coef3c).real();
                          } else {
                            soap_vector_by_type(nn, l0) += \
                              w3j*mult*(coef1c*coef2c*coef3c).imag();
                          }
                          count++;
                          } // m3
                          } // m2
                          } // m1
                          l0++;
                        } // l3
                      } // l2
                    } // l1
                    nn++;
                  } // n3
                } // n2
              } // n1
            } // if count triplet
          } // coef3
        } // coef2
      } // coef1

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }

    } // center
  }  // end function

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<
      Mngr>::initialize_percenter_bispectrum_soap_vectors() {
    size_t n_row{math::pow(this->max_radial, 3_z)};
    size_t n_col{0};
    double max_ang{static_cast<double>(this->max_angular)};
    if (this->inversion_symmetry == false) {
      n_col = static_cast<size_t>(1.0 + 2.0*max_ang + \
              1.5*math::pow(max_ang, 2_z) + \
              math::pow(max_ang, 3_z)*0.5);
    } else {
      n_col = static_cast<size_t>(std::floor(((math::pow(max_ang + 1.0, 2_z) + 1)* (2*(max_ang + 1.0) + 3))/8.0));
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

      std::vector<internal::SortedKey<Key_t>> triplet_list{};
      auto & center_type{center.get_atom_type()};
      Key_t triplet_type{center_type, center_type, center_type};
      // TODO(felix) optimize this loop
      for (const auto& el1: coefficients) {
        triplet_type[0] = el1.first[0];
        for (const auto& el2: coefficients) {
          triplet_type[1] = el2.first[0];
          for (const auto& el3: coefficients) {
            triplet_type[2] = el3.first[0];
            triplet_list.emplace_back(is_not_sorted, triplet_type);
          }
        }
      }
      // initialize the power spectrum with the proper dimension
      soap_vector.resize(triplet_list, n_row, n_col);
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<
      Mngr>::initialize_percenter_powerspectrum_soap_vectors() {
    size_t n_row{math::pow(this->max_radial, 2_z)};
    size_t n_col{this->max_angular + 1};

    // clear the data container and resize it
    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    // identify the species in each environment and initialize soap_vectors
    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      internal::Sorted<true> is_sorted{};

      std::vector<internal::SortedKey<Key_t>> pair_list{};
      auto & center_type{center.get_atom_type()};
      Key_t pair_type{center_type, center_type};
      // avoid checking the order in pair_type by ensuring it has already been
      // done
      pair_list.emplace_back(is_sorted, pair_type);
      for (const auto & el1 : coefficients) {
        auto && neigh1_type{el1.first[0]};
        if (center_type <= neigh1_type) {
          pair_type[0] = center_type;
          pair_type[1] = neigh1_type;
        } else {
          pair_type[1] = center_type;
          pair_type[0] = neigh1_type;
        }

        pair_list.emplace_back(is_sorted, pair_type);

        for (const auto & el2 : coefficients) {
          auto && neigh2_type{el2.first[0]};
          if (neigh1_type <= neigh2_type) {
            pair_type[0] = neigh1_type;
            pair_type[1] = neigh2_type;
            pair_list.emplace_back(is_sorted, pair_type);
          }
        }
      }
      // initialize the power spectrum with the proper dimension
      soap_vector.resize(pair_list, n_row, n_col);
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAPInvariant<
      Mngr>::initialize_percenter_radialspectrum_soap_vectors() {
    size_t n_row{this->max_radial};
    size_t n_col{1};

    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    for (auto center : this->structure_manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{this->soap_vectors[center]};
      Key_t element_type{0};

      std::unordered_set<Key_t, internal::Hash<Key_t>> keys{};
      for (const auto & el1 : coefficients) {
        keys.insert({el1.first[0]});
      }
      keys.insert({center.get_atom_type()});
      // initialize the radial spectrum to 0 and the proper size
      soap_vector.resize(keys, n_row, n_col, 0);
    }
  }

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
