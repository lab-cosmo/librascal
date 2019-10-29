/**
 * @file   calculator_spherical_covariants.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Michael Willatt <michael.willatt@epfl.ch>
 * @author Andrea Grisafi <andrea.grisafi@epfl.ch>
 *
 * @date   12 March 2019
 *
 * @brief  Compute covariant features of the spherical expansion
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

#ifndef SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_COVARIANTS_HH_
#define SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_COVARIANTS_HH_

#include "math/math_utils.hh"
#include "rascal_utility.hh"
#include "representations/calculator_base.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "representations/calculator_spherical_invariants.hh"
#include "structure_managers/property.hh"
#include "structure_managers/property_block_sparse.hh"
#include "structure_managers/structure_manager.hh"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>

namespace rascal {

  namespace internal {
    enum class SphericalCovariantsType { LambdaSpectrum, End_ };

    template <SphericalCovariantsType SpectrumType>
    struct SphericalCovariantsPrecomputation {};

    template <>
    struct SphericalCovariantsPrecomputation<
        SphericalCovariantsType::LambdaSpectrum>
        : SphericalInvariantsPrecomputationBase {
      using Parent = SphericalInvariantsPrecomputationBase;
      using Hypers_t = typename SphericalInvariantsPrecomputationBase::Hypers_t;

      explicit SphericalCovariantsPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular").get<size_t>();
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
        this->lambda = hypers.at("lam").get<size_t>();

        // get the number of non zero elements in the w3j
        // Use the well-known selection rules for the Wigner 3j symbols
        int n_elements{0};
        size_t l3{this->lambda};
        for (size_t l1{0}; l1 < this->max_angular + 1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular + 1; ++l2) {
            // Triangle rule
            if (l1 < static_cast<size_t>(std::abs<int>(l2 - l3)) ||
                l1 > l2 + l3) {
              continue;
            }
            if (this->inversion_symmetry == true) {
              if ((l1 + l2 + l3) % 2 == 1) {
                continue;
              }
            }
            for (size_t m1{0}; m1 < 2 * l1 + 1; m1++) {
              int m1s{static_cast<int>(m1 - l1)};
              for (size_t m2{0}; m2 < 2 * l2 + 1; m2++) {
                int m2s{static_cast<int>(m2 - l2)};
                for (size_t m3{0}; m3 < 2 * l3 + 1; m3++) {
                  int m3s{static_cast<int>(m3 - l3)};
                  // Selection rule for m
                  if ((m1s + m2s + m3s != 0) && (m1s + m2s - m3s != 0)) {
                    continue;
                  }
                  ++n_elements;
                }
              }
            }
          }
        }

        this->wigner_3js.resize(n_elements);
        n_elements = 0;
        wig_table_init(2 * (this->max_angular + 1), 3);
        wig_temp_init(2 * (this->max_angular + 1));
        for (size_t l1{0}; l1 < this->max_angular + 1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular + 1; ++l2) {
            if ((l1 < static_cast<size_t>(std::abs<int>(l2 - l3))) ||
                (l1 > l2 + l3)) {
              continue;
            }
            if (this->inversion_symmetry == true) {
              if ((l1 + l2 + l3) % 2 == 1) {
                continue;
              }
            }
            for (size_t m1{0}; m1 < 2 * l1 + 1; m1++) {
              int m1s{static_cast<int>(m1 - l1)};
              for (size_t m2{0}; m2 < 2 * l2 + 1; m2++) {
                int m2s{static_cast<int>(m2 - l2)};
                for (size_t m3{0}; m3 < 2 * l3 + 1; m3++) {
                  int m3s{static_cast<int>(m3 - l3)};
                  if ((m1s + m2s + m3s != 0) && (m1s + m2s - m3s != 0)) {
                    continue;
                  }
                  double w3j1{wig3jj(2 * l1, 2 * l2, 2 * l3, 2 * m1s, 2 * m2s,
                                     2 * m3s)};
                  double w3j2{wig3jj(2 * l1, 2 * l2, 2 * l3, 2 * m1s, 2 * m2s,
                                     -2 * m3s)};
                  if (m3s > 0) {
                    this->wigner_3js(n_elements) =
                        (w3j2 + math::pow(-1, m3s) * w3j1) * math::INV_SQRT_TWO;
                  } else if (m3s < 0) {
                    this->wigner_3js(n_elements) =
                        ((w3j1 - math::pow(-1, m3s) * w3j2)) *
                        math::INV_SQRT_TWO;
                  } else if (m3s == 0) {
                    this->wigner_3js(n_elements) = w3j1;
                  }
                  //*/
                  // change to the following for agreement with SOAPFAST
                  //(different definition of the real spherical harmonics)
                  /*
                  if (m3s > 0) { this->wigner_3js.push_back((w3j1 + pow(-1,
                  m3s)*w3j2)/sqrt(2.0)); } else if (m3s == 0) {
                  this->wigner_3js.push_back(w3j1); } else if (m3s < 0) {
                  this->wigner_3js.push_back(((w3j2 - pow(-1,
                  m3s)*w3j1))/sqrt(2.0));
                  }
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
      Eigen::ArrayXd wigner_3js{};
      size_t lambda{2};
    };

  }  // namespace internal

  template <internal::SphericalCovariantsType Type, class Hypers>
  auto make_spherical_covariants_precompute(const Hypers & hypers) {
    return std::static_pointer_cast<
        internal::SphericalInvariantsPrecomputationBase>(
        std::make_shared<internal::SphericalCovariantsPrecomputation<Type>>(
            hypers));
  }

  template <internal::SphericalCovariantsType Type>
  auto downcast_spherical_covariants_precompute(
      const std::shared_ptr<internal::SphericalInvariantsPrecomputationBase> &
          spherical_covariants_precompute) {
    return std::static_pointer_cast<
        internal::SphericalCovariantsPrecomputation<Type>>(
        spherical_covariants_precompute);
  }

  class CalculatorSphericalCovariants : public CalculatorBase {
   public:
    using Hypers_t = typename CalculatorBase::Hypers_t;
    using Key_t = typename CalculatorBase::Key_t;

    template <class StructureManager>
    using Property_t =
        BlockSparseProperty<double, 1, 0, StructureManager, Key_t>;
    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, 0, StructureManager, Key_t>;

    template <class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;
    template <class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    explicit CalculatorSphericalCovariants(const Hypers_t & hyper)
        : CalculatorBase{}, rep_expansion{hyper} {
      this->set_default_prefix("spherical_covariants_");
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    CalculatorSphericalCovariants(const CalculatorSphericalCovariants & other) =
        delete;

    //! Move constructor
    CalculatorSphericalCovariants(CalculatorSphericalCovariants && other) =
        default;

    //! Destructor
    virtual ~CalculatorSphericalCovariants() = default;

    //! Copy assignment operator
    CalculatorSphericalCovariants &
    operator=(const CalculatorSphericalCovariants & other) = delete;

    //! Move assignment operator
    CalculatorSphericalCovariants &
    operator=(CalculatorSphericalCovariants && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::SphericalCovariantsType;
      this->max_radial = hypers.at("max_radial").get<size_t>();
      this->max_angular = hypers.at("max_angular").get<size_t>();
      this->spherical_covariants_type_str =
          hypers.at("soap_type").get<std::string>();
      this->lambda = hypers.at("lam").get<size_t>();
      this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
      this->normalize = hypers.at("normalize").get<bool>();

      if (this->spherical_covariants_type_str == "LambdaSpectrum") {
        this->spherical_covariants_type =
            SphericalCovariantsType::LambdaSpectrum;
        this->precompute_spherical_covariants[enumValue(
            SphericalCovariantsType::LambdaSpectrum)] =
            make_spherical_covariants_precompute<
                SphericalCovariantsType::LambdaSpectrum>(hypers);

      } else {
        throw std::logic_error("Requested Spherical Covariants type \'" +
                               this->spherical_covariants_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'LambdaSpectrum\'.");
      }

      this->set_name(hypers);
    }

    /**
     * Compute representation for a given structure manager.
     *
     * @tparam StructureManager a (single or collection)
     * of structure manager(s) (in an iterator) held in shared_ptr
     *
     * TODO(felix) add mechanism to check if the StructureManager is
     * compatible with the representation
     */
    template <class StructureManager>
    void compute(StructureManager & managers);

    /**
     * loop over a collection of manangers if it is an iterator.
     * Or just call compute_impl
     */
    template <
        internal::SphericalCovariantsType Type, class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<Type>(manager);
      }
    }

    //! single manager case
    template <internal::SphericalCovariantsType Type, class StructureManager,
              std::enable_if_t<
                  not(internal::is_proper_iterator<StructureManager>::value),
                  int> = 0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<Type>(manager);
    }

    //! compute representation lambda spectrum
    template <
        internal::SphericalCovariantsType Type,
        std::enable_if_t<
            Type == internal::SphericalCovariantsType::LambdaSpectrum, int> = 0,
        class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

    template <class Invariants, class ExpansionCoeff, class StructureManager>
    void initialize_per_center_lambda_soap_vectors(
        Invariants & soap_vector, ExpansionCoeff & expansions_coefficients,
        std::shared_ptr<StructureManager> manager);

   protected:
    size_t max_radial{};
    size_t max_angular{};

    CalculatorSphericalExpansion rep_expansion;
    internal::SphericalCovariantsType spherical_covariants_type{};
    std::array<std::shared_ptr<internal::SphericalInvariantsPrecomputationBase>,
               internal::enumSize<internal::SphericalCovariantsType>()>
        precompute_spherical_covariants{};
    std::string spherical_covariants_type_str{};

    bool inversion_symmetry{false};
    size_t lambda{0};
    bool normalize{true};
  };

  template <class StructureManager>
  void CalculatorSphericalCovariants::compute(StructureManager & managers) {
    using internal::SphericalCovariantsType;
    switch (this->spherical_covariants_type) {
    case SphericalCovariantsType::LambdaSpectrum:
      this->compute_loop<SphericalCovariantsType::LambdaSpectrum>(managers);
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  template <class Invariants, class ExpansionCoeff, class StructureManager>
  void CalculatorSphericalCovariants::initialize_per_center_lambda_soap_vectors(
      Invariants & soap_vectors, ExpansionCoeff & expansions_coefficients,
      std::shared_ptr<StructureManager> manager) {
    using math::pow;

    size_t n_row{pow(this->max_radial, 2_n)};
    // number of combinations of l1 and l2 satisfying the triangle constraint
    size_t n_col{0};
    if (this->inversion_symmetry == false) {
      n_col = static_cast<size_t>(
          (2 + this->lambda - 3 * pow(this->lambda, 2_n) +
           2 * this->max_angular + 4 * this->lambda * this->max_angular) /
          2 * (2 * this->lambda + 1));
    } else {
      n_col =
          static_cast<size_t>(
              std::ceil(pow(this->max_angular + 1, 2_n) / 2.0) -
              pow(1.0 + std::floor((this->lambda - 1) / 2.0), 2) -
              std::floor(pow(this->max_angular + 1 - this->lambda, 2_n) / 2.0) *
                  (this->lambda % 2) -
              (std::ceil(pow(this->max_angular + 1 - this->lambda, 2_n) / 2.0) -
               (this->max_angular - this->lambda + 1)) *
                  (1.0 - this->lambda % 2)) *
          (2 * this->lambda + 1);
      if (this->lambda % 2 == 1) {
        n_col = static_cast<size_t>(
            0.5 *
                (2 + this->lambda - 3 * pow(this->lambda, 2_n) +
                 2 * this->max_angular + 4 * this->lambda * this->max_angular) *
                (2 * this->lambda + 1) -
            static_cast<int>(n_col));
      }
    }

    // clear the data container and resize it
    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);
    soap_vectors.resize();

    // identify the species in each environment and initialize soap_vectors
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
      internal::Sorted<false> is_not_sorted{};

      std::vector<internal::SortedKey<Key_t>> pair_list{};
      auto center_type{center.get_atom_type()};
      Key_t pair_type{center_type, center_type};
      // TODO(felix) optimize this loop
      // the species are not sorted by construction so there are sorted
      // explicitly and many redundant combinations are present in pair_list
      for (const auto & el1 : coefficients) {
        pair_type[0] = el1.first[0];
        for (const auto & el2 : coefficients) {
          pair_type[1] = el2.first[0];
          pair_list.emplace_back(is_not_sorted, pair_type);
        }
      }
      // initialize the power spectrum with the proper dimension
      soap_vector.resize(pair_list, n_row, n_col, 0);
    }
  }

  template <internal::SphericalCovariantsType Type,
            std::enable_if_t<
                Type == internal::SphericalCovariantsType::LambdaSpectrum, int>,
            class StructureManager>
  void CalculatorSphericalCovariants::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using PropExp_t =
        typename CalculatorSphericalExpansion::Property_t<StructureManager>;
    // using PropGradExp_t = typename
    // CalculatorSphericalExpansion::Property_t<StructureManager>;
    using Prop_t = Property_t<StructureManager>;
    // using PropGrad_t = PropertyGradient_t<StructureManager>;
    using internal::enumValue;
    using internal::SphericalCovariantsType;
    using math::pow;
    using complex = std::complex<double>;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_spherical_covariants_precompute<
        SphericalCovariantsType::LambdaSpectrum>(
        this->precompute_spherical_covariants[enumValue(
            SphericalCovariantsType::LambdaSpectrum)])};
    auto & wigner_3js{precomputation->wigner_3js};

    // Compute the spherical expansions of the current structure
    rep_expansion.compute(manager);
    auto && expansions_coefficients{
        manager->template get_property_ref<PropExp_t>(
            rep_expansion.get_name())};

    auto && soap_vectors{
        manager->template get_property_ref<Prop_t>(this->get_name())};

    // if the representation has already been computed for the current
    // structure then do nothing
    if (soap_vectors.is_updated()) {
      return;
    }

    this->initialize_per_center_lambda_soap_vectors(
        soap_vectors, expansions_coefficients, manager);

    Key_t p_type{0, 0};
    internal::SortedKey<Key_t> pair_type{p_type};

    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};

      for (const auto & el1 : coefficients) {
        pair_type[0] = el1.first[0];
        auto & coef1{el1.second};
        for (const auto & el2 : coefficients) {
          pair_type[1] = el2.first[0];
          auto & coef2{el2.second};

          size_t & l3{this->lambda};
          if (soap_vector.count(pair_type) == 1) {
            auto && soap_vector_by_type{soap_vector[pair_type]};

            size_t nn{0};
            for (size_t n1{0}; n1 < this->max_radial; n1++) {
              for (size_t n2{0}; n2 < this->max_radial; n2++) {
                size_t l0{0};
                int wigner_count{0};
                for (size_t l1{0}; l1 < this->max_angular + 1; l1++) {
                  for (size_t l2{0}; l2 < this->max_angular + 1; l2++) {
                    if ((l1 < static_cast<size_t>(std::abs<int>(l2 - l3))) ||
                        (l1 > l2 + l3)) {
                      continue;
                    }
                    if (this->inversion_symmetry == true) {
                      if ((l1 + l2 + l3) % 2 == 1) {
                        continue;
                      }
                    }
                    for (size_t m3{0}; m3 < 2 * l3 + 1; m3++) {
                      int m3s{static_cast<int>(m3 - l3)};
                      for (size_t m1{0}; m1 < 2 * l1 + 1; m1++) {
                        int m1s{static_cast<int>(m1 - l1)};
                        size_t lm1{pow(l1, 2_n) + m1};
                        for (size_t m2{0}; m2 < 2 * l2 + 1; m2++) {
                          int m2s{static_cast<int>(m2 - l2)};
                          if ((m1s + m2s + m3s != 0) &&
                              (m1s + m2s - m3s != 0)) {
                            continue;
                          }
                          size_t lm2{pow(l2, 2_n) + m2};
                          complex coef1c, coef2c;
                          double w3j = wigner_3js[wigner_count];
                          // usual formulae for converting from real to complex
                          // These are the inverses of the formulae we used to
                          // define the real spherical harmonics in
                          // src/math/spherical_harmonics.hh
                          if (m1s > 0) {
                            coef1c = pow(-1.0, m1s) *
                                     complex(coef1(n1, lm1),
                                             coef1(n1, lm1 - 2 * m1s));
                          } else if (m1s == 0) {
                            coef1c =
                                complex(coef1(n1, lm1), 0.0) * math::SQRT_TWO;
                          } else if (m1s < 0) {
                            coef1c = complex(coef1(n1, lm1 - 2 * m1s),
                                             -coef1(n1, lm1));
                          }
                          if (m2s > 0) {
                            coef2c = pow(-1.0, m2s) *
                                     complex(coef2(n2, lm2),
                                             coef2(n2, lm2 - 2 * m2s));
                          } else if (m2s == 0) {
                            coef2c =
                                complex(coef2(n2, lm2), 0.0) * math::SQRT_TWO;
                          } else if (m2s < 0) {
                            coef2c = complex{coef2(n2, lm2 - 2 * m2s),
                                             -coef2(n2, lm2)};
                          }
                          coef1c *= math::INV_SQRT_TWO;
                          coef2c *= math::INV_SQRT_TWO;

                          // combine the coefficients with Wigner 3j symbols
                          // TODO(andrea, michael) avoid complex arithmetic
                          complex i_unit{0.0, 1.0};
                          if ((l1 + l2 + l3) % 2 == 0) {
                            if (m3s < 0) {
                              soap_vector_by_type(nn, l0) +=
                                  w3j * (i_unit * coef1c * coef2c).real();
                            } else {
                              soap_vector_by_type(nn, l0) +=
                                  w3j * (coef1c * coef2c).real();
                            }
                          } else if (this->inversion_symmetry == false) {
                            if (m3s < 0) {
                              soap_vector_by_type(nn, l0) +=
                                  w3j * (i_unit * coef1c * coef2c).imag();
                            } else {
                              soap_vector_by_type(nn, l0) +=
                                  w3j * (coef1c * coef2c).imag();
                            }
                          }
                          wigner_count++;
                        }  // m2
                      }    // m1
                      l0++;
                    }  // m3
                  }    // l2
                }      // l1
                nn++;
              }  // n2
            }    // n1
          }      // if pair_type count
        }        // coef1
      }          // coef1

      // the SQRT_TWO factor comes from the fact that
      // the upper diagonal of the species is not considered
      soap_vector.multiply_off_diagonal_elements_by(math::SQRT_TWO);

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }
    }  // center
  }    // compute_lambdaspectrum

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_COVARIANTS_HH_
