/**
 * file   calculator_spherical_invariants.hh
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

#ifndef SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
#define SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_

#include "representations/calculator_base.hh"
#include "representations/calculator_spherical_expansion.hh"
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
    enum class SOAPType { RadialSpectrum, PowerSpectrum, End_ };

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

      using Hypers_t = CalculatorBase::Hypers_t;
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
        this->l_factors.resize(math::pow(this->max_angular + 1, 2));

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

  class CalculatorSphericalInvariants : public CalculatorBase {
   public:

    using Hypers_t = typename CalculatorBase::Hypers_t;
    using Key_t = typename CalculatorBase::Key_t;

    template<class StructureManager>
    using Property_t =
        BlockSparseProperty<double, 1, 0, StructureManager, Key_t>;

    template<class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;
    template<class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    using internal::enumValue;
    using internal::SOAPType;
    using internal::SOAPPrecomputationBase;

    CalculatorSphericalInvariants(const Hypers_t & hyper)
      : rep_expansion{hyper} {
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    CalculatorSphericalInvariants(const CalculatorSphericalInvariants & other) =
        delete;

    //! Move constructor
    CalculatorSphericalInvariants(CalculatorSphericalInvariants && other) =
        default;

    //! Destructor
    virtual ~CalculatorSphericalInvariants() = default;

    //! Copy assignment operator
    CalculatorSphericalInvariants &
    operator=(const CalculatorSphericalInvariants & other) = delete;

    //! Move assignment operator
    CalculatorSphericalInvariants &
    operator=(CalculatorSphericalInvariants && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
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
      } else {
        throw std::logic_error("Requested SOAP type \'" + this->soap_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'PowerSpectrum or RadialSpectrum\'.");
      }
    }


    /**
     * Compute representation for a given structure manager.
     *
     * @tparam StructureManager a (single or collection in an iterator)
     * of structure manager(s) held in shared_ptr
     *
     * TODO(felix) add mechanism to check if the StructureManager is
     * compatible with the representation
     */
    template <class StructureManager>
    void compute(StructureManager& managers);

    /**
     * loop over a collection of manangers if it is an iterator,
     * i.e. std::iterator_traits<T>::value_type is defined.
     * Or just call compute_impl
     */
    template <class StructureManager,
              SOAPType BodyOrder,
              std::enable_if_t<internal::is_iterator<StructureManager>::value, int> = 0>
    inline void compute_loop(StructureManager& managers) {
      for (auto& manager : managers) {
        ComputeHelper<StructureManager, BodyOrder>::run(*this, manager);
      }
    }

    //! single manager case
    template <class StructureManager,
              SOAPType BodyOrder,
              std::enable_if_t<(not internal::is_iterator<StructureManager>::value), int> = 0>
    inline void compute_loop(StructureManager& managers) {
      ComputeHelper<StructureManager, BodyOrder>::run(*this, manager);
    }

    /**
     *
     */
    template <class StructureManager, SOAPType BodyOrder>
    struct ComputeHelper {};

    template <class StructureManager>
    struct ComputeHelper<StructureManager, SOAPType::RadialSpectrum> {
      static void run(const CalculatorSphericalInvariants& calculator,
                      std::shared_ptr<StructureManager>& manager) {
        calculator.compute_radialspectrum(manager);
      }
    };

    template <class StructureManager>
    struct ComputeHelper<StructureManager, SOAPType::PowerSpectrum> {
      static void run(const CalculatorSphericalInvariants& calculator,
                      std::shared_ptr<StructureManager>& manager) {
        calculator.compute_powerspectrum(manager);
      }
    };

    //! compute representation \nu == 1
    template <class StructureManager>
    void compute_radialspectrum(std::shared_ptr<StructureManager>& manager);

    //! compute representation \nu == 2
    template <class StructureManager>
    void compute_powerspectrum(std::shared_ptr<StructureManager>& manager);


    //! initialize the soap vectors with only the keys needed for each center
    template <class StructureManager, class Invariants, class ExpansionCoeff>
    void initialize_percenter_powerspectrum_soap_vectors(Invariants& soap_vector, ExpansionCoeff& expansions_coefficients, std::shared_ptr<StructureManager>& manager);

    template <class StructureManager, class Invariants, class ExpansionCoeff>
    void initialize_percenter_radialspectrum_soap_vectors(Invariants& soap_vector, ExpansionCoeff& expansions_coefficients, std::shared_ptr<StructureManager>& manager);

    inline std::string get_name() {
      return calculator_name;
    }

   protected:
    size_t max_radial{};
    size_t max_angular{};
    bool normalize{};

    CalculatorSphericalExpansion rep_expansion;
    
    SOAPType soap_type{};
    //! collection of precomputation for the different body order
    std::array<std::shared_ptr<SOAPPrecomputationBase>,
               enumSize<SOAPType>()>
        precompute_soap{};
    std::string soap_type_str{};

    constexpr char calculator_name[] = "spherical_invariant";

  };

  template <class StructureManager>
  void CalculatorSphericalInvariants::compute(StructureManager& managers) {
    using internal::SOAPType;
    switch (this->soap_type) {
    case SOAPType::RadialSpectrum:
      this->compute_loop<SOAPType::RadialSpectrum>(managers);
      break;
    case SOAPType::PowerSpectrum:
      this->compute_loop<SOAPType::PowerSpectrum>(managers);
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  template <class StructureManager>
  void CalculatorSphericalInvariants::compute_powerspectrum(std::shared_ptr<StructureManager>& manager) {
    using internal::enumValue;
    using internal::SOAPType;
    using math::pow;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_soap_precompute<SOAPType::PowerSpectrum>(
        this->precompute_soap[enumValue(SOAPType::PowerSpectrum)])};
    auto & l_factors{precomputation->l_factors};

    // TODO(felix) use the updated mech of the prop to avoid recomputing
    // if the prop already exists and is uptodate
    // Compute the spherical expansions of the current structure
    rep_expansion.compute(manager);

    auto&& expansions_coefficients{this->get_property<CalculatorSphericalExpansion::Property_t>(manager, rep_expansion.get_name())};

    auto&& soap_vectors{this->get_property<Property_t>(manager, this->get_name())};
    this->initialize_percenter_powerspectrum_soap_vectors(soap_vectors, expansions_coefficients, manager);

    Key_t pair_type{0, 0};
    // use special container to tell that there is not need to sort when
    // using operator[] of soap_vector
    internal::SortedKey<Key_t> spair_type{pair_type};

    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};

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

  template <class StructureManager>
  void CalculatorSphericalInvariants::compute_radialspectrum(std::shared_ptr<StructureManager>& manager) {
    using math::pow;

    rep_expansion.compute(manager);

    auto&& expansions_coefficients{this->get_property<CalculatorSphericalExpansion::Property_t>(manager, rep_expansion.get_name())};

    auto&& soap_vectors{this->get_property<Property_t>(manager, this->get_name())};
    this->initialize_percenter_radialspectrum_soap_vectors(soap_vectors, expansions_coefficients, manager);
    Key_t element_type{0};

    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};

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

  template <class StructureManager, class Invariants, class ExpansionCoeff>
  void CalculatorSphericalInvariants::initialize_percenter_powerspectrum_soap_vectors(Invariants& soap_vector, ExpansionCoeff& expansions_coefficients, std::shared_ptr<StructureManager>& manager) {
    size_t n_row{static_cast<size_t>(pow(this->max_radial, 2))};
    size_t n_col{this->max_angular + 1};

    // clear the data container and resize it
    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);
    soap_vectors.resize();

    // identify the species in each environment and initialize soap_vectors
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
      internal::Sorted<true> is_sorted{};

      std::vector<internal::SortedKey<Key_t>> pair_list{};
      auto & center_type{center.get_atom_type()};
      Key_t pair_type{center_type, center_type};
      // avoid checking the order in pair_type by ensuring it has already been
      // done
      internal::SortedKey<Key_t> spair_type{is_sorted, pair_type};

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

  template <class StructureManager, class Invariants, class ExpansionCoeff>
  void CalculatorSphericalInvariants::initialize_percenter_radialspectrum_soap_vectors(Invariants& soap_vector, ExpansionCoeff& expansions_coefficients, std::shared_ptr<StructureManager>& manager) {
    size_t n_row{this->max_radial};
    size_t n_col{1};

    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);
    soap_vectors.resize();

    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
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

#endif  // SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
