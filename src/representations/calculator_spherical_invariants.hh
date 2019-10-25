/**
 * @file   calculator_spherical_invariants.hh
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

#include "math/math_utils.hh"
#include "rascal_utility.hh"
#include "representations/calculator_base.hh"
#include "representations/calculator_spherical_expansion.hh"
#include "structure_managers/property_block_sparse.hh"
#include "structure_managers/structure_manager.hh"

#include <wigxjpf.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <exception>
#include <unordered_set>
#include <vector>

namespace rascal {

  namespace internal {
    enum class SphericalInvariantsType {
      RadialSpectrum,
      PowerSpectrum,
      BiSpectrum,
      End_
    };

    /**
     * Base class for the specification of the atomic smearing.
     */
    struct SphericalInvariantsPrecomputationBase {
      //! Constructor
      SphericalInvariantsPrecomputationBase() = default;
      //! Destructor
      virtual ~SphericalInvariantsPrecomputationBase() = default;
      //! Copy constructor
      SphericalInvariantsPrecomputationBase(
          const SphericalInvariantsPrecomputationBase & other) = delete;
      //! Move constructor
      SphericalInvariantsPrecomputationBase(
          SphericalInvariantsPrecomputationBase && other) = default;
      //! Copy assignment operator
      SphericalInvariantsPrecomputationBase &
      operator=(const SphericalInvariantsPrecomputationBase & other) = delete;
      //! Move assignment operator
      SphericalInvariantsPrecomputationBase &
      operator=(SphericalInvariantsPrecomputationBase && other) = default;

      using Hypers_t = CalculatorBase::Hypers_t;
    };

    template <SphericalInvariantsType SpectrumType>
    struct SphericalInvariantsPrecomputation {};

    template <>
    struct SphericalInvariantsPrecomputation<
        SphericalInvariantsType::RadialSpectrum>
        : SphericalInvariantsPrecomputationBase {
      using Hypers_t = typename SphericalInvariantsPrecomputationBase::Hypers_t;
      explicit SphericalInvariantsPrecomputation(const Hypers_t &) {}
    };

    template <>
    struct SphericalInvariantsPrecomputation<
        SphericalInvariantsType::PowerSpectrum>
        : SphericalInvariantsPrecomputationBase {
      using Hypers_t = typename SphericalInvariantsPrecomputationBase::Hypers_t;
      explicit SphericalInvariantsPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular");
        this->l_factors.resize(this->max_angular + 1);
        for (size_t l{0}; l < this->max_angular + 1; ++l) {
          this->l_factors(l) = math::pow(std::sqrt(2 * l + 1), -1);
        }
      }

      size_t max_angular{0};
      //! factor of 1 / sqrt(2*l+1) in front of the powerspectrum
      Eigen::VectorXd l_factors{};
    };

    template <>
    struct SphericalInvariantsPrecomputation<
        SphericalInvariantsType::BiSpectrum>
        : SphericalInvariantsPrecomputationBase {
      using Hypers_t = typename SphericalInvariantsPrecomputationBase::Hypers_t;
      explicit SphericalInvariantsPrecomputation(const Hypers_t & hypers) {
        this->max_angular = hypers.at("max_angular").get<size_t>();
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
        // get the number of non zero elements in the w3j
        int n_elements{0};
        for (size_t l1{0}; l1 < this->max_angular + 1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular + 1; ++l2) {
            for (size_t l3{0}; l3 < this->max_angular + 1; ++l3) {
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
                    if (m1s + m2s + m3s != 0) {
                      continue;
                    }
                    ++n_elements;
                  }
                }
              }
            }
          }
        }

        this->wigner_w3js.resize(n_elements);
        n_elements = 0;
        wig_table_init(2 * (this->max_angular + 1), 3);
        wig_temp_init(2 * (this->max_angular + 1));
        for (size_t l1{0}; l1 < this->max_angular + 1; ++l1) {
          for (size_t l2{0}; l2 < this->max_angular + 1; ++l2) {
            for (size_t l3{0}; l3 < this->max_angular + 1; ++l3) {
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
                    if (m1s + m2s + m3s != 0) {
                      continue;
                    }
                    this->wigner_w3js(n_elements) = wig3jj(
                        2 * l1, 2 * l2, 2 * l3, 2 * m1s, 2 * m2s, 2 * m3s);
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
      Eigen::ArrayXd wigner_w3js{};
    };
  }  // namespace internal

  template <internal::SphericalInvariantsType Type, class Hypers>
  auto make_spherical_invariants_precompute(const Hypers & hypers) {
    return std::static_pointer_cast<
        internal::SphericalInvariantsPrecomputationBase>(
        std::make_shared<internal::SphericalInvariantsPrecomputation<Type>>(
            hypers));
  }

  template <internal::SphericalInvariantsType Type>
  auto downcast_spherical_invariants_precompute(
      const std::shared_ptr<internal::SphericalInvariantsPrecomputationBase> &
          spherical_invariants_precompute) {
    return std::static_pointer_cast<
        internal::SphericalInvariantsPrecomputation<Type>>(
        spherical_invariants_precompute);
  }

  class CalculatorSphericalInvariants : public CalculatorBase {
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

    explicit CalculatorSphericalInvariants(const Hypers_t & hyper)
        : CalculatorBase{}, rep_expansion{hyper} {
      this->set_default_prefix("spherical_invariants_");
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
      using internal::enumValue;
      using internal::SphericalInvariantsPrecomputationBase;
      using internal::SphericalInvariantsType;

      this->max_radial = hypers.at("max_radial").get<size_t>();
      this->max_angular = hypers.at("max_angular").get<size_t>();
      this->normalize = hypers.at("normalize").get<bool>();
      this->spherical_invariants_type_str =
          hypers.at("soap_type").get<std::string>();

      if (hypers.find("compute_gradients") != hypers.end()) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }

      if (this->spherical_invariants_type_str == "PowerSpectrum") {
        this->spherical_invariants_type =
            SphericalInvariantsType::PowerSpectrum;
        this->precompute_spherical_invariants[enumValue(
            SphericalInvariantsType::PowerSpectrum)] =
            make_spherical_invariants_precompute<
                SphericalInvariantsType::PowerSpectrum>(hypers);
      } else if (this->spherical_invariants_type_str == "RadialSpectrum") {
        this->spherical_invariants_type =
            SphericalInvariantsType::RadialSpectrum;
        this->precompute_spherical_invariants[enumValue(
            SphericalInvariantsType::RadialSpectrum)] =
            make_spherical_invariants_precompute<
                SphericalInvariantsType::RadialSpectrum>(hypers);
        if (this->max_angular > 0) {
          throw std::logic_error("max_angular should be 0 with RadialSpectrum");
        }
      } else if (this->spherical_invariants_type_str == "BiSpectrum") {
        this->spherical_invariants_type =
            internal::SphericalInvariantsType::BiSpectrum;
        this->precompute_spherical_invariants[enumValue(
            SphericalInvariantsType::BiSpectrum)] =
            make_spherical_invariants_precompute<
                SphericalInvariantsType::BiSpectrum>(hypers);
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();

      } else {
        throw std::logic_error("Requested SphericalInvariants type \'" +
                               this->spherical_invariants_type_str +
                               "\' has not been implemented.  Must be one of" +
                               ": \'PowerSpectrum or RadialSpectrum\', " +
                               "\'BiSpectrum\'.");
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
        internal::SphericalInvariantsType BodyOrder, class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<BodyOrder>(manager);
      }
    }

    //! single manager case
    template <
        internal::SphericalInvariantsType BodyOrder, class StructureManager,
        std::enable_if_t<
            not(internal::is_proper_iterator<StructureManager>::value), int> =
            0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<BodyOrder>(manager);
    }

    //! compute representation @f$ \nu == 1 @f$
    template <
        internal::SphericalInvariantsType BodyOrder,
        std::enable_if_t<BodyOrder ==
                             internal::SphericalInvariantsType::RadialSpectrum,
                         int> = 0,
        class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

    //! compute representation @f$ \nu == 2 @f$
    template <internal::SphericalInvariantsType BodyOrder,
              std::enable_if_t<
                  BodyOrder == internal::SphericalInvariantsType::PowerSpectrum,
                  int> = 0,
              class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

    //! compute representation @f$ \nu == 3 @f$
    template <internal::SphericalInvariantsType BodyOrder,
              std::enable_if_t<
                  BodyOrder == internal::SphericalInvariantsType::BiSpectrum,
                  int> = 0,
              class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

    //! initialize the soap vectors with only the keys needed for each center
    template <class StructureManager, class Invariants,
              class InvariantsDerivative, class ExpansionCoeff>
    void initialize_per_center_powerspectrum_soap_vectors(
        Invariants & soap_vector, InvariantsDerivative & soap_vector_gradients,
        ExpansionCoeff & expansions_coefficients,
        std::shared_ptr<StructureManager> manager);

    template <class StructureManager, class Invariants,
              class InvariantsDerivative, class ExpansionCoeff>
    void initialize_per_center_radialspectrum_soap_vectors(
        Invariants & soap_vector, InvariantsDerivative & soap_vector_gradients,
        ExpansionCoeff & expansions_coefficients,
        std::shared_ptr<StructureManager> manager);

    template <class StructureManager, class Invariants, class ExpansionCoeff>
    void initialize_per_center_bispectrum_soap_vectors(
        Invariants & soap_vector, ExpansionCoeff & expansions_coefficients,
        std::shared_ptr<StructureManager> manager);

   protected:
    size_t max_radial{};
    size_t max_angular{};
    bool normalize{};
    bool compute_gradients{};
    bool inversion_symmetry{false};

    CalculatorSphericalExpansion rep_expansion;

    internal::SphericalInvariantsType spherical_invariants_type{};
    //! collection of precomputation for the different body order
    std::array<std::shared_ptr<internal::SphericalInvariantsPrecomputationBase>,
               internal::enumSize<internal::SphericalInvariantsType>()>
        precompute_spherical_invariants{};
    std::string spherical_invariants_type_str{};
  };

  template <class StructureManager>
  void CalculatorSphericalInvariants::compute(StructureManager & managers) {
    using internal::SphericalInvariantsType;
    switch (this->spherical_invariants_type) {
    case SphericalInvariantsType::RadialSpectrum:
      this->compute_loop<SphericalInvariantsType::RadialSpectrum>(managers);
      break;
    case SphericalInvariantsType::PowerSpectrum:
      this->compute_loop<SphericalInvariantsType::PowerSpectrum>(managers);
      break;
    case SphericalInvariantsType::BiSpectrum:
      this->compute_loop<SphericalInvariantsType::BiSpectrum>(managers);
      break;
    default:
      // Will never reach here (it's an enum...)
      break;
    }
  }

  template <
      internal::SphericalInvariantsType BodyOrder,
      std::enable_if_t<
          BodyOrder == internal::SphericalInvariantsType::PowerSpectrum, int>,
      class StructureManager>
  void CalculatorSphericalInvariants::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using PropExp_t =
        typename CalculatorSphericalExpansion::Property_t<StructureManager>;
    using PropGradExp_t =
        typename CalculatorSphericalExpansion::PropertyGradient_t<
            StructureManager>;
    using Prop_t = Property_t<StructureManager>;
    using PropGrad_t = PropertyGradient_t<StructureManager>;
    using internal::enumValue;
    constexpr static int n_spatial_dimensions = StructureManager::dim();
    using internal::SphericalInvariantsType;
    using math::pow;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_spherical_invariants_precompute<
        SphericalInvariantsType::PowerSpectrum>(
        this->precompute_spherical_invariants[enumValue(
            SphericalInvariantsType::PowerSpectrum)])};
    auto & l_factors{precomputation->l_factors};

    // Compute the spherical expansions of the current structure
    rep_expansion.compute(manager);

    auto && expansions_coefficients{
        manager->template get_property_ref<PropExp_t>(
            rep_expansion.get_name())};

    // No error if gradients not computed; just an empty array in that case
    auto && expansions_coefficients_gradient{
        manager->template get_property_ref<PropGradExp_t>(
            rep_expansion.get_gradient_name())};

    auto && soap_vectors{
        manager->template get_property_ref<Prop_t>(this->get_name())};

    auto && soap_vector_gradients{
        manager->template get_property_ref<PropGrad_t>(
            this->get_gradient_name())};

    // if the representation has already been computed for the current
    // structure then do nothing
    if (soap_vectors.is_updated()) {
      return;
    }

    this->initialize_per_center_powerspectrum_soap_vectors(
        soap_vectors, soap_vector_gradients, expansions_coefficients, manager);

    Key_t pair_type{0, 0};
    // use special container to tell that there is not need to sort when
    // using operator[] of soap_vector
    internal::SortedKey<Key_t> spair_type{pair_type};

    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
      // Compute the Powerspectrum coefficients
      for (const auto & el1 : coefficients) {
        spair_type[0] = el1.first[0];
        auto & coef1{el1.second};

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
          // multiply with the precomputed factors
          soap_vector_by_pair *= l_factors.asDiagonal();
        }  // for el1 : coefficients
      }    // for el2 : coefficients

      // the SQRT_TWO factor comes from the fact that
      // the upper diagonal of the species is not considered
      soap_vector.multiply_off_diagonal_elements_by(math::SQRT_TWO);

      // normalize the soap vector
      double soap_vector_norm{1.0};
      if (this->normalize) {
        soap_vector_norm = soap_vector.norm();
        soap_vector.normalize();
      }
      if (this->compute_gradients) {
        auto ii_pair = center.get_atom_ii();
        auto & grad_center_coefficients{
            expansions_coefficients_gradient[ii_pair]};
        auto & soap_center_gradient{soap_vector_gradients[ii_pair]};
        for (const auto & grad_species_1 : grad_center_coefficients) {
          spair_type[0] = grad_species_1.first[0];
          const auto & expansion_coefficients_1{
              coefficients[grad_species_1.first]};
          const auto & grad_center_coefficients_1{grad_species_1.second};
          for (const auto & grad_species_2 : grad_center_coefficients) {
            spair_type[1] = grad_species_2.first[0];
            // Half-iteration over species, but not over radial basis index 'n'
            if (spair_type[0] > spair_type[1]) {
              continue;
            }
            const auto & expansion_coefficients_2{
                coefficients[grad_species_2.first]};
            const auto & grad_center_coefficients_2{grad_species_2.second};
            auto && soap_center_gradient_by_species_pair{
                soap_center_gradient[spair_type]};

            // Sum the gradients wrt the central atom position
            size_t n1n2{0};
            size_t l_block_idx{0};
            // TODO(max) looks like this nested loop could be nicely
            // encapsulated in a separate function
            for (size_t cartesian_idx{0}; cartesian_idx < 3; ++cartesian_idx) {
              size_t cartesian_offset_n{cartesian_idx * this->max_radial};
              size_t cartesian_offset_n1n2{
                  cartesian_idx *
                  static_cast<size_t>(math::pow(this->max_radial, 2))};
              n1n2 = 0;
              for (size_t n1{0}; n1 < this->max_radial; ++n1) {
                for (size_t n2{0}; n2 < this->max_radial; ++n2) {
                  // NOTE(max) this is included in the l=0 case, no?
                  // soap_center_gradient_by_species_pair(n1n2, 0) =
                  // coef1(n1, 0) * coef2(n2, 0);
                  // pos = 1;
                  l_block_idx = 0;
                  for (size_t l{0}; l < this->max_angular + 1; ++l) {
                    size_t l_block_size{2 * l + 1};
                    // do the reduction over m (with vectorization)
                    // Leibniz rule for the expansion coefficients
                    // clang-format off
                    soap_center_gradient_by_species_pair(
                            n1n2 + cartesian_offset_n1n2, l) =
                      ((expansion_coefficients_1.block(
                                n1, l_block_idx,
                                1,  l_block_size).array() *
                        grad_center_coefficients_2.block(
                                n2 + cartesian_offset_n, l_block_idx,
                                1, l_block_size).array()).sum()) +
                      ((grad_center_coefficients_1.block(
                                n1 + cartesian_offset_n, l_block_idx,
                                1, l_block_size).array() *
                        expansion_coefficients_2.block(
                                n2, l_block_idx,
                                1,  l_block_size).array()).sum());
                    // clang-format on
                    l_block_idx += l_block_size;
                  }
                  ++n1n2;
                }  // for n2
              }    // for n1
            }      // for cartesian_idx

            // The gradients also need the 1/sqrt(2l + 1) factors
            soap_center_gradient_by_species_pair *= l_factors.asDiagonal();
            // Scale off-(species-diagonal) elements so that the dot product
            // comes out right (carried over from soap vectors to gradients)
            if (spair_type[0] != spair_type[1]) {
              soap_center_gradient_by_species_pair *= math::SQRT_TWO;
            }

            // Sum the gradients wrt the neighbour atom position
            for (auto neigh : center) {
              auto & grad_neigh_coefficients{
                  expansions_coefficients_gradient[neigh]};
              auto & soap_neigh_gradient{soap_vector_gradients[neigh]};
              auto && soap_neigh_gradient_by_species_pair{
                  soap_neigh_gradient[spair_type]};
              soap_neigh_gradient_by_species_pair.setZero();

              auto neigh_type = neigh.get_atom_type();
              if ((neigh_type != spair_type[0]) and
                  (neigh_type != spair_type[1])) {
                // Save the cost of iteration
                // TODO(max) eliminate the zeroing once we figure out how to
                // avoid storing these empty species blocks
                continue;
              }

              // TODO(max) is there a symmetry here we can exploit?
              if (neigh_type == spair_type[0]) {
                const auto & grad_neigh_coefficients_1{
                    grad_neigh_coefficients[grad_species_1.first]};
                for (size_t cartesian_idx{0}; cartesian_idx < 3;
                     ++cartesian_idx) {
                  size_t cartesian_offset_n{cartesian_idx * this->max_radial};
                  size_t cartesian_offset_n1n2{
                      cartesian_idx * math::pow(this->max_radial, 2_n)};
                  n1n2 = 0;
                  for (size_t n1{0}; n1 < this->max_radial; ++n1) {
                    for (size_t n2{0}; n2 < this->max_radial; ++n2) {
                      l_block_idx = 0;
                      for (size_t l{0}; l < this->max_angular + 1; ++l) {
                        size_t l_block_size{2 * l + 1};
                        // clang-format off
                        soap_neigh_gradient_by_species_pair(
                                n1n2 + cartesian_offset_n1n2, l) +=
                          (grad_neigh_coefficients_1.block(
                                n1 + cartesian_offset_n, l_block_idx,
                                1,                       l_block_size).array() *
                           expansion_coefficients_2.block(
                                n2, l_block_idx,
                                1,  l_block_size).array()).sum();
                        // clang-format on
                        l_block_idx += l_block_size;
                      }
                      ++n1n2;
                    }  // for n2
                  }    // for n1
                }      // for cartesian_idx
              }        // if (neigh_type == spair_type[0])

              // Same as above, only gradient wrt neighbour of type 2
              // Necessary because we're only doing a half-iteration over
              // species pairs
              // TODO(max) consider changing to full iteration and just one of
              // these if-nested-horrible-for-loop blocks
              if (neigh_type == spair_type[1]) {
                const auto & grad_neigh_coefficients_2{
                    grad_neigh_coefficients[grad_species_2.first]};
                for (size_t cartesian_idx{0}; cartesian_idx < 3;
                     ++cartesian_idx) {
                  size_t cartesian_offset_n{cartesian_idx * this->max_radial};
                  size_t cartesian_offset_n1n2{
                      cartesian_idx * math::pow(this->max_radial, 2_n)};
                  n1n2 = 0;
                  for (size_t n1{0}; n1 < this->max_radial; ++n1) {
                    for (size_t n2{0}; n2 < this->max_radial; ++n2) {
                      l_block_idx = 0;
                      for (size_t l{0}; l < this->max_angular + 1; ++l) {
                        size_t l_block_size{2 * l + 1};
                        // clang-format off
                        soap_neigh_gradient_by_species_pair(
                                n1n2 + cartesian_offset_n1n2, l) +=
                          (expansion_coefficients_1.block(
                                n1, l_block_idx,
                                1,  l_block_size).array() *
                           grad_neigh_coefficients_2.block(
                                n2 + cartesian_offset_n, l_block_idx,
                                1, l_block_size).array()).sum();
                        // clang-format on
                        l_block_idx += l_block_size;
                      }
                      ++n1n2;
                    }  // for n2
                  }    // for n1
                }      // for cartesian_idx
              }        // if (neigh_type == spair_type[1])

              // Same factors as for the gradient wrt center
              soap_neigh_gradient_by_species_pair *= l_factors.asDiagonal();
              if (spair_type[0] != spair_type[1]) {
                soap_neigh_gradient_by_species_pair *= math::SQRT_TWO;
              }
            }  // for neigh : center
          }    // for grad_species_2 : grad_coefficients
        }      // for grad_species_1 : grad_coefficients

        // NOTE(max) the multiplications below have already been done within the
        // species pair loop above -- not sure which way is more efficient

        // soap_center_gradient.multiply_off_diagonal_elements_by(
        //         math::SQRT_TWO);
        // for (auto neigh : center) {
        //   auto & soap_neigh_gradient{this->soap_vector_gradients[neigh]};
        //   soap_neigh_gradient.multiply_off_diagonal_elements_by(
        //           math::SQRT_TWO);
        // }

        // Update the gradients to include SOAP vector normalization
        // Note that this expects the soap vectors to be normalized already, and
        // the norm stored separately
        if (this->normalize) {
          // Note that the normalization _must_ be done in a separate loop to
          // the loop over species pairs above, because it includes a sum over
          // all species pairs of the already-computed gradient vectors
          //
          // First, perform a dot product of each gradient _component_ with the
          // corresponding normalized soap vector
          Eigen::Vector3d soap_vector_dot_center_gradient{};
          soap_vector_dot_center_gradient.setZero();
          // (n1, n2) * l
          size_t grad_component_size{
              static_cast<size_t>(math::pow(this->max_radial, 2)) *
              (this->max_angular + 1)};
          for (auto soap_grad_spair : soap_center_gradient) {
            auto & soap_gradient_by_species_pair{soap_grad_spair.second};
            const auto & soap_vector_by_species_pair{
                soap_vector[soap_grad_spair.first]};
            Eigen::Map<Eigen::Matrix<double, n_spatial_dimensions,
                                     Eigen::Dynamic, Eigen::RowMajor>>
                soap_gradient_dim_N(soap_gradient_by_species_pair.data(),
                                    n_spatial_dimensions, grad_component_size);
            const Eigen::Map<const Eigen::VectorXd> soap_vector_N(
                soap_vector_by_species_pair.data(), grad_component_size);
            soap_vector_dot_center_gradient +=
                (soap_gradient_dim_N * soap_vector_N);
          }
          // Now update each species-pair-block using the dot-product just
          // computed
          for (auto soap_grad_spair : soap_center_gradient) {
            auto & soap_gradient_by_species_pair{soap_grad_spair.second};
            const auto & soap_vector_by_species_pair{
                soap_vector[soap_grad_spair.first]};
            Eigen::Map<Eigen::Matrix<double, n_spatial_dimensions,
                                     Eigen::Dynamic, Eigen::RowMajor>>
                soap_gradient_dim_N(soap_gradient_by_species_pair.data(),
                                    n_spatial_dimensions, grad_component_size);
            const Eigen::Map<const Eigen::VectorXd> soap_vector_N(
                soap_vector_by_species_pair.data(), grad_component_size);
            soap_gradient_dim_N =
                ((soap_gradient_dim_N -
                  soap_vector_dot_center_gradient * soap_vector_N.transpose()) /
                 soap_vector_norm);
          }

          // Aaand do the same thing for the gradients wrt neighbouring atoms
          for (auto neigh : center) {
            auto & soap_neigh_gradient{soap_vector_gradients[neigh]};
            Eigen::Vector3d soap_vector_dot_neigh_gradient{};
            soap_vector_dot_neigh_gradient.setZero();
            // First, dot product between soap vector and _neighbour_ gradient
            for (auto soap_grad_spair : soap_neigh_gradient) {
              auto & soap_gradient_by_species_pair{soap_grad_spair.second};
              const auto & soap_vector_by_species_pair{
                  soap_vector[soap_grad_spair.first]};
              Eigen::Map<Eigen::Matrix<double, n_spatial_dimensions,
                                       Eigen::Dynamic, Eigen::RowMajor>>
                  soap_gradient_dim_N(soap_gradient_by_species_pair.data(),
                                      n_spatial_dimensions,
                                      grad_component_size);
              const Eigen::Map<const Eigen::VectorXd> soap_vector_N(
                  soap_vector_by_species_pair.data(), grad_component_size);
              soap_vector_dot_neigh_gradient +=
                  (soap_gradient_dim_N * soap_vector_N);
            }
            // Then, update each species-pair-block
            for (auto soap_grad_spair : soap_neigh_gradient) {
              auto & soap_gradient_by_species_pair{soap_grad_spair.second};
              const auto & soap_vector_by_species_pair{
                  soap_vector[soap_grad_spair.first]};
              Eigen::Map<Eigen::Matrix<double, n_spatial_dimensions,
                                       Eigen::Dynamic, Eigen::RowMajor>>
                  soap_gradient_dim_N(soap_gradient_by_species_pair.data(),
                                      n_spatial_dimensions,
                                      grad_component_size);
              const Eigen::Map<const Eigen::VectorXd> soap_vector_N(
                  soap_vector_by_species_pair.data(), grad_component_size);
              soap_gradient_dim_N =
                  ((soap_gradient_dim_N - soap_vector_dot_neigh_gradient *
                                              soap_vector_N.transpose()) /
                   soap_vector_norm);
            }  // for soap_grad_spair : soap_center_gradient
          }    // for neigh : center
        }      // if normalize
      }        // if compute gradients
    }          // for center : manager
  }            // compute_powerspectrum()

  template <
      internal::SphericalInvariantsType BodyOrder,
      std::enable_if_t<
          BodyOrder == internal::SphericalInvariantsType::RadialSpectrum, int>,
      class StructureManager>
  void CalculatorSphericalInvariants::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using PropExp_t =
        typename CalculatorSphericalExpansion::Property_t<StructureManager>;
    using PropGradExp_t =
        typename CalculatorSphericalExpansion::PropertyGradient_t<
            StructureManager>;
    using Prop_t = Property_t<StructureManager>;
    using PropGrad_t = PropertyGradient_t<StructureManager>;
    using math::pow;

    rep_expansion.compute(manager);

    auto && expansions_coefficients{
        manager->template get_property_ref<PropExp_t>(
            rep_expansion.get_name())};

    auto & expansions_coefficients_gradient{
        manager->template get_property_ref<PropGradExp_t>(
            rep_expansion.get_gradient_name())};

    auto && soap_vectors{
        manager->template get_property_ref<Prop_t>(this->get_name())};

    auto && soap_vector_gradients{
        manager->template get_property_ref<PropGrad_t>(
            this->get_gradient_name())};
    // if the representation has already been computed for the current
    // structure then do nothing
    if (soap_vectors.is_updated()) {
      return;
    }

    this->initialize_per_center_radialspectrum_soap_vectors(
        soap_vectors, soap_vector_gradients, expansions_coefficients, manager);
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

      if (this->compute_gradients) {
        auto ii_pair = center.get_atom_ii();
        auto & grad_center_coefficients{
            expansions_coefficients_gradient[ii_pair]};
        auto & soap_center_gradient{soap_vector_gradients[ii_pair]};

        for (const auto & el : grad_center_coefficients) {
          element_type[0] = el.first[0];
          auto && coef_grad_center{el.second};
          soap_center_gradient[element_type] += coef_grad_center;
        }

        for (auto neigh : center) {
          auto && grad_neigh_coefficients{
              expansions_coefficients_gradient[neigh]};
          auto && soap_neigh_gradient{soap_vector_gradients[neigh]};
          for (const auto & el : grad_neigh_coefficients) {
            element_type[0] = el.first[0];
            auto && coef_grad_neigh{el.second};
            soap_neigh_gradient[element_type] += coef_grad_neigh;
          }
        }  // for (auto neigh : center)

        if (this->normalize) {
          double coefficients_norm_inv{1. / coefficients.norm()};
          double coefficients_norm_inv3{math::pow(coefficients_norm_inv, 3_n)};

          soap_center_gradient.multiply_elements_by(coefficients_norm_inv);

          Eigen::Array3d norm_grad_center = Eigen::Array3d::Zero();
          for (const auto & el : grad_center_coefficients) {
            element_type[0] = el.first[0];
            auto && coef_grad_center{el.second};
            auto && coef_by_type{coefficients[element_type]};
            for (size_t cartesian_idx{0}; cartesian_idx < 3; ++cartesian_idx) {
              size_t cartesian_offset_n{cartesian_idx * this->max_radial};
              norm_grad_center[cartesian_idx] +=
                  (coef_grad_center
                       .block(cartesian_offset_n, 0, this->max_radial, 1)
                       .array() *
                   coef_by_type.array())
                      .sum();
            }
          }
          norm_grad_center *= coefficients_norm_inv3;

          for (const auto & el : coefficients) {
            element_type[0] = el.first[0];
            auto && coef{el.second};
            auto && soap_center_gradient_by_type{
                soap_center_gradient[element_type]};
            for (size_t cartesian_idx{0}; cartesian_idx < 3; ++cartesian_idx) {
              size_t cartesian_offset_n{cartesian_idx * this->max_radial};
              soap_center_gradient_by_type.block(cartesian_offset_n, 0,
                                                 this->max_radial, 1) -=
                  coef * norm_grad_center[cartesian_idx];
            }
          }

          for (auto neigh : center) {
            auto && grad_neigh_coefficients{
                expansions_coefficients_gradient[neigh]};
            auto && soap_neigh_gradient{soap_vector_gradients[neigh]};

            soap_neigh_gradient.multiply_elements_by(coefficients_norm_inv);

            Eigen::Array3d norm_grad_neigh = Eigen::Array3d::Zero();
            for (const auto & el : grad_neigh_coefficients) {
              element_type[0] = el.first[0];
              auto && coef_grad_neigh{el.second};
              auto && coef_by_type{coefficients[element_type]};
              for (size_t cartesian_idx{0}; cartesian_idx < 3;
                   ++cartesian_idx) {
                size_t cartesian_offset_n{cartesian_idx * this->max_radial};
                norm_grad_neigh[cartesian_idx] +=
                    (coef_grad_neigh
                         .block(cartesian_offset_n, 0, this->max_radial, 1)
                         .array() *
                     coef_by_type.array())
                        .sum();
              }
            }
            norm_grad_neigh *= coefficients_norm_inv3;

            for (const auto & el : coefficients) {
              element_type[0] = el.first[0];
              auto && coef{el.second};
              auto && soap_neigh_gradient_by_type{
                  soap_neigh_gradient[element_type]};
              for (size_t cartesian_idx{0}; cartesian_idx < 3;
                   ++cartesian_idx) {
                size_t cartesian_offset_n{cartesian_idx * this->max_radial};
                soap_neigh_gradient_by_type.block(cartesian_offset_n, 0,
                                                  this->max_radial, 1) -=
                    coef * norm_grad_neigh[cartesian_idx];
              }
            }
          }  // for (auto neigh : center)
        }    // if (this->normalize)
      }      // if (this->compute_gradients)
    }        // for (auto center : manager)
  }

  template <
      internal::SphericalInvariantsType BodyOrder,
      std::enable_if_t<
          BodyOrder == internal::SphericalInvariantsType::BiSpectrum, int>,
      class StructureManager>
  void CalculatorSphericalInvariants::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using PropExp_t =
        typename CalculatorSphericalExpansion::Property_t<StructureManager>;
    // using PropGradExp_t = typename
    // CalculatorSphericalExpansion::PropertyGradient_t<StructureManager>;
    using Prop_t = Property_t<StructureManager>;
    // using PropGrad_t = PropertyGradient_t<StructureManager>;
    using internal::enumValue;
    using internal::SphericalInvariantsType;
    using math::pow;

    auto precomputation{downcast_spherical_invariants_precompute<
        SphericalInvariantsType::BiSpectrum>(
        this->precompute_spherical_invariants[enumValue(
            SphericalInvariantsType::BiSpectrum)])};
    auto & wigner_w3js{precomputation->wigner_w3js};

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

    this->initialize_per_center_bispectrum_soap_vectors(
        soap_vectors, expansions_coefficients, manager);

    using complex = std::complex<double>;

    // factor that takes into acount the missing equivalent off diagonal
    // element with respect to the key (or species) index
    double mult{1.0};
    Key_t trip_type{0, 0, 0};
    internal::SortedKey<Key_t> triplet_type{trip_type};
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};

      for (const auto & el1 : coefficients) {
        triplet_type[0] = el1.first[0];
        auto & coef1{el1.second};
        for (const auto & el2 : coefficients) {
          triplet_type[1] = el2.first[0];
          auto & coef2{el2.second};
          for (const auto & el3 : coefficients) {
            triplet_type[2] = el3.first[0];
            auto & coef3{el3.second};
            // triplet multiplicity
            // all the same
            if (triplet_type[0] == triplet_type[1] &&
                triplet_type[1] == triplet_type[2]) {
              mult = 1.0;
            } else if (triplet_type[0] == triplet_type[1] ||
                       triplet_type[0] == triplet_type[2] ||
                       triplet_type[1] == triplet_type[2]) {
              // two the same
              mult = std::sqrt(3.0);
            } else {  // all different
              mult = std::sqrt(6.0);
            }

            if (soap_vector.count(triplet_type) == 1) {
              auto && soap_vector_by_type{soap_vector[triplet_type]};

              size_t nn{0};
              for (size_t n1{0}; n1 < this->max_radial; n1++) {
                for (size_t n2{0}; n2 < this->max_radial; n2++) {
                  for (size_t n3{0}; n3 < this->max_radial; n3++) {
                    size_t l0{0};
                    int wigner_count{0};
                    for (size_t l1{0}; l1 < this->max_angular + 1; l1++) {
                      for (size_t l2{0}; l2 < this->max_angular + 1; l2++) {
                        for (size_t l3{0}; l3 < this->max_angular + 1; l3++) {
                          if (this->inversion_symmetry == true) {
                            if ((l1 + l2 + l3) % 2 == 1) {
                              continue;
                            }
                          }

                          if ((l1 <
                               static_cast<size_t>(std::abs<int>(l2 - l3))) ||
                              (l1 > l2 + l3)) {
                            continue;
                          }

                          for (size_t m1{0}; m1 < 2 * l1 + 1; m1++) {
                            int m1s{static_cast<int>(m1 - l1)};
                            size_t lm1{math::pow(l1, 2_n) + m1};
                            for (size_t m2{0}; m2 < 2 * l2 + 1; m2++) {
                              int m2s{static_cast<int>(m2 - l2)};
                              size_t lm2{math::pow(l2, 2_n) + m2};
                              for (size_t m3{0}; m3 < 2 * l3 + 1; m3++) {
                                int m3s{static_cast<int>(m3 - l3)};
                                if (m1s + m2s + m3s != 0) {
                                  continue;
                                }
                                size_t lm3{math::pow(l3, 2_n) + m3};
                                double w3j = wigner_w3js[wigner_count];
                                complex coef1c, coef2c, coef3c;
                                // usual formulae for converting from real to
                                // complex
                                // see src/math/spherical_harmonics.hh for the
                                // inverse transformation
                                // TODO(andrea, michael) avoid complex
                                // arithmetic
                                if (m1s > 0) {
                                  coef1c = math::pow(-1, m1s) *
                                           complex{coef1(n1, lm1),
                                                   coef1(n1, lm1 - 2 * m1s)};
                                } else if (m1s == 0) {
                                  coef1c = complex(coef1(n1, lm1), 0.0) *
                                           math::SQRT_TWO;
                                } else if (m1s < 0) {
                                  coef1c = complex{coef1(n1, lm1 - 2 * m1s),
                                                   -coef1(n1, lm1)};
                                }
                                if (m2s > 0) {
                                  coef2c = math::pow(-1.0, m2s) *
                                           complex{coef2(n2, lm2),
                                                   coef2(n2, lm2 - 2 * m2s)};
                                } else if (m2s == 0) {
                                  coef2c = complex(coef2(n2, lm2), 0.0) *
                                           math::SQRT_TWO;
                                } else if (m2s < 0) {
                                  coef2c = complex{coef2(n2, lm2 - 2 * m2s),
                                                   -coef2(n2, lm2)};
                                }
                                if (m3s > 0) {
                                  coef3c = math::pow(-1.0, m3s) *
                                           complex{coef3(n3, lm3),
                                                   coef3(n3, lm3 - 2 * m3s)};
                                } else if (m3s == 0) {
                                  coef3c = complex(coef3(n3, lm3), 0.0) *
                                           math::SQRT_TWO;
                                } else if (m3s < 0) {
                                  coef3c = complex{coef3(n3, lm3 - 2 * m3s),
                                                   -coef3(n3, lm3)};
                                }
                                coef1c *= math::INV_SQRT_TWO;
                                coef2c *= math::INV_SQRT_TWO;
                                coef3c *= math::INV_SQRT_TWO;

                                // The descriptor components are purely real or
                                // imaginary
                                // depending on the divisibility of l1 + l2 +l3
                                // by 2.
                                if ((l1 + l2 + l3) % 2 == 0) {
                                  soap_vector_by_type(nn, l0) +=
                                      w3j * mult *
                                      (coef1c * coef2c * coef3c).real();
                                } else {
                                  soap_vector_by_type(nn, l0) +=
                                      w3j * mult *
                                      (coef1c * coef2c * coef3c).imag();
                                }
                                wigner_count++;
                              }  // m3
                            }    // m2
                          }      // m1
                          l0++;
                        }  // l3
                      }    // l2
                    }      // l1
                    nn++;
                  }  // n3
                }    // n2
              }      // n1
            }        // if count triplet
          }          // coef3
        }            // coef2
      }              // coef1

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }
    }  // center
  }    // end function

  template <class StructureManager, class Invariants, class ExpansionCoeff>
  void
  CalculatorSphericalInvariants::initialize_per_center_bispectrum_soap_vectors(
      Invariants & soap_vectors, ExpansionCoeff & expansions_coefficients,
      std::shared_ptr<StructureManager> manager) {
    size_t n_row{math::pow(this->max_radial, 3_n)};
    size_t n_col{0};
    double max_ang{static_cast<double>(this->max_angular)};
    if (this->inversion_symmetry == false) {
      n_col = static_cast<size_t>(1.0 + 2.0 * max_ang +
                                  1.5 * math::pow(max_ang, 2_n) +
                                  math::pow(max_ang, 3_n) * 0.5);
    } else {
      n_col = static_cast<size_t>(std::floor(
          ((math::pow(max_ang + 1.0, 2_n) + 1) * (2 * (max_ang + 1.0) + 3)) /
          8.0));
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

      std::vector<internal::SortedKey<Key_t>> triplet_list{};
      auto center_type{center.get_atom_type()};
      Key_t triplet_type{center_type, center_type, center_type};
      // TODO(felix) optimize this loop
      for (const auto & el1 : coefficients) {
        triplet_type[0] = el1.first[0];
        for (const auto & el2 : coefficients) {
          triplet_type[1] = el2.first[0];
          for (const auto & el3 : coefficients) {
            triplet_type[2] = el3.first[0];
            triplet_list.emplace_back(is_not_sorted, triplet_type);
          }
        }
      }
      // initialize the power spectrum with the proper dimension
      soap_vector.resize(triplet_list, n_row, n_col, 0);
    }
  }

  template <class StructureManager, class Invariants,
            class InvariantsDerivative, class ExpansionCoeff>
  void CalculatorSphericalInvariants::
      initialize_per_center_powerspectrum_soap_vectors(
          Invariants & soap_vectors,
          InvariantsDerivative & soap_vector_gradients,
          ExpansionCoeff & expansions_coefficients,
          std::shared_ptr<StructureManager> manager) {
    constexpr static int n_spatial_dimensions = StructureManager::dim();
    size_t n_row{math::pow(this->max_radial, 2_n)};
    size_t n_col{this->max_angular + 1};

    // clear the data container and resize it
    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);
    soap_vectors.resize();

    if (this->compute_gradients) {
      soap_vector_gradients.clear();
      soap_vector_gradients.set_shape(n_spatial_dimensions * n_row, n_col);
      soap_vector_gradients.resize();
    }

    // identify the species in each environment and initialize soap_vectors
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
      internal::Sorted<true> is_sorted{};

      std::vector<internal::SortedKey<Key_t>> pair_list{};
      auto center_type{center.get_atom_type()};
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

      if (this->compute_gradients) {
        auto ii_pair = center.get_atom_ii();
        // The gradient wrt center is nonzero for all species pairs
        soap_vector_gradients[ii_pair].resize(
            pair_list, n_spatial_dimensions * n_row, n_col, 0.);

        // TODO(max,felix) needs work
        /*
        // Neighbour gradients need a separate pair list because if the species
        // of j is not the same as either of the species for that SOAP entry,
        // the gradient is zero.
        for (auto neigh : center) {
          std::vector<internal::SortedKey<Key_t>> grad_pair_list{};
          for (const auto & el1 : coefficients) {
            auto && neigh_1_type{el1.first[0]};
            for (const auto & el2 : coefficients) {
              auto && neigh_2_type{el2.first[0]};
              auto neigh_type = neigh.get_atom_type();
              if (neigh_1_type <= neigh_2_type) {
                pair_type[0] = neigh_1_type;
                pair_type[1] = neigh_2_type;
                if ((neigh_type == pair_type[0]) or
                    (neigh_type == pair_type[1])) {
                  grad_pair_list.emplace_back(is_sorted, pair_type);
                }
              }
            }
          }
          */
        for (auto neigh : center) {
          soap_vector_gradients[neigh].resize(
              pair_list, n_spatial_dimensions * n_row, n_col, 0.);
        }
      }  // if compute gradients
    }    // for center : manager
  }

  template <class StructureManager, class Invariants,
            class InvariantsDerivative, class ExpansionCoeff>
  void CalculatorSphericalInvariants::
      initialize_per_center_radialspectrum_soap_vectors(
          Invariants & soap_vectors,
          InvariantsDerivative & soap_vector_gradients,
          ExpansionCoeff & expansions_coefficients,
          std::shared_ptr<StructureManager> manager) {
    constexpr static int n_spatial_dimensions = StructureManager::dim();
    size_t n_row{this->max_radial};
    size_t n_col{1};

    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);
    soap_vectors.resize();

    if (this->compute_gradients) {
      soap_vector_gradients.clear();
      soap_vector_gradients.set_shape(n_spatial_dimensions * n_row, n_col);
      soap_vector_gradients.resize();
    }

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

      if (this->compute_gradients) {
        auto ii_pair = center.get_atom_ii();
        // The gradient wrt center is nonzero for all species pairs
        soap_vector_gradients[ii_pair].resize(
            keys, n_spatial_dimensions * n_row, n_col, 0);
        for (auto neigh : center) {
          soap_vector_gradients[neigh].resize(
              keys, n_spatial_dimensions * n_row, n_col, 0);
        }
      }  // if compute gradients
    }
  }

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
