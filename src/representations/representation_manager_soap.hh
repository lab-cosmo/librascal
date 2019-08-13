/**
 * @file   representation_manager_soap.hh
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
    enum class SOAPType { RadialSpectrum, PowerSpectrum, End_ };

    /**
     * Base class for the specification of the SOAP precomputations
     * (it's different for radial vs. full spectrum)
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
        this->l_factors.resize(this->max_angular + 1);

        for (size_t l{0}; l < this->max_angular + 1; ++l) {
          double l_factor{math::pow(std::sqrt(2 * l + 1), -1)};
          this->l_factors(l) = l_factor;
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

  template <class StructureManager>
  class RepresentationManagerSOAP : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using Hypers_t = RepresentationManagerBase::Hypers_t;
    using Key_t = std::vector<int>;
    using SparseProperty_t =
        BlockSparseProperty<double, 1, 0, Manager_t, Key_t>;
    using SparsePropertyGradient_t =
        BlockSparseProperty<double, 2, 0, Manager_t, Key_t>;
    using Data_t = typename SparseProperty_t::Data_t;

    RepresentationManagerSOAP(ManagerPtr_t sm, const Hypers_t & hyper)
        : soap_vectors{*sm}, soap_vector_gradients{*sm}, structure_manager{sm},
          rep_expansion{std::move(sm), hyper} {
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
      using internal::enumValue;
      using internal::SOAPType;

      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      this->normalize = hypers.at("normalize").get<bool>();
      this->soap_type_str = hypers.at("soap_type").get<std::string>();
      if (hypers.find("compute_gradients") != hypers.end()) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }

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

    std::vector<Precision_t> & get_representation_raw_data() {
      return this->dummy;
    }

    Data_t & get_representation_sparse_raw_data() {
      return this->soap_vectors.get_raw_data();
    }

    /**
     * Return a reference to the internal sparse data storage
     *
     * @todo(max) this should really be a const reference, but that screws
     * things up further down the line (when indexing the sparse property)
     */
    SparseProperty_t & get_representation_sparse() {
      return this->soap_vectors;
    }

    /**
     * Return a reference to the internal sparse storage of the gradients
     */
    SparsePropertyGradient_t & get_gradient_sparse() {
      return this->soap_vector_gradients;
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
    SparsePropertyGradient_t soap_vector_gradients;

    //! initialize the soap vectors with only the keys needed for each center
    void initialize_percenter_powerspectrum_soap_vectors();

    void initialize_percenter_radialspectrum_soap_vectors();

   protected:
    size_t max_radial{};
    size_t max_angular{};
    bool normalize{};
    bool compute_gradients{};
    ManagerPtr_t structure_manager;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;
    internal::SOAPType soap_type{};
    //! collection of precomputation for the different body order
    std::array<std::shared_ptr<internal::SOAPPrecomputationBase>,
               internal::enumSize<internal::SOAPType>()>
        precompute_soap{};
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

  // template <class Mngr>
  // void RepresentationManagerSOAP<Mngr>::compute_bispectrum() {
  //   rep_expansion.compute();
  //   auto& expansions_coefficients{rep_expansion.expansions_coefficients};

  //   size_t n_row{pow(this->max_radial, 3)};
  //   size_t n_col{*pow(this->max_angular, 3)*pow((2*this->max_angular + 1),
  // 3)};

  //   this->soap_vectors.clear();
  //   this->soap_vectors.set_shape(n_row, n_col);
  //   this->soap_vectors.resize();

  //   for (auto center : this->structure_manager) {
  //     auto& coefficients{expansions_coefficients[center]};
  //     auto& soap_vector{this->soap_vectors[center]};
  //     Key_t triplet_type{0, 0, 0};
  //     for (const auto& el1 : coefficients) {
  //       triplet_type[0] = el1.first[0];
  //       auto& coef1{el1.second};
  //       for (const auto& el2 : coefficients) {
  //         triplet_type[1] = el2.first[0];
  //         auto& coef2{el2.second};
  //         for (const auto& el2 : coefficients) {
  //           triplet_type[2] = el3.first[0];
  //           auto& coef3{el3.second};

  //           if (soap_vector.count(triplet_type) == 0) {
  //             soap_vector[triplet_type] = Dense_t::Zero(n_row, n_col);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_powerspectrum() {
    using internal::enumValue;
    using internal::n_spatial_dimensions;
    using internal::SOAPType;
    using math::pow;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_soap_precompute<SOAPType::PowerSpectrum>(
        this->precompute_soap[enumValue(SOAPType::PowerSpectrum)])};
    auto & l_factors{precomputation->l_factors};

    // Compute the spherical expansions of the current structure
    rep_expansion.compute();
    auto & expansions_coefficients{rep_expansion.expansions_coefficients};
    // No error if gradients not computed; just an empty array in that case
    auto & expansions_coefficients_gradient{
        rep_expansion.expansions_coefficients_gradient};

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
      soap_vector.multiply_offdiagonal_elements_by(math::SQRT_TWO);

      // normalize the soap vector
      double soap_vector_norm{1.0};
      if (this->normalize) {
        soap_vector_norm = soap_vector.norm();
        soap_vector.normalize();
      }

      if (this->compute_gradients) {
        auto & grad_center_coefficients{
            expansions_coefficients_gradient[center]};
        auto & soap_center_gradient{this->soap_vector_gradients[center]};
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
              auto & soap_neigh_gradient{this->soap_vector_gradients[neigh]};
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
                      cartesian_idx *
                      static_cast<size_t>(math::pow(this->max_radial, 2))};
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
                      cartesian_idx *
                      static_cast<size_t>(math::pow(this->max_radial, 2))};
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
        /*
        soap_center_gradient.multiply_offdiagonal_elements_by(
                math::SQRT_TWO);
        for (auto neigh : center) {
          auto & soap_neigh_gradient{this->soap_vector_gradients[neigh]};
          soap_neigh_gradient.multiply_offdiagonal_elements_by(
                  math::SQRT_TWO);
        }
        */

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
          // TODO(max,felix) is it appropriate to use FeatureManagerBlockSparse
          // in this situation?
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
            auto & soap_neigh_gradient{this->soap_vector_gradients[neigh]};
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

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_radialspectrum() {
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
  void RepresentationManagerSOAP<
      Mngr>::initialize_percenter_powerspectrum_soap_vectors() {
    using internal::n_spatial_dimensions;
    size_t n_row{static_cast<size_t>(pow(this->max_radial, 2))};
    size_t n_col{this->max_angular + 1};

    // clear the data container and resize it
    this->soap_vectors.clear();
    this->soap_vectors.set_shape(n_row, n_col);
    this->soap_vectors.resize();

    if (this->compute_gradients) {
      this->soap_vector_gradients.clear();
      this->soap_vector_gradients.set_shape(n_spatial_dimensions * n_row,
                                            n_col);
      this->soap_vector_gradients.resize();
    }

    auto & expansions_coefficients{rep_expansion.expansions_coefficients};

    // TODO(max) can we use the same pair_list for the gradients as for the
    // expansion coefficients?

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
        // The gradient wrt center is nonzero for all species pairs
        soap_vector_gradients[center].resize(
            pair_list, n_spatial_dimensions * n_row, n_col);

        //TODO(max,felix) needs work
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
              pair_list, n_spatial_dimensions * n_row, n_col);
        }
      }  // if compute gradients
    }    // for center : manager
  }

  template <class Mngr>
  void RepresentationManagerSOAP<
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
