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
    using internal::SOAPType;
    using math::pow;

    // get the relevant precomputation object and unpack the useful infos
    auto precomputation{downcast_soap_precompute<SOAPType::PowerSpectrum>(
        this->precompute_soap[enumValue(SOAPType::PowerSpectrum)])};
    auto & l_factors{precomputation->l_factors};

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

        // multiply with the precomputed factors
        auto coef1{el1.second * l_factors.asDiagonal()};

        for (const auto & el2 : coefficients) {
          pair_type[1] = el2.first[0];
          // avoid computing p^{ab} and p^{ba} since p^{ab} = p^{ba}^T
          if (pair_type[0] > pair_type[1]) {
            continue;
          }

          auto & coef2{el2.second};

          auto && soap_vector_by_pair{soap_vector[pair_type]};
          // TODO(felix) understand why this is slower than below
          // size_t n1n2{0};
          // auto& n_max{this->max_radial};
          // auto l_max{this->max_angular + 1};
          // for (size_t n1{0}; n1 < n_max; ++n1) {
          //   for (size_t l{0}; l < l_max; ++l) {
          //     auto& pos{lm_blocks[l][0]};
          //     auto& size{lm_blocks[l][1]};
          //     // do the reduction over m and iteration over n2
          //     // (with vectorization)
          //     soap_vector_by_pair.block(n1n2, l, n_max, 1).noalias() =
          //     (coef2.block(0, pos, n_max, size) * coef1.row(n1).segment(pos,
          //     size).asDiagonal()).rowwise().sum();
          //   }
          //   n1n2 += n_max;
          // }
          size_t n1n2{0};
          size_t pos, size;
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
      for (const auto & el : soap_vector) {
        auto && pair_type{el.first};
        if (pair_type[0] != pair_type[1]) {
          soap_vector[pair_type] *= math::SQRT_TWO;
        }
      }

      // normalize the soap vector
      if (this->normalize) {
        soap_vector.normalize();
      }
    }
  }

  template <class Mngr>
  void RepresentationManagerSOAP<Mngr>::compute_radialspectrum() {
    rep_expansion.compute();
    using math::pow;

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
      soap_vector.resize(keys, n_row, n_col, 0);

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
