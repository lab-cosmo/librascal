/**
 * @file   rascal/representations/calculator_spherical_invariants.hh
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

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_

#include "rascal/math/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/structure_managers/property_block_sparse.hh"
#include "rascal/structure_managers/structure_manager.hh"
#include "rascal/utils/utils.hh"

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
    };

    //! factor of 1 / sqrt(2*l+1) in front of the powerspectrum
    inline Eigen::VectorXd precompute_l_factors(size_t max_angular) {
      Eigen::VectorXd l_factors{};
      l_factors.resize(max_angular + 1);
      for (size_t l{0}; l < max_angular + 1; ++l) {
        l_factors(l) = math::pow(std::sqrt(2 * l + 1), -1);
      }
      return l_factors;
    }

    //! Wigner pre-factors
    inline Eigen::ArrayXd precompute_wigner_w3js(size_t max_angular,
                                                 bool inversion_symmetry) {
      Eigen::ArrayXd wigner_w3js{};
      // get the number of non zero elements in the w3j
      int n_elements{0};
      for (size_t l1{0}; l1 < max_angular + 1; ++l1) {
        for (size_t l2{0}; l2 < max_angular + 1; ++l2) {
          for (size_t l3{0}; l3 < max_angular + 1; ++l3) {
            if ((l1 < static_cast<size_t>(std::abs<int>(l2 - l3))) ||
                (l1 > l2 + l3)) {
              continue;
            }
            if (inversion_symmetry) {
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

      wigner_w3js.resize(n_elements);
      n_elements = 0;
      wig_table_init(2 * (max_angular + 1), 3);
      wig_temp_init(2 * (max_angular + 1));
      for (size_t l1{0}; l1 < max_angular + 1; ++l1) {
        for (size_t l2{0}; l2 < max_angular + 1; ++l2) {
          for (size_t l3{0}; l3 < max_angular + 1; ++l3) {
            if ((l1 < static_cast<size_t>(std::abs<int>(l2 - l3))) ||
                (l1 > l2 + l3)) {
              continue;
            }
            if (inversion_symmetry) {
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
                  wigner_w3js(n_elements) =
                      wig3jj(2 * l1, 2 * l2, 2 * l3, 2 * m1s, 2 * m2s, 2 * m3s);
                  ++n_elements;
                }
              }
            }
          }
        }
      }
      wig_temp_free();
      wig_table_free();

      return wigner_w3js;
    }
  }  // namespace internal

  class CalculatorSphericalInvariants : public CalculatorBase {
   public:
    using Hypers_t = typename CalculatorBase::Hypers_t;
    using Key_t = typename CalculatorBase::Key_t;

    template <class StructureManager>
    using Property_t = BlockSparseProperty<double, 1, StructureManager, Key_t>;

    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, StructureManager, Key_t>;

    template <class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;

    template <class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    template <class StructureManager>
    using SpectrumNorm_t = Property<double, 1, StructureManager, 1>;

    explicit CalculatorSphericalInvariants(const Hypers_t & hypers)
        : CalculatorBase{}, rep_expansion{hypers} {
      this->set_default_prefix("spherical_invariants_");
      this->set_hyperparameters(hypers);
      this->hypers = hypers;
    }

    //! Copy constructor
    CalculatorSphericalInvariants(const CalculatorSphericalInvariants & other) =
        delete;

    //! Move constructor
    CalculatorSphericalInvariants(
        CalculatorSphericalInvariants && other) noexcept
        : CalculatorBase{std::move(other)},
          unique_pair_list{other.unique_pair_list},
          coeff_indices_map{other.coeff_indices_map}, key_map{other.key_map},
          coeff_indices{other.coeff_indices},
          is_sparsified{other.is_sparsified},
          max_radial{std::move(other.max_radial)}, max_angular{std::move(
                                                       other.max_angular)},
          shape{std::move(other.shape)}, normalize{std::move(other.normalize)},
          compute_gradients{std::move(other.compute_gradients)},
          inversion_symmetry{std::move(other.inversion_symmetry)},
          rep_expansion{std::move(other.rep_expansion)},
          type{std::move(other.type)}, l_factors{std::move(other.l_factors)},
          wigner_w3js{std::move(other.wigner_w3js)} {}
    //! Destructor
    virtual ~CalculatorSphericalInvariants() = default;

    //! Copy assignment operator
    CalculatorSphericalInvariants &
    operator=(const CalculatorSphericalInvariants & other) = delete;

    //! Move assignment operator
    CalculatorSphericalInvariants &
    operator=(CalculatorSphericalInvariants && other) = default;

    //!
    struct PowerSpectrumCoeffIndex {
      PowerSpectrumCoeffIndex(std::uint32_t n1, std::uint32_t n2,
                              std::uint32_t l, std::uint32_t n1n2,
                              std::uint32_t l_block_size,
                              std::uint32_t l_block_idx, double l_factor)
          : n1{n1}, n2{n2}, l_block_size{l_block_size},
            l_block_idx{l_block_idx}, l{l}, n1n2{n1n2}, l_factor{l_factor} {}
      // indices refering to coefficients of the spherical expansion
      std::uint32_t n1;
      std::uint32_t n2;
      std::uint32_t l_block_size;
      std::uint32_t l_block_idx;
      // indices refering to the power spectrum
      std::uint32_t l;
      std::uint32_t n1n2;
      double l_factor;

      bool operator==(const PowerSpectrumCoeffIndex & other) const {
        bool ret{this->n1 == other.n1 and this->n2 == other.n2 and
                 this->l_block_size == other.l_block_size and
                 this->l_block_idx == other.l_block_idx and
                 this->l == other.l and this->n1n2 == other.n1n2 and
                 math::relative_error(this->l_factor, other.l_factor) < 1e-10};
        return ret;
      }
    };
    //! list of possible input pairs when using coefficient_subselection
    std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
        unique_pair_list{};
    //! list of coefficient that will be computed, mapping keys to a list
    //! invariant coefficients that will be computed
    std::map<internal::SortedKey<Key_t>, std::vector<PowerSpectrumCoeffIndex>,
             internal::CompareSortedKeyLess>
        coeff_indices_map{};
    //! map the possible keys of the coefficient to the keys of the
    //! powerspectrum. If sparsification then it maps to (0,0) always otherwise
    //! it maps onto itself
    std::map<internal::SortedKey<Key_t>, internal::SortedKey<Key_t>,
             internal::CompareSortedKeyLess>
        key_map{};
    //! list of coefficient that will be computed if there is no sparsification
    std::vector<PowerSpectrumCoeffIndex> coeff_indices{};
    //! tell if the calculator has been sparsified using
    //! the coefficient_subselection input
    bool is_sparsified{false};

    void set_hyperparameters(const Hypers_t & hypers) override {
      using internal::SphericalInvariantsType;

      this->normalize = hypers.at("normalize").get<bool>();
      auto soap_type = hypers.at("soap_type").get<std::string>();

      if (hypers.find("compute_gradients") != hypers.end()) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }

      this->max_radial = hypers.at("max_radial").get<size_t>();
      this->max_angular = hypers.at("max_angular").get<size_t>();

      if (hypers.find("coefficient_subselection") != hypers.end() and
          (soap_type == "PowerSpectrum")) {  // NOLINT
        this->is_sparsified = true;
        auto coefficient_subselection =
            hypers.at("coefficient_subselection").get<json>();
        // get the indices of the subselected PowerSpectrum coefficients:
        // p_{abn_1n_2l} where a and b refer to atomic species and should be
        // lexicographically sorted so that a <= b, n_1 and n_2 refer to
        // radial basis index and l refer to the angular index of the
        // spherical harmonic
        auto sp_a = coefficient_subselection.at("a")
                        .get<std::vector<typename Key_t::value_type>>();
        auto sp_b = coefficient_subselection.at("b")
                        .get<std::vector<typename Key_t::value_type>>();
        auto radial_n1 =
            coefficient_subselection.at("n1").get<std::vector<std::uint32_t>>();
        auto radial_n2 =
            coefficient_subselection.at("n2").get<std::vector<std::uint32_t>>();
        auto angular_l =
            coefficient_subselection.at("l").get<std::vector<std::uint32_t>>();

        std::uint32_t angular_max{
            *std::max_element(angular_l.begin(), angular_l.end())};
        std::uint32_t radial1_max{
            *std::max_element(radial_n1.begin(), radial_n1.end())};
        std::uint32_t radial2_max{
            *std::max_element(radial_n2.begin(), radial_n2.end())};
        // check if the inputs are sane
        {
          if (sp_a.size() != sp_b.size()) {
            throw std::logic_error(
                R"(coefficient_subselection should have elements of the same size but: a.size() != b.size())");
          } else if (sp_a.size() != radial_n1.size()) {
            throw std::logic_error(
                R"(coefficient_subselection should have elements of the same size but: a.size() != n1.size())");
          } else if (sp_a.size() != radial_n2.size()) {
            throw std::logic_error(
                R"(coefficient_subselection should have elements of the same size but: a.size() != n2.size())");
          } else if (sp_a.size() != angular_l.size()) {
            throw std::logic_error(
                R"(coefficient_subselection should have elements of the same size but: a.size() != l.size())");
          }

          if (angular_max >= this->max_angular + 1) {
            std::stringstream err_str{};
            err_str << "'max(l)' in coefficient_subselection is larger"
                    << " than 'max_angular+1'" << angular_max
                    << " >= " << this->max_angular + 1 << std::endl;
            throw std::runtime_error(err_str.str());
          }
          if (radial1_max >= this->max_radial) {
            std::stringstream err_str{};
            err_str << "'max(n1)' in coefficient_subselection is larger"
                    << " than 'max_radial'" << radial1_max
                    << " >= " << this->max_radial << std::endl;
            throw std::runtime_error(err_str.str());
          }
          if (radial2_max >= this->max_radial) {
            std::stringstream err_str{};
            err_str << "'max(n2)' in coefficient_subselection is larger"
                    << " than 'max_radial'" << radial2_max
                    << " >= " << this->max_radial << std::endl;
            throw std::runtime_error(err_str.str());
          }
        }
        this->shape[0] = 1;
        this->shape[1] = sp_a.size();

        internal::Sorted<true> is_sorted{};
        // collect all unique sorted pairs
        for (size_t i_pair{0}; i_pair < sp_a.size(); ++i_pair) {
          if (sp_a[i_pair] > sp_b[i_pair]) {
            throw std::logic_error(
                R"(coefficient_subselection should have only lexicographically
              sorted pairs of atomic species in the selection but: a > b)");
          }
          Key_t pair{{sp_a[i_pair], sp_b[i_pair]}};
          this->unique_pair_list.insert({is_sorted, pair});
        }
        // initialize coeff_indices_map
        std::vector<PowerSpectrumCoeffIndex> init_vec{};
        for (const auto & key : this->unique_pair_list) {
          this->coeff_indices_map[key] = init_vec;
        }

        // precompute the positions and sizes of the angular momentum channels
        // of the linear storage of the spherical expansion coefficients
        // related to the spherical harmonics
        std::vector<std::uint32_t> l_block_sizes{}, l_block_ids{};
        std::uint32_t pos{0};
        for (std::uint32_t l{0}; l < this->max_angular + 1; ++l) {
          std::uint32_t size{2 * l + 1};
          l_block_sizes.push_back(size);
          l_block_ids.push_back(pos);
          pos += size;
        }

        this->l_factors = internal::precompute_l_factors(this->max_angular);

        // fill coeff_indices_map
        for (size_t i_pair{0}; i_pair < sp_a.size(); ++i_pair) {
          Key_t pair{{sp_a[i_pair], sp_b[i_pair]}};
          internal::SortedKey<Key_t> spair{is_sorted, pair};
          // insertion order n1, n2, l, n1n2, l_block_size, l_block_idx
          // n1n2 is set to zero because we store the sparsified features as
          // one row
          this->coeff_indices_map[spair].emplace_back(
              radial_n1[i_pair], radial_n2[i_pair],
              static_cast<std::uint32_t>(i_pair), 0,
              l_block_sizes[angular_l[i_pair]], l_block_ids[angular_l[i_pair]],
              this->l_factors[angular_l[i_pair]]);
        }

        Key_t sparsified_type{0, 0};
        // there are 118 atomic elements atm
        for (int sp1{1}; sp1 < 130; sp1++) {
          for (int sp2{1}; sp2 < 130; sp2++) {
            if (sp1 <= sp2) {
              Key_t pair_type{sp1, sp2};
              internal::SortedKey<Key_t> spair_type{is_sorted, pair_type};
              internal::SortedKey<Key_t> ssparsified_type{is_sorted,
                                                          sparsified_type};
              this->key_map[spair_type] = ssparsified_type;
            }
          }
        }

      } else {  // Default false (compute all coefficents)
        this->l_factors = internal::precompute_l_factors(this->max_angular);
        this->is_sparsified = false;
        this->shape[0] = math::pow(this->max_radial, 2_size_t);
        this->shape[1] = this->max_angular + 1;
        for (size_t n1{0}; n1 < this->max_radial; ++n1) {
          for (size_t n2{0}; n2 < this->max_radial; ++n2) {
            size_t pos{0};
            for (size_t l{0}; l < this->max_angular + 1; ++l) {
              size_t size{2 * l + 1};
              this->coeff_indices.emplace_back(n1, n2, l,
                                               n1 * this->max_radial + n2, size,
                                               pos, this->l_factors[l]);
              pos += size;
            }
          }
        }
        internal::Sorted<true> is_sorted{};
        // there are 118 atomic elements atm
        for (int sp1{1}; sp1 < 130; sp1++) {
          for (int sp2{1}; sp2 < 130; sp2++) {
            if (sp1 <= sp2) {
              Key_t pair_type{sp1, sp2};
              internal::SortedKey<Key_t> spair_type{is_sorted, pair_type};
              this->key_map[spair_type] = spair_type;
            }
          }
        }
      }

      if (soap_type == "PowerSpectrum") {
        this->type = SphericalInvariantsType::PowerSpectrum;
      } else if (soap_type == "RadialSpectrum") {
        this->type = SphericalInvariantsType::RadialSpectrum;
        if (this->max_angular > 0) {
          throw std::logic_error("max_angular should be 0 with RadialSpectrum");
        }
      } else if (soap_type == "BiSpectrum") {
        this->type = internal::SphericalInvariantsType::BiSpectrum;
        this->inversion_symmetry = hypers.at("inversion_symmetry").get<bool>();
        this->wigner_w3js = internal::precompute_wigner_w3js(
            this->max_angular, this->inversion_symmetry);
      } else {
        throw std::logic_error(
            "Requested SphericalInvariants type \'" + soap_type +
            "\' has not been implemented.  Must be one of" +
            ": 'PowerSpectrum', 'RadialSpectrum', or 'BiSpectrum'.");
      }

      this->set_name(hypers);
    }

    bool operator==(const CalculatorSphericalInvariants & other) const {
      bool grad_match{this->does_gradients() == other.does_gradients()};
      bool main_hypers_match{
          this->max_radial == other.max_radial and
          this->max_angular == other.max_angular and
          this->normalize == other.normalize and
          this->inversion_symmetry == other.inversion_symmetry and
          this->type == other.type and
          (this->l_factors.array() == other.l_factors.array()).all() and
          (this->wigner_w3js.array() == other.wigner_w3js.array()).all()};
      bool rep_expansion_match{this->rep_expansion == other.rep_expansion};
      bool sparsification_match{
          this->unique_pair_list == other.unique_pair_list and
          this->coeff_indices_map == other.coeff_indices_map and
          this->key_map == other.key_map and
          this->coeff_indices == other.coeff_indices and
          this->is_sparsified == other.is_sparsified};
      return (grad_match and main_hypers_match and rep_expansion_match and
              sparsification_match);
    }

    /**
     * Returns if the calculator is able to compute gradients of the
     * representation w.r.t. atomic positions ?
     */
    bool does_gradients() const override { return this->compute_gradients; }

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

    /**
     * Update the gradients \grad_k p^{i} to include normalization, N_i,
     * resulting in \grad_k \tilde{p}^{i}.
     * We have:
     * \grad_k \tilde{p}^{i} = \grad_k p^{i} / N_i
             - \tilde{p}^{i} [\tilde{p}^{i} \cdot \grad_k p^{i} / N_i],
     * where $\cdot$ is a dot product between vectors.
     * Note that this expects the soap vectors to be normalized already, and
     * the norm stored separately.
     */
    template <class StructureManager, class SpectrumNorm>
    void update_gradients_for_normalization(
        Property_t<StructureManager> & soap_vectors,
        PropertyGradient_t<StructureManager> & soap_vector_gradients,
        std::shared_ptr<StructureManager> manager, SpectrumNorm & inv_norms,
        const size_t & grad_component_size) {
      using MapSoapGradFlat_t = Eigen::Map<
          Eigen::Matrix<double, ThreeD, Eigen::Dynamic, Eigen::RowMajor>>;
      using ConstMapSoapFlat_t = const Eigen::Map<const Eigen::VectorXd>;
      // divide all gradients with the normalization factor N_i
      for (auto center : manager) {
        for (auto neigh : center.pairs_with_self_pair()) {
          soap_vector_gradients[neigh].multiply_elements_by(inv_norms[center]);
        }
      }

      // \tilde{p}^{i} \cdot \grad_k p^{i} / N_i
      Eigen::Vector3d soap_vector_dot_gradient{};

      // compute the dot product and update the gradients to be normalized
      for (auto center : manager) {
        const int atom_i_tag{center.get_atom_tag()};
        for (auto neigh : center.pairs_with_self_pair()) {
          auto && atom_j = neigh.get_atom_j();
          const int atom_j_tag{atom_j.get_atom_tag()};
          // const bool is_center_atom{manager->is_center_atom(neigh)};
          // compute grad contribution only if the neighbour is _not_ an image
          // of the center (because then it moves with the center)
          if (atom_j_tag == atom_i_tag and manager->is_ghost_atom(neigh)) {
            continue;
          }
          const auto & soap_vector = soap_vectors[center];
          auto & soap_vector_gradients_by_neigh = soap_vector_gradients[neigh];
          soap_vector_dot_gradient.setZero();
          // make sure to iterate over keys that are present in both soap_vector
          // and soap_vector_gradients_by_neigh
          const auto keys_grad = soap_vector_gradients_by_neigh.get_keys();
          const auto keys_intersect = soap_vector.intersection(keys_grad);
          // compute \tilde{p}^{i} \cdot \grad_k p^{i} / N_i
          for (const auto & key : keys_intersect) {
            auto soap_gradient_by_species_pair =
                soap_vector_gradients_by_neigh[key];
            const auto & soap_vector_by_species_pair = soap_vector[key];
            // reshape for easy dot prod
            MapSoapGradFlat_t soap_gradient_dim_N(
                soap_gradient_by_species_pair.data(), ThreeD,
                grad_component_size);
            ConstMapSoapFlat_t soap_vector_N(soap_vector_by_species_pair.data(),
                                             grad_component_size);
            // dot
            soap_vector_dot_gradient += (soap_gradient_dim_N * soap_vector_N);
          }  // for (const auto& key : keys_intersect)

          // Now update each species-pair-block using the dot-product just
          // computed
          // for (const auto & key : keys_grad) {
          for (const auto & key : soap_vector.get_keys()) {
            auto soap_gradient_by_species_pair =
                soap_vector_gradients_by_neigh[key];
            const auto & soap_vector_by_species_pair = soap_vector[key];
            // reshape for easy dot prod
            MapSoapGradFlat_t soap_gradient_dim_N(
                soap_gradient_by_species_pair.data(), ThreeD,
                grad_component_size);
            ConstMapSoapFlat_t soap_vector_N(soap_vector_by_species_pair.data(),
                                             grad_component_size);
            // compute \tilde{p}^{i} [\tilde{p}^{i} \cdot \grad_k p^{i} / N_i]
            // as an outer product
            soap_gradient_dim_N -=
                soap_vector_dot_gradient * soap_vector_N.transpose();
          }  // for (const auto& key : keys_intersect)
        }    // (auto neigh : center.pairs_with_self_pair())
      }      // (auto center : manager)
    }

   protected:
    size_t max_radial{};
    size_t max_angular{};
    // shape of the inner dense section of the computed invariant coefficients
    std::array<size_t, 2> shape{{0, 0}};
    bool normalize{};
    bool compute_gradients{};
    bool inversion_symmetry{false};

    CalculatorSphericalExpansion rep_expansion;

    internal::SphericalInvariantsType type{};

    //! precomputed l-factors the PowerSpectrum
    Eigen::VectorXd l_factors{};

    //! precomputed wigner symbols for the BiSpectrum
    Eigen::ArrayXd wigner_w3js{};
  };

  template <class StructureManager>
  void CalculatorSphericalInvariants::compute(StructureManager & managers) {
    using internal::SphericalInvariantsType;
    switch (this->type) {
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
    using internal::SphericalInvariantsType;
    using math::pow;

    // Compute the spherical expansions of the current structure
    rep_expansion.compute(manager);

    constexpr bool ExcludeGhosts{true};
    auto && expansions_coefficients{*manager->template get_property<PropExp_t>(
        rep_expansion.get_name(), true, true, ExcludeGhosts)};

    // No error if gradients not computed; just an empty array in that case
    auto && expansions_coefficients_gradient{
        *manager->template get_property<PropGradExp_t>(
            rep_expansion.get_gradient_name(), true, true)};

    auto && soap_vectors{*manager->template get_property<Prop_t>(
        this->get_name(), true, true, ExcludeGhosts)};

    auto && soap_vector_gradients{*manager->template get_property<PropGrad_t>(
        this->get_gradient_name(), true, true)};

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

    // to store the norm of the soap vectors
    SpectrumNorm_t<StructureManager> soap_vector_norm_inv{
        *manager, "power spectrums inverse norms", true};
    soap_vector_norm_inv.resize();

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
          auto && soap_vector_by_pair{soap_vector[this->key_map[spair_type]]};
          const auto & coef_ids{coeff_indices_map[spair_type]};
          for (const auto & coef_idx : coef_ids) {
            // multiply with the constant 1 / \sqrt(2l+1)
            soap_vector_by_pair(coef_idx.n1n2, coef_idx.l) =
                (coef1
                     .block(coef_idx.n1, coef_idx.l_block_idx, 1,
                            coef_idx.l_block_size)
                     .array() *
                 coef2
                     .block(coef_idx.n2, coef_idx.l_block_idx, 1,
                            coef_idx.l_block_size)
                     .array())
                    .sum() *
                coef_idx.l_factor;
          }

          // the \sqrt(2) factor to account for the missing (b,a) components
          if (spair_type[0] < spair_type[1]) {
            for (const auto & coef_idx : coef_ids) {
              soap_vector_by_pair(coef_idx.n1n2, coef_idx.l) *= math::SQRT_TWO;
            }
          }
        }  // for el1 : coefficients
      }    // for el2 : coefficients

      // normalize the soap vector
      if (this->normalize) {
        double norm_inv{1. / soap_vector.normalize_and_get_norm()};
        soap_vector_norm_inv[center] = norm_inv;
      }

      if (this->compute_gradients) {
        const int atom_i_tag{center.get_atom_tag()};
        // c^{i}
        auto & coefficients{expansions_coefficients[center]};
        std::vector<Key_t> keys_coef{coefficients.get_keys()};

        // Sum the gradients wrt the neighbour atom position
        // compute the \grad_k p^{i} coeffs where k is either i or j
        for (auto neigh : center.pairs_with_self_pair()) {
          auto && atom_j = neigh.get_atom_j();
          const int atom_j_tag{atom_j.get_atom_tag()};
          // const bool is_center_atom{manager->is_center_atom(neigh)};
          // compute grad contribution only if the neighbour is _not_ an image
          // of the center (because then it moves with the center)
          if (atom_j_tag == atom_i_tag and manager->is_ghost_atom(neigh)) {
            continue;
          }

          // \grad_k c^{i}
          auto & grad_neigh_coefficients{
              expansions_coefficients_gradient[neigh]};
          std::vector<Key_t> keys_coef_grad_neigh{
              grad_neigh_coefficients.get_keys()};
          // \grad_k p^{i}
          auto & soap_neigh_gradient{soap_vector_gradients[neigh]};

          std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
              grad_neigh_keys{};
          // \grad_k p^{iab} = \grad_k c^{i a} c^{i b} + c^{i a} \grad_k c^{i b}
          // by definition \grad_k c^{i a} is non zero for one key 'a' for k!=i
          // so either a == b and we compute one term with a factor of 2 or only
          // one of the two terms is non zero hence the swap of entry when
          // spair_type[0] > spair_type[1] == true
          for (const auto & coef_key_1 : keys_coef_grad_neigh) {
            // \grad_k c^{i a}
            const auto & grad_neigh_coefficients_1{
                grad_neigh_coefficients[coef_key_1]};

            for (const auto & coef_key_2 : keys_coef) {
              // c^{i b}
              const auto & expansion_coefficients_2{coefficients[coef_key_2]};
              bool sorted{true}, equal{false};
              // make sure spair_type has sorted entries
              if (coef_key_1[0] > coef_key_2[0]) {
                sorted = false;
                spair_type[0] = coef_key_2[0];
                spair_type[1] = coef_key_1[0];
              } else if (coef_key_1[0] == coef_key_2[0]) {
                equal = true;
                spair_type[0] = coef_key_1[0];
                spair_type[1] = coef_key_2[0];
              } else {
                spair_type[0] = coef_key_1[0];
                spair_type[1] = coef_key_2[0];
              }

              grad_neigh_keys.insert(spair_type);
              // \grad_k p^{i ab}
              auto soap_neigh_gradient_by_species_pair{
                  soap_neigh_gradient[this->key_map[spair_type]]};
              const auto & coef_ids{coeff_indices_map[spair_type]};

              // computes  \grad_k c^{i a}_{n_1} c^{i b}_{n_2}
              if (sorted or equal) {
                for (size_t cartesian_idx{0}; cartesian_idx < 3;
                     ++cartesian_idx) {
                  const size_t cartesian_offset_n{cartesian_idx *
                                                  this->max_radial};
                  const size_t cartesian_offset_n1n2{cartesian_idx *
                                                     this->shape[0]};
                  for (const auto & coef_idx : coef_ids) {
                    // clang-format off
                    soap_neigh_gradient_by_species_pair(
                      coef_idx.n1n2 + cartesian_offset_n1n2, coef_idx.l) +=
                      (grad_neigh_coefficients_1.block(
                          coef_idx.n1 + cartesian_offset_n,
                          coef_idx.l_block_idx, 1,
                          coef_idx.l_block_size).array() *
                       expansion_coefficients_2.block(
                          coef_idx.n2, coef_idx.l_block_idx,
                          1,  coef_idx.l_block_size).array()).sum()
                      * coef_idx.l_factor;
                    // clang-format on
                  }  // for const auto& coef_idx : coef_ids
                }    // for cartesian_idx
              }      // if (sorted or equal)

              // computes c^{i a}_{n_1} \grad_k c^{i b}_{n_2}
              if (not sorted or equal) {
                for (size_t cartesian_idx{0}; cartesian_idx < 3;
                     ++cartesian_idx) {
                  const size_t cartesian_offset_n{cartesian_idx *
                                                  this->max_radial};
                  const size_t cartesian_offset_n1n2{cartesian_idx *
                                                     this->shape[0]};
                  for (const auto & coef_idx : coef_ids) {
                    // clang-format off
                    soap_neigh_gradient_by_species_pair(
                      coef_idx.n1n2 + cartesian_offset_n1n2, coef_idx.l) +=
                      (grad_neigh_coefficients_1.block(
                            coef_idx.n2 + cartesian_offset_n,
                            coef_idx.l_block_idx,
                            1, coef_idx.l_block_size).array() *
                       expansion_coefficients_2.block(
                            coef_idx.n1, coef_idx.l_block_idx,
                            1, coef_idx.l_block_size).array()).sum()
                      * coef_idx.l_factor;
                    // clang-format on
                  }  // for const auto& coef_idx : coef_ids
                }    // for cartesian_idx
              }      // if (not sorted or equal)
            }        // keys_coef
          }          // keys_coef_grad_neigh

          // multiply with \sqrt(2) factor to account
          // for the missing (b,a) components
          for (const auto & key : grad_neigh_keys) {
            if (key[0] != key[1]) {
              auto soap_neigh_gradient_by_species_pair{
                  soap_neigh_gradient[this->key_map[key]]};
              const auto & coef_ids{coeff_indices_map[key]};
              for (size_t cartesian_idx{0}; cartesian_idx < 3;
                   ++cartesian_idx) {
                const size_t cartesian_offset_n1n2{cartesian_idx *
                                                   this->shape[0]};
                for (const auto & coef_idx : coef_ids) {
                  soap_neigh_gradient_by_species_pair(
                      coef_idx.n1n2 + cartesian_offset_n1n2, coef_idx.l) *=
                      math::SQRT_TWO;
                }
              }
            }
          }
        }  // for neigh : center
      }    // if compute gradients
    }      // for center : manager

    if (this->normalize and this->compute_gradients) {
      const size_t grad_component_size{this->shape[0] * this->shape[1]};
      this->update_gradients_for_normalization(
          soap_vectors, soap_vector_gradients, manager, soap_vector_norm_inv,
          grad_component_size);
    }  // if normalize and compute_gradients
  }    // compute_powerspectrum()

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
    constexpr bool ExcludeGhosts{true};

    rep_expansion.compute(manager);

    auto && expansions_coefficients{*manager->template get_property<PropExp_t>(
        rep_expansion.get_name(), true, true, ExcludeGhosts)};

    auto & expansions_coefficients_gradient{
        *manager->template get_property<PropGradExp_t>(
            rep_expansion.get_gradient_name(), true, true)};

    auto && soap_vectors{*manager->template get_property<Prop_t>(
        this->get_name(), true, true, ExcludeGhosts)};

    auto && soap_vector_gradients{*manager->template get_property<PropGrad_t>(
        this->get_gradient_name(), true, true)};
    // if the representation has already been computed for the current
    // structure then do nothing
    if (soap_vectors.is_updated()) {
      return;
    }

    this->initialize_per_center_radialspectrum_soap_vectors(
        soap_vectors, soap_vector_gradients, expansions_coefficients, manager);
    Key_t element_type{0};

    // to store the norm of the soap vectors
    SpectrumNorm_t<StructureManager> soap_vector_norm_inv{
        *manager, "radial spectrum inverse norms", ExcludeGhosts};
    soap_vector_norm_inv.resize();

    for (auto center : manager) {
      const auto & coefficients{expansions_coefficients[center]};
      auto & soap_vector{soap_vectors[center]};
      auto keys = coefficients.get_keys();
      for (const auto & key : keys) {
        soap_vector[key] = coefficients[key];
      }

      // normalize the soap vector
      if (this->normalize) {
        double norm_inv{1. / soap_vector.normalize_and_get_norm()};
        soap_vector_norm_inv[center] = norm_inv;
      }

      if (this->compute_gradients) {
        for (auto neigh : center.pairs_with_self_pair()) {
          auto && grad_neigh_coefficients{
              expansions_coefficients_gradient[neigh]};
          auto && soap_neigh_gradient{soap_vector_gradients[neigh]};
          auto keys_grad = grad_neigh_coefficients.get_keys();
          for (const auto & key : keys_grad) {
            soap_neigh_gradient[key] += grad_neigh_coefficients[key];
          }
        }  // for (auto neigh : center.pairs_with_self_pair())
      }    // if (this->compute_gradients)
    }      // for (auto center : manager)

    if (this->normalize and this->compute_gradients) {
      this->update_gradients_for_normalization(
          soap_vectors, soap_vector_gradients, manager, soap_vector_norm_inv,
          this->max_radial);
    }  // if (this->normalize and this->compute_gradients)
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
    using Prop_t = Property_t<StructureManager>;
    using internal::SphericalInvariantsType;
    using math::pow;

    rep_expansion.compute(manager);

    constexpr bool ExcludeGhosts{true};
    auto && expansions_coefficients{*manager->template get_property<PropExp_t>(
        rep_expansion.get_name(), true, true, ExcludeGhosts)};

    auto && soap_vectors{*manager->template get_property<Prop_t>(
        this->get_name(), true, true, ExcludeGhosts)};

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
                            size_t lm1{math::pow(l1, 2_size_t) + m1};
                            for (size_t m2{0}; m2 < 2 * l2 + 1; m2++) {
                              int m2s{static_cast<int>(m2 - l2)};
                              size_t lm2{math::pow(l2, 2_size_t) + m2};
                              for (size_t m3{0}; m3 < 2 * l3 + 1; m3++) {
                                int m3s{static_cast<int>(m3 - l3)};
                                if (m1s + m2s + m3s != 0) {
                                  continue;
                                }
                                size_t lm3{math::pow(l3, 2_size_t) + m3};
                                double w3j = this->wigner_w3js[wigner_count];
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
        soap_vector.normalize_and_get_norm();
      }
    }  // center
  }    // end function

  template <class StructureManager, class Invariants, class ExpansionCoeff>
  void
  CalculatorSphericalInvariants::initialize_per_center_bispectrum_soap_vectors(
      Invariants & soap_vectors, ExpansionCoeff & expansions_coefficients,
      std::shared_ptr<StructureManager> manager) {
    size_t n_row{math::pow(this->max_radial, 3_size_t)};
    size_t n_col{0};
    double max_ang{static_cast<double>(this->max_angular)};
    if (this->inversion_symmetry == false) {
      n_col = static_cast<size_t>(1.0 + 2.0 * max_ang +
                                  1.5 * math::pow(max_ang, 2_size_t) +
                                  math::pow(max_ang, 3_size_t) * 0.5);
    } else {
      n_col = static_cast<size_t>(
          std::floor(((math::pow(max_ang + 1.0, 2_size_t) + 1) *
                      (2 * (max_ang + 1.0) + 3)) /
                     8.0));
    }

    // clear the data container and resize it
    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);

    std::vector<
        std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>>
        keys_list{};
    // identify the species in each environment and initialize soap_vectors
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      internal::Sorted<false> is_not_sorted{};

      std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
          triplet_list{};
      auto center_type{center.get_atom_type()};
      Key_t triplet_type{center_type, center_type, center_type};
      for (const auto & el1 : coefficients) {
        triplet_type[0] = el1.first[0];
        for (const auto & el2 : coefficients) {
          triplet_type[1] = el2.first[0];
          for (const auto & el3 : coefficients) {
            triplet_type[2] = el3.first[0];
            triplet_list.insert({is_not_sorted, triplet_type});
          }
        }
      }
      // initialize the power spectrum with the proper dimension
      keys_list.emplace_back(triplet_list);
    }
    soap_vectors.resize(keys_list);
    soap_vectors.setZero();
  }

  template <class StructureManager, class Invariants,
            class InvariantsDerivative, class ExpansionCoeff>
  void CalculatorSphericalInvariants::
      initialize_per_center_powerspectrum_soap_vectors(
          Invariants & soap_vectors,
          InvariantsDerivative & soap_vector_gradients,
          ExpansionCoeff & expansions_coefficients,
          std::shared_ptr<StructureManager> manager) {
    if (this->is_sparsified) {
      size_t n_row{this->shape[0]};
      size_t n_col{this->shape[1]};
      std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
          keys_list{};
      internal::Sorted<true> is_sorted{};
      Key_t pair_type{0, 0};
      // avoid checking the order in pair_type by ensuring it has already been
      // done
      internal::SortedKey<Key_t> spair_type{is_sorted, pair_type};
      keys_list.insert(spair_type);
      // clear the data container and resize it
      soap_vectors.clear();
      soap_vectors.set_shape(n_row, n_col);
      soap_vectors.resize(keys_list);
      soap_vectors.setZero();
      if (this->compute_gradients) {
        soap_vector_gradients.clear();
        soap_vector_gradients.set_shape(ThreeD * n_row, n_col);
        soap_vector_gradients.resize(keys_list);
        soap_vector_gradients.setZero();
      } else {
        soap_vector_gradients.clear();
        soap_vector_gradients.resize();
      }
    } else {
      size_t n_row{this->shape[0]};
      size_t n_col{this->shape[1]};
      this->coeff_indices_map.clear();
      // clear the data container and resize it
      soap_vectors.clear();
      soap_vectors.set_shape(n_row, n_col);

      if (this->compute_gradients) {
        soap_vector_gradients.clear();
        soap_vector_gradients.set_shape(ThreeD * n_row, n_col);
      }

      std::vector<
          std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>>
          keys_list{};
      std::vector<
          std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>>
          keys_list_grad{};

      // identify the species in each environment and initialize soap_vectors
      for (auto center : manager) {
        auto atom_i_tag = center.get_atom_tag();
        auto & coefficients{expansions_coefficients[center]};
        internal::Sorted<true> is_sorted{};

        std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
            pair_list{};
        int center_type{center.get_atom_type()};
        Key_t pair_type{center_type, center_type};

        pair_list.insert({is_sorted, pair_type});
        for (const auto & el1 : coefficients) {
          auto && neigh1_type{el1.first[0]};
          if (center_type <= neigh1_type) {
            pair_type[0] = center_type;
            pair_type[1] = neigh1_type;
          } else {
            pair_type[1] = center_type;
            pair_type[0] = neigh1_type;
          }

          pair_list.insert({is_sorted, pair_type});

          for (const auto & el2 : coefficients) {
            auto && neigh2_type{el2.first[0]};
            if (neigh1_type <= neigh2_type) {
              pair_type[0] = neigh1_type;
              pair_type[1] = neigh2_type;
              pair_list.insert({is_sorted, pair_type});
            }
          }
        }
        keys_list.emplace_back(pair_list);
        if (this->compute_gradients) {
          keys_list_grad.emplace_back(pair_list);
          auto & coef = expansions_coefficients[center];
          // Neighbour gradients need a separate pair list because if the
          // species of j are not the same as either of the species for that
          // SOAP entry, the gradient is zero. since we compute \grad_j p{i ab}
          // we need the species present in the environment of c^{i}
          for (auto neigh : center.pairs()) {
            int neigh_type{neigh.get_atom_type()};
            auto atom_j = neigh.get_atom_j();
            auto atom_j_tag = atom_j.get_atom_tag();
            std::set<internal::SortedKey<Key_t>, internal::CompareSortedKeyLess>
                grad_pair_list{};
            // grad contribution is not zero if the neighbour is _not_ an
            // image of the center
            if (atom_j_tag != atom_i_tag) {
              // list of keys present in the neighbor environment (contains
              // center_type by definition)
              std::vector<Key_t> keys_j{coef.get_keys()};

              for (const auto & neigh_1_type : keys_j) {
                for (const auto & neigh_2_type : keys_j) {
                  if (neigh_1_type[0] <= neigh_2_type[0]) {
                    if ((neigh_type == neigh_1_type[0]) or
                        (neigh_type == neigh_2_type[0])) {
                      pair_type[0] = neigh_1_type[0];
                      pair_type[1] = neigh_2_type[0];
                      grad_pair_list.insert({is_sorted, pair_type});
                    }
                  }
                }
              }
              if (not this->normalize) {
                keys_list_grad.emplace_back(grad_pair_list);
              } else {
                // the normalization: \grad_k \tilde{p}^{i} ~ \tilde{p}^{i}
                // the keys for the gradient are the same as for the
                // representation
                keys_list_grad.emplace_back(pair_list);
              }
            } else {
              keys_list_grad.emplace_back(grad_pair_list);
            }
          }  // auto neigh : center.pairs()
        }    // if compute_gradients

        // initialize coeff_indices_map if the representation is not sparsified
        for (const auto & key : pair_list) {
          if (not this->coeff_indices_map.count(key)) {
            this->coeff_indices_map[key] = coeff_indices;
          }
        }
      }  // for center : manager

      soap_vectors.resize(keys_list);
      soap_vectors.setZero();
      if (this->compute_gradients) {
        soap_vector_gradients.resize(keys_list_grad);
        soap_vector_gradients.setZero();
      } else {
        soap_vector_gradients.resize();
      }
    }
  }

  template <class StructureManager, class Invariants,
            class InvariantsDerivative, class ExpansionCoeff>
  void CalculatorSphericalInvariants::
      initialize_per_center_radialspectrum_soap_vectors(
          Invariants & soap_vectors,
          InvariantsDerivative & soap_vector_gradients,
          ExpansionCoeff & expansions_coefficients,
          std::shared_ptr<StructureManager> manager) {
    size_t n_row{this->max_radial};
    size_t n_col{1};

    soap_vectors.clear();
    soap_vectors.set_shape(n_row, n_col);

    if (this->compute_gradients) {
      soap_vector_gradients.clear();
      soap_vector_gradients.set_shape(ThreeD * n_row, n_col);
    }

    std::vector<std::set<Key_t>> keys_list{};
    std::vector<std::set<Key_t>> keys_list_grad{};
    for (auto center : manager) {
      auto & coefficients{expansions_coefficients[center]};
      Key_t element_type{0};

      std::set<Key_t> keys{};
      for (const auto & el1 : coefficients) {
        keys.insert({el1.first[0]});
      }
      keys.insert({center.get_atom_type()});
      // initialize the radial spectrum to 0 and the proper size
      keys_list.emplace_back(keys);
      keys_list_grad.emplace_back(keys);
      for (auto neigh : center.pairs()) {
        (void)neigh;  // to avoid compiler warning
        keys_list_grad.emplace_back(keys);
      }
    }

    soap_vectors.resize(keys_list);
    soap_vectors.setZero();
    if (this->compute_gradients) {
      soap_vector_gradients.resize(keys_list_grad);
      soap_vector_gradients.setZero();
    } else {
      soap_vector_gradients.resize();
    }
  }
}  // namespace rascal

namespace nlohmann {
  /**
   * Special specialization of the json serialization for non default
   * constructible type.
   */
  template <>
  struct adl_serializer<rascal::CalculatorSphericalInvariants> {
    static rascal::CalculatorSphericalInvariants from_json(const json & j) {
      return rascal::CalculatorSphericalInvariants{j};
    }

    static void to_json(json & j,
                        const rascal::CalculatorSphericalInvariants & t) {
      j = t.hypers;
    }
  };
}  // namespace nlohmann

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
