/**
 * @file   calculator_pair_distances.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Sahil Shah <sahil.shah@epfl.ch>
 *
 * @date   13 Feb 2019
 *
 * @brief  compute pair distances
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

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_PAIR_DISTANCES_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_PAIR_DISTANCES_HH_

#include "rascal/math/utils.hh"
#include "rascal/utils/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/representations/cutoff_functions.hh"
#include "rascal/structure_managers/property_block_sparse.hh"
#include "rascal/structure_managers/structure_manager.hh"

#include <wigxjpf.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <exception>
#include <unordered_set>
#include <vector>

namespace rascal {

  class CalculatorPairDistances : public CalculatorBase {
   public:
    using Hypers_t = typename CalculatorBase::Hypers_t;
    using Key_t = typename CalculatorBase::Key_t;

    template <class StructureManager>
    using Property_t =
        BlockSparseProperty<double, 2, StructureManager, Key_t>;

    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, StructureManager, Key_t>;

    template <class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;

    template <class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    explicit CalculatorPairDistances(const Hypers_t & hyper)
        : CalculatorBase{} {
      this->set_default_prefix("pair_distances_");
      this->set_hyperparameters(hyper);
    }

    //! Copy constructor
    CalculatorPairDistances(const CalculatorPairDistances & other) =
        delete;

    //! Move constructor
    CalculatorPairDistances(CalculatorPairDistances && other) =
        default;

    //! Destructor
    virtual ~CalculatorPairDistances() = default;

    //! Copy assignment operator
    CalculatorPairDistances &
    operator=(const CalculatorPairDistances & other) = delete;

    //! Move assignment operator
    CalculatorPairDistances &
    operator=(CalculatorPairDistances && other) = default;

    void set_hyperparameters(const Hypers_t & hypers) {
      using internal::CutoffFunctionType;

      if (hypers.find("cutoff_per_pair_type") != hypers.end()) {
          this->cutoff_per_pair_type = hypers.at("cutoff_per_pair_type").get<bool>();
      }
      if (this->cutoff_per_pair_type) {
          throw std::logic_error("Cutoff per pair type is not yet supported" );
      }
      auto fc_hypers = hypers.at("cutoff_function").get<json>();
      auto fc_type = fc_hypers.at("type").get<std::string>();
      this->interaction_cutoff = fc_hypers.at("cutoff").at("value");
      this->cutoff_smooth_width = fc_hypers.at("smooth_width").at("value");
      if (fc_type == "ShiftedCosine") {
        this->cutoff_function_type = CutoffFunctionType::ShiftedCosine;
        this->cutoff_function =
            make_cutoff_function<CutoffFunctionType::ShiftedCosine>(fc_hypers);
      } else if (fc_type == "RadialScaling") {
        this->cutoff_function_type = CutoffFunctionType::RadialScaling;
        this->cutoff_function =
            make_cutoff_function<CutoffFunctionType::RadialScaling>(fc_hypers);
      } else {
        throw std::logic_error("Requested cutoff function type \'" + fc_type +
                               "\' has not been implemented.  Must be one of" +
                               ": \'ShiftedCosine\' or 'RadialScaling'.");
      }

      if (hypers.find("compute_gradients") != hypers.end()) {
        this->compute_gradients = hypers.at("compute_gradients").get<bool>();
      } else {  // Default false (don't compute gradients)
        this->compute_gradients = false;
      }

      // TODO(max) rethink this in the context of more general rescaling
      //           functions
      // Also: Round power to int if close? Or only allow integer powers?
      if (hypers.find("distance_powers") != hypers.end()) {
        this->distance_powers = hypers.at("distance_powers").get<
          std::vector<double>>();
      }
      this->descriptor_dimension = this->distance_powers.size();
      this->set_name(hypers);
    }

    /**
     * Compute representation for a given structure manager.
     *
     * @tparam StructureManager a (single or collection)
     * of structure manager(s) (in an iterator) held in shared_ptr
     *
     */
    template <class StructureManager>
    void compute(StructureManager & managers);

    /**
     * loop over a collection of manangers if it is an iterator.
     * Or just call compute_impl
     */
    template <class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl(manager);
      }
    }

    //! single manager case
    template <class StructureManager,
        std::enable_if_t<
            not(internal::is_proper_iterator<StructureManager>::value), int> =
            0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl(manager);
    }

    //! compute representation @f$ \nu == 1 @f$
    template <class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

   protected:
    double interaction_cutoff{};
    double cutoff_smooth_width{};
    bool compute_gradients{};
    bool cutoff_per_pair_type{false};
    std::vector<double> distance_powers{1.0,};
    int descriptor_dimension{1};

    std::shared_ptr<internal::CutoffFunctionBase> cutoff_function{};
    internal::CutoffFunctionType cutoff_function_type{};
  };

  template <class StructureManager>
  void CalculatorPairDistances::compute(StructureManager & managers) {
    this->compute_loop(managers);
  }

  template <class StructureManager>
  void CalculatorPairDistances::compute_impl(
      std::shared_ptr<StructureManager> manager) {
    using Prop_t = Property_t<StructureManager>;
    // using PropGrad_t = PropertyGradient_t<StructureManager>;
    // constexpr static int n_spatial_dimensions = StructureManager::dim();
    using math::pow;

    auto && pair_distances{
        *manager->template get_property<Prop_t>(this->get_name(), true, true)};

    /* Gradients not yet implemented
     * auto && pair_distance_gradients{
        manager->template get_property_ptr<PropGrad_t>(
            this->get_gradient_name())};
    */

    // if the representation has already been computed for the current
    // structure then do nothing
    if (pair_distances.is_updated()) {
      return;
    }

    pair_distances.clear();
    pair_distances.set_shape(this->descriptor_dimension, 1);

    std::vector<std::set<Key_t>> keys_list{};
    for (auto center : manager) {
      for (auto neigh : center.pairs()) {
        Key_t pair_type{0, 0};
        if (center.get_atom_type() <= neigh.get_atom_type()) {
            pair_type[0] = center.get_atom_type();
            pair_type[1] =  neigh.get_atom_type();
        } else {
            pair_type[1] = center.get_atom_type();
            pair_type[0] = neigh.get_atom_type();
        }
        std::set<Key_t> keys{{pair_type}};
        keys_list.emplace_back(keys);
      }
    }
    pair_distances.resize(keys_list);

    for (auto center : manager) {
      for (auto neigh : center.pairs()) {
        Key_t pair_type{0, 0};
        if (center.get_atom_type() <= neigh.get_atom_type()) {
            pair_type[0] = center.get_atom_type();
            pair_type[1] =  neigh.get_atom_type();
        } else {
            pair_type[1] = center.get_atom_type();
            pair_type[0] = neigh.get_atom_type();
        }
        internal::SortedKey<Key_t> spair_type{pair_type};
        std::vector<internal::SortedKey<Key_t>> pair_list{spair_type};
        // pair_distances[neigh].resize(pair_list);
        for (int pow_idx{0}; pow_idx < this->descriptor_dimension; ++pow_idx) {
          pair_distances[neigh][spair_type](pow_idx) = std::pow(
              manager->get_distance(neigh),
              this->distance_powers[pow_idx]);
        }
        // TODO(max) need to store the cutoff function values separately
      }        // for neigh : center
    }          // for center : manager
  }            // compute_impl()

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_PAIR_DISTANCES_HH_
