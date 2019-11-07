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
// #include "representations/calculator_spherical_invariants.hh"
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

  class CalculatorSphericalInvariants : public CalculatorBase {
   public:
    using Hypers_t = typename CalculatorBase::Hypers_t;
    using Key_t = typename CalculatorBase::Key_t;

    template <class StructureManager>
    using Property_t =
        BlockSparseProperty<double, 2, 0, StructureManager, Key_t>;

    template <class StructureManager>
    using PropertyGradient_t =
        BlockSparseProperty<double, 2, 0, StructureManager, Key_t>;

    template <class StructureManager>
    using Dense_t = typename Property_t<StructureManager>::Dense_t;

    template <class StructureManager>
    using Data_t = typename Property_t<StructureManager>::Data_t;

    template <class StructureManager, size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    explicit CalculatorPairDistances(const Hypers_t & hyper)
        : CalculatorBase{}, rep_expansion{hyper} {
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
      
      
	  // Implement scaling function for power laws here


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
    template <class StructureManager,
        std::enable_if_t<internal::is_proper_iterator<StructureManager>::value,
                         int> = 0>
    void compute_loop(StructureManager & managers) {
      for (auto & manager : managers) {
        this->compute_impl<BodyOrder>(manager);
      }
    }

    //! single manager case
    template <class StructureManager,
        std::enable_if_t<
            not(internal::is_proper_iterator<StructureManager>::value), int> =
            0>
    void compute_loop(StructureManager & manager) {
      this->compute_impl<BodyOrder>(manager);
    }

    //! compute representation @f$ \nu == 1 @f$
    template <class StructureManager>
    void compute_impl(std::shared_ptr<StructureManager> manager);

   protected:
    double interaction_cutoff{};
    double cutoff_smooth_width{};
    bool compute_gradients{};
    bool cutoff_per_pair_type{false};
    
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
    using PropGrad_t = PropertyGradient_t<StructureManager>;
    constexpr static int n_spatial_dimensions = StructureManager::dim();
    using math::pow;

    auto && pair_distances{
        manager->template get_property_ref<Prop_t>(this->get_name())};

    auto && pair_distance_gradients{
        manager->template get_property_ref<PropGrad_t>(
            this->get_gradient_name())};

    // if the representation has already been computed for the current
    // structure then do nothing
    if (pair_distances.is_updated()) {
      return;
    }

    Key_t pair_type{0, 0};
    // use special container to tell that there is not need to sort when
    // using operator[] of soap_vector
    internal::SortedKey<Key_t> spair_type{pair_type};
	
	size_t number_of_pairs = 0;
	for (auto center : manager){
	  number_of_pairs += center.get_size();
	}
	
	pair_distances.resize(number_of_pairs);
	
    for (auto center : manager) {
	  for (auto neigh : center) {
		  pair_distances[spair_type][neigh] = neigh.get_distance();
      // auto & soap_vector{soap_vectors[center]};
      /*for (const auto & el1 : coefficients) {
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
          

        }  // for el1 : coefficients
      }    // for el2 : coefficients
      */
      


      if (this->compute_gradients) {
		// Compute Gradients
      }        // if compute gradients
    }          // for center : manager
  }            // compute_powerspectrum()

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_SPHERICAL_INVARIANTS_HH_
