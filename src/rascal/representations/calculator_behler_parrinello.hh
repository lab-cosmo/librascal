/**
 * @file   calculator_behler_parrinello.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   10 Sep 2019
 *
 * @brief  Behler-Parrinello implementation
 *
 * Copyright © 2019 Till Junge
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_HH_

#include "rascal/representations/behler_feature.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/representations/cutoff_functions.hh"
#include "rascal/structure_managers/adaptor_filter.hh"
#include "rascal/structure_managers/property.hh"
#include "rascal/structure_managers/species_manager.hh"
#include "rascal/utils/tuple_standardisation.hh"

namespace rascal {

  class CalculatorBehlerParrinello : public CalculatorBase {
   public:
    using Parent = CalculatorBase;
    using Hypers_t = typename Parent::Hypers_t;
    // type of the data structure for the representation feaures
    // currently only here to make it compile
    template <class StructureManager>
    using Property_t =
        Property<double, 1, 1, StructureManager, Eigen::Dynamic, 1>;
    // type of the datastructure used to register the list of valid
    // hyperparameters
    using ReferenceHypers_t = Parent::ReferenceHypers_t;

    //! Default constructor
    CalculatorBehlerParrinello() = delete;

    //! Construct for a input parameter json
    inline explicit CalculatorBehlerParrinello(const Hypers_t & parameters);

    //! Copy constructor
    CalculatorBehlerParrinello(const CalculatorBehlerParrinello & other) = delete;

    //! Move constructor
    CalculatorBehlerParrinello(CalculatorBehlerParrinello && other) = default;

    //! Destructor
    virtual ~CalculatorBehlerParrinello() = default;

    //! Copy assignment operator
    CalculatorBehlerParrinello &
    operator=(const CalculatorBehlerParrinello & other) = delete;

    //! Move assignment operator
    CalculatorBehlerParrinello &
    operator=(CalculatorBehlerParrinello && other) = default;

    //! Pure Virtual Function to set hyperparameters of the representation
    void set_hyperparameters(const Hypers_t & hyper) { this->set_name(hyper); }

    template <class StructureManager>
    inline void compute(std::shared_ptr<StructureManager> manager) {
      // temporary meaningless definition
      manager->template get_property_ref<Property_t<StructureManager>>(
          this->get_name());
    }

   protected:
    //! unique cutoff function used for all input nodes
    internal::CutoffFunctionType cutoff_fun{};
    //! set of all cutoff values for optimisation
    std::set<double> cutoffs{};
    //! reference the requiered hypers
    const ReferenceHypers_t reference_hypers{
        {"bla", {}},
        {"bla list", {"it", "em"}},
    };
  };
}  // namespace rascal

#include "calculator_behler_parrinello_impl.hh"
#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_HH_
