/**
 * @file   calculator_behler_parrinello_dense.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   12 Sep 2019
 *
 * @brief  Implementation of a dense (i.e. non-sparse) Behler-Parrinello-type
 * descriptor calculater (e.g., as would be normal in a single-species
 * calculation)
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

#include "calculator_base.hh"
#include "symmetry_functions.hh"

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_

namespace rascal {
  //! forward declaration
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  class BehlerFeatureBase;

  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  class CalculatorBehlerParrinelloDense : public CalculatorBase {
   public:
    using Parent = CalculatorBase;
    using Hypers_t = typename Parent::Hypers_t;
    using BehlerFeatureBase_t =
        BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>;

    // type of the datastructure used to register the list of valid
    // hyperparameters
    using ReferenceHypers_t = Parent::ReferenceHypers_t;

    //! Default constructor
    CalculatorBehlerParrinelloDense() = delete;

    //! Constructor with input json
    explicit CalculatorBehlerParrinelloDense(const Hypers_t & parameters);

    //! Copy constructor
    CalculatorBehlerParrinelloDense(
        const CalculatorBehlerParrinelloDense & other) = delete;

    //! Move constructor
    CalculatorBehlerParrinelloDense(CalculatorBehlerParrinelloDense && other) =
        default;

    //! Destructor
    virtual ~CalculatorBehlerParrinelloDense() = default;

    //! Copy assignment operator
    CalculatorBehlerParrinelloDense &
    operator=(const CalculatorBehlerParrinelloDense & other) = delete;

    //! Move assignment operator
    CalculatorBehlerParrinelloDense &
    operator=(CalculatorBehlerParrinelloDense && other) = default;

    template <class StructureManager>
    inline void compute(StructureManager & manager);

   protected:
    /**
     * stores base class refs to the Symmetry functions to be evaluated. The
     * main loop iterates over this container and applies each function to
     * all clusters of an input structure manager, storing the results to a
     * provided input property (which is not necessarily attached to the
     * same structure manager)
     */
    std::vector<std::unique_ptr<BehlerFeatureBase_t>> behler_features{};
    /**
     * reference of the requiered hypers, these are used to check the parameters
     * for building the behler_features vector at construction
     */
    const ReferenceHypers_t reference_hypers{{"unit_style", {}},
                                             {"scaling", {}},
                                             {"cutoff_function", {}},
                                             {"symmetry_functions", {}}};
  };
}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_HH_
