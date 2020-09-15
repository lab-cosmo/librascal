/**
 * @file   calculator_behler_parrinello_dense.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   13 Sep 2019
 *
 * @brief  Implementation of CalculatorBehlerParrinelloDense
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

#include "calculator_behler_parrinello_dense.hh"

#include "behler_feature.hh"
#include "cutoff_functions_inlineable.hh"

namespace rascal {

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  CalculatorBehlerParrinelloDense<CompatibilityMode_, SymFunTypes...>::
      CalculatorBehlerParrinelloDense(const Hypers_t & parameters)
      : CalculatorBase{} {
    this->set_default_prefix("dense_behler_parrinello_");
    // simple check (just existence of keys)
    this->check_hyperparameters(this->reference_hypers, parameters);
    // true parameter checks
    auto unit_style{UnitStyle::make(parameters.at("unit_style"))};
    // 1)  create shared_ptr to cutoff_function
    auto cut_fun{CutoffFunctionBase::make_shared(
        unit_style, parameters.at("cutoff_function"))};

    // 2) iterate through sym_function_params and fill vector of
    // behlerfeatures
    for (auto && sym_fun_params : parameters.at("symmetry_functions")) {
      this->behler_features.push_back(std::move(
          BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>::make_unique(
              cut_fun, unit_style, sym_fun_params)));
    }
  }

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  template <class StructureManager>
  void
  CalculatorBehlerParrinelloDense<CompatibilityMode_, SymFunTypes...>::compute(
      StructureManager & manager) {
    for (auto && behler_feature : this->behler_features) {
      behler_feature->compute(manager);
    }
  }

  /* ---------------------------------------------------------------------- */
  template class CalculatorBehlerParrinelloDense<
      true, SymmetryFunctionType::Gaussian, SymmetryFunctionType::AngularNarrow,
      SymmetryFunctionType::AngularWide>;

}  // namespace rascal
