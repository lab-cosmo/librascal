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
    auto unit_style{[&parameters]() {
      auto unit_style_label{parameters.at("unit_style").get<std::string>()};
      if (unit_style_label == "metal") {
        return units::metal;
      } else {
        throw std::runtime_error("unable to handle unit style '" +
                                 unit_style_label + "'.");
      }
    }()};
    // 1)  create shared_ptr to cutoff_function
    auto cut_fun{[&parameters, &unit_style]() {
      auto cut_fun_label{
          parameters.at("cutoff_function").at(name).get<std::string>()};
      if (cut_fun_label == "Cosine") {
        return std::make_shared<CutoffFunction<InlCutoffFunctionType::Cosine>>(
            unit_style, parameters.at("cutoff_function"));
      } else if (false /*there is no other cutoff function for the moment*/) {
      } else {
        throw std::runtime_error("unable to handle "
                                 "cutoff_function style '" +
                                 cut_fun_label + "'.");
      }
    }()};

    // 2) iterate through sym_function_params and fill vector of
    // behlerfeatures
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
