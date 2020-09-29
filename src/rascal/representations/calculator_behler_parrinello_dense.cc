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

#include <map>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  template <bool CompatibilityMode_, SymmetryFunctionType... SymFunTypes>
  CalculatorBehlerParrinelloDense<CompatibilityMode_, SymFunTypes...>::
      CalculatorBehlerParrinelloDense(const Hypers_t & parameters,
                                      const UnitStyle & unit_style)

        : CalculatorBase{}, cutoff_function_type{json_io::get<std::string>(
                              parameters, "cutoff_function_type")} {
    this->set_default_prefix("dense_behler_parrinello_");
    // simple check (just existence of keys)
    this->check_hyperparameters(this->reference_hypers, parameters);
    // 1)  create shared_ptr to cutoff_function
    std::string cut_fun_type{json_io::get<std::string>(
                              parameters, "cutoff_function_type")};

    // 2) iterate through sym_function_params and fill vector of
    // behlerfeatures
    std::map<std::string, std::shared_ptr<CutoffFunctionBase>> cut_funs;
    for (auto && sym_fun_params :
         json_io::get(parameters, "symmetry_functions")) {
      auto && cut_fun_params{json_io::get(sym_fun_params, "cutoff_function")};
      auto && cut_fun{cut_funs[CutoffFunctionBase::identifier(
          unit_style, cut_fun_params)]};

      if (cut_fun == nullptr) {
        cut_fun = CutoffFunctionBase::make_shared(unit_style, cut_fun_params);
      }
      this->behler_features.push_back(std::move(
          BehlerFeatureBase<CompatibilityMode_, SymFunTypes...>::make_unique(
              cut_fun, unit_style, sym_fun_params)));
    }
  }

  /* ---------------------------------------------------------------------- */
  template class CalculatorBehlerParrinelloDense<
      true, SymmetryFunctionType::Gaussian, SymmetryFunctionType::AngularNarrow,
      SymmetryFunctionType::AngularWide>;

  template class CalculatorBehlerParrinelloDense<
      true, SymmetryFunctionType::Gaussian,
      SymmetryFunctionType::AngularNarrow>;

}  // namespace rascal
