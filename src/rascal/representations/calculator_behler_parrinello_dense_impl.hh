/**
 * @file   calculator_behler_parrinello_dense_impl.hh
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

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_IMPL_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_IMPL_HH_

namespace rascal {

  /* ---------------------------------------------------------------------- */
  CalculatorBehlerParrinelloDense::CalculatorBehlerParrinelloDense(
      const Hypers_t & parameters)
      : CalculatorBase{} {
    this->set_default_prefix("dense_behler_parrinello_");
    // simple check (just existence of keys)
    this->check_hyperparameters(this->reference_hypers, parameters);
    // true parameter checks
  }

  /* ---------------------------------------------------------------------- */
  template <class StructureManager>
  void CalculatorBehlerParrinelloDense::compute(StructureManager & manager) {
    for (auto && behler_feature : this->behler_features) {
      behler_feature->compute(manager);
    }
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_DENSE_IMPL_HH_
