/**
 * @file   calculator_behler_parrinello_impl.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   10 Sep 2019
 *
 * @brief  Implementation of CalculatorBehlerParrinello
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

#ifndef SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_IMPL_HH_
#define SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_IMPL_HH_

namespace rascal {

  /* ---------------------------------------------------------------------- */
  CalculatorBehlerParrinello::CalculatorBehlerParrinello(
      const Hypers_t & parameters) {
    // Extract the options and hyperparameters
    this->set_hyperparameters(parameters);
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CALCULATOR_BEHLER_PARRINELLO_IMPL_HH_
