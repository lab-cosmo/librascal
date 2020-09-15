
/**
 * @file   cutoff_functions_inlineable.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   09 Dec 2019
 *
 * @brief  bundle inlineable cutoff functions
 *
 * Copyright © 2019 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

#include "cutoff_functions_inlineable.hh"

namespace rascal {

  /* ---------------------------------------------------------------------- */
  std::shared_ptr<CutoffFunctionBase>
  CutoffFunctionBase::make_shared(const units::UnitStyle & unit_style,
                                  Hypers_t hypers) {
    auto cut_fun_label{hypers.at("name").get<std::string>()};
    if (cut_fun_label == "Cosine") {
      return std::make_shared<CutoffFunction<InlCutoffFunctionType::Cosine>>(
          unit_style, hypers);
    } else if (false /*there is no other cutoff function for the moment*/) {
    } else {
      throw std::runtime_error("unable to handle "
                               "cutoff_function style '" +
                               cut_fun_label + "'.");
    }
  }

}  // namespace rascal
