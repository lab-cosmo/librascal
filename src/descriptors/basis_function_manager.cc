/**
 * file   basis_function_manager.cc
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   08 May 2018
 *
 * @brief  manager for basis functions, which are the input of the
 * neural network. e.g. symmetry functions or spherical harmonics
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "descriptors/basis_function_manager.cc"

using unit = BasisFunManager::unit;
using BasisFunType = BasisFunManager::BasisFunType;

namespace rascal {
  /* ---------------------------------------------------------------------- */
  enum class BasisFunManager::BasisFunType: int
  {
    One=0; // yields only cutoff function
    Gaussian=1;
    cosine=2;
    angular1=3; // depends on 3 distances
    angular2=4; // depends on 2 distances

  }
  enum class BasisFunManager::CutoffFunType: int
  {
    cosine=1;
    cosine_with_shift=2; // as used by Kobayashi et al 2017
    tanh=3;
  }

  constexpr unint BasisFunManager::get_nhyper(const BasisFunType& fun_type) {
    switch(fun_type) {
    case BasisFunType::One:
      return 0;
    case BasisFunType::Gaussian:
      return 2;
    case BasisFunType::cosine:
      return 1;
    case BasisFunType::angular1:
      return 3;
    case BasisFunType::angular2:
      return 3;
    default:
      return 0;
    }
  }

}
