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

#include <cmath>

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

  /* ---------------------------------------------------------------------- */
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

  /* ---------------------------------------------------------------------- */
  template<BasisFunType fun_type>
  inline double BasisFunManager::comp_fun(const double * const param,
					  const double * rij) {
    switch (fun_type) {
    case BasisFunType::One: {
      return 1.;
    }
    case BasisFunType::Gaussian: {
      const auto & eta = param[0]; // Gaussian width
      const auto & Rs = param[1];    // position
      const auto & Rij = rij[0];
      const auto & dR = Rij - Rs;
      return exp(-eta * dR * dR);

    }
    case BasisFunType::cosine: {
      const auto & kappa = param [0];
      const auto & Rij = rij[0];
      return cos(kappa * Rij);
    }
      // remark: possible direct use of cos
    case BasisFunType::angular1: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      const auto & Rjk = rij[3];
      return pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
	* exp(-eta * (Rij*Rij + Rik*Rik + Rjk*Rjk));
    }
      // remark: possible direct use of cos
    case BasisFunType::angular2: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      return pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
	* exp(-eta * (Rij*Rij + Rik*Rik));
    }
    default:
      throw std::runtime_error("Basis function not specified")
    }
  }

  template<BasisFunType fun_type>
  inline double BasisFunManager::comp_Dfun(const double * const param,
					   const double * rij) {
    switch (fun_type) {
    case BasisFunType::One: {
      return 1.; // need to keep a '1.' -> keep derivative of cutoff
    }
    case BasisFunType::Gaussian: {
      const auto & eta = param[0]; // Gaussian width
      const auto & Rs = param[1];    // position
      const auto & Rij = rij[0];
      const auto & dR = Rij - Rs;
      return -2. * eta * dR * exp(-eta * dR * dR);

    }
    case BasisFunType::cosine: {
      const auto & kappa = param [0];
      const auto & Rij = rij[0];
      return -kappa * sin(kappa * Rij);
    }
    case BasisFunType::angular1: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      const auto & Rjk = rij[3];
      return /* pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
	      * exp(-eta * (Rij*Rij + Rik*Rik + Rjk*Rjk)); */
    }
    case BasisFunType::angular2: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      return /* pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
	      * exp(-eta * (Rij*Rij + Rik*Rik)) */;
    }
    default:
      throw std::runtime_error("Basis function not specified")
    }
  }


}
