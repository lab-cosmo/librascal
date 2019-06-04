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
 * Copyright  2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "representations/basis_function_manager.hh"

#include <cmath>

using unit = rascal::BasisFunManager::unit;
using BasisFunType = rascal::BasisFunManager::BasisFunType;

namespace rascal {  // NOLINT
  /* ---------------------------------------------------------------------- */
  enum class BasisFunManager::BasisFunType : int {  // Symmetry Functions
    One = 0;  // yields only cutoff function
    Gaussian = 1;
    cosine = 2;
    angular1 = 3;  // depends on 3 distances
    angular2 = 4;  // depends on 2 distances
  };

  enum class BasisFunManager::CutoffFunType : int {
    cosine = 1; cosine_with_shift = 2;  // as used by Kobayashi et al 2017
    tanh = 3;
  };

  /* ------------------------------------------------------------------ */
  constexpr unint BasisFunManager::get_nhyper(const BasisFunType & fun_type) {
    switch (fun_type) {
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

  /* -------------------------------------------------------------------- */
  template <BasisFunType fun_type>
  inline double BasisFunManager::comp_fun(const double * const param,
                                          const double * const rij) {
    switch (fun_type) {
    case BasisFunType::One: {
      return 1.;
    }
    case BasisFunType::Gaussian: {
      const auto & eta = param[0];  // Gaussian width
      const auto & Rs = param[1];   // position
      const auto & Rij = rij[0];
      const auto & dR = Rij - Rs;
      return exp(-eta * dR * dR);
    }
    case BasisFunType::cosine: {
      const auto & kappa = param[0];
      const auto & Rij = rij[0];
      return cos(kappa * Rij);
    }
    // remark: possible direct use of cosine
    case BasisFunType::angular1: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      const auto & Rjk = rij[3];
      return pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta) *
             exp(-eta * (Rij * Rij + Rik * Rik + Rjk * Rjk));
    }
    // remark: possible direct use of cosine?
    case BasisFunType::angular2: {
      const auto & zeta = param[0];
      const auto & lambda = param[1];
      const auto & eta = param[2];
      const auto & thetaijk = rij[0];
      const auto & Rij = rij[1];
      const auto & Rik = rij[2];
      return pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta) *
             exp(-eta * (Rij * Rij + Rik * Rik));
    }
    default:
      throw std::runtime_error("Basis function not specified")
    }
  }

  /* -------------------------------------------------------------------- */
  template <BasisFunType fun_type, T>
  // Return type should be templated - as input type * input type
  decltype(auto) BasisFunManager::comp_Dfun(const double * const param,
                                            const double * const rij) {
    switch (fun_type) {
    case BasisFunType::One: {
      return 1.;  // need to keep a '1.' -> keep derivative of cutoff
    }
    case BasisFunType::Gaussian: {
      const auto & eta = param[0];  // Gaussian width
      const auto & Rs = param[1];   // position
      const auto & Rij = rij[0];
      const auto & dR = Rij - Rs;
      return -2. * eta * dR * exp(-eta * dR * dR);  // returns double
    }
    case BasisFunType::cosine: {
      const auto & kappa = param[0];
      const auto & Rij = rij[0];
      return -kappa * sin(kappa * Rij);  // returns double
    }
    // case BasisFunType::angular1: {
    //   const auto & zeta = param[0];
    //   const auto & lambda = param[1];
    //   const auto & eta = param[2];
    //   const auto & thetaijk = rij[0];
    //   const auto & Rij = rij[1];
    //   const auto & Rik = rij[2];
    //   const auto & Rjk = rij[3];
    //   return 3 * pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
    //        * exp(-eta * (Rij*Rij + Rik*Rik + Rjk*Rjk));
    // }
    // case BasisFunType::angular2: {
    //   const auto & zeta = param[0];
    //   const auto & lambda = param[1];
    //   const auto & eta = param[2];
    //   const auto & thetaijk = rij[0];
    //   const auto & Rij = rij[1];
    //   const auto & Rik = rij[2];
    //   return  // returns double * 2
    // /* pow(2., 1 - zeta) * pow(1. + lambda * cos(thetaijk), zeta)
    //      * exp(-eta * (Rij*Rij + Rik*Rik)) */;
    // }
    default:
      throw std::runtime_error("Basis function not specified");
    }
  }

  /* -------------------------------------------------------------------- */
  template <CutoffFunType cfun_type>
  inline double BasisFunManager::comp_fc(const double * const param,
                                         const double * const rij) {
    switch (cfun_type) {
    case CutoffFunType::cosine: {
      const auto & Rc = param[0];
      if (rij < Rc) {
        return 0.5 * (cos(rij * M_PI / Rc) + 1.);

      } else {
        return 0.;
      }
    }
    case CutoffFunType::cosine_with_shift: {
      const auto & Rc = param[0];
      const auto & alpha = param[1];
      const auto Rstar = alpha * Rc;
      if (rij < Rstar) {
        return 1.;
      } else if (rij < Rc) {
        0.5 * (cos((rij - Rstar) * M_PI / (Rc - Rstar)) + 1.);
      } else {
        return 0.;
      }
    }
    case CutoffFunType::tanh: {
      const auto & Rc = param[0];
      if (rij < Rc) {
        return pow(tanh(1. - rij / Rc), 3.);
      } else {
        return 0.;
      }
    }
    }
  }

  /* -------------------------------------------------------------------- */
  template <CutoffFunType cfun_type>
  inline double
  BasisFunManager::comp_Dfc(const double * const param,  // Sigmoid functions
                            const double * const rij) {
    switch (cfun_type) {
    case CutoffFunType::cosine: {
      const auto & Rc = param[0];
      if (rij < Rc) {
        return -0.5 * sin(rij * M_PI / Rc);
      } else {
        return 0.;
      }
    }
    case CutoffFunType::cosine_with_shift: {
      const auto & Rc = param[0];
      const auto & alpha = param[1];
      const auto Rstar = alpha * Rc;
      if (rij < Rstar) {
        return 0.;
      } else if (rij < Rc) {
        return -0.5 * sin((rij - Rstar) * M_PI / (Rc - Rstar)) * M_PI /
               (Rc - Rstar);
      } else {
        return 0.;
      }
    }
    case CutoffFunType::tanh: {
      const auto & Rc = param[0];
      if (rij < Rc) {
        -3. / Rc * pow(sinh(1. - rij / Rc), 2.) / pow(cosh(1. - rij / Rc), 4.);
      } else {
        return 0;
      }
    }
    }
  }
}  // namespace rascal
