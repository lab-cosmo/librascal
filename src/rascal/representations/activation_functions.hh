/**
 * file   activation_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   14 Jan 2019
 *
 * @brief  implementation of activation functions for neural nets
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with librascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_REPRESENTATIONS_ACTIVATION_FUNCTIONS_HH_
#define SRC_RASCAL_REPRESENTATIONS_ACTIVATION_FUNCTIONS_HH_

#include <cmath>

#include "Eigen/Dense"

namespace rascal {

  enum ActivationFunType { Identity, Logistic, CentredLogistic };

  /* ---------------------------------------------------------------------- */
  template <ActivationFunType FunType>
  struct ActivationFun {};

  /* ---------------------------------------------------------------------- */
  template <>
  struct ActivationFun<ActivationFunType::Identity> {
    static constexpr size_t NbParams{0};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 0>;

    static Arr_t eval_function(const Eigen::ArrayBase<Derived> & x) {
      return x;
    }
    static Arr_t eval_derivative(const Eigen::ArrayBase<Derived> & x) {
      return 1.;
    }
  }

  /* ---------------------------------------------------------------------- */
  template <>
  struct ActivationFun<ActivationFunType::Logistic> {
    static constexpr size_t NbParams{0};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 0>;

    static Arr_t eval_function(const Eigen::ArrayBase<Derived> & x) {
      return 1. / (1. + exp(-x));
    }
    static Arr_t eval_derivative(const Eigen::ArrayBase<Derived> & x) {
      auto && f_val{eval_function(x)};
      return f_val * (1. - f_val);
    }
  }

  /* ---------------------------------------------------------------------- */
  template <>
  struct ActivationFun<ActivationFunType::CentredLogistic>
      : public ActivationFun<ActivationFunType::Logistic> {
    static constexpr size_t NbParams{0};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 0>;
    using Parent = ActivationFun<ActivationFunType::Logistic>;

    template <class Derived>
    static Arr_t eval_function(const Eigen::ArrayBase<Derived> & x) {
      return Parent::eval_function(x) - .5;
    }
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_ACTIVATION_FUNCTIONS_HH_
