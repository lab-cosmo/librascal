/**
 * file   cutoff_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   17 Dec 2018
 *
 * @brief  implementation of cutoff functions for neural nets
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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


#ifndef SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
#define SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
#include <cmath>

namespace rascal {

  /* ---------------------------------------------------------------------- */
  enum class CutoffFunType { Cosine, CosineShifted, Tanh };

  /* ---------------------------------------------------------------------- */
  template <CutoffFunType FunType>
  struct CutoffFun {};

  /* ---------------------------------------------------------------------- */
  template <>
  struct CutoffFun<CutoffFunType::Cosine> {
    static constexpr size_t NbParams{1};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 1>;

    static double eval_function(const Eigen::MatrixBase<ParamShape> & params,
                                const double & r_ij) {
      auto && r_c{params(0)};
      if (r_ij < r_c) {
        return 0.5 * (cos(r_ij * M_PI / r_c) + 1.);
      } else {
        return 0.;
      }
    }

    static double eval_derivative(const Eigen::MatrixBase<ParamShape> & params,
                                  const double & r_ij) {
      auto && r_c{params(0)};
      if (r_ij < r_c) {
        return -0.5 * sin(r_ij * M_PI / r_c);
      } else {
        return 0.;
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <>
  struct CutoffFun<CutoffFunType::CosineShifted> {
    static constexpr size_t NbParams{2};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 1>;

    static double eval_function(const Eigen::MatrixBase<ParamShape> & params,
                                const double & r_ij) {
      auto && r_c{params(0)};
      auto && alpha{params(1)};
      const auto Rstar = alpha * r_c;
      if (r_ij < Rstar) {
        return 1.;
      } else if (r_ij < r_c) {
        0.5 * (cos((r_ij - Rstar) * M_PI / (r_c - Rstar)) + 1.);
      } else {
        return 0.;
      }
    }

    static double eval_derivative(const Eigen::MatrixBase<ParamShape> & params,
                                  const double & alpha, const double & r_ij) {
      auto && r_c{params(0)};
      auto && alpha{params(1)};
      const auto Rstar = alpha * r_c;
      if (r_ij < r_star) {
        return 0.;
      } else if (r_ij < r_c) {
        return -0.5 * sin((r_ij - r_star) * M_PI / (r_c - r_star)) * M_PI /
               (r_c - r_star);
      } else {
        return 0.;
      }
    }
  };

  /* ---------------------------------------------------------------------- */
  template <>
  struct CutoffFun<CutoffFunType::Tanh> {
    static constexpr size_t NbParams{1};
    using ParamShape = Eigen::MatrixBase<double, NbParams, 1>;

    static double eval_function(const Eigen::MatrixBase<ParamShape> & params,
                                const double & r_ij) {
      auto && r_c{params(0)};
      if (r_ij < r_c) {
        return pow(tanh(1. - r_ij / r_c), 3.);
      } else {
        return 0.;
      }
    }

    static double eval_derivative(const Eigen::MatrixBase<ParamShape> & params,
                                  const double & r_ij) {
      auto && r_c{params(0)};
      if (r_ij < r_c) {
        -3. / r_c * pow(sinh(1. - r_ij / r_c), 2.) /
            pow(cosh(1. - r_ij / r_c), 4.);
      } else {
        return 0;
      }
    }
  };

}  // namespace rascal


#endif  // SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
