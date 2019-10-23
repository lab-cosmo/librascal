/**
 * @file   bessel.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   27 May 2019
 *
 * @brief Implementation of the modified spherical bessel of the 1st kind
 *
 * Copyright  2019  Felix Musil, Max Veit COSMO (EPFL), LAMMM (EPFL)
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

#include "rascal/math/bessel.hh"
using namespace rascal::math;

void ModifiedSphericalBessel::precompute(
    size_t l_max, const Eigen::Ref<const Eigen::VectorXd> & x_v) {
  this->x_v = x_v.array();
  this->n_max = x_v.size();
  this->bessel_values.resize(this->n_max, l_max + 1);
  this->bessel_arg.resize(this->n_max);
  this->bessel_arg_i.resize(this->n_max);
  this->exp_bessel_arg.resize(this->n_max);
  this->bessel_arg_pow.resize(this->n_max);
  this->efac.resize(this->n_max);
  this->l_max = l_max;
  this->order_max = static_cast<int>(this->l_max + 1);
  this->igammas.resize(2);
  int ii{0};
  for (int order{this->order_max - 2}; order < this->order_max; ++order) {
    this->hyp1f1s.emplace_back(static_cast<double>(order + 1),
                               static_cast<double>(2 * order + 2));
    this->igammas[ii] = 1. / std::tgamma(1.5 + order);
    ii++;
  }
}

void ModifiedSphericalBessel::upward_recursion(double distance, double fac_a,
                                               int n_rows) {
  auto vals = this->bessel_values.bottomRows(n_rows);
  // i_0(z) = sinh(z) / z
  vals.col(0) = (Eigen::exp(-fac_a * (x_v.tail(n_rows) - distance).square()) -
                 Eigen::exp(-fac_a * (x_v.tail(n_rows) + distance).square())) *
                0.5 * this->bessel_arg_i.tail(n_rows);
  // i_1(z) = cosh(z)/z - i_0(z)/z
  vals.col(1) = ((Eigen::exp(-fac_a * (x_v.tail(n_rows) - distance).square()) +
                  Eigen::exp(-fac_a * (x_v.tail(n_rows) + distance).square())) *
                 0.5 * this->bessel_arg_i.tail(n_rows)) -
                vals.col(0) * this->bessel_arg_i.tail(n_rows);

  for (int order{2}; order < this->order_max; ++order) {
    vals.col(order) = vals.col(order - 2) - vals.col(order - 1) *
                                                (2. * order - 1.) *
                                                this->bessel_arg_i.tail(n_rows);
  }
}

void ModifiedSphericalBessel::downward_recursion(double distance, double fac_a,
                                                 int n_rows) {
  auto vals = this->bessel_values.topRows(n_rows);
  this->exp_bessel_arg = Eigen::exp(-this->bessel_arg.head(n_rows));
  this->efac = std::exp(-fac_a * distance * distance) *
               Eigen::exp(-fac_a * this->x_v.head(n_rows).square());
  for (int i_order{0}; i_order < 2; ++i_order) {
    int order{this->order_max - 2 + i_order};
    auto & hyp1f1{this->hyp1f1s[i_order]};
    for (int ii{0}; ii < n_rows; ++ii) {
      vals(ii, order) = this->exp_bessel_arg[ii] * this->igammas[i_order] *
                        math::pow(this->bessel_arg[ii] * 0.5, order) * 0.5 *
                        math::SQRT_PI * hyp1f1.calc(2. * this->bessel_arg[ii]);
    }
    vals.col(order) *= this->efac;
  }

  for (int order{this->order_max - 3}; order >= 0; --order) {
    vals.col(order) = vals.col(order + 2) + vals.col(order + 1) *
                                                (2. * order + 3.) *
                                                this->bessel_arg_i.head(n_rows);
  }
}

void ModifiedSphericalBessel::calc(double distance, double fac_a) {
  this->bessel_arg = (2. * fac_a * distance) * this->x_v;
  this->bessel_arg_i = this->bessel_arg.inverse();

  if (this->l_max == 0) {
    // recursions are not valid for l_max==0 so direct computation
    // i_0(z) = sinh(z) / z
    this->bessel_values.col(0) =
        (Eigen::exp(-fac_a * (x_v - distance).square()) -
         Eigen::exp(-fac_a * (x_v + distance).square())) *
        0.5 * this->bessel_arg_i;
  } else {
    // for l_max > 0 downward/upward_recursion functions are applicable
    // find the index where bessel_arg is larger than 50
    // (bessel_arg is sorted by increasing order)
    int n_down{0};
    for (; n_down < this->n_max; ++n_down) {
      if (this->bessel_arg[n_down] > 50) {
        ++n_down;
        break;
      }
    }
    // apply downward recurence where bessel_arg < 50
    if (n_down > 0) {
      this->downward_recursion(distance, fac_a, n_down);
    }

    // apply upward recurence where bessel_arg > 50
    int n_up{this->n_max - n_down};
    if (n_up > 0) {
      this->upward_recursion(distance, fac_a, n_up);
    }
  }

  // Set small values to 0 because the recursion looses accuracy for very
  // small values. Also on the python side it avoids some unexpected
  // interpretation of values that are strictly speaking outside of the
  // range of double precision
  bessel_values = bessel_values.unaryExpr([](double d) {
    if (d < 1e-100) {
      return 0.;
    } else {
      return d;
    }
  });
}
