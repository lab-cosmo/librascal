/**
 * @file   rascal/math/hyp1f1.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   11 May 2019
 *
 * @brief Implementation of the confluent hypergeometric function
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "rascal/math/hyp1f1.hh"

using namespace rascal::math;            // NOLINT
using namespace rascal::math::internal;  // NOLINT

/**
 * Utilities to apply the recurrence relations of 1F1 from maximal
 * values of a and b to minimal values.
 * The recurrence path is chosen so that derivatives of 1F1 are also
 * a result.
 *
 * G refers to the modified 1F1 defined below
 */
static double recurence_G_to_val_downward(double a, double b, double z,
                                          double M1p2p, double M1p1p) {
  return (z * (a - b) * M1p2p + M1p1p * b) / a;
}

static double recurence_G_to_der_downward(double, double b, double z,
                                          double M2p3p, double M1p2p) {
  return z * M2p3p + M1p2p * (b + 1);
}

Hyp1f1Series::Hyp1f1Series(double a, double b, size_t mmax, double tolerance)
    : a{a}, b{b}, mmax{mmax}, prefac{std::tgamma(a) / std::tgamma(b)},
      tolerance{tolerance} {
  // when a == b, 1F1 is an exponential
  if (std::abs(1 - this->b / this->a) < DBL_FTOL) {
    this->is_exp = true;
  } else {
    coeff.resize(mmax);
    coeff_derivative.resize(mmax);
    double u{0.};
    // precomputes the (a)_i / (b)_i / i! coefficients
    coeff(0) = a / b;
    coeff_derivative(0) = (a + 1) / (b + 1);
    for (size_t i{1}; i < mmax; ++i) {
      u = (a + i) / ((i + 1) * (b + i));
      coeff(i) = coeff(i - 1) * u;
      u = (a + i + 1) / ((i + 1) * (b + i + 1));
      coeff_derivative(i) = coeff_derivative(i - 1) * u;
    }
  }
}

double Hyp1f1Series::calc(double z, double z2, double ez2, bool derivative,
                          int n_terms) {
  using math::pow;
  double result{0.};
  if (not this->is_exp) {
    result = this->prefac * this->calc(z, derivative, n_terms) * ez2;
  } else {
    result = this->prefac * std::exp(z + z2);
  }
  return result;
}

double Hyp1f1Series::calc(double z, bool derivative, int n_terms) {
  double result{0.};
  if (not this->is_exp) {
    result = this->hyp1f1(z, derivative, n_terms);
  } else {
    result = std::exp(z);
  }
  return result;
}

double Hyp1f1Series::hyp1f1(double z, bool derivative, int n_terms) {
  using math::pow;
  size_t mmax{0};
  if (n_terms == -1) {
    // use adaptive number of terms in the series expansion
    mmax = this->mmax;
  } else if (n_terms > -1) {
    // use a fixed number of terms (useful when computing numerical
    // derivatives)
    mmax = static_cast<size_t>(n_terms);
  } else {
    throw std::runtime_error("n_terms should be >= -1");
  }
  double res{0.};
  if (not derivative) {
    res = this->sum(z, this->coeff, mmax, n_terms);
  } else {
    res =
        this->sum(z, this->coeff_derivative, mmax, n_terms) * this->a / this->b;
  }
  return res;
}

double Hyp1f1Series::sum(double z, const Eigen::VectorXd & coefficient,
                         size_t mmax, int n_terms) {
  // perform the sum
  double res{1.0}, a1{1.0}, zpow{z}, z4{z * z};
  z4 *= z4;
  if (n_terms == -1) {
    // adaptive sum. computes several terms at a time to save on
    // and on bailout tests (typical n. of terms needed is ~20)
    for (size_t i{0}; i < mmax - 3; i += 4) {
      a1 = zpow * (coefficient(i) +
                   z * (coefficient(i + 1) +
                        z * (coefficient(i + 2) + z * coefficient(i + 3))));
      if (a1 < this->tolerance * res) {
        res += a1;
        break;
      }
      zpow *= z4;
      res += a1;
    }
  } else {
    // TODO(mc) this could be done with telescopic sums - and might be
    // faster than checking for bailout condition
    for (int i{0}; i < n_terms; ++i) {
      a1 = coefficient(i) * zpow;
      zpow *= z;
      res += a1;
    }
  }
  if (res > DOVERFLOW) {
    std::stringstream error{};
    error << "Hyp1f1Series series expansion: a=" << std::to_string(this->a)
          << " b=" << std::to_string(this->b) << " z=" << std::to_string(z)
          << std::endl;
    throw std::overflow_error(error.str());
  }
  return res;
}



Hyp1f1Asymptotic::Hyp1f1Asymptotic(double a, double b, size_t mmax,
                                   double tolerance)
    : a{a}, b{b}, prefac{std::tgamma(b) / std::tgamma(a)}, tolerance{tolerance},
      mmax{mmax} {
  double intpart;
  double f2{std::modf(2 * (a - b), &intpart)};
  if (std::abs(f2) < 1e-14) {
    this->is_n_and_l = true;
  }

  if (std::abs(1 - this->b / this->a) < DBL_FTOL) {
    this->is_exp = true;
  } else {
    coeff.resize(mmax);
    coeff_derivative.resize(mmax);
    // precomputes the (1-a)_i / (b-a)_i / i! coefficients for the
    // computation of hyp2f0
    double bma{b - a}, oma{1 - a}, u{0.}, oma1{1 - a - 1};
    coeff(0) = bma * oma;
    coeff_derivative(0) = bma * oma1;
    for (size_t i{1}; i < mmax; ++i) {
      u = (bma + i) * (oma + i) / (i + 1);
      coeff(i) = coeff(i - 1) * u;
      u = (bma + i) * (oma1 + i) / (i + 1);
      coeff_derivative(i) = coeff_derivative(i - 1) * u;
    }
  }
}

double Hyp1f1Asymptotic::z_power_a_b(double z) {
  using math::pow;
  double fac{0.};
  if (this->is_n_and_l) {
    fac = std::sqrt(pow(z, static_cast<int>(2 * (this->a - this->b))));
  } else {
    fac = pow(z, this->a - this->b);
  }
  return fac;
}

double Hyp1f1Asymptotic::calc(double z, double z2, bool derivative,
                              int n_terms) {
  using math::pow;
  double result{0.};
  if (not this->is_exp) {
    auto fac{this->z_power_a_b(z)};
    result = this->hyp2f0(z, derivative, n_terms) * std::exp(z + z2) * fac;
  } else {
    result = std::exp(z + z2) / this->prefac;
  }
  return result;
}

double Hyp1f1Asymptotic::calc(double z, bool derivative, int n_terms) {
  if (not this->is_exp) {
    auto fac{this->z_power_a_b(z)};
    double result{this->prefac * std::exp(z) * fac *
                  this->hyp2f0(z, derivative, n_terms)};
    if (std::isnan(result)) {
      result = DOVERFLOW;
    }
    return result;
  } else {
    return std::exp(z);
  }
}

double Hyp1f1Asymptotic::calc_neg(double z, bool derivative, int n_terms) {
  auto fac{this->z_power_a_b(z)};
  double result{this->prefac * fac *
                this->hyp2f0(z, derivative, n_terms)};
  if (std::isnan(result)) {
    result = DOVERFLOW;
  }
  return result;
}

double Hyp1f1Asymptotic::hyp2f0(double z, bool derivative, int n_terms) {
  using math::pow;

  size_t mmax{0};
  if (n_terms == -1) {
    // use adaptive number of terms in the series expansion
    mmax = this->mmax;
  } else if (n_terms > -1) {
    // use a fixed number of terms (usefull when computing numerical
    // derivatives)
    mmax = static_cast<size_t>(n_terms);
  } else {
    throw std::runtime_error("n_terms should be >= -1");
  }

  double res{0.};
  if (not derivative) {
    res = this->sum(z, this->coeff, mmax, n_terms);
  } else {
    res = this->sum(z, this->coeff_derivative, mmax, n_terms);
  }
  return res;
}

double Hyp1f1Asymptotic::sum(double z, const Eigen::VectorXd & coefficient,
                             size_t mmax, int n_terms) {
  double iz{1.0 / z};
  double res{1.}, izpow{1.}, s_i{1.};
  // perform the sum
  for (size_t i{0}; i < mmax; ++i) {
    izpow *= iz;
    s_i = coefficient(i) * izpow;
    if (res > 0 and std::fabs(s_i) < this->tolerance * res and n_terms == -1) {
      break;
    }
    res += s_i;
  }

  if (res > DOVERFLOW) {
    std::stringstream error{};
    error << "Hyp1f1Asymptotic expansion: a=" << std::to_string(this->a)
          << " b=" << std::to_string(this->b) << " z=" << std::to_string(z)
          << std::endl;
    throw std::overflow_error(error.str());
  }

  return res;
}

Hyp1f1::Hyp1f1(double a, double b, size_t mmax, double tolerance)
    : hyp1f1_series{a, b, mmax, tolerance},
      hyp1f1_asymptotic{a, b, mmax, tolerance},
      hyp1f1_series_neg{b-a, b, mmax, tolerance},
      hyp1f1_asymptotic_neg{b-a, b, mmax, tolerance}, a{a}, b{b}, tolerance{
                                                                tolerance} {
  

  // now we try to determine what is the switching point between
  // power series and asymptotic expansion. basically we choose
  // the method that requires fewer terms for a chosen target accuracy.
  // the asymptotic expansion tends to blow up at the switching point.
  if (std::abs(1 - b / a) > DBL_FTOL) {
    this->z_asympt = this->find_switching_point(this->z_asympt, tolerance,
                              this->hyp1f1_series,
                              this->hyp1f1_asymptotic);
    // fix the number of terms needed for the numerical derivative
    // with nterms_s and nterms_a
    this->hyp1f1_series.calc(this->z_asympt);
    this->hyp1f1_asymptotic.calc(this->z_asympt);
  }

  if (std::abs(1 - b / (b-a)) > DBL_FTOL) {
    this->z_asympt_neg = this->find_switching_point(this->z_asympt_neg, tolerance,
                              this->hyp1f1_series_neg,
                              this->hyp1f1_asymptotic_neg);
    // fix the number of terms needed for the numerical derivative
    // with nterms_s and nterms_a
    this->hyp1f1_series_neg.calc(this->z_asympt_neg);
    this->hyp1f1_asymptotic_neg.calc(this->z_asympt_neg);
  }
}

double Hyp1f1::find_switching_point(double z_asympt, double tolerance, internal::Hyp1f1Series & hyp1f1_series, internal::Hyp1f1Asymptotic & hyp1f1_asymptotic) {
  int max_it{100};
  int i_it{0};
  double z_below{z_asympt}, z_above{z_asympt};
  double h1f1_s{0.}, h1f1_a{0.};

  h1f1_s = hyp1f1_series.calc(z_asympt);
  h1f1_a = hyp1f1_asymptotic.calc(z_asympt);
  // brackets the switching point
  if (std::fabs(1 - h1f1_s / h1f1_a) > tolerance) {
    i_it = 0;
    while (std::fabs(1 - h1f1_s / h1f1_a) >
               100 * tolerance and
           i_it < max_it) {
      z_asympt *= 1.5;
      h1f1_s = hyp1f1_series.calc(z_asympt);
      h1f1_a = hyp1f1_asymptotic.calc(z_asympt);
      i_it++;
    }
    z_above = z_asympt;
  } else {
    i_it = 0;
    while (std::fabs(1 - h1f1_s / h1f1_a) <
               100 * tolerance and
           i_it < max_it) {
      z_asympt *= 0.5;
      h1f1_s = hyp1f1_series.calc(z_asympt);
      h1f1_a = hyp1f1_asymptotic.calc(z_asympt);
      i_it++;
    }
    z_below = z_asympt;
  }
  // and now bisects until we are reasonably close to an accurate
  // determination
  z_asympt = (z_above + z_below) * 0.5;
  h1f1_s = hyp1f1_series.calc(z_asympt);
  h1f1_a = hyp1f1_asymptotic.calc(z_asympt);
  i_it = 0;
  while (z_above - z_below > tolerance and i_it < max_it) {
    if (std::abs(1 - h1f1_s / h1f1_a) > tolerance) {
      z_below = z_asympt;
    } else {
      z_above = z_asympt;
    }
    z_asympt = (z_above + z_below) * 0.5;
    h1f1_s = hyp1f1_series.calc(z_asympt);
    h1f1_a = hyp1f1_asymptotic.calc(z_asympt);
    i_it++;
  }
  return z_asympt;
}

double Hyp1f1::calc_numerical_derivative(double z, double h) {
  if (z < 0.) {
    double fzp{hyp1f1(this->a, this->b, z + h)};
    double fzm{hyp1f1(this->a, this->b,z - h)};
    return (fzp - fzm) / (2 * h);
  } else if (z > this->z_asympt) {
    double fzp{this->hyp1f1_asymptotic.calc(z + h)};
    int n_terms{this->hyp1f1_asymptotic.n_terms};
    double fzm{this->hyp1f1_asymptotic.calc(z - h, false, n_terms)};
    return (fzp - fzm) / (2 * h);
  } else {
    double fzp{this->hyp1f1_series.calc(z + h)};
    int n_terms{this->hyp1f1_series.n_terms};
    double fzm{this->hyp1f1_series.calc(z - h, false, n_terms)};
    return (fzp - fzm) / (2 * h);
  }
}

double Hyp1f1::calc(double z, double z2, double ez2, bool derivative) {
  double res{};
  if (z > this->z_asympt) {
    res = this->hyp1f1_asymptotic.calc(z, z2, derivative);
  } else {
    res = this->hyp1f1_series.calc(z, z2, ez2, derivative);
  }
  assert(std::isfinite(res));
  return res;
}

void Hyp1f1SphericalExpansion::precompute(size_t max_radial,
                                          size_t max_angular) {
  this->max_angular = max_angular;
  this->max_radial = max_radial;
  this->values.resize(max_radial, max_angular + 1);
  this->derivatives.resize(max_radial, max_angular + 1);
  this->z.resize(max_radial);
  this->dz_dr.resize(max_radial);

  for (size_t n_radial{0}; n_radial < max_radial; n_radial++) {
    for (size_t l_angular{0}; l_angular < max_angular + 1; l_angular++) {
      auto a{this->get_a(n_radial, l_angular)};
      auto b{this->get_b(l_angular)};
      hyp1f1.emplace_back(a, b, precomputation_size, tolerance);
    }
  }
}

double Hyp1f1SphericalExpansion::calc(size_t n_radial, size_t l_angular,
                                      double r_ij, double alpha, double beta) {
  int ipos{this->get_pos(n_radial, l_angular)};
  double z{math::pow(alpha * r_ij, 2) / (alpha + beta)};
  double z2{-alpha * r_ij * r_ij};
  double ez2{std::exp(z2)};
  return this->hyp1f1[ipos].calc(z, z2, ez2);
}

void Hyp1f1SphericalExpansion::calc(double r_ij, double alpha,
                                    const Vector_Ref & fac_b, bool derivative) {
  if (not this->recursion or this->max_angular < 3) {
    // recursion needs 4 evaluations of 1F1 so not worth it if l_max < 3
    this->calc_direct(r_ij, alpha, fac_b, derivative);
  } else {
    this->calc_recursion(r_ij, alpha, fac_b);
  }
}

void Hyp1f1SphericalExpansion::calc_recursion(double r_ij, double alpha,
                                              const Vector_Ref & fac_b) {
  double M1p2p{0.}, M2p3p{0.}, MP1p2p{0.}, MP2p3p{0.}, M1p1p{0.}, Moo{0.},
      MP1p1p{0.}, MPoo{0.};

  double alpha_rij{alpha * r_ij};
  double z2{-r_ij * alpha_rij};
  double ez2{std::exp(z2)};
  this->z = (alpha_rij * alpha) * (alpha + fac_b.array()).inverse();
  this->dz_dr = this->z;
  this->z *= r_ij;
  this->dz_dr *= 2;

  for (size_t n_radial{0}; n_radial < this->max_radial; ++n_radial) {
    // get the starting points for the recursion
    auto & z{this->z(n_radial)};
    int l_angular{static_cast<int>(this->max_angular)};
    int ipos{this->get_pos(n_radial, l_angular)};
    M1p2p = this->hyp1f1[ipos].calc(z, z2, ez2);
    this->values(n_radial, l_angular) = M1p2p;
    M2p3p = this->hyp1f1[ipos].calc(z, z2, ez2, true);
    this->derivatives(n_radial, l_angular) = M2p3p;

    ipos = this->get_pos(n_radial, l_angular - 1);
    MP1p2p = this->hyp1f1[ipos].calc(z, z2, ez2);
    this->values(n_radial, l_angular - 1) = MP1p2p;
    MP2p3p = this->hyp1f1[ipos].calc(z, z2, ez2, true);
    this->derivatives(n_radial, l_angular - 1) = MP2p3p;
    l_angular -= 2;
    for (; l_angular > 0; l_angular -= 2) {
      auto a{this->get_a(n_radial, l_angular)};
      auto b{this->get_b(l_angular)};
      M1p1p = recurence_G_to_der_downward(a, b, z, M2p3p, M1p2p);
      Moo = recurence_G_to_val_downward(a, b, z, M1p2p, M1p1p);
      this->values(n_radial, l_angular) = Moo;
      this->derivatives(n_radial, l_angular) = M1p1p;
      M2p3p = M1p1p;
      M1p2p = Moo;

      a = this->get_a(n_radial, l_angular - 1);
      b = this->get_b(l_angular - 1);
      MP1p1p = recurence_G_to_der_downward(a, b, z, MP2p3p, MP1p2p);
      MPoo = recurence_G_to_val_downward(a, b, z, MP1p2p, MP1p1p);
      this->values(n_radial, l_angular - 1) = MPoo;
      this->derivatives(n_radial, l_angular - 1) = MP1p1p;
      MP2p3p = MP1p1p;
      MP1p2p = MPoo;
    }
    // makes sure l == 0 is taken care of
    if (this->max_angular % 2 == 0) {
      auto a{this->get_a(n_radial, 0)};
      auto b{this->get_b(0)};
      M1p1p = recurence_G_to_der_downward(a, b, z, M2p3p, M1p2p);
      Moo = recurence_G_to_val_downward(a, b, z, M1p2p, M1p1p);
      this->values(n_radial, 0) = Moo;
      this->derivatives(n_radial, 0) = M1p1p;
    }

    this->derivatives.row(n_radial) *= this->dz_dr(n_radial);
  }
  // here is where dG/dz*dz/dr is computed
  this->derivatives -= (2 * alpha_rij) * this->values;
}

void Hyp1f1SphericalExpansion::calc_direct(double r_ij, double alpha,
                                           const Vector_Ref & fac_b,
                                           bool derivative) {
  // computes some intermediates that accelerate calculations further
  // down
  double alpha_rij{alpha * r_ij};
  double z2{-r_ij * alpha_rij};
  double ez2{std::exp(z2)};
  this->z = (alpha_rij * alpha) * (alpha + fac_b.array()).inverse();
  this->dz_dr = this->z;
  this->z *= r_ij;
  this->dz_dr *= 2;

  for (size_t n_radial{0}; n_radial < this->max_radial; n_radial++) {
    auto & z{this->z(n_radial)};
    for (size_t l_angular{0}; l_angular < this->max_angular + 1; l_angular++) {
      int ipos{this->get_pos(n_radial, l_angular)};
      this->values(n_radial, l_angular) = this->hyp1f1[ipos].calc(z, z2, ez2);
    }
  }
  if (derivative) {
    for (size_t n_radial{0}; n_radial < this->max_radial; n_radial++) {
      auto & z{this->z(n_radial)};
      for (size_t l_angular{0}; l_angular < this->max_angular + 1;
           l_angular++) {
        int ipos{this->get_pos(n_radial, l_angular)};
        this->derivatives(n_radial, l_angular) =
            this->hyp1f1[ipos].calc(z, z2, ez2, true);
      }
      this->derivatives.row(n_radial) *= this->dz_dr(n_radial);
    }
    this->derivatives -= (2 * alpha * r_ij) * this->values;
  }
}
