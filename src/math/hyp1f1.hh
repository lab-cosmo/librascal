/**
 * @file   hyp1f1.hh
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

#ifndef SRC_MATH_HYP1F1_HH_
#define SRC_MATH_HYP1F1_HH_

#include "math/math_utils.hh"

#include <map>
#include <vector>

namespace rascal {
  namespace math {

    // This is kept as a header only for a small but not marginal performance
    // gain

    /**
     * Utilities to apply the recurrence relations of 1F1 from maximal
     * values of a and b to minimal values.
     * The recurrence path is chosen so that derivatives of 1F1 are also
     * a result.
     *
     * G refers to the modified 1F1 defined below
     */
    inline double recurence_to_val_downward(double a, double b, double z,
                                            double M1p2p, double M1p1p) {
      return z * (a - b) * M1p2p / (b * (b + 1)) + M1p1p;
    }

    inline double recurence_to_der_downward(double a, double b, double z,
                                            double M2p3p, double M1p2p) {
      return z * (a + 1) * M2p3p / ((b + 2) * (b + 1)) + M1p2p;
    }

    inline double recurence_G_to_val_downward(double a, double b, double z,
                                              double M1p2p, double M1p1p) {
      return (z * (a - b) * M1p2p + M1p1p * b) / a;
    }

    inline double recurence_G_to_der_downward(double, double b, double z,
                                              double M2p3p, double M1p2p) {
      return z * M2p3p + M1p2p * (b + 1);
    }

    /**
     * Computes the 1F1 with the direct sum
     *  @f[
     *      1F1(a,b,z) = \sum_{j=0}^{\infty} \frac{(a)_j}{(b)_jj!} z^{j}
     *  @f]
     *
     *  @f[
     *      G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
     *                    * 1F1(a,b,z)
     *  @f]
     */
    class Hyp1f1Series {
     protected:
      double a, b;
      size_t mmax;
      double prefac;
      double tolerance;
      bool is_exp{false};

      Eigen::VectorXd coeff{};
      Eigen::VectorXd coeff_derivative{};

     public:
      size_t n_terms{0};

      Hyp1f1Series(double a, double b, size_t mmax, double tolerance = 1e-14)
          : a{a}, b{b}, mmax{mmax}, prefac{std::tgamma(a) / std::tgamma(b)},
            tolerance{tolerance} {
        // when a == b, 1F1 is an exponential
        if (std::abs(1 - this->b / this->a) < dbl_ftol) {
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
      /**
       * Computes G(a,b,z)
       * @param z -> pow(alpha * r_ij, 2) / (alpha + beta)
       * @param z2 -> argument of exp(-alpha*r_ij^2)
       * @param ez2 -> exp(-alpha*r_ij^2)
       *
       * @warning the derivative outputed from this function is not dG/dz
       * but @f$ d1F1/dz * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
       * @f$. We do this to avoid computing both d1F1/dz and 1F1 when asking for
       * gradients and perform this step in Hyp1f1SphericalExpansion.
       */
      double calc(double z, double z2, double ez2, bool derivative = false,
                  int n_terms = -1) {
        using math::pow;
        double result{0.};
        if (not this->is_exp) {
          result = this->prefac * this->calc(z, derivative, n_terms) * ez2;
        } else {
          result = this->prefac * std::exp(z + z2);
        }
        return result;
      }

      double calc(double z, bool derivative = false, int n_terms = -1) {
        double result{0.};
        if (not this->is_exp) {
          result = this->hyp1f1(z, derivative, n_terms);
        } else {
          result = std::exp(z);
        }
        return result;
      }
      //! Computes 1F1

      double hyp1f1(double z, bool derivative, int n_terms) {
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
        this->n_terms = this->mmax;
        double res{0.};
        if (not derivative) {
          res = this->sum(z, this->coeff, mmax, n_terms);
        } else {
          res = this->sum(z, this->coeff_derivative, mmax, n_terms) * this->a /
                this->b;
        }
        return res;
      }

      double sum(double z, const Eigen::VectorXd & coefficient, size_t mmax,
                 int n_terms) {
        // perform the sum
        double res{1.0}, a1{1.0}, zpow{z}, z4{z * z};
        z4 *= z4;
        if (n_terms == -1) {
          // adaptive sum. computes several terms at a time to save on
          // and on bailout tests (typical n. of terms needed is ~20)
          for (size_t i{0}; i < mmax - 3; i += 4) {
            a1 = zpow * (coefficient(i) + z * (coefficient(i + 1) +
                                               z * (coefficient(i + 2) +
                                                    z * coefficient(i + 3))));
            if (a1 < this->tolerance * res) {
              this->n_terms = i;
              res += a1;
              break;
            }
            zpow *= z4;
            res += a1;
          }
        } else {
          // TODO(mc) this could be done with telescopic sums - and might be
          // faster than checking for bailout condition
          for (size_t i{0}; i < static_cast<size_t>(n_terms); ++i) {
            a1 = coefficient(i) * zpow;
            zpow *= z;
            res += a1;
          }
        }
        if (res > DOVERFLOW) {
          std::stringstream error{};
          error << "Hyp1f1Series series expansion: a="
                << std::to_string(this->a) << " b=" << std::to_string(this->b)
                << " z=" << std::to_string(z) << std::endl;
          throw std::overflow_error(error.str());
        }
        return res;
      }
    };

    /**
     * Computes the 1F1 with the asymptotic limit
     * @f[
     *      1F1(a,b,z) \sim \exp{z} z^{a-b} \frac{\Gamma{b}}{\Gamma{a}}
     *                    \sum_{j=0}^{\infty} \frac{(b-a)_j(1-a)_j}{j!} z^{-j}
     * @f]
     *
     * @f[
     *  G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
     *                    * 1F1(a,b,z)
     * @f]
     */
    class Hyp1f1Asymptotic {
     private:
      double a, b, prefac;
      bool is_n_and_l{false};
      bool is_exp{false};
      double tolerance;
      size_t mmax;
      Eigen::VectorXd coeff{};
      Eigen::VectorXd coeff_derivative{};

      double z_power_a_b(double z) {
        using math::pow;
        double fac{0.};
        if (this->is_n_and_l) {
          fac = std::sqrt(pow(z, static_cast<int>(2 * (this->a - this->b))));
        } else {
          fac = pow(z, this->a - this->b);
        }
        return fac;
      }

     public:
      size_t n_terms{0};

      Hyp1f1Asymptotic(double a, double b, size_t mmax,
                       double tolerance = 1e-14)
          : a{a}, b{b}, prefac{std::tgamma(b) / std::tgamma(a)},
            tolerance{tolerance}, mmax{mmax} {
        double intpart;
        double f2{std::modf(2 * (a - b), &intpart)};
        if (std::abs(f2) < 1e-14) {
          this->is_n_and_l = true;
        }

        if (std::abs(1 - this->b / this->a) < dbl_ftol) {
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

      /**
       * Computes G(a,b,z)
       * @param z -> @f$ pow(alpha * r_ij, 2) / (alpha + beta) @f$
       * @param z2 -> argument of @f$ exp(-alpha*r_ij^2) @f$
       *
       * @warning the derivative outputed from this function is not dG/dz
       * but @f$ d1F1/dz * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
       * @f$. We do this to avoid computing both d1F1/dz and 1F1 when asking for
       * gradients and perform this step in Hyp1f1SphericalExpansion.
       */
      double calc(double z, double z2, bool derivative = false,
                  int n_terms = -1) {
        using math::pow;
        double result{0.};
        if (not this->is_exp) {
          auto fac{this->z_power_a_b(z)};
          result =
              this->hyp2f0(z, derivative, n_terms) * std::exp(z + z2) * fac;
        } else {
          result = std::exp(z + z2) / this->prefac;
        }
        return result;
      }

      //! Computes 1F1
      double calc(double z, bool derivative = false, int n_terms = -1) {
        using math::pow;
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

      //! computes hyp2f0 with arg1 = b-a and arg2 = 1-a arg3 = 1 / z
      double hyp2f0(double z, bool derivative, int n_terms) {
        using math::pow;
        this->n_terms = this->mmax;

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

      double sum(double z, const Eigen::VectorXd & coefficient, size_t mmax,
                 int n_terms) {
        double iz{1.0 / z};
        double res{1.}, izpow{1.}, s_i{1.};
        // perform the sum
        for (size_t i{0}; i < mmax; ++i) {
          izpow *= iz;
          s_i = coefficient(i) * izpow;
          if (res > 0 and std::fabs(s_i) < this->tolerance * res and
              n_terms == -1) {
            this->n_terms = i;
            break;
          }
          res += s_i;
        }

        if (res > DOVERFLOW) {
          std::stringstream error{};
          error << "Hyp1f1Asymptotic expansion: a=" << std::to_string(this->a)
                << " b=" << std::to_string(this->b)
                << " z=" << std::to_string(z) << std::endl;
          throw std::overflow_error(error.str());
        }

        return res;
      }
    };

    /**
     * Computes the confluent hypergeometric function \f${}_1F_1(a,b,z)\f$ for
     * a given a and b and variable argument.
     *
     * The class to bundles the two definitions of the 1F1 (Hyp1f1Asymptotic
     * and Hyp1f1Series) and implements a switch rational (with bisection) to
     * go between the two at construction.
     * It works because we probe test at construction the domain of
     * applicability of the two definitions.
     */
    class Hyp1f1 {
     private:
      Hyp1f1Series hyp1f1_series;
      Hyp1f1Asymptotic hyp1f1_asymptotic;
      double a, b;
      double tolerance;

      double z_asympt{1.};

      size_t nterms_a{0}, nterms_s{0};
      double h1f1_a{0}, h1f1_s{0};
      double z_above{0.}, z_below{0.};

      void update_switching_point() {
        this->h1f1_s = this->hyp1f1_series.calc(this->z_asympt);
        this->h1f1_a = this->hyp1f1_asymptotic.calc(this->z_asympt);

        this->nterms_s = this->hyp1f1_series.n_terms;
        this->nterms_a = this->hyp1f1_asymptotic.n_terms;
      }

      /**
       * Find the largest z for which the 2 definitions agree within tolerance
       * using bisection.
       */
      void find_switching_point() {
        int max_it{100};
        int i_it{0};
        this->update_switching_point();
        // brackets the switching point
        this->z_below = this->z_above = this->z_asympt;
        if (std::fabs(1 - this->h1f1_s / this->h1f1_a) > this->tolerance) {
          i_it = 0;
          while (std::fabs(1 - this->h1f1_s / this->h1f1_a) >
                     100 * this->tolerance and
                 i_it < max_it) {
            this->z_asympt *= 1.5;
            this->update_switching_point();
            i_it++;
          }
          this->z_above = this->z_asympt;
        } else {
          i_it = 0;
          while (std::fabs(1 - this->h1f1_s / this->h1f1_a) <
                     100 * this->tolerance and
                 i_it < max_it) {
            this->z_asympt *= 0.5;
            this->update_switching_point();
            i_it++;
          }
          this->z_below = this->z_asympt;
        }
        // and now bisects until we are reasonably close to an accurate
        // determination
        this->z_asympt = (this->z_above + this->z_below) * 0.5;
        this->update_switching_point();
        i_it = 0;
        while (this->z_above - this->z_below > this->tolerance and
               i_it < max_it) {
          if (std::abs(1 - this->h1f1_s / this->h1f1_a) > this->tolerance) {
            this->z_below = this->z_asympt;
          } else {
            this->z_above = this->z_asympt;
          }
          this->z_asympt = (this->z_above + this->z_below) * 0.5;
          this->update_switching_point();
          i_it++;
        }
        // to be safe take a slightly larger switching point
        this->z_asympt = this->z_asympt * 1.1;
      }

     public:
      Hyp1f1(double a, double b, size_t mmax = 500, double tolerance = 1e-13)
          : hyp1f1_series{a, b, mmax, tolerance}, hyp1f1_asymptotic{a, b, mmax,
                                                                    tolerance},
            a{a}, b{b}, tolerance{tolerance} {
        // now we try to determine what is the switching point between
        // power series and asymptotic expansion. basically we choose
        // the method that requires fewer terms for a chosen target accuracy.
        // the asymptotic expansion tends to blow up at the switching point.
        if (std::abs(1 - this->b / this->a) > dbl_ftol) {
          this->find_switching_point();
          // fix the number of terms needed for the numerical derivative
          // with nterms_s and nterms_a
          this->update_switching_point();
        }
      }

      double get_z_switch() { return this->z_asympt; }

      //! Compute @f${}_1F_1(a,b,z)@f$
      double calc(double z, bool derivative = false) {
        if (z > this->z_asympt) {
          return this->hyp1f1_asymptotic.calc(z, derivative);
        } else {
          return this->hyp1f1_series.calc(z, derivative);
        }
      }

      double calc_numerical_derivative(double z, double h) {
        if (z > this->z_asympt) {
          double fzp{this->hyp1f1_asymptotic.calc(z + h)};
          size_t n_terms{this->hyp1f1_asymptotic.n_terms};
          double fzm{this->hyp1f1_asymptotic.calc(z - h, false, n_terms)};
          return (fzp - fzm) / (2 * h);
        } else {
          double fzp{this->hyp1f1_series.calc(z + h)};
          size_t n_terms{this->hyp1f1_series.n_terms};
          double fzm{this->hyp1f1_series.calc(z - h, false, n_terms)};
          return (fzp - fzm) / (2 * h);
        }
      }

      /**
       * Computes G(a,b,z)
       * @param z -> @f$ pow(alpha * r_ij, 2) / (alpha + beta) @f$
       * @param z2 -> argument of @f$ exp(-alpha*r_ij^2) @f$
       * @param ez2 -> @f$ exp(-alpha*r_ij^2) @f$
       *
       * @warning the derivative outputed from this function is not @f$dG/dz@f$
       * but @f$ \frac{\mathrm{d}\,{}_1F_1}{\mathrm{d}\,z} *
       * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2} @f$.
       * We do this to avoid computing both @f$d{}_1F_1/dz@f$ and @f${}_1F_1@f$
       * when asking for gradients and perform this step in
       * Hyp1f1SphericalExpansion.
       */
      double calc(double z, double z2, double ez2, bool derivative = false) {
        if (z > this->z_asympt) {
          return this->hyp1f1_asymptotic.calc(z, z2, derivative);
        } else {
          return this->hyp1f1_series.calc(z, z2, ez2, derivative);
        }
      }
    };

    /**
     * Computes 1F1 and its derivative with respect to z for a range of a and b
     * for @f$l < l_\mathrm{max} + 1@f$ and @f$n < n_\mathrm{max}@f$ where:
     * @f$ a = 0.5 * (n + l + 3) @f$
     * @f$ b = l + 1.5 @f$
     *
     * For efficiency the function computed is:
     * @f[
     *   G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
     *                     1F1(a,b,z)
     * @f]
     *
     * It can use the recurence relationships of the 1F1 to speed things up.
     *
     * This class is tailored to work with the GTO basis in
     * CalculatorSphericalExpansion.
     */
    class Hyp1f1SphericalExpansion {
     protected:
      using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;
      using Vector_t = Eigen::VectorXd;
      using Vector_Ref = typename Eigen::Ref<const Vector_t>;

      size_t max_angular{0}, max_radial{0};
      std::vector<Hyp1f1> hyp1f1{};

      Matrix_t values{};
      Matrix_t derivatives{};
      double tolerance;
      size_t precomputation_size;
      bool recursion;
      Eigen::ArrayXd z{};
      Eigen::ArrayXd dz_dr{};

      int get_pos(int n_radial, int l_angular) {
        return l_angular + (this->max_angular + 1) * n_radial;
      }

      double get_a(int n_radial, int l_angular) {
        return 0.5 * (n_radial + l_angular + 3);
      }

      double get_b(int l_angular) { return l_angular + 1.5; }

     public:
      Hyp1f1SphericalExpansion(bool recursion = false, double tolerance = 1e-14,
                               size_t precomputation_size = 200)
          : tolerance{tolerance},
            precomputation_size{precomputation_size}, recursion{recursion} {}

      //! initialize the 1F1 computers
      void precompute(size_t max_radial, size_t max_angular) {
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

      //! helper function to compute z for one set of n, l, z
      double calc(size_t n_radial, size_t l_angular, double z) {
        int ipos{this->get_pos(n_radial, l_angular)};
        return this->hyp1f1[ipos].calc(z);
      }

      //! helper function to compute z for one set of n, l, custom z
      double calc(size_t n_radial, size_t l_angular, double r_ij, double alpha,
                  double beta) {
        int ipos{this->get_pos(n_radial, l_angular)};
        double z{math::pow(alpha * r_ij, 2) / (alpha + beta)};
        double z2{-alpha * r_ij * r_ij};
        double ez2{std::exp(z2)};
        return this->hyp1f1[ipos].calc(z, z2, ez2);
      }

      //! the work horse that computes G for all possible n, l values
      void calc(double r_ij, double alpha, const Vector_Ref & fac_b,
                bool derivative = false) {
        if (not this->recursion or this->max_angular < 3) {
          // recursion needs 4 evaluations of 1F1 so not worth it if l_max < 3
          this->calc_direct(r_ij, alpha, fac_b, derivative);
        } else {
          this->calc_recursion(r_ij, alpha, fac_b);
        }
      }

      /**
       *  Computes G and dG/dz*dz/dr using recursion relations with
       *  z = (alpha * r_ij)**2 / (alpha + fac_b)
       */
      void calc_recursion(double r_ij, double alpha, const Vector_Ref & fac_b) {
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

      //! computes G by direct evaluation
      void calc_direct(double r_ij, double alpha, const Vector_Ref & fac_b,
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
          for (size_t l_angular{0}; l_angular < this->max_angular + 1;
               l_angular++) {
            int ipos{this->get_pos(n_radial, l_angular)};
            this->values(n_radial, l_angular) =
                this->hyp1f1[ipos].calc(z, z2, ez2);
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

      //! get a reference to the computed G values
      Matrix_Ref get_values() { return Matrix_Ref(this->values); }

      //! get a reference to the computed G derivatives
      Matrix_Ref get_derivatives() { return Matrix_Ref(this->derivatives); }
    };

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_HYP1F1_HH_
