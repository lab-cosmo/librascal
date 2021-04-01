/**
 * @file   rascal/math/hyp1f1.hh
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

#ifndef SRC_RASCAL_MATH_HYP1F1_HH_
#define SRC_RASCAL_MATH_HYP1F1_HH_

#include "rascal/math/utils.hh"

#include <vector>

namespace rascal {
  namespace math {
    namespace internal {
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
       private:
        double a, b;
        size_t mmax;
        double prefac;
        double tolerance;
        bool is_exp{false};

        Eigen::VectorXd coeff{};
        Eigen::VectorXd coeff_derivative{};

       public:
        // Parameter used to fix the sum when evaluating the numerical
        // derivatives.
        int n_terms{100};

        Hyp1f1Series(double a, double b, size_t mmax, double tolerance = 1e-14);

        /**
         * Computes G(a,b,z)
         * @param z -> @f$ pow(alpha * r_ij, 2) / (alpha + beta) @f$
         * @param z2 -> argument of @f$ exp(-alpha*r_ij^2) @f$
         * @param ez2 -> @f$ exp(-alpha*r_ij^2) @f$
         *
         * @warning the derivative outputed from this function is not dG/dz
         * but @f$ d1F1/dz * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha
         * r_{ij}^2}
         * @f$. We do this to avoid computing both d1F1/dz and 1F1 when asking
         * for gradients and perform this step in Hyp1f1SphericalExpansion.
         */
        double calc(double z, double z2, double ez2, bool derivative = false,
                    int n_terms = -1);

        double calc(double z, bool derivative = false, int n_terms = -1);

        //! Computes 1F1
        double hyp1f1(double z, bool derivative, int n_terms);

        double sum(double z, const Eigen::VectorXd & coefficient, size_t mmax,
                   int n_terms);
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

        double z_power_a_b(double z);

       public:
        // Parameter used to fix the sum when evaluating the numerical
        // derivatives.
        int n_terms{20};

        Hyp1f1Asymptotic(double a, double b, size_t mmax,
                         double tolerance = 1e-14);
        /**
         * Computes G(a,b,z)
         * @param z -> @f$ pow(alpha * r_ij, 2) / (alpha + beta) @f$
         * @param z2 -> argument of @f$ exp(-alpha*r_ij^2) @f$
         *
         * @warning the derivative outputed from this function is not dG/dz
         * but @f$ d1F1/dz * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha
         * r_{ij}^2}
         * @f$. We do this to avoid computing both d1F1/dz and 1F1 when asking
         * for gradients and perform this step in Hyp1f1SphericalExpansion.
         */
        double calc(double z, double z2, bool derivative = false,
                    int n_terms = -1);

        //! Computes 1F1
        double calc(double z, bool derivative = false, int n_terms = -1);

        double calc_neg(double z, bool derivative = false, int n_terms = -1);

        //! computes hyp2f0 with arg1 = b-a and arg2 = 1-a arg3 = 1 / z
        double hyp2f0(double z, bool derivative, int n_terms);

        double sum(double z, const Eigen::VectorXd & coefficient, size_t mmax,
                   int n_terms);
      };
    }  // namespace internal


    inline double hyp1f1_sum(double a, double b, double z) {
      double epsilon{1e-6};
      int k{0};
      double term{1.};
      double result{1.};
      double abssum{result};

      for (; k < 1000; k++) {
        term *= (a + k) * z / (b + k) / (k + 1);
        abssum += std::abs(term);
        result += term;
        if (std::abs(term) <= math::DBL_FTOL * std::abs(result)) {
          break;
        }
      }
      assert(std::isfinite(result));
      if (k * math::DBL_FTOL * abssum <= epsilon * std::abs(result)) {
        return result;
      } else {
        std::stringstream error{};
        error << "hyp1f1 series expansion: a=" << std::to_string(a)
              << " b=" << std::to_string(b) << " z=" << std::to_string(z)
              << std::endl;
        throw std::overflow_error(error.str());
      }
    }

    inline double hyp1f1(double a, double b, double x) {
      double res{0.};
      double z{x};
      if (x < 0.) {
        a = b-a;
        z = -x;
      }
      if (std::abs(1 - b / a) < math::DBL_FTOL) {
        res = std::exp(z);
      } else if (std::abs(a - b - 1) < math::DBL_FTOL) {
        res = (1 + z / b) * std::exp(z);
      } else {
        res = hyp1f1_sum(a,b,z);
      }
      
      if (x < 0.) {
        res *= std::exp(-z);
      }
      assert(std::isfinite(res));
      return res;
    }

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
      internal::Hyp1f1Series hyp1f1_series;
      internal::Hyp1f1Asymptotic hyp1f1_asymptotic;
      internal::Hyp1f1Series hyp1f1_series_neg;
      internal::Hyp1f1Asymptotic hyp1f1_asymptotic_neg;
      double a, b;
      double tolerance;

      double z_asympt{1.}, z_asympt_neg{1.};

      /**
       * Find the largest z for which the 2 definitions agree within tolerance
       * using bisection.
       */
      double
      find_switching_point(double z_init, double tolerance,
                           internal::Hyp1f1Series & hyp1f1_series,
                           internal::Hyp1f1Asymptotic & hyp1f1_asymptotic);

     public:
      Hyp1f1(double a, double b, size_t mmax = 500, double tolerance = 1e-13);

      double get_z_switch() { return this->z_asympt; }

      //! Compute @f${}_1F_1(a,b,z)@f$
      double calc(double z, bool derivative = false) {
        double res{};
        if (z >= 0.) {
          if (z > this->z_asympt) {
            res = this->hyp1f1_asymptotic.calc(z, derivative);
          } else {
            res = this->hyp1f1_series.calc(z, derivative);
          }
        } else {
          if (not derivative) {
            res = hyp1f1(this->a, this->b, z);
          } else {
            res = this->a / this->b * hyp1f1(this->a + 1, this->b + 1, z);
          }
        }

        assert(std::isfinite(res));
        return res;
      }

      double calc_numerical_derivative(double z, double h);

      /**
       * Computes G(a,b,z)
       * @param z -> @f$ pow(alpha * r_ij, 2) / (alpha + beta) @f$
       * @param z2 -> argument of @f$ exp(-alpha*r_ij^2) @f$
       * @param ez2 -> @f$ exp(-alpha*r_ij^2) @f$
       *
       * @warning the derivative outputed from this function is not
       * @f$dG/dz@f$ but @f$ \frac{\mathrm{d}\,{}_1F_1}{\mathrm{d}\,z} *
       * \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2} @f$.
       * We do this to avoid computing both @f$d{}_1F_1/dz@f$ and
       * @f${}_1F_1@f$ when asking for gradients and perform this step in
       * Hyp1f1SphericalExpansion.
       */
      double calc(double z, double z2, double ez2, bool derivative = false);
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
     private:
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
      void precompute(size_t max_radial, size_t max_angular);

      //! helper function to compute z for one set of n, l, z
      double calc(size_t n_radial, size_t l_angular, double z) {
        int ipos{this->get_pos(n_radial, l_angular)};
        return this->hyp1f1[ipos].calc(z);
      }

      //! helper function to compute z for one set of n, l, custom z
      double calc(size_t n_radial, size_t l_angular, double r_ij, double alpha,
                  double beta);

      //! the work horse that computes G for all possible n, l values
      void calc(double r_ij, double alpha, const Vector_Ref & fac_b,
                bool derivative = false);

      /**
       *  Computes G and dG/dz*dz/dr using recursion relations with
       *  @f$ z = (alpha * r_ij)**2 / (alpha + fac_b) @f$
       */
      void calc_recursion(double r_ij, double alpha, const Vector_Ref & fac_b);

      //! computes G by direct evaluation
      void calc_direct(double r_ij, double alpha, const Vector_Ref & fac_b,
                       bool derivative);

      //! get a reference to the computed G values
      Matrix_Ref get_values() { return Matrix_Ref(this->values); }

      //! get a reference to the computed G derivatives
      Matrix_Ref get_derivatives() { return Matrix_Ref(this->derivatives); }
    };

  }  // namespace math
}  // namespace rascal

#endif  // SRC_RASCAL_MATH_HYP1F1_HH_
