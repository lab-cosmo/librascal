/**
 * @file   rascal/math/bessel.hh
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

#ifndef SRC_RASCAL_MATH_BESSEL_HH_
#define SRC_RASCAL_MATH_BESSEL_HH_

#include "rascal/math/hyp1f1.hh"
#include "rascal/math/utils.hh"

#include <Eigen/Core>

#include <vector>

namespace rascal {
  namespace math {

    /**
     * Computes the modified spherical bessel function of the first kind (MBSF)
     * as
     * \f[
     *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2*a*r*x_n)
     * \f]
     * Just call precompute() once, calc(), and use get_values() to access
     * the values
     *
     * Note that this uses the "exp-exp complete-square" optimized function,
     * where the Bessel function is multiplied by two exponentials that keep
     * it from blowing up.
     *
     * The recursion relation used here is:
     * \f[
     *    i_0(x) = sinh(x) / x
     *    i_1(x) = (x*cosh(x) + sinh(x)) / x^2
     *    i_n(x) = i_{n-2}(x) - (2n - 1)/x * i_{n-1}(x)
     * \f]
     * from
     * http://mathworld.wolfram.com/ModifiedSphericalBesselFunctionoftheFirstKind.html
     */
    class ModifiedSphericalBessel {
     public:
      ModifiedSphericalBessel() = default;

      /**
       * Compute all the MBSFs for the given x-values up to the given order
       *
       * The MBSFs are accurate when the expected value is > 1e-100. Below this
       * threshold the MBSFs are set to 0 because of the numerical noise
       * arrising below 1e-150.
       */
      void calc(double distance, double fac_a);

      /**
       * Initialize arrays for the computation of
       * @param x_v      Eigen::Array of x-values (part of the argument of the
       *                 first exponential; see equation above)
       */
      void precompute(size_t l_max,
                      const Eigen::Ref<const Eigen::VectorXd> & x_v);

      /**
       * Return a reference to the precomputed Bessel function values
       *
       * @return Eigen::Array (2-D) of Bessel function values.  The different
       *         arguments (xs) go along the rows, while the column indexes
       *         the orders (n-values).
       *         Note that a reference is returned to avoid unnecessary
       *         copies.
       */
      Eigen::Ref<Eigen::ArrayXXd> get_values() { return bessel_values; }

     private:
      /**
       * Compute the MBSFs times two exponentials that complete the square
       * using upward recursion. This is stable when 2*a*r*x_n > 50, so
       * it is useful to avoid overflow/underflow in the individual terms of f.
       *
       * @param distance (single) value of distance
       *                (second exponential argument) in the equation above
       * @param fac_a  Scaling factor for the exponential arguments
       * @param n_rows number of rows where the recursion is applicable
       *               from the bottom
       */
      void upward_recursion(double distance, double fac_a, int n_rows);

      /**
       * Compute the MBSFs times two exponentials using downward recurence
       * which is stable in general but when a, r and/or x_n becomes too large
       * one of the terms of f might overflow/underflow while f is finite.
       * The recurence relation is initialized using the confluent
       * hypergeometric function.
       *
       * \f[
       *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2*a*r*x_n)
       * \f]
       * using the representation of Modified Bessel function as 1F1
       * \f[
       *  i_l(x) = \exp{-x} \frac{\sqrt{\pi}}{4\Gamma{1.5 + n}}
       *            (\frac{x}{2})^{l} 1F1(l+1, 2l+2, 2x)
       * \f]
       *
       * @param distance (single) value of distance
       *                (second exponential argument) in the equation above
       * @param fac_a  Scaling factor for the exponential arguments
       * @param n_rows number of rows where the recursion is applicable
       *               from the top
       */
      void downward_recursion(double distance, double fac_a, int n_rows);

      Eigen::ArrayXXd bessel_values{};
      Eigen::ArrayXd bessel_arg{};
      Eigen::ArrayXd bessel_arg_i{};
      Eigen::ArrayXd bessel_arg_pow{};
      Eigen::ArrayXd exp_bessel_arg{};
      Eigen::ArrayXd x_v{};
      Eigen::ArrayXd efac{};
      std::vector<Hyp1f1> hyp1f1s{};
      Eigen::ArrayXd igammas{};

      int order_max{};
      size_t l_max{};
      int n_max{};
    };

  }  // namespace math
}  // namespace rascal

#endif  // SRC_RASCAL_MATH_BESSEL_HH_
