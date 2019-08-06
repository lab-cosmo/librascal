/**
 * file   bessel.hh
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

#ifndef SRC_MATH_BESSEL_HH_
#define SRC_MATH_BESSEL_HH_

#include "math/math_utils.hh"
#include "math/hyp1f1.hh"

namespace rascal {
  namespace math {
    /**
     * Compute the modified spherical Bessel function of the first kind
     *
     * This function, hereafter denoted MBSF, is notated as i_n(x).  Here,
     * n must be an integer and x must be a real number > 0.  All orders n
     * are computed and returned, up to a specified maximum order.
     *
     * @param x     Argument of the MBSF
     * @param l_max Maximum order to compute
     *
     * @return Eigen::Array (1-D) of Bessel function values
     *
     * Since the MSBF blows up exponentially with increasing argument (x),
     * what is returned here is e^{-x}*i_n(x), which is bounded (decays to
     * leading order as 1/x).  If you want i_n(x) itself, multiply the
     * result by e^x (at your own peril).
     *
     * The recursion relation used here is:
     * \[
     *    i_0(x) = sinh(x) / x
     *    i_1(x) = (x*cosh(x) + sinh(x)) / x^2
     *    i_n(x) = i_{n-2}(x) - (2n - 1)/x * i_{n-1}(x)
     * \]
     * from http://mathworld.wolfram.com/ModifiedSphericalBesselFunctionoftheFirstKind.html // NOLINT
     *
     */
    inline Eigen::ArrayXd bessel_i_exp_allorders(double x, size_t order_max) {
      if (order_max < 1) {
        order_max = 1;
      }
      Eigen::ArrayXd function_values(order_max);
      // Note: Since x is not likely to be close to zero, there is no
      // signficant improvement in accuracy from using std::expm1 instead
      function_values(0) = (1. - std::exp(-2. * x)) / (2. * x);
      function_values(1) =
          ((x - 1.) + std::exp(-2. * x) * (x + 1.)) / (2. * x * x);
      for (size_t order{2}; order < order_max; ++order) {
        function_values(order) =
            - function_values(order - 1) * (2. * order - 1.) / x
            + function_values(order - 2);
      }
      return function_values;
    }

    /**
     * Compute the modified spherical Bessel function of the first kind
     *
     * Vectorized version for multiple arguments (x_v)
     *
     * @param x_v   Eigen::Array of MBSF arguments
     * @param order_max Maximum order to compute
     *
     * @return Eigen::Array (2-D) of Bessel function values.  The different
     *         arguments (x_v) go along the rows, while the column indexes
     *         the orders (n-values).
     */
    inline Eigen::ArrayXXd
    bessel_i_exp_allorders(const Eigen::Ref<const Eigen::ArrayXd> & x_v,
                           size_t order_max) {
      if (order_max < 1) {
        order_max = 1;
      }
      Eigen::ArrayXXd function_values(x_v.size(), order_max);
      function_values.col(0) = (1. - Eigen::exp(-2. * x_v)) / (2. * x_v);
      function_values.col(1) =
          ((x_v - 1.) + Eigen::exp(-2. * x_v) * (x_v + 1.)) /
          (2. * x_v.square());
      for (size_t order{2}; order < order_max; ++order) {
        function_values.col(order) =
            - function_values.col(order - 1) * (2. * order - 1.) / x_v
            + function_values.col(order - 2);
      }
      return function_values;
    }

    inline Eigen::ArrayXXd
    bessel_i_allorders(const Eigen::Ref<const Eigen::ArrayXd> & x_v,
                           int order_max) {
      if (order_max < 1) {
        order_max = 1;
      }

      Eigen::ArrayXXd function_values(x_v.size(), order_max);
      Eigen::ArrayXd epx_x_v = Eigen::exp(-x_v);
      for (int order{0}; order < order_max; ++order) {
        Hyp1f1 hyp1f1{static_cast<double>(order+1), static_cast<double>(2*order+2)};
        double igamma{1/std::tgamma(1.5+order)};
        for (int ii{0}; ii < x_v.size(); ++ii) {
          function_values(ii, order) = epx_x_v[ii] * igamma * pow(x_v[ii] * 0.5, order) * 0.5 * std::sqrt(PI) * hyp1f1.calc(2*x_v[ii]);
        }

      }
      return function_values;
    }
    /**
     * Optimized MBSFs times two exponentials that complete the square
     *
     * This is for the lucky case that we're computing something of the form
     * \[
     *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2ar*x_n)
     * \]
     * where the exponentials complete the square of the cosh and sinh
     * arguments that are used to build the (modified) Bessel functions.
     *
     * @param x_v      Eigen::Array of x-values (part of the argument of the
     *                 first exponential; see equation above)
     * @param r        (single) value of r (second exponential argument)
     *                 in the equation above
     * @param a_scale  Scaling factor for the exponential arguments
     * @param l_max    Maximum order of MBSF to compute
     *
     * @return     Eigen::Array (2-D) of function values.  The x_v
     *             (ultimately, radial indices) correspond to rows, while the
     *             columns correspond to MBSF orders (ultimately l-channels).
     */
    inline Eigen::ArrayXXd bessel_i_exp_exp_complete_square(
        const Eigen::Ref<const Eigen::ArrayXd> & x_v, const double r,
        const double a_scale, size_t order_max) {
      using Prec_t = long double;
      if (order_max < 1) {
        order_max = 1;
      }
      Eigen::Array<Prec_t, Eigen::Dynamic, Eigen::Dynamic> function_values(x_v.size(), order_max);
      Eigen::Array<Prec_t, Eigen::Dynamic, 1> bessel_arg(x_v.size());
      Eigen::Array<Prec_t, Eigen::Dynamic, 1> epx_p_arg(x_v.size());
      Eigen::Array<Prec_t, Eigen::Dynamic, 1> epx_m_arg(x_v.size());
      Prec_t fff{static_cast<Prec_t>(2. * a_scale * r)};
      for (int ii{0}; ii < x_v.size(); ++ii) {
        bessel_arg[ii] = static_cast<Prec_t>(1.) / (fff*static_cast<Prec_t>(x_v[ii]));
        epx_p_arg[ii] = std::exp(static_cast<Prec_t>(-a_scale * (x_v[ii] + r) * (x_v[ii] + r)));
        epx_m_arg[ii] = std::exp(static_cast<Prec_t>(-a_scale * (x_v[ii] - r) * (x_v[ii] - r)));
      }
      // bessel_arg = static_cast<Prec_t>(2. * a_scale * r) * x_v.cast<Prec_t>();
      // bessel_arg = bessel_arg.inverse();
      // TODO(max) is it faster to allocate arrays for the exp results, or
      //           does the cost of allocating memory outweigh the savings
      //           of evaluating it again _once_?
      // i_0(z) = sinh(z)/z
      Prec_t fac{0.5};
      // function_values.col(0) = (Eigen::exp(-a_scale * (x_v - r).square()) -
      //                           Eigen::exp(-a_scale * (x_v + r).square())).cast<Prec_t>() *
      //                          (fac*bessel_arg);
      function_values.col(0) = (epx_m_arg - epx_p_arg) * fac * bessel_arg;
      // i_1(z) = cosh(z)/z - i_0(z)/z
      // function_values.col(1) = ((Eigen::exp(-a_scale * (x_v - r).square()) +
      //                            Eigen::exp(-a_scale * (x_v + r).square())).cast<Prec_t>() *
      //                           (fac*bessel_arg)) -
      //                          (function_values.col(0) * bessel_arg);
      function_values.col(1) = (epx_m_arg + epx_p_arg) * fac * bessel_arg -
                               (function_values.col(0) * bessel_arg);
      for (size_t order{2}; order < order_max; ++order) {
        function_values.col(order) = - function_values.col(order - 1) *
                              static_cast<Prec_t>(2. * order - 1.) * bessel_arg
                                       + function_values.col(order - 2);
      }
      return function_values.cast<double>();
    }

    /**
     * Computes of the MSBFs as
     * \[
     *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2*a*r*x_n)
     * \]
     * Just precompte() once, calc(), and use get_values() to access
     * the values
     *
     * Note that this uses the "exp-exp complete-square" optimized function,
     * where the Bessel function is multiplied by two exponentials that keep
     * it from blowing up.
     */
    class ModifiedSphericalBessel {
     public:
      using ref = Eigen::Ref<const Eigen::ArrayXXd>;

      ModifiedSphericalBessel() = default;

      /**
       * Initialize arrays for the computation of
       * @param x_v      Eigen::Array of x-values (part of the argument of the
       *                 first exponential; see equation above)
       */
      void precompute(const size_t & l_max, const Eigen::Ref<const Eigen::VectorXd> & x_v) {
        this->x_v = x_v.array();
        this->n_max = x_v.size();
        this->bessel_values.resize(this->n_max, l_max + 1);
        this->bessel_arg.resize(this->n_max);
        this->exp_bessel_arg.resize(this->n_max);
        this->bessel_arg_pow.resize(this->n_max);
        this->efac.resize(this->n_max);
        this->l_max = l_max;
        this->order_max = static_cast<int>(this->l_max+1);
        this->igammas.resize(this->order_max);
        for (int order{0}; order < this->order_max; ++order) {
          this->hyp1f1s.emplace_back(static_cast<double>(order+1), static_cast<double>(2*order+2));
          this->igammas[order] = 1./std::tgamma(1.5+order);
        }
      }

      /**
       * Optimized MBSFs times two exponentials that complete the square
       *
       * \[
       *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2*a*r*x_n)
       * \]
       * using the representation of Modified Bessel function as 1F1
       * \[
       *  i_l(x) = \exp{-x} \frac{\sqrt{\pi}}{4\Gamma{1.5 + n}}
       *            (\frac{x}{2})^{l} 1F1(l+1, 2l+2, 2x)
       * \]
       *
       * @param r        (single) value of r (second exponential argument)
       *                 in the equation above
       * @param a_scale  Scaling factor for the exponential arguments
       */
      inline void
      bessel_i_allorders_complete_square(const double & r, const double & a_scale) {
        this->bessel_arg = 2. * a_scale * r * this->x_v;
        this->exp_bessel_arg = Eigen::exp(-this->bessel_arg);
        this->bessel_arg_pow = Eigen::ArrayXd::Constant(this->n_max, 1.);
        for (int order{0}; order < this->order_max; ++order) {
          auto & hyp1f1{this->hyp1f1s[order]};
          for (int ii{0}; ii < this->n_max; ++ii) {
            this->bessel_values(ii, order) = this->exp_bessel_arg[ii] * this->igammas[order] * this->bessel_arg_pow[ii] * 0.5 * math::SQRT_PI * hyp1f1.calc(2.*this->bessel_arg[ii]);
          }
          this->bessel_arg_pow *= 0.5 * this->bessel_arg;
        }

        this->efac = std::exp(-a_scale*r*r) * Eigen::exp(-a_scale*this->x_v.square());
        for (int order{0}; order < this->order_max; ++order) {
          this->bessel_values.col(order) *= this->efac;
        }
      }

      /**
       * Optimized MBSFs times two exponentials that complete the square
       *
       * This is for the lucky case that we're computing something of the form
       * \[
       *    f(r; x_n, a) = e^{-ar^2} e^{-ax_n^2} i_l(2*a*r*x_n)
       * \]
       * where the exponentials complete the square of the cosh and sinh
       * arguments that are used to build the (modified) Bessel functions.
       *
       * @param r        (single) value of r (second exponential argument)
       *                 in the equation above
       * @param a_scale  Scaling factor for the exponential arguments
       */
      void bessel_i_recursive_complete_square(const double & r, const double & a_scale) {

        this->bessel_arg = 2. * a_scale * r * this->x_v;
        this->bessel_arg = this->bessel_arg.inverse();
        // i_0(z) = sinh(z) / z
        this->bessel_values.col(0) =
                                (Eigen::exp(-a_scale * (x_v - r).square()) -
                                Eigen::exp(-a_scale * (x_v + r).square())) *
                                0.5 * this->bessel_arg;
        // i_1(z) = cosh(z)/z - i_0(z)/z
        this->bessel_values.col(1) =
                                ((Eigen::exp(-a_scale * (x_v - r).square()) +
                                  Eigen::exp(-a_scale * (x_v + r).square())) *
                                  0.5*bessel_arg)
                                - this->bessel_values.col(0) * bessel_arg;

        for (int order{2}; order < this->order_max; ++order) {
          this->bessel_values.col(order) = this->bessel_values.col(order - 2)
                                - this->bessel_values.col(order - 1) *
                                  (2. * order - 1.) * bessel_arg;
        }
      }

      /**
       * Compute all the MBSFs for the given x-values up to the given order
       *
       * The smallest x-value should not be too small (x > 0.05), the results
       * become inacurate between very small values (<~1e-30). For our purpose
       * such errors are irelevant since the radial integral does not have
       * additional factors >~1e10.
       */
      inline void calc(const double & distance, const double & fac_a) {
        double fac{distance*distance*fac_a};
        // the criteria for fac has be found empirically
        if (fac > 40) {
          this->bessel_i_recursive_complete_square(distance, fac_a);
        } else {
          this->bessel_i_allorders_complete_square(distance, fac_a);
        }
      }

      /**
       * Return a reference to the precomputed Bessel function values
       *
       * @return Eigen::Array (2-D) of Bessel function values.  The different
       *         arguments (xs) go along the rows, while the column indexes
       *         the orders (n-values).
       *         Note that a reference is returned to avoid unnecessary
       *         copies.
       */
      auto get_values() { return ref(bessel_values); }

      Eigen::ArrayXXd bessel_values{};

      Eigen::ArrayXd bessel_arg{};
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

#endif  // SRC_MATH_BESSEL_HH_
