/**
 * file   math_utils.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 * @author  Max Veit <max.veit@epfl.ch>
 *
 * @date   14 October 2018
 *
 * @brief contains the implementation of miscellaneous math functions
 *
 * Copyright  2018  Felix Musil, Max Veit, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_MATH_MATH_UTILS_HH_
#define SRC_MATH_MATH_UTILS_HH_

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <cstdint>

namespace rascal {
  /**
   * User defined literal operator to initialize size_t literals
   *
   * linter wants std::uint64_t to avoid the C style type but litterals are
   * of type unsigned long long while std::uint64_t on 64bits machines stands
   * for unsigned long. Both have equivalent underlying storage but not the
   * same type...
   */
  constexpr std::size_t operator"" _n(unsigned long long int n) {  // NOLINT
    return n;
  }

  namespace math {

    // Reminder: C++ floating-point literals are automatically of type double
    /// Pi to more digits than anyone could possibly need
    const double PI = 3.14159265358979323846264338327950288419716939937510;
    const double SQRT_TWO = std::sqrt(2.0);
    const double INV_SQRT_TWO = std::sqrt(0.5);
    const double SQRT_THREE = std::sqrt(3.0);

    /// How small a number must be to be considered effectively zero
    const double dbl_ftol = 100.0 * std::numeric_limits<double>::epsilon();

    /// How large a number must be to be considered infinity
    const double DOVERFLOW = std::numeric_limits<double>::infinity() / 100.;

    // define some usefull matrix type
    using Matrix_t =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Vector_t = Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>;

    using MatrixX2_t = Eigen::Matrix<double, Eigen::Dynamic, 2>;
    using Matrix_Ref = typename Eigen::Ref<const Matrix_t>;
    using MatrixX2_Ref = typename Eigen::Ref<const MatrixX2_t>;
    using Vector_Ref = typename Eigen::Ref<const Vector_t>;
    /**
     * Define integer powers and wrap the different cases under the same name
     */
    namespace details {
      //! unsingned integer power
      template <typename Scalar_>
      inline Scalar_ pow_u(Scalar_ x, size_t n) {
        Scalar_ value{1};

        /* repeated squaring method
         * returns 0.0^0 = 1.0, so continuous in x
         * (from GSL)
         */
        do {
          if (n & 1)
            value *= x; /* for n odd */
          n >>= 1;
          x *= x;
        } while (n);

        return value;
      }

      //! integer power
      template <typename Scalar_>
      inline double pow_i(const Scalar_ & x, const int & n) {
        size_t un{0};
        double value{static_cast<double>(x)};

        if (n < 0) {
          value = 1.0 / x;
          un = static_cast<size_t>(-n);
        } else {
          un = static_cast<size_t>(n);
        }

        return pow_u(value, un);
      }
    }  // namespace details

    //! integer power
    inline double pow(const double & x, const int & n) {
      return details::pow_i(x, n);
    }

    //! integer power
    inline double pow(const int & x, const int & n) {
      return details::pow_i(x, n);
    }

    //! integer power
    inline double pow(const size_t & x, const int & n) {
      return details::pow_i(x, n);
    }

    //! unsingned integer power
    inline double pow(const double & x, const std::size_t & n) {
      return details::pow_u(x, n);
    }

    //! unsingned integer power
    inline int pow(const int & x, const std::size_t & n) {
      return details::pow_u(x, n);
    }

    //! unsingned integer power
    inline size_t pow(const size_t & x, const std::size_t & n) {
      return details::pow_u(x, n);
    }

    //! general power
    inline double pow(const double & x, const double & n) {
      return std::pow(x, n);
    }

    /**
     * Defines an integer power functor, which is used as CustomUnaryOp object
     * to apply elementwise integer power operations on an Eigen::MatrixBase
     * object e.g. :
     *
     *   Eigen::VectorXd vec(3);
     *   vec << 5,10,20;
     *   auto functor{MakePositiveIntegerPower(2)};
     *   vec.unaryExpr(functor);
     *   // result is vec = [25,100,400]
     *
     * @tparam Scalar can be of type: double, int and uint
     */
    template <typename Scalar>
    struct MakePositiveIntegerPower {
      typedef Scalar result_type;
      size_t b;
      explicit MakePositiveIntegerPower(const size_t & b) : b{b} {}

      Scalar operator()(const Scalar & a) const { return pow(a, this->b); }
    };

    /**
     * Compute a cosine-type switching function for smooth cutoffs
     *
     * @param cutoff Outer (strict) cutoff, beyond which this function becomes
     *               zero
     *
     * @param smooth_width Width over which the smoothing function extends;
     *                     the function becomes one for r less than
     *                     cutoff - smooth_width
     *
     * @param r Distance at which to evaluate the switching function
     *
     * The functional form is:
     *
     * sw(r) = 1/2 + 1/2 cos(pi * (r - cutoff + smooth_width) / smooth_width)
     *
     * if r is within the cutoff region (cutoff - smooth_width < r <= cutoff);
     * if r is outside (> cutoff) the function is zero; if r is inside, the
     * function is 1.
     *
     * Specifying smooth_width less than cutoff is not an error.
     * If smooth_width is equal to zero the result will just be a step
     * function.
     *
     */
    inline double switching_function_cosine(const double & r,
                                            const double & cutoff,
                                            const double & smooth_width) {
      if (r <= (cutoff - smooth_width)) {
        return 1.0;
      } else if (r > cutoff) {
        return 0.0;
      }
      double r_scaled{PI * (r - cutoff + smooth_width) / smooth_width};
      return (0.5 * (1. + std::cos(r_scaled)));
    }

    /**
     * Compute the derivative of the cosine-type switching function
     *
     * @param cutoff Outer (strict) cutoff, beyond which this function becomes
     *               zero
     *
     * @param smooth_width Width over which the smoothing function extends;
     *                     the function becomes one for r less than
     *                     cutoff - smooth_width
     *
     * @param r Distance at which to evaluate the derivative
     *
     * The functional form is:
     *
     * dsw/dr (r) = -pi/(2*smooth_width) * sin(pi * (r - cutoff + smooth_width)
     *                                              / smooth_width)
     */
    inline double
    derivative_switching_funtion_cosine(const double & r, const double & cutoff,
                                        const double & smooth_width) {
      if (r <= (cutoff - smooth_width)) {
        return 0.0;
      } else if (r > cutoff) {
        return 0.0;
      }
      double r_scaled{PI * (r - cutoff + smooth_width) / smooth_width};
      return (-0.5 * PI / smooth_width * std::sin(r_scaled));
    }

  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_MATH_UTILS_HH_
