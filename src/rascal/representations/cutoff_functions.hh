/**
 * @file   rascal/representations/cutoff_functions.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   10 May 2019
 *
 * @brief  bundle the cutoff functions for simple use in the representations
 *
 * Copyright © 2019 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
#define SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_

#include "rascal/math/utils.hh"
#include "rascal/representations/calculator_base.hh"
#include "rascal/utils/utils.hh"

#include <Eigen/Dense>

#include <memory>
#include <vector>

namespace rascal {

  namespace internal {
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
    inline double switching_function_cosine(double r, double cutoff,
                                            double smooth_width) {
      if (r <= (cutoff - smooth_width)) {
        return 1.0;
      } else if (r > cutoff) {
        return 0.0;
      }
      double r_scaled{math::PI * (r - cutoff + smooth_width) / smooth_width};
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
    inline double derivative_switching_funtion_cosine(double r, double cutoff,
                                                      double smooth_width) {
      if (r <= (cutoff - smooth_width)) {
        return 0.0;
      } else if (r > cutoff) {
        return 0.0;
      }
      double r_scaled{math::PI * (r - cutoff + smooth_width) / smooth_width};
      return (-0.5 * math::PI / smooth_width * std::sin(r_scaled));
    }

    /**
     * List of implemented cutoff function
     */
    enum class CutoffFunctionType { ShiftedCosine, RadialScaling };

    struct CutoffFunctionBase {
      //! Constructor
      CutoffFunctionBase() = default;
      //! Destructor
      virtual ~CutoffFunctionBase() = default;
      //! Copy constructor
      CutoffFunctionBase(const CutoffFunctionBase & other) = delete;
      //! Move constructor
      CutoffFunctionBase(CutoffFunctionBase && other) = default;
      //! Copy assignment operator
      CutoffFunctionBase & operator=(const CutoffFunctionBase & other) = delete;
      //! Move assignment operator
      CutoffFunctionBase & operator=(CutoffFunctionBase && other) = default;

      using Hypers_t = CalculatorBase::Hypers_t;

      //! Pure Virtual Function to set hyperparameters of the cutoff function
      virtual void set_hyperparameters(const Hypers_t &) = 0;

      // TODO(felix) having these as pure virtual changes the performance of the
      // tests
      //! Pure Virtual Function to evaluate the cutoff function
      // virtual double f_c(double distance) = 0;
      //! Pure Virtual Function to evaluate the derivative of the cutoff
      //! function with respect to the distance
      // virtual double df_c(double distance) = 0;
    };

    // Empty general template class implementing the cutoff functions
    // It should never be instantiated.
    template <CutoffFunctionType Type>
    struct CutoffFunction : CutoffFunctionBase {};

    template <>
    struct CutoffFunction<internal::CutoffFunctionType::ShiftedCosine>
        : CutoffFunctionBase {
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
      }

      void set_hyperparameters(const Hypers_t & hypers) override {
        this->cutoff = hypers.at("cutoff").at("value").get<double>();
        this->smooth_width =
            hypers.at("smooth_width").at("value").get<double>();
      }

      double f_c(double distance) {
        return switching_function_cosine(distance, this->cutoff,
                                         this->smooth_width);
      }

      double df_c(double distance) {
        return derivative_switching_funtion_cosine(distance, this->cutoff,
                                                   this->smooth_width);
      }

      //! keep the hypers
      Hypers_t hypers{};
      //! cutoff radii
      double cutoff{0.};
      //! interval into which the smoothing happens [cutoff-smooth_width,cutoff]
      double smooth_width{0.};
    };

    /**
     * Computes the RadialScaling switching function as expressed in equation
     * 21 of Willatt, M. J., Musil, F., & Ceriotti, M. (2018).
     * https://doi.org/10.1039/c8cp05921g
     *
     *        ╭ 1 / (r/r_0)^m, if c == 0
     *        |
     * u(r) = ┤ 1, if m == 0
     *        |
     *        ╰ c / (c + (r/r_0)^m), else
     *
     * c -> rate
     * r_0 -> scale
     * m -> exponent
     *
     * multiplied by the cosine switching function defined in
     * math::switching_function_cosine() (which comes with additional parameters
     * `cutoff` and `smooth_width`).
     *
     * Typically c == 1, r_0 > 0 and m is a positive integer.
     *
     * Derivatives for the radial scaling component are:
     *
     *         ╭ -m /( (r/r_0)^m * r), if c == 0
     *         |
     * u'(r) = ┤ 0, if m == 0
     *         |
     *         ╰ -m c (r/r_0)^m / (r * (c + (r/r_0)^m)^2), otherwise
     *
     * These are combined with the cosine switching function derivatives using
     * the Leibniz product rule to get the derivative of the final radial
     * modulation function.
     */
    template <>
    struct CutoffFunction<internal::CutoffFunctionType::RadialScaling>
        : CutoffFunctionBase {
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
      }

      void set_hyperparameters(const Hypers_t & hypers) override {
        this->cutoff = hypers.at("cutoff").at("value").get<double>();
        this->smooth_width =
            hypers.at("smooth_width").at("value").get<double>();
        this->rate = hypers.at("rate").at("value").get<double>();
        if (this->rate < 0) {
          throw std::runtime_error("RadialScaling's rate should be positive");
        }
        this->exponent = hypers.at("exponent").at("value").get<int>();
        this->scale = hypers.at("scale").at("value").get<double>();
      }

      double value(double distance) {
        double factor{0.};
        if (std::abs(this->rate) <= math::DBL_FTOL) {
          factor = 1. / math::pow(distance / this->scale, this->exponent);
        } else if (this->exponent == 0) {
          factor = 1.;
        } else {
          factor = this->rate / (this->rate + math::pow(distance / this->scale,
                                                        this->exponent));
        }
        return factor;
      }

      double grad(double distance) {
        double factor{0.};
        if (std::abs(this->rate) <= math::DBL_FTOL) {
          factor = -this->exponent / distance /
                   math::pow(distance / this->scale, this->exponent);
        } else if (this->exponent == 0) {
          factor = 0.;
        } else {
          double ff{math::pow(distance / this->scale, this->exponent)};
          factor = -this->rate * this->exponent * ff / distance /
                   math::pow(this->rate + ff, 2_size_t);
        }
        return factor;
      }

      double f_c(double distance) {
        return this->value(distance) *
               switching_function_cosine(distance, this->cutoff,
                                         this->smooth_width);
      }

      double df_c(double distance) {
        double df_c1{this->grad(distance) *
                     switching_function_cosine(distance, this->cutoff,
                                               this->smooth_width)};
        double df_c2{this->value(distance) *
                     derivative_switching_funtion_cosine(distance, this->cutoff,
                                                         this->smooth_width)};
        return df_c1 + df_c2;
      }

      //! keep the hypers
      Hypers_t hypers{};
      //! cutoff radii
      double cutoff{0.};
      //! interval into which the smoothing happens [cutoff-smooth_width,cutoff]
      double smooth_width{0.};
      //! rate c
      double rate{0.};
      //! exponent m
      int exponent{0};
      //! scale r_0
      double scale{1.};
    };

  }  // namespace internal

  template <internal::CutoffFunctionType Type, class Hypers>
  auto make_cutoff_function(const Hypers & fc_hypers) {
    return std::static_pointer_cast<internal::CutoffFunctionBase>(
        std::make_shared<internal::CutoffFunction<Type>>(fc_hypers));
  }

  template <internal::CutoffFunctionType Type>
  auto downcast_cutoff_function(
      std::shared_ptr<internal::CutoffFunctionBase> & cutoff_function) {
    return std::static_pointer_cast<internal::CutoffFunction<Type>>(
        cutoff_function);
  }

}  // namespace rascal

#endif  // SRC_RASCAL_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
