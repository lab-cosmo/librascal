/**
 * @file   rascal/representations/cutoff_functions.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   17 Dec 2018
 *
 * @brief  bundle the cutoff functions for simple use in the representations
 *
 * Copyright © 2018 Till Junge, Markus Stricker, Max Veit, Felix Musil, COSMO
 * (EPFL), LAMMM (EPFL)
 *
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
#include "rascal/utils.hh"

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
    enum class CutoffFunctionType {
      Cosine,
      CosineShifted,
      RadialScaling,
      End_
    };

    /**
     * forward declaration for compute dispatch
     */
    template <CutoffFunctionType Type>
    class CutoffFunction;

    class CutoffFunctionBase {
     public:
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

      //! Main worker (raison d'être)
      template <class StructureManager>
      inline void compute(StructureManager & manager) const;

      /**
       * The identifier string must provide a unique name for a property to
       * store precomputed cutoff function values. For this, the naming scheme
       * must guarantee that two different parameter sets (e.g., cut-off radii)
       * generate diffent names, and should guarantee that the same function
       * with the same  parameters generates the same name.
       */
      virtual const std::string & get_identifier() const = 0;

     protected:
      //! Main worker (raison d'être)
      template <CutoffFunctionType CutFunType, class StructureManager>
      inline void compute_helper(StructureManager & manager) const;
    };

    /* ---------------------------------------------------------------------- */
    template <class CutoffFunctionType>
    class CutoffFunctionComputer : public CutoffFunctionBase {
     public:
      template <class StructureManager>
      void compute(StructureManager & manager) const {
        auto & typed_this{static_cast<CutoffFunctionType &>(*this)};
        if (not manager.has_property()) {
          // still to be filled
        }
      }
    };

    // Empty general template class implementing the cutoff functions
    // It should never be instantiated.
    template <CutoffFunctionType Type>
    class CutoffFunction : public CutoffFunctionComputer<CutoffFunction<Type>> {
    };

    /**
     * Cosine cutoff function as in Behler, can only be used with strict
     * managers
     */
    template <>
    class CutoffFunction<CutoffFunctionType::Cosine>
        : public CutoffFunctionComputer<
              CutoffFunction<CutoffFunctionType::Cosine>> {
     public:
      using Hypers_t = typename CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers)
          : hypers{hypers},
            cutoff{hypers.at("cutoff").at("value").get<double>()},
            identifier{this->make_identifier()} {}

      inline double f_c(const double & distance) const {
        assert(distance <= this->cutoff);
        return .5 * (std::cos(math::PI * distance / this->cutoff) + 1.);
      }

      inline double df_c(const double & distance) {
        assert(distance <= this->cutoff);
        auto && scaled_dist{math::PI * distance / this->cutoff};
        return -.5 * scaled_dist * std::sin(scaled_dist);
      }

      const std::string & get_identifier() const { return this->identifier; }

     protected:
      std::string make_identifier() {
        std::stringstream id{};
        id.precision(14);
        id << "Cosine_" << this->cutoff;
        return id.str();
      }
      //! keep the hypers
      const Hypers_t hypers;
      //! cutoff radii
      const double cutoff;
      const std::string identifier;
    };

    template <>
    class CutoffFunction<internal::CutoffFunctionType::CosineShifted>
        : public CutoffFunctionBase {
     public:
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers)
          : hypers{hypers},
            cutoff{hypers.at("cutoff").at("value").get<double>()},
            smooth_width{hypers.at("smooth_width").at("value").get<double>()},
            identifier{this->make_identifier()} {}

      double f_c(double distance) {
        return switching_function_cosine(distance, this->cutoff,
                                         this->smooth_width);
      }

      double df_c(double distance) {
        return derivative_switching_funtion_cosine(distance, this->cutoff,
                                                   this->smooth_width);
      }

      const std::string & get_identifier() const { return this->identifier; }

     private:
      std::string make_identifier() const {
        std::stringstream id{};
        id.precision(14);
        id << "CosineShifted_" << this->cutoff << "_" << this->smooth_width;
        return id.str();
      }
      //! keep the hypers
      Hypers_t hypers;
      //! cutoff radii
      double cutoff;
      //! interval into which the smoothing happens [cutoff-smooth_width,cutoff]
      double smooth_width;
      std::string identifier;
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
    class CutoffFunction<CutoffFunctionType::RadialScaling>
        : public CutoffFunctionComputer<
              CutoffFunction<CutoffFunctionType::RadialScaling>> {
     public:
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers)
          : hypers{hypers},
            cutoff{hypers.at("cutoff").at("value").get<double>()},
            smooth_width{hypers.at("smooth_width").at("value").get<double>()},
            rate{hypers.at("rate").at("value").get<double>()},
            exponent{hypers.at("exponent").at("value").get<int>()},
            scale{hypers.at("scale").at("value").get<double>()},
            identifier{this->make_identifier()} {
        if (this->rate < 0) {
          throw std::runtime_error("RadialScaling's rate should be positive");
        }
      }

      double value(double distance) {
        double factor{0.};
        if (this->rate > math::DBL_FTOL) {
          factor = this->rate / (this->rate + math::pow(distance / this->scale,
                                                        this->exponent));
        } else if (this->exponent == 0) {
          factor = 1.;
        } else {
          factor = 1. / math::pow(distance / this->scale, this->exponent);
        }
        return factor;
      }

      double grad(double distance) {
        double factor{0.};
        if (this->rate > math::DBL_FTOL) {
          double ff{math::pow(distance / this->scale, this->exponent)};
          factor = -this->rate * this->exponent * ff / distance /
                   math::pow(this->rate + ff, 2_size_t);
        } else if (this->exponent == 0) {
          factor = 0.;
        } else {
          factor = -this->exponent / distance /
                   math::pow(distance / this->scale, this->exponent);
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

      template <class StructureManager>
      void compute(StructureManager & manager) const;

      const std::string & get_identifier() const { return this->identifier; }

     protected:
      std::string make_identifier() {
        std::stringstream id{};
        id.precision(14);
        id << "RadialScaling_" << this->cutoff << "_" << this->smooth_width
           << "_" << this->rate << "_" << this->exponent << "_" << this->scale;
        return id.str();
      }

      //! keep the hypers
      Hypers_t hypers;
      //! cutoff radii
      double cutoff;
      //! interval into which the smoothing happens [cutoff-smooth_width,cutoff]
      double smooth_width;
      //! rate c
      double rate;
      //! exponent m
      int exponent;
      //! scale r_0
      double scale;
      std::string identifier;
    };

    template <CutoffFunctionType CutFunType, class StructureManager>
    inline void
    CutoffFunctionBase::compute_helper(StructureManager & manager) const {
      switch (CutFunType) {
      case CutoffFunctionType::CosineShifted: {
        auto & cutoff_function{static_cast<
            const CutoffFunction<CutoffFunctionType::CosineShifted> &>(*this)};
        cutoff_function.compute(manager);
        break;
      }
      case CutoffFunctionType::RadialScaling: {
        auto & cutoff_function{static_cast<
            const CutoffFunction<CutoffFunctionType::RadialScaling> &>(*this)};
        cutoff_function.compute(manager);
        break;
      }
      default:
        throw std::runtime_error("Unknown cutoff function type");
        break;
      }
    }
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
