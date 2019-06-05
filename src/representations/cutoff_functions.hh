/**
 * file   cutoff_functions.hh
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

#ifndef SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
#define SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_

#include "rascal_utility.hh"
#include "math/math_utils.hh"

#include <vector>
#include <memory>

#include <Eigen/Dense>

namespace rascal {

  namespace internal {

    /**
     * List of implemented cutoff function
     */
    enum class CutoffFunctionType { Cosine, RadialScaling, End_ };

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

      using Hypers_t = RepresentationManagerBase::Hypers_t;

      //! Pure Virtual Function to set hyperparameters of the cutoff function
      virtual void set_hyperparameters(const Hypers_t &) = 0;

      // TODO(felix) test is having these pure virtual changes performance
      //! Pure Virtual Function to evaluate the cutoff function
      // virtual double f_c(const double& distance) = 0;
      //! Pure Virtual Function to evaluate the derivative of the cutoff
      //! function with respect to the distance
      // virtual double df_c(const double& distance) = 0;
    };

    // Empty general template class implementing the cutoff functions
    // It should never be instantiated.
    template <CutoffFunctionType Type>
    struct CutoffFunction : CutoffFunctionBase {};

    template <>
    struct CutoffFunction<internal::CutoffFunctionType::Cosine>
        : CutoffFunctionBase {
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
      }

      void set_hyperparameters(const Hypers_t & hypers) {
        this->cutoff = hypers.at("cutoff").at("value").get<double>();
        this->smooth_width =
            hypers.at("smooth_width").at("value").get<double>();
      }

      inline double f_c(const double & distance) {
        return math::switching_function_cosine(distance, this->cutoff,
                                               this->smooth_width);
      }

      inline double df_c(const double & distance) {
        return math::derivative_switching_funtion_cosine(distance, this->cutoff,
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
     * with the cosine switching function.
     */
    template <>
    struct CutoffFunction<internal::CutoffFunctionType::RadialScaling>
        : CutoffFunctionBase {
      using Hypers_t = CutoffFunctionBase::Hypers_t;
      explicit CutoffFunction(const Hypers_t & hypers) {
        this->set_hyperparameters(hypers);
      }

      void set_hyperparameters(const Hypers_t & hypers) {
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

      inline double f_c(const double & distance) {
        double factor{0.};
        if (this->rate < math::dbl_ftol) {
          factor = math::pow(distance / this->scale, -this->exponent);
        } else if (this->exponent == 0) {
          factor = 1.;
        } else {
          factor = this->rate / (this->rate + math::pow(distance / this->scale,
                                                        this->exponent));
        }
        return factor * math::switching_function_cosine(distance, this->cutoff,
                                                        this->smooth_width);
      }

      inline double df_c(const double & distance) {
        double factor{0.};
        if (this->rate < math::dbl_ftol) {
          factor = -this->exponent / distance *
                   math::pow(distance / this->scale, -this->exponent);
        } else if (this->exponent == 0) {
          factor = 0.;
        } else {
          double ff{math::pow(distance / this->scale, this->exponent)};
          factor = this->rate * this->exponent * ff / distance *
                   math::pow(this->rate + ff, -2);
        }
        return factor * math::derivative_switching_funtion_cosine(
                            distance, this->cutoff, this->smooth_width);
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
  decltype(auto) make_cutoff_function(const Hypers & fc_hypers) {
    return std::static_pointer_cast<internal::CutoffFunctionBase>(
        std::make_shared<internal::CutoffFunction<Type>>(fc_hypers));
  }

  template <internal::CutoffFunctionType Type>
  decltype(auto) downcast_cutoff_function(
      std::shared_ptr<internal::CutoffFunctionBase> & cutoff_function) {
    return std::static_pointer_cast<internal::CutoffFunction<Type>>(
        cutoff_function);
  }

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CUTOFF_FUNCTIONS_HH_
