/**
 * file   hyp1f1.hh
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

#include "math/math_interface.hh"
#include "math/math_utils.hh"

#include  <map>

namespace rascal {
  namespace math {
    // using MatrixX2_t = Eigen::Matrix<double, Eigen::Dynamic, 2>;

    /**
     * Utilities to apply the recurrence relations of 1F1 from maximal
     * values of a and b to minimal values.
     * The recurrence path is chosen so that derivatives of 1F1 are also
     * a result.
     */
    inline double recurence_to_val_downward(const double& a, const double& b, const double& z, const double& M1p2p, const double& M1p1p) {
      return z * (a - b) * M1p2p / b /(b + 1) + M1p1p;
    }

    inline double recurence_to_der_downward(const double& a, const double& b, const double& z, const double& M2p3p, const double& M1p2p) {
      return z * (a + 1) * M2p3p / (b + 2) / (b + 1) + M1p2p;
    }

    // inline MatrixX2_t hyp1f1_recurrence(const int& max_angular, const int& n, const double& z, const double& f) {
    //   // first
    //   MatrixX2_t results(max_angular, 2);

    // }

    class Hyp1f1Series {
     protected:
      double a,b;
      size_t mmax;
      double prefac;
      double tolerance;
      bool is_exp{false};

      Eigen::VectorXd coeff{};
      Eigen::VectorXd coeff_derivative{};

     public:
      size_t n_terms{0};

      Hyp1f1Series(const double& a, const double& b, const size_t& mmax, const double& tolerance = 1e-14)
          : a{a}, b{b}, mmax{mmax},
            prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance} {

          if (std::abs(1 - this->b/this->a) < dbl_ftol) {
            this->is_exp = true;
          } else {
            coeff.resize(mmax);
            coeff_derivative.resize(mmax);
            double u{0.};
            // precomputes the (a)_i / (b)_i / i! coefficients
            coeff(0)=a/b;
            coeff_derivative(0) = (a+1) / (b+1);
            for (size_t i{1}; i < mmax; ++i) {
                u = (a + i) / ((i + 1) * (b + i));
                coeff(i) = coeff(i-1) * u;
                u = (a + i + 1) / ((i + 1) * (b + i + 1));
                coeff_derivative(i) = coeff_derivative(i-1) * u;
            }
          }
      }
      //! Computes G(a,b,z)
      inline double calc(const double& r_ij, const double& alpha, const double& beta, const bool& derivative) {
        using math::pow;
        double z{pow(alpha*r_ij, 2) / (alpha + beta)};
        return this->prefac*this->calc(z, derivative)*std::exp(-alpha*r_ij*r_ij);
      }

      inline double calc(const double& z, const bool& derivative = false) {
        if (not this->is_exp) {
          return this->hyp1f1(z, derivative);
        } else {
          return std::exp(z);
        }
      }
      //! Computes 1F1
      inline double hyp1f1(const double& z, const bool& derivative) {
        using math::pow;
        this->n_terms = this->mmax;
        Eigen::Ref<Eigen::VectorXd> coefficient(this->coeff);
        if (derivative) {
          coefficient = Eigen::Ref<Eigen::VectorXd>(this->coeff_derivative);
        }
        // perform the sum
        double res{1.0}, a1{1.0}, zpow{1.0};
        for (size_t i{0}; i < this->mmax; ++i) {
          zpow *= z;
          a1 = coefficient(i) * zpow;
          if (a1 < this->tolerance*res) {
            this->n_terms = i;
            break;
          }
          res += a1;

          if (res > DOVERFLOW){
            std::stringstream error{};
            error << "Hyp1f1Series series expansion: a="
                  << std::to_string(this->a) << " b="
                  << std::to_string(this->b) <<" z=" << std::to_string(z)
                  << std::endl;
            throw std::overflow_error(error.str());
          }
        }

        if (derivative) {
          res *= this->a / this->b;
        }

        return res;
      }
    };

    /**
     * Computes the 1F1 with the asymptotic limit
     *  1F1(a,b,z) \sim \frac{\exp{z} z^{a-b} \Gamma{b}}{\Gamma{a}}
     *                    \sum_{j=0}^{\infty} \frac{(b-a)_j(1-a)_j}{j!} z^{-j}
     *
     *  G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
     *                    * 1F1(a,b,z)
     */
    class Hyp1f1Asymptotic {
     private:
      double a,b,prefac;
      bool is_n_and_l{false};
      bool is_exp{false};
      double tolerance;
      size_t mmax;
      Eigen::VectorXd coeff{};
      Eigen::VectorXd coeff_derivative{};

      inline double z_power_a_b(const double& z) {
        using math::pow;
        double fac{0.};
        if (this->is_n_and_l){
          fac = std::sqrt(pow(z,static_cast<int>(2*(this->a-this->b))));
        } else {
          fac = pow(z,this->a-this->b);
        }
        return fac;
      }

     public:
      size_t n_terms{0};

      Hyp1f1Asymptotic(const double& a, const double& b, const size_t& mmax, const double& tolerance = 1e-14 )
        : a{a}, b{b}, prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance}, mmax{mmax} {
        double intpart;
        double f2{std::modf(2*(a-b), &intpart)};
        if (std::abs(f2) < 1e-14) {
          this->is_n_and_l = true;
        }

        if (std::abs(1 - this->b/this->a) < dbl_ftol) {
          this->is_exp = true;
        } else {
          coeff.resize(mmax);
          coeff_derivative.resize(mmax);
          // precomputes the (1-a)_i / (b-a)_i / i! coefficients for the
          // computation of hyp2f0
          double bma{b-a}, oma{1-a}, u{0.},oma1{1-a-1};
          coeff(0) = bma * oma;
          coeff_derivative(0) = bma * oma1;
          for (size_t i{1}; i < mmax; ++i) {
              u = (bma+i)*(oma+i)/(i+1);
              coeff(i) = coeff(i-1) * u;
              u = (bma+i)*(oma1+i)/(i+1);
              coeff_derivative(i) = coeff_derivative(i-1) * u;
          }
        }
      }

      //! Computes G(a,b,z)
      inline double calc(const double& r_ij, const double& alpha, const double& beta, const bool& derivative = false) {
        using math::pow;
        // argument of 1F1
        double z{pow(alpha*r_ij, 2) / (alpha + beta)};
        // simplification of the argument with exp(-alpha*r_ij^2)
        double z2{-alpha*beta*pow(r_ij, 2) / (alpha + beta)};

        auto fac{this->z_power_a_b(z)};

        return this->hyp2f0(z, derivative)*std::exp(z2)*fac;
      }

      //! Computes 1F1
      inline double calc(const double& z, const bool& derivative = false) {
        using math::pow;
        if (not this->is_exp) {
          auto fac{this->z_power_a_b(z)};
          double result{this->prefac*std::exp(std::log(this->hyp2f0(z, derivative))+z)*fac};
          if (std::isnan(result)) {result = DOVERFLOW;}
          return result;
        } else {
          return std::exp(z);
        }
      }

      //! computes hyp2f0 with arg1 = b-a and arg2 = 1-a arg3 = 1 / z
      inline double hyp2f0(const double& z, const bool& derivative) {
        using math::pow;
        this->n_terms = this->mmax;

        Eigen::Ref<Eigen::VectorXd> coefficient(this->coeff);
        if (derivative) {
          coefficient = Eigen::Ref<Eigen::VectorXd>(this->coeff_derivative);
        }

        double iz{1.0/z};
        double res{1.}, izpow{1.}, s_i{1.};
        // perform the sum
        // std::cout << "######" << std::endl;
        for (size_t i{0}; i < this->mmax; ++i) {
          izpow *= iz;
          s_i = coefficient(i) * izpow;
          if (res > 0 and std::fabs(s_i) < this->tolerance*res) {
            this->n_terms = i;
            break;
          }
          // std::cout << s_i << std::endl;
          res += s_i;
          if (res > DOVERFLOW){
            std::stringstream error{};
            error << "Hyp1f1Asymptotic expansion: a="
                  << std::to_string(this->a) << " b="
                  << std::to_string(this->b) <<" z=" << std::to_string(z)
                  << std::endl;
            throw std::overflow_error(error.str());
          }
        }
        // std::cout << "## "<<res << std::endl;
        return res;
      }
    };

    class Hyp1f1 {
     private:
      Hyp1f1Series hyp1f1_series;
      Hyp1f1Asymptotic hyp1f1_asymptotic;
      double a, b;
      double tolerance;
      size_t mmax;

      size_t nterms_a{0}, nterms_s{0};
      double h1f1_a{0}, h1f1_s{0};
      double z_above{0.}, z_below{0.};

      void update_switching_point() {
        this->h1f1_s = this->hyp1f1_series.calc(this->z_asympt);
        this->h1f1_a = this->hyp1f1_asymptotic.calc(this->z_asympt);

        this->nterms_s = this->hyp1f1_series.n_terms;
        this->nterms_a = this->hyp1f1_asymptotic.n_terms;
      }

      void find_switching_point() {
        int max_it{100};
        int i_it{0};
        this->update_switching_point();
        // brackets the switching point
        this->z_below = this->z_above = this->z_asympt;
        // std::cout << "# "<< std::fabs(1-this->h1f1_s/this->h1f1_a) << std::endl;
        if(std::fabs(1-this->h1f1_s/this->h1f1_a) > this->tolerance) {
          i_it = 0;
          while(std::fabs(1-this->h1f1_s/this->h1f1_a) > this->tolerance and i_it < max_it) {
            this->z_asympt *= 2.0;
            this->update_switching_point();
            i_it++;
            // std::cout << "### "<< std::fabs(1-this->h1f1_s/this->h1f1_a) << std::endl;
          }
          this->z_above = this->z_asympt;
        }
        else {
          i_it = 0;
          while (std::fabs(1-this->h1f1_s/this->h1f1_a) < this->tolerance and i_it < max_it) {
            this->z_asympt *= 0.5;
            this->update_switching_point();
            i_it++;
          }
          this->z_below = this->z_asympt;
        }
        /* and now bisects until we are reasonably close to an accurate
        determination */
        this->z_asympt = (this->z_above + this->z_below) * 0.5;
        this->update_switching_point();
        i_it = 0;
        while (this->z_above - this->z_below > this->tolerance and i_it < max_it) {
          if (std::abs(1-this->h1f1_s/this->h1f1_a) > this->tolerance) {
              this->z_below = this->z_asympt;
          } else {
              this->z_above = this->z_asympt;
          }
          this->z_asympt = (this->z_above + this->z_below) * 0.5;
          this->update_switching_point();
          i_it++;
        }
        this->z_asympt = this->z_asympt * 1.1;
        // std::cout << this->z_asympt << std::endl;
      }

     public:
      double z_asympt{1.};
      Hyp1f1(const double& a, const double& b, const size_t& mmax = 500, const double& tolerance = 1e-13)
        : hyp1f1_series{a,b,mmax,tolerance},hyp1f1_asymptotic{a,b,mmax,tolerance}, a{a}, b{b},tolerance{tolerance}, mmax{mmax} {

        /*now we try to determine what is the switching point between
        power series and asymptotic expansion. basically we choose
        the method that requires fewer terms for a chosen target accuracy.
        the asymptotic expansion tends to blow up at the switching point.*/
        if (std::abs(1 - this->b/this->a) > dbl_ftol) {
          this->find_switching_point();
        }
      }

      //! Compute 1F1(a,b,z)
      inline double calc(const double& z, const bool& derivative = false) {
        if (z > this->z_asympt) {
          return this->hyp1f1_asymptotic.calc(z, derivative);
        } else {
          return this->hyp1f1_series.calc(z, derivative);
        }
      }

      //! Compute G(a,b,z)
      inline double calc(const double& r_ij, const double& alpha, const double& beta, const bool& derivative = false) {
        double z{math::pow(alpha*r_ij, 2) / (alpha + beta)};
        if (z > this->z_asympt) {
          return this->hyp1f1_asymptotic.calc(r_ij, alpha, beta, derivative);
        } else {
          return this->hyp1f1_series.calc(r_ij, alpha, beta, derivative);
        }
      }
    };

    class Hyp1f1SphericalExpansion {
     protected:
      int max_angular, max_radial;
      std::vector<Hyp1f1> hyp1f1{};
      std::map<std::array<int,2>,int> map{};
      // Eigen::
     public:
      Hyp1f1SphericalExpansion(const int& max_angular, const int& max_radial, const double& tolerance = 1e-13, const size_t& precomputation_size = 200 )
        :max_angular{max_angular},max_radial{max_radial}  {
        int ipos{0};
        for(int l{0}; l < max_angular+1; l++) {
          for(int n{0}; n < max_radial; n++) {
            std::array<int,2> key{l,n};
            this->map[key] = ipos;
            ipos++;
          }
        }
      }
    };


  }  // namespace math
}  // namespace rascal

#endif  // SRC_MATH_HYP1F1_HH_
