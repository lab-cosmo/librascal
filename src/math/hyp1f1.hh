// /**
//  * file   math_utils.hh
//  *
//  * @author  Felix Musil <felix.musil@epfl.ch>
//  * @author  Max Veit <max.veit@epfl.ch>
//  *
//  * @date   14 October 2018
//  *
//  * @brief contains the implementation of miscellaneous math functions
//  *
//  * Copyright  2018  Felix Musil, Max Veit, COSMO (EPFL), LAMMM (EPFL)
//  *
//  * Rascal is free software; you can redistribute it and/or
//  * modify it under the terms of the GNU Lesser General Public License as
//  * published by the Free Software Foundation, either version 3, or (at
//  * your option) any later version.
//  *
//  * Rascal is distributed in the hope that it will be useful, but
//  * WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  * Lesser General Public License for more details.
//  *
//  * You should have received a copy of the GNU Lesser General Public License
//  * along with this software; see the file LICENSE. If not, write to the
//  * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
//  * Boston, MA 02111-1307, USA.
//  */

// #ifndef SRC_MATH_MATH_UTILS_HH_
// #define SRC_MATH_MATH_UTILS_HH_

// #include "math_interface.hh"

// #include <Eigen/Dense>
// #include <cmath>
// #include <limits>

// namespace rascal {
//   namespace math {
//     using MatrixX2_t = Eigen::Matrix<double, Eigen::Dynamic, 2>;

//     /**
//      * Utilities to apply the recurrence relations of 1F1 from maximal
//      * values of a and b to minimal values.
//      * The recurrence path is chosen so that derivatives of 1F1 are also
//      * a result.
//      */
//     inline double recurence_to_val_downward(const double& a, const double& b, const double& z, const double& M1p2p, const double& M1p1p) {
//       return z * (a - b) * M1p2p / b /(b + 1) + M1p1p;
//     }

//     inline double recurence_to_der_downward(const double& a, const double& b, const double& z, const double& M2p3p, const double& M1p2p) {
//       return z * (a + 1) * M2p3p / (b + 2) / (b + 1) + M1p2p;
//     }

//     inline MatrixX2_t hyp1f1_recurrence(const int& l_max, const int& n, const double& z, const double& f) {
//       // first
//       MatrixX2_t results(l_max, 2);

//     }


//     /**
//      * Implementation of the 1F1 function optimized for its use in
//      * represenation_manager_spherical_expansion.hh
//      *
//      * The function that is computed is actually
//      *  G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
//      *                    * 1F1(a,b,z)
//      *
//      * where 1F1(a,b,z) = \sum_{i=0}^{\infty} \frac{(a)_i}{i! (b)_i} z^{i}
//      * and (a)_i is the Pochhammer symbol
//      */

//     /**
//      * Computes the 1F1 with the series expansion
//      *
//      */
//     class Hyp1f1Series {
//      protected:
//       double a,b;
//       size_t mmax;
//       double prefac;
//       double tolerance;

//       Eigen::VectorXd coeff{};

//      public:
//       Hyp1f1Series(const double& a, const double& b, const size_t& mmax, const double& tolerance)
//           : a{a}, b{b}, mmax{mmax},
//             prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance} {
//         coeff.resize(mmax);
//         double u{0.};
//         // precomputes the (a)_i / (b)_i / i! coefficients
//         coeff(0)=a/b;
//         for (size_t i{1}; i < mmax; ++i) {
//             u = (a + i) / ((i + 1) * (b + i));
//             coeff(i) = coeff(i-1) * u;
//         }
//       }
//       //! Computes G(a,b,z)
//       inline double calc(const double& r_ij, const double& alpha, const double& beta) {
//         using math::pow;
//         double z{pow(alpha*r_ij, 2) / (alpha + beta)};
//         return this->prefac*this->calc(z)*std::exp(-alpha*r_ij*r_ij);
//       }
//       //! Computes 1F1
//       inline double calc(const double& z) {
//         using math::pow;
//         double res{1.}, a1{1.}, a2{1.};
//         auto&& a{this->a};
//         auto&& b{this->b};

//         // perform the sum
//         for (size_t i{1}; i < this->mmax; ++i) {
//           a1 *= coeff(i) * z;
//           res += a1;

//           if (res > 1e100){
//             break; // overflow
//           }

//           // check if two successive terms are within the tolerance
//           // note that a,b,z>0 so no absolute value is needed
//           if (a1 < this->tolerance*res or a2 < this->tolerance*res) {
//             break;
//           }
//           a2 = a1;
//         }

//         return res;
//       }
//     };

//     /**
//      * Computes the 1F1 with the asymptotic limit
//      *  1F1(a,b,z) \sim \frac{\exp{z} z^{a-b} \Gamma{b}}{\Gamma{a}}
//      *                    \sum_{j=0}^{\infty} \frac{(b-a)_j(1-a)_j}{j!} z^{-j}
//      *
//      *  G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
//      *                    * 1F1(a,b,z)
//      */
//     class Hyp1f1Asymptotic {
//      private:
//       double a,b,prefac;
//       bool is_n_and_l{false};
//       double tolerance;
//       size_t mmax;
//       Eigen::VectorXd coeff{};

//      public:
//       Hyp1f1Asymptotic(const double& a, const double& b, const size_t& mmax, const double& tolerance)
//         : a{a}, b{b}, mmax{mmax}, prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance} {
//           double intpart;
//           double f2{std::modf(2*(a-b), &intpart)};
//           if (std::abs(f2) < 1e-14) {
//             this->is_n_and_l = true;
//           }


//           coeff.resize(mmax);
//           // precomputes the (1-a)_i / (b-a)_i / i! coefficients for the
//           // computation of hyp2f0
//           double bma{b-a}, oma{1-a}, u{0.};
//           coeff(0) = bma * oma;
//           for (size_t i{1}; i < mmax; ++i) {
//               u = (bma+i)*(oma+i)/(i+1);
//               coeff(i) = coeff(i-1) * u;
//           }
//       }

//       //! Computes G(a,b,z)
//       inline double calc(const double& r_ij, const double& alpha, const double& beta) {
//         using math::pow;
//         size_t imax{this->mmax};
//         auto&& a{this->a};
//         auto&& b{this->b};

//         // argument of 1F1
//         double z{pow(alpha*r_ij, 2) / (alpha + beta)};
//         // simplification of the argument with exp(-alpha*r_ij^2)
//         double z2{-alpha*beta*pow(r_ij, 2) / (alpha + beta)};

//         double fac{0.};
//         if (this->is_n_and_l){
//           fac = std::sqrt(pow(z,static_cast<int>(2*(a-b))));
//         } else {
//           fac = pow(z,a-b);
//         }
//         return this->hyp2f0(z)*std::exp(z2)*fac;
//       }

//       //! Computes 1F1
//       inline double calc(const double& z) {
//         using math::pow;
//         size_t imax{this->mmax};
//         auto&& a{this->a};
//         auto&& b{this->b};

//         return this->prefac*std::exp(std::log(this->hyp2f0(z))+z)*pow(z,a-b);
//       }

//       //! computes hyp2f0 with arg1 = b-a and arg2 = 1-a arg3 = 1 / z
//       inline double hyp2f0(const double& z) {
//         using math::pow;
//         size_t imax{this->mmax};
//         auto&& a{this->a};
//         auto&& b{this->b};

//         double iz{1.0/z};
//         double res{1.}, a1{1.}, a2{1.};
//         // perform the sum
//         for (size_t i{1}; i < this->mmax; ++i) {
//           a1 *= coeff(i) * iz;
//           res += a1;

//           if (res > 1e100){
//             break; // overflow
//           }
//           // check if two successive terms are within the tolerance
//           // note that a,b,z>0 so no absolute value is needed
//           if (a1 < this->tolerance*res or a2 < this->tolerance*res) {
//             break;
//           }
//           a2 = a1;
//         }
//         return res;
//       }
//     };

//     /**
//      * Implements the switch between the series and asymptotic expansion
//      * methods to compute 1F1
//      *
//      * Note that:
//      * G(a,b,z) = \frac{\Gamma(a)}{\Gamma(b)} * \exp{-\alpha r_{ij}^2}
//      *                    * 1F1(a,b,z)
//      */
//     class Hyp1f1 {
//      private:
//       Hyp1f1Series hyp1f1_series;
//       Hyp1f1Asymptotic hyp1f1_asymptotic;
//       double a, b;
//       double tolerance;
//       size_t mmax;
//       double z_asympt{0.5};

//      public:
//       Hyp1f1(const double& a, const double& b, size_t mmax = 200, double tolerance = 1e-14)
//         : hyp1f1_series{a,b,mmax,tolerance},hyp1f1_asymptotic{a,b,mmax,tolerance}, a{a}, b{b},tolerance{tolerance}, mmax{mmax} {
//         // now we try to determine what is the switching point between
//         // power series and asymptotic expansion
//         double fs{hyp1f1_series.calc(z_asympt)};
//         double fa{hyp1f1_asymptotic.calc(z_asympt)};
//         while(std::abs(fs-fa) > 10*tolerance) {
//           z_asympt = z_asympt + 0.5;
//           fs = hyp1f1_series.calc(z_asympt);
//           fa = hyp1f1_asymptotic.calc(z_asympt);
//           if (z_asympt > 100) {
//             throw std::runtime_error("Could not find the switch value");
//           }
//         }
//       }

//       inline double calc(const double& z) {
//         if (z < this->z_asympt) {
//             return this->hyp1f1_series.calc(z);
//         } else {
//             return this->hyp1f1_asymptotic.calc(z);
//         }
//       }

//       inline double calc(const double& r_ij, const double& alpha, const double& beta) {
//         double z{math::pow(alpha*r_ij, 2) / (alpha + beta)};
//         if (z < this->z_asympt) {
//             return this->hyp1f1_series.calc(r_ij,alpha,beta);
//         } else {
//             return this->hyp1f1_asymptotic.calc(r_ij,alpha,beta);
//         }
//       }
//     };


//   }  // namespace math
// }  // namespace rascal

// #endif  // SRC_MATH_MATH_UTILS_HH_
