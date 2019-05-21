#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <sstream>

#include "math/math_interface.hh"
#include "math/math_utils.hh"
// #include "math/hyp1f1.hh"


namespace rascal {

  const double tol = 1e-14 * 100;  // it's in percent

  class Hyp1f1Series {
     protected:
      double a,b;
      size_t mmax;
      double prefac;
      double tolerance;

      Eigen::VectorXd coeff{};

     public:
      Hyp1f1Series(const double& a, const double& b, const size_t& mmax, const double& tolerance = 1e-14)
          : a{a}, b{b}, mmax{mmax},
            prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance} {
          coeff.resize(mmax);
          double u{0.};
          // precomputes the (a)_i / (b)_i / i! coefficients
          coeff(0)=a/b;
          for (size_t i{1}; i < mmax; ++i) {
              u = (a + i) / ((i + 1) * (b + i));
              coeff(i) = coeff(i-1) * u;
          }
      }
      //! Computes G(a,b,z)
      inline double calc(const double& r_ij, const double& alpha, const double& beta) {
        using math::pow;
        double z{pow(alpha*r_ij, 2) / (alpha + beta)};
        return this->prefac*this->calc(z)*std::exp(-alpha*r_ij*r_ij);
      }
      //! Computes 1F1
      inline double calc(const double& z) {
        using math::pow;

        auto&& a{this->a};
        auto&& b{this->b};

        // perform the sum
        double res{1.}, a1{1.};
        for (size_t i{0}; i < this->mmax; ++i) {
          a1 *= coeff(i) * z;
          if (a1 < this->tolerance*res) {
            break;
          }

          res += a1;

          if (res > 1e100){
            break; // overflow
          }
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
      double tolerance;
      size_t mmax;
      Eigen::VectorXd coeff{};

     public:
      Hyp1f1Asymptotic(const double& a, const double& b, const size_t& mmax, const double& tolerance = 1e-14 )
        : a{a}, b{b}, mmax{mmax}, prefac{std::tgamma(b)/std::tgamma(a)}, tolerance{tolerance} {
          double intpart;
          double f2{std::modf(2*(a-b), &intpart)};
          if (std::abs(f2) < 1e-14) {
            this->is_n_and_l = true;
          }


          coeff.resize(mmax);
          // precomputes the (1-a)_i / (b-a)_i / i! coefficients for the
          // computation of hyp2f0
          double bma{b-a}, oma{1-a}, u{0.};
          coeff(0) = bma * oma;
          for (size_t i{1}; i < mmax; ++i) {
              u = (bma+i)*(oma+i)/(i+1);
              coeff(i) = coeff(i-1) * u;
          }
      }

      //! Computes G(a,b,z)
      inline double calc(const double& r_ij, const double& alpha, const double& beta) {
        using math::pow;
        size_t imax{this->mmax};
        auto&& a{this->a};
        auto&& b{this->b};

        // argument of 1F1
        double z{pow(alpha*r_ij, 2) / (alpha + beta)};
        // simplification of the argument with exp(-alpha*r_ij^2)
        double z2{-alpha*beta*pow(r_ij, 2) / (alpha + beta)};

        double fac{0.};
        if (this->is_n_and_l){
          fac = std::sqrt(pow(z,static_cast<int>(2*(a-b))));
        } else {
          fac = pow(z,a-b);
        }
        return this->hyp2f0(z)*std::exp(z2)*fac;
      }

      //! Computes 1F1
      inline double calc(const double& z) {
        using math::pow;
        size_t imax{this->mmax};
        auto&& a{this->a};
        auto&& b{this->b};

        return this->prefac*std::exp(std::log(this->hyp2f0(z))+z)*pow(z,a-b);
      }

      //! computes hyp2f0 with arg1 = b-a and arg2 = 1-a arg3 = 1 / z
      inline double hyp2f0(const double& z) {
        using math::pow;
        size_t imax{this->mmax};
        auto&& a{this->a};
        auto&& b{this->b};

        double iz{1.0/z};
        double res{1.}, a1{1.}, a2{1.};
        // perform the sum
        for (size_t i{1}; i < this->mmax; ++i) {
          a1 *= coeff(i) * iz;
          res += a1;

          if (res > 1e100){
            break; // overflow
          }
          // check if two successive terms are within the tolerance
          // note that a,b,z>0 so no absolute value is needed
          if (a1 < this->tolerance*res or a2 < this->tolerance*res) {
            break;
          }
          a2 = a1;
        }
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
      double z_asympt{0.5};

     public:
      Hyp1f1(const double& a, const double& b, size_t mmax = 500, double tolerance = 1e-14)
        : hyp1f1_series{a,b,mmax,tolerance},hyp1f1_asymptotic{a,b,mmax,tolerance}, a{a}, b{b},tolerance{tolerance}, mmax{mmax} {
        // now we try to determine what is the switching point between
        // power series and asymptotic expansion
        double fs{hyp1f1_series.calc(z_asympt)};
        double fa{hyp1f1_asymptotic.calc(z_asympt)};
        std::cout <<fs << ", "<< fa<< std::endl;
        while(std::abs(fs-fa) > 10*tolerance) {

          z_asympt = z_asympt + 0.5;
          fs = hyp1f1_series.calc(z_asympt);
          fa = hyp1f1_asymptotic.calc(z_asympt);
          if (z_asympt > 100) {
            throw std::runtime_error("Could not find the switch value");
          }
        }
        std::cout << this->z_asympt<< ", "<< tolerance << std::endl;
      }

      inline double calc(const double& z) {
        if (z < this->z_asympt) {
            return this->hyp1f1_series.calc(z);
        } else {
            return this->hyp1f1_asymptotic.calc(z);
        }
      }

      inline double calc(const double& r_ij, const double& alpha, const double& beta) {
        double z{math::pow(alpha*r_ij, 2) / (alpha + beta)};
        if (z < this->z_asympt) {
            return this->hyp1f1_series.calc(r_ij,alpha,beta);
        } else {
            return this->hyp1f1_asymptotic.calc(r_ij,alpha,beta);
        }
      }
    };

}  // namespace rascal

class h1f1_series {
    private:
    double a,b;
    size_t mmax;
    Eigen::VectorXd coeff, dcoeff;

    public:
    h1f1_series(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax}
    {
        coeff.resize(mmax);
        dcoeff.resize(mmax);
        coeff(0)=a/b; dcoeff(0)=1;
        for (int i=1; i<mmax; ++i) {
            dcoeff(i)=(a+i)/((i+1)*(b+i));
            coeff(i)=coeff(i-1)*dcoeff(i);
        }
    }

    inline double calc(double z) {
        double res;
        int imax = mmax;
        auto a = this->a; auto b = this->b;
        imax = int((z-b+sqrt((z-b)*(z-b)+4*a*z))/2 *1.5)  + 8;
        if (imax>mmax) imax = mmax;
        res=coeff(imax-1)*z;
        for (int i=imax-2; i>=0; --i) {
            res=z*(coeff(i)+res);
        }
        return res+1;
    }

};

class h1f1_asympt {
    private:
    double a,b,prefac;
    size_t mmax;
    Eigen::VectorXd coeff, dcoeff;

    public:
    h1f1_asympt(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax}
    {
        coeff.resize(mmax);
        dcoeff.resize(mmax);
        double bma=b-a, oma=1-a;
        this->prefac = std::tgamma(b)/std::tgamma(a);
        coeff(0)=bma*oma; dcoeff(0)=1;
        for (int i=1; i<mmax; ++i) {
            dcoeff(i)=(bma+i)*(oma+i)/(i+1);
            coeff(i)=coeff(i-1)*dcoeff(i);
        }
    }

    inline double calc(double z) {
        double res;
        int imax = mmax;
        auto a = this->a; auto b = this->b;
        double iz=1.0/z;
        res=coeff(imax-1)*iz;
        for (int i=imax-2; i>=0; --i) {
            res=iz*(coeff(i)+res);
        }
        return this->prefac*std::exp(std::log(res+1)+z)*std::pow(z,a-b);
    }
};

#define H1F1_MINM 10
class h1f1_hybrid {
    private:
    double a, b, a_prefactor, z_asympt;
    double tolerance{1e-10};
    size_t mmax;
    // expansion coefficients for power series and asymptotic expansion
    Eigen::VectorXd p_coeff, a_coeff;


    inline double calc_power_series(double z, size_t& nterms) {
        double res;
        int imax = this->mmax;
        double zpow{z}, p_term;
        p_term = p_coeff(0)*zpow;
        res=1.0 + p_term;
        int i;
        for (i=1; i<mmax; ++i) {
            zpow *= z; // calc up so we can stop when we're OK
            p_term = zpow * p_coeff(i);
            if (res*tolerance > p_term) {
                break;
            }
            res += p_term;
        }
        nterms = i;
        return res;
    }

    inline double calc_asymptotic(double z, size_t& nterms) {
        double res;
        int imax = this->mmax;
        double iz{1.0/z}, izpow{iz}, a_term;
        a_term = a_coeff(0)*izpow;

        res=1.0 + a_term;
        int i;
        for (i=1; i<mmax; ++i) {
            izpow *= iz; // calc up so we can stop when we're OK
            a_term = izpow * a_coeff(i);
            //std::cout<<"res "<<res <<" < " <<a_term<<"\n";
            if (res > 0 and res*tolerance > std::fabs(a_term)) {
                break;
            }
            res += a_term;
            if (res>1e100) { i=mmax;  // exit but signal it's not stable
               }
        }
        nterms = i;
        return this->a_prefactor*std::exp(std::log(res)+z)*std::pow(z,a-b);
    }

    public:
    h1f1_hybrid(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax} {
        // initializes the expansion, computing enough terms to compute h1f1 for z=zmax
        p_coeff.resize(mmax);
        a_coeff.resize(mmax);

        double bma=b-a, oma=1-a;
        this->a_prefactor = std::tgamma(b)/std::tgamma(a);

        // compute expansion coefficients
        p_coeff(0)=a/b;
        a_coeff(0)=bma*oma;
        for (int i=1; i<mmax; ++i) {
            p_coeff(i)=p_coeff(i-1)*(a+i)/((b+i)*(i+1));
            a_coeff(i)=a_coeff(i-1)*(bma+i)*(oma+i)/(i+1);
        }

        // now we try to determine what is the switching point between
        // power series and asymptotic expansion
        size_t nterms_a, nterms_p;
        double z_above, z_below;
        z_asympt = 1;
        z_above = calc_power_series(z_asympt, nterms_p);
        z_below = calc_asymptotic(z_asympt, nterms_a);
        //std::cout<<" ps: "<<z_above << " n: "<<nterms_p<<"\n";
        //std::cout<<" as: "<<z_below << " n: "<<nterms_a<<"\n";
        if (nterms_p < nterms_a) {
            z_below = z_asympt;
            while (nterms_p<nterms_a) {
                z_asympt *= 2.0;
                calc_power_series(z_asympt, nterms_p);
                calc_asymptotic(z_asympt, nterms_a);
            }
            z_above = z_asympt;
        }
        else {
            z_above = z_asympt;
            while (nterms_p > nterms_a) {
                z_asympt *= 0.5;
                calc_power_series(z_asympt, nterms_p);
                calc_asymptotic(z_asympt, nterms_a);
            }
            z_below = z_asympt;
        }

        double z_new;
        z_new = (z_above + z_below)*0.5;
        calc_power_series(z_new, nterms_p);
        calc_asymptotic(z_new, nterms_a);
        while (z_above - z_below > 1e-8) {
            if (nterms_p > nterms_a) {
                z_above = z_new;
            } else {
                z_below = z_new;
            }
            z_new = (z_above + z_below)*0.5;
            calc_power_series(z_new, nterms_p);
            calc_asymptotic(z_new, nterms_a);
        }
        z_asympt = (z_above + z_below) * 0.5;
    }

    inline double calc(double z) {
        size_t nterms;
        if (z<z_asympt) {
            return calc_power_series(z, nterms);
        } else {
            return calc_asymptotic(z, nterms);
        }
    }

};

class h1f1_cephes {
    private:
    double a,b; size_t mmax;

    public:
    h1f1_cephes(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax} {}

    inline double calc(double z) { return rascal::math::hyp1f1(this->a,this->b,z); }
};

#define NITER 10000
template<class T>
void timeit(double a, double b, double z, double m)
{

    double t{0};
    auto calculator = T(a,b,m);

    auto start = std::chrono::high_resolution_clock::now();

    for (int i=0; i<NITER; ++i)
    {   t+=i+calculator.calc(z); }

    auto finish = std::chrono::high_resolution_clock::now();
    t-=(NITER)*(NITER-1.0)*0.5;
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "a: "<<a <<" b: "<<b<<" z: "<<z<<" elapsed: "<<elapsed.count()/NITER<< " value: "<<t/NITER<<"\n";

}

using rascal::Hyp1f1;
using rascal::Hyp1f1Series;
using rascal::Hyp1f1Asymptotic;

int main(int argc, char * argv[])
{

    std::vector<int> ls{{8}};
    std::vector<int> ns{{3}};
    for (auto& l : ls) {
        for (auto& n : ns) {
            double a{0.5 * (3+n+l)};
            double b{l+1.5};
            std::cout<< "h1f1_michele" << std::endl;
            // timeit<h1f1_hybrid>(a,b, 0.01, 200);
            timeit<h1f1_hybrid>(a,b, 0.1, 200);
            // timeit<h1f1_hybrid>(a,b, 1, 200);
            // timeit<h1f1_hybrid>(a,b, 10, 200);
            // timeit<h1f1_hybrid>(a,b, 100, 200);
            // timeit<h1f1_hybrid>(a,b, 200, 200);
            std::cout<<"\n";

            std::cout<< "h1f1" << std::endl;
            // timeit<Hyp1f1Series>(a,b, 0.01, 200);
            timeit<Hyp1f1Series>(a,b, 0.1, 200);
            // timeit<Hyp1f1Series>(a,b, 1, 200);
            // timeit<Hyp1f1Series>(a,b, 10, 200);
            // timeit<Hyp1f1Series>(a,b, 100, 200);
            // timeit<Hyp1f1Series>(a,b, 200, 200);
            std::cout<<"\n";
        }
    }
    // std::cout<< "h1f1_hybrid" << std::endl;



}
