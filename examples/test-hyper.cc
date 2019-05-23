#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <sstream>

#include "math/hyp1f1.hh"

namespace rascal {

  const double tol = 1e-14 * 100;  // it's in percent

}  // namespace rascal

class h1f1_cephes {
 private:
  double a, b;
  size_t mmax;

 public:
  h1f1_cephes(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax} {}

  inline double calc(double z) {
    return rascal::math::hyp1f1(this->a, this->b, z);
  }
};

#define NITER 1
template <class T>
void timeit(double a, double b, double z, double m) {
  double t{0};
  auto calculator = T(a, b, m);

  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < NITER; ++i) {
    t += i + calculator.calc(z);
  }

  auto finish = std::chrono::high_resolution_clock::now();
  t -= (NITER) * (NITER - 1.0) * 0.5;
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "a: " << a << " b: " << b << " z: " << z
            << " elapsed: " << elapsed.count() / NITER
            << " value: " << t / NITER << "\n";
}

using rascal::math::Hyp1f1;
using rascal::math::Hyp1f1Asymptotic;
using rascal::math::Hyp1f1Series;

// template<class T>
// double dfridr(T &func, const double x, const double h, double &err) {
//   const int ntab=10;
//   const double con=1.4, con2=(con*con);
//   const double big=numeric_limits<double>::max();
//   const double safe=2.0;
//   int i,j;
//   double errt,fac,hh,ans;
//   MatDoub a(ntab,ntab);
//   if (h == 0.0) throw("h must be nonzero in dfridr.");
//   hh=h;
//   a[0][0]=(func(x+hh)-func(x-hh))/(2.0*hh);
//   err=big;
//   for (i=1;i<ntab;i++) {
//     hh /= con;
//     a[0][i]=(func(x+hh)-func(x-hh))/(2.0*hh);
//     fac=con2;
//     for (j=1;j<=i;j++) {
//       a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
//       fac=con2*fac;
//       errt=std::max(std::abs(a[j][i]-a[j-1][i]),std::abs(a[j][i]-a[j-1][i-1]));
//       if (errt <= err) {
//         err=errt;
//         ans=a[j][i];
//       }
//     }
//     if (std::abs(a[i][i]-a[i-1][i-1]) >= safe*err) break;
//     }
//   return ans;
// }

int main() {
  std::cout.setf(std::ios_base::scientific);
  std::cout.precision(14);
  // double err{};
  // auto res{dfridr([]std::exp , 20, 1e-6, err)};
  // double rel_err{std::abs(1-res/std::exp(20))};
  // std::cout<< rel_err << std::endl;
  // std::vector<int> ls{{4,16}};
  // std::vector<int> ns{{0,2,4,5,7,10,13,16}};
  // std::vector<int> ls{{16}};
  // std::vector<int> ns{{8}};
  // for (auto& l : ls) {
  //     for (auto& n : ns) {
  // double a{0.5 * (3+n+l)};
  // double b{l+1.5};
  // double a{13};
  // double b{17.5};
  // // std::cout<< "h1f1_cephes" << std::endl;
  // // timeit<h1f1_cephes>(a,b, 0.01, 200);
  // // timeit<h1f1_cephes>(a,b, 0.1, 200);
  // // timeit<h1f1_cephes>(a,b, 1, 200);
  // // timeit<h1f1_cephes>(a,b, 10, 200);
  // // timeit<h1f1_cephes>(a,b, 100, 200);
  // // timeit<h1f1_cephes>(a,b, 300, 200);
  // // std::cout<<"\n";

  // std::cout<< "h1f1" << std::endl;
  // // timeit<Hyp1f1>(a,b, 0.01, 200);
  // // timeit<Hyp1f1>(a,b, 0.1, 200);
  // timeit<Hyp1f1>(a,b, 10, 200);
  // timeit<Hyp1f1Series>(a,b, 10, 200);
  // // timeit<Hyp1f1>(a,b, 10, 200);
  // // timeit<Hyp1f1>(a,b, 100, 200);
  // // timeit<Hyp1f1>(a,b, 300, 200);
  // std::cout<<"\n";
  //     }
  // }
  // std::cout<< "h1f1_hybrid" << std::endl;
  double a{3.5};
  double b{4.5};
  Hyp1f1Series hyp{a, b, 300};
  Hyp1f1Series hyp_d{a + 1, b + 1, 300};
  Hyp1f1Asymptotic hypa{a, b, 300};
  Hyp1f1Asymptotic hypa_d{a + 1, b + 1, 300};
  Hyp1f1 hyph{a, b, 200, 1e-12};
  h1f1_cephes ccc{a, b, 300};
  h1f1_cephes dccc{a + 1, b + 1, 300};
  double z{80};
  double h{1e-6};
  // std::cout << hyp.calc(z, false) << std::endl;
  // std::cout << hypa.calc(z, false) << std::endl;
  // std::cout << hyph.calc(z, false) << std::endl;
  // std::cout << hyp.calc(z, true) << std::endl;
  // std::cout << a/b*hyp_d.calc(z, false) << std::endl;
  // std::cout << hypa.calc(z, true) << std::endl;
  // std::cout << a/b*hypa_d.calc(z, false) << std::endl;
  // std::cout << hyph.calc(z, true) << std::endl;
  //    std::cout << (ccc.calc(z+h) - ccc.calc(z-h)) / (2*h)<< std::endl;
  //    std::cout << (hypa.calc(z+h) - hypa.calc(z-h)) / (2*h)<< std::endl;
  //    std::cout << (hyph.calc(z+h) - hyph.calc(z-h)) / (2*h)<< std::endl;
  //    std::cout << z+h << ", " << z-h << std::endl;
  //    std::cout << hyph.calc(z+h) << ", " << hyph.calc(z-h) << std::endl;
  //    // std::cout << (hyph.hyp1f1_asymptotic.calc(z+h) -
  //    hyph.hyp1f1_asymptotic.calc(z-h)) / (2*h)<< std::endl;
  std::cout << hyph.calc(z) << std::endl;
  std::cout << hyph.calc(z, true) << std::endl;
  std::cout << hyph.calc_numerical_derivative(z, h) << std::endl;
}
