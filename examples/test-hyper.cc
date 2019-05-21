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
    double a,b; size_t mmax;

    public:
    h1f1_cephes(double a, double b, size_t mmax) : a{a}, b{b}, mmax{mmax} {}

    inline double calc(double z) { return rascal::math::hyp1f1(this->a,this->b,z); }
};

#define NITER 100
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

using rascal::math::Hyp1f1;

int main()
{
    std::cout.setf(std::ios_base::scientific);
    std::cout.precision(14);
    std::vector<int> ls{{4,16}};
    std::vector<int> ns{{0,2,4,5,7,10,13,16}};
    for (auto& l : ls) {
        for (auto& n : ns) {
            double a{0.5 * (3+n+l)};
            double b{l+1.5};
            std::cout<< "h1f1_cephes" << std::endl;
            timeit<h1f1_cephes>(a,b, 0.01, 200);
            timeit<h1f1_cephes>(a,b, 0.1, 200);
            timeit<h1f1_cephes>(a,b, 1, 200);
            timeit<h1f1_cephes>(a,b, 10, 200);
            timeit<h1f1_cephes>(a,b, 100, 200);
            timeit<h1f1_cephes>(a,b, 300, 200);
            std::cout<<"\n";

            std::cout<< "h1f1" << std::endl;
            timeit<Hyp1f1>(a,b, 0.01, 200);
            timeit<Hyp1f1>(a,b, 0.1, 200);
            timeit<Hyp1f1>(a,b, 1, 200);
            timeit<Hyp1f1>(a,b, 10, 200);
            timeit<Hyp1f1>(a,b, 100, 200);
            timeit<Hyp1f1>(a,b, 300, 200);
            std::cout<<"\n";
        }
    }
    // std::cout<< "h1f1_hybrid" << std::endl;



}
