#include <Eigen/Dense>

int main() {
  Eigen::ArrayXd x(10);
  Eigen::ArrayXd res(10);
  Eigen::ArrayXi exponents(10);
  x = Eigen::ArrayXd::Random(10);
  exponents = Eigen::ArrayXi::LinSpaced(10, 0, 9);
  // res = Eigen::pow(x, exponents);

  return (0);
}

