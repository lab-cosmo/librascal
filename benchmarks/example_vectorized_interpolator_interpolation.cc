#include <chrono>
#include "benchmark_interpolator.hh"

#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

using namespace rascal::math;
using namespace rascal::internal;


namespace Eigen{
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(typename Matrix::Index));
    out.write((char*) (&cols), sizeof(typename Matrix::Index));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
}
template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    typename Matrix::Index rows=0, cols=0;
    in.read((char*) (&rows),sizeof(typename Matrix::Index));
    in.read((char*) (&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
    in.close();
}
} // Eigen::

inline bool file_exists(const char* name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 
}

static constexpr int ITERATIONS = 200;

int main(){
  auto intp{InterpolatorVectorized <
    InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
    GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
    SearchMethod<SearchMethod_t::Uniform>
      >()};
  int max_radial{3};
  int max_angular{max_radial-1};
  json fc_hypers{
       {"type", "Constant"},
       {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
      };
  json hypers{{"gaussian_density", fc_hypers},
            {"max_radial", max_radial},
            {"max_angular", max_angular},
            {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
  };
  auto radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
  std::function<Matrix_t(double)> func = [&radial_contr](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x, 0.5);};

  double x1{0};
  double x2{8};
  double mean_error_bound{1e-10};
  size_t nb_points = 1e6;
  size_t nb_iterations = 100000;
  const char* filename{"interpolator_vectorized_grid.dat"};

  if (not(file_exists(filename))) {
    std::cout << "Grid file does not exists, has to be computed." << std::endl;
    intp.initialize(func, x1, x2, mean_error_bound); 
    Eigen::write_binary(filename, intp.grid);
  } else {
    std::cout << "Grid file exists, is read." << std::endl;
    Vector_t grid;
    Eigen::read_binary(filename, grid);
    intp.initialize(func, x1, x2, grid);
  }
  std::cout << "grid size=" << intp.grid.size() << std::endl;

  Vector_t points_tmp  = Vector_t::LinSpaced(nb_points,x1,x2);
  Vector_t points = Vector_t::Zero(nb_points);
  srand(50000000);
  for (size_t i{0};i<nb_points;i++) {
    points(i) = points_tmp(rand() % nb_points);
  }

  Matrix_t mat_tmp = Matrix_t::Zero(max_radial, max_angular+1);
  std::chrono::duration<double> elapsed{};
  auto start = std::chrono::high_resolution_clock::now();
  for (int j{0}; j < ITERATIONS; j++) {
    for (size_t i{0}; i<nb_iterations;i++) {
      mat_tmp = intp.interpolate(points(i % nb_points));
    }
  }
  auto finish = std::chrono::high_resolution_clock::now();

  elapsed = finish - start;
  std::cout << "Interpolation of " << nb_iterations << " points"
            << " elapsed: " << elapsed.count()/ITERATIONS << " seconds"            
            << std::endl;

  return 0;
}
