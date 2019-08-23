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

//inline bool file_exists(const char* name) {
//  struct stat buffer;   
//  return (stat (name, &buffer) == 0); 
//}

static constexpr int ITERATIONS = 100;

int main(){
  auto intp{Interpolator <
    InterpolationMethod<InterpolationMethod_t::CubicSpline>,
    GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
    SearchMethod<SearchMethod_t::Uniform>
      >()};
  //auto func = [](double x) {return std::exp(-std::pow((x-1)/0.5,2)/2);};

  double n = 5;
  double l = 5;
  double a = 0.5*(n+l+3);
  double b = l+1.5;
  auto hyp1f1 = Hyp1f1(a, b, 200, 1e-15);

  double x1{0};
  double x2{8};
  std::function<double(double)> func = [&hyp1f1](double x) {return hyp1f1.calc(x);};
  double mean_error_bound{1e-5};
  size_t nb_points = 1e6;
  //size_t nb_iterations = 100000;
  std::vector<size_t> nbs_iterations = {1000,10000,100000};//,1000000};
  const char* filename{"interpolator_hyp1f1_grid.dat"};

  if (not(file_exists(filename))) {
    std::cout << "Grid file does not exists, has to be computed." << std::endl;
    intp.initialize(func, x1, x2, mean_error_bound); 
    Eigen::write_binary(filename, intp.grid);
  } else {
    std::cout << "Grid file exists, is read." << std::endl;
    Vector_t grid;
    Eigen::read_binary(filename, grid);
    intp.initialize(func, grid);
  }
  std::cout << "grid size=" << intp.grid.size() << std::endl;

  Vector_t points_tmp  = Vector_t::LinSpaced(nb_points,x1,x2);
  Vector_t points = Vector_t::Zero(nb_points);
  int SEED =  1597463007; //1597463007;
  srand(SEED); //50000000
  for (size_t i{0};i<nb_points;i++) {
    points(i) = points_tmp(rand() % nb_points);
  }



  std::chrono::duration<double> elapsed{};


  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {
        points_tmp(i % nb_points) =  hyp1f1.calc(points(i % points.size()));
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout << std::fixed;
    std::cout 
      << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "Hyp1f1 of " << nb_iterations << " points"
              << std::endl;
  }

  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {
        points_tmp(i % nb_points) = intp.interpolate(points(i % nb_points));
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout 
              << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "Interpolation of " << nb_iterations << " points"
              << std::endl;
  }

  // RADIAL CONTRIBUTION
  int max_radial{1};
  int max_angular{max_radial};
  json fc_hypers{
       {"type", "Constant"},
       {"gaussian_sigma", {{"value", 0.5}, {"unit", "A"}}}
      };
  json hypers{{"gaussian_density", fc_hypers},
            {"max_radial", max_radial},
            {"max_angular", max_angular},
            {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
  };
  auto radial_contr{RadialContribution<RadialBasisType::GTO>(hypers)};
  func = [&radial_contr](double x) {
    radial_contr.compute_neighbour_contribution(x, 0.5);
    return radial_contr.radial_integral_neighbour(0,0);
  };

  const char* filename2{"interpolator_radial_grid.dat"};
  if (not(file_exists(filename2))) {
    std::cout << "Grid file does not exists, has to be computed." << std::endl;
    intp.initialize(func, x1, x2, mean_error_bound); 
    Eigen::write_binary(filename2, intp.grid);
  } else {
    std::cout << "Grid file exists, is read." << std::endl;
    Vector_t grid;
    Eigen::read_binary(filename2, grid);
    intp.initialize(func, grid);
  }
  std::cout << "grid size=" << intp.grid.size() << std::endl;


  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {
        points_tmp(i % points.size()) = radial_contr.compute_contribution<AtomicSmearingType::Constant>(points(i % points.size()), 0.5)(0,0);
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout << std::fixed;
    std::cout 
      << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "Radial of " << nb_iterations << " points"
              << std::endl;
  }

  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {
        points_tmp(i % nb_points) = intp.interpolate(points(i % nb_points));
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout 
              << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "Interpolation of " << nb_iterations << " points"
              << std::endl;
  }

  // RADIAL CONTRIBUTION VECTORIZED
  max_radial = 5;
  max_angular = max_radial;
  hypers = {{"gaussian_density", fc_hypers},
              {"max_radial", max_radial},
              {"max_angular", max_angular},
              {"cutoff_function", {{"cutoff",{{"value", 2.0}, {"unit", "A"}}}}}
    };
  radial_contr = RadialContribution<RadialBasisType::GTO>(hypers);
  std::function<Matrix_t(double)> func_vec = [&radial_contr](double x) {return radial_contr.compute_contribution<AtomicSmearingType::Constant>(x, 0.5);};

  auto intp_vec{InterpolatorVectorized <
    InterpolationMethod<InterpolationMethod_t::CubicSplineVectorized>,
    GridRational<GridType_t::Uniform, RefinementMethod_t::Exponential>,
    SearchMethod<SearchMethod_t::Uniform>
      >()};
  const char* filename_radial_vec{"interpolator_radial_vec_grid.dat"};
  if (not(file_exists(filename_radial_vec))) {
    std::cout << "Grid file does not exists, has to be computed." << std::endl;
    intp_vec.initialize(func_vec, x1, x2, mean_error_bound); 
    Eigen::write_binary(filename_radial_vec, intp.grid);
  } else {
    std::cout << "Grid file exists, is read." << std::endl;
    Vector_t grid;
    Eigen::read_binary(filename_radial_vec, grid);
    intp_vec.initialize(func_vec, grid);
  }
  std::cout << "grid size=" << intp_vec.grid.size() << std::endl;

  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {        
        Matrix_t points_vec_tmp = radial_contr.compute_contribution<AtomicSmearingType::Constant>(points(i % points.size()), 0.5);
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout << std::fixed;
    std::cout 
      << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "RadialVec of " << nb_iterations << " points"
              << std::endl;
  }

  std::cout << std::endl;

  for (size_t nb_iterations : nbs_iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j{0}; j < ITERATIONS; j++) {
      for (size_t i{0}; i<nb_iterations;i++) {
        Matrix_t points_vec_tmp = intp_vec.interpolate(points(i % nb_points)); // ~x41
        //Matrix_t points_vec_tmp = intp_vec.interpolate_optimal(points(i % nb_points)); // ~x41
        //Matrix_t points_vec_tmp = points(i % nb_points) * Matrix_t::Ones(max_radial, max_angular+1); // ~x65
      }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    elapsed = finish - start;
    std::cout 
              << " elapsed: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()/double(ITERATIONS) << " ns "
      << "Interpolation of " << nb_iterations << " points"
              << std::endl;
  }



  return 0;
}
