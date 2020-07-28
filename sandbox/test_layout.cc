#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <iomanip>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "rascal/utils/json_io.hh"


using Vector_t = Eigen::VectorXd;
using Vector_CRef = const typename Eigen::Ref<const Vector_t>;
using Vector_CMap = const typename Eigen::Map<const Vector_t>;

using VectorG_t = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using VectorG_CRef = const typename Eigen::Ref<const VectorG_t>;
using VectorG_CMap = const typename Eigen::Map<const VectorG_t>;

using Matrix_t =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Matrix_CRef = const typename Eigen::Ref<const Matrix_t>;
using Matrix_CMap = const typename Eigen::Map<const Matrix_t>;

using Tj3nm_t = typename Eigen::Tensor<double, 4, Eigen::RowMajor>;

using Tnm_t = typename Eigen::Tensor<double, 2, Eigen::RowMajor>;
using Tnm_CMap_t = const typename Eigen::TensorMap<const Tnm_t>;

// TensorFixedSize<double, Sizes<>>

using namespace rascal;  // NOLINT


class Timer {
 private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

 public:
	Timer() : m_beg(clock_t::now()) { }

	void reset() {
		m_beg = clock_t::now();
	}

	double elapsed() const {
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

double std_dev(const Vector_t& vec) {
  double std_dev = std::sqrt((vec.array() - vec.mean()).array().square().sum()/(vec.size()-1));
  return std_dev;
}

void compute_0(const std::vector<Matrix_t> & dIljn, const std::vector<Matrix_t> & Yljm,
                const VectorG_t& rij,  std::vector<Tj3nm_t> & Rl3nm) {
  for (size_t l{0}; l <dIljn.size(); l++) {
    auto & dIjn = dIljn[l];
    auto & Yjm = Yljm[l];
    auto & R3nm = Rl3nm[l];
    const auto n_max = dIjn.cols();
    const int m_max = 2*l+1;
    auto inter = Matrix_t(n_max, m_max);
    for (int i_neigh{0}; i_neigh < dIjn.rows(); i_neigh++) {
      for (int i_der{0}; i_neigh < 3; i_der++) {
        Eigen::array<int, 4> offsets = {i_neigh, i_der, 0, 0};
        Eigen::array<int, 4> extents = {1, 1, n_max , m_max};
        inter = rij(i_neigh, i_der) * dIjn.row(i_neigh).transpose() * Yjm.row(i_neigh);
        auto inter_ten = Tnm_CMap_t(inter.data(), n_max, m_max);
        R3nm.slice(offsets, extents) = inter_ten;
      }
    }
  }
}

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide setup json filename as argument";
    std::cerr << std::endl;
    return -1;
  }

  json input = json_io::load(argv[1]);

  auto n_max = input["n"].get<int>();
  auto l_max = input["l"].get<int>();
  auto n_neighs = input["n_neigh"].get<std::vector<int>>();
  auto N_ITERATIONS = input["N_ITERATIONS"].get<int>();
  int n_neighbors{0};

  for (auto & n_neigh : n_neighs) {
    n_neighbors += n_neigh;
  }

  std::vector<Matrix_t> dIljn{};
  std::vector<Matrix_t> Yljm{};
  std::vector<Tj3nm_t> Rl3nm{};
  for (int l{0}; l <l_max; l++) {
    dIljn.push_back(Matrix_t::Random(n_neighbors, n_max));
    Yljm.push_back(Matrix_t::Random(n_neighbors, 2*l+1));
    Rl3nm.push_back(Tj3nm_t(n_neighbors, 3, n_max ,2*l+1));
  }
  auto rij = VectorG_t::Random(n_neighbors, 3);


  Vector_t elapsed{N_ITERATIONS};
  Timer timer{};

  for (int looper{0}; looper < N_ITERATIONS; looper++) {
    timer.reset();
    compute_0(dIljn, Yljm, rij, Rl3nm);
    elapsed[looper] = timer.elapsed();
  }



  json cmp0{};
  cmp0["elapsed_mean"] = elapsed.mean();
  cmp0["elapsed_std"] = std_dev(elapsed);
  cmp0["elapsed"] = elapsed;

  json results{};
  results["compute_0"] = cmp0;
  std::ofstream o(argv[2]);
  o << std::setw(2) << results << std::endl;
}