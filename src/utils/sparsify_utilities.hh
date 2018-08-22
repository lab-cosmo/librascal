
#include "basic_types.hh"

namespace rascal {
namespace utils {
  using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::tuple<Eigen::ArrayXi, Eigen::ArrayXd> fps(const Eigen::Ref<const RowMatrixXd>& feature_matrix, int n_sparse=0, int i_first_point=0);
  
    
} //namespace utils
} // namespace rascal
