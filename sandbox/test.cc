#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <list>
#include <random>
#include <string>
#include <algorithm>
#include <iterator>

#include <Eigen/Core>
#include <Eigen/CXX11/Tensor>


using Vector_t = Eigen::VectorXd;
using Vector_CRef = typename const Eigen::Ref<const Vector_t>;
using Vector_CMap = typename const Eigen::Map<const Vector_t>;

using VectorG_t = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
using VectorG_CRef = typename const Eigen::Ref<const VectorG_t>;
using VectorG_CMap = typename const Eigen::Map<const VectorG_t>;

using Matrix_t =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Matrix_CRef = typename const Eigen::Ref<const Matrix_t>;
using Matrix_CMap = typename const Eigen::Map<const Matrix_t>;

using T3nm_t = Tensor<data_type, rank>