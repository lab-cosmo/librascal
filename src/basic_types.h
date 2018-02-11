#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Dense>

namespace py = pybind11;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>  Matrix_f;
typedef Matrix::Scalar Scalar;
constexpr bool rowMajor = Matrix::Flags & Eigen::RowMajorBit;

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::RowMajor>  Positions_c;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 3, Eigen::ColMajor>  Positions_f;

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

Matrix cdist(Matrix &X,Matrix &Y);

Matrix cdist_pf(Positions_f &X,Positions_f &Y);

void export_proteus_basic_types(py::module& m);



/*
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/buffer_info.h>
namespace py = pybind11;


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  Matrix;
*/
