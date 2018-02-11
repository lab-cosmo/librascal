#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;



/* Bind MatrixXd (or some other Eigen type) to Python */
typedef Eigen::MatrixXd Matrix;

typedef Matrix::Scalar Scalar;
constexpr bool rowMajor = Matrix::Flags & Eigen::RowMajorBit;

py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
    .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                                /* Pointer to buffer */
            sizeof(Scalar),                          /* Size of one scalar */
            py::format_descriptor<Scalar>::format(), /* Python struct-style format descriptor */
            2,                                       /* Number of dimensions */
            { m.rows(), m.cols() },                  /* Buffer dimensions */
            { sizeof(Scalar) * (rowMajor ? m.cols() : 1),
              sizeof(Scalar) * (rowMajor ? 1 : m.rows()) }
                                                     /* Strides (in bytes) for each index */
        );
     })
    .def("__init__", [](Matrix &m, py::buffer b) {
        typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

        /* Request a buffer descriptor from Python */
        py::buffer_info info = b.request();

        /* Some sanity checks ... */
        if (info.format != py::format_descriptor<Scalar>::format())
            throw std::runtime_error("Incompatible format: expected a double array!");

        if (info.ndim != 2)
            throw std::runtime_error("Incompatible buffer dimension!");

        auto strides = Strides(
            info.strides[rowMajor ? 0 : 1] / (py::ssize_t)sizeof(Scalar),
            info.strides[rowMajor ? 1 : 0] / (py::ssize_t)sizeof(Scalar));

        auto map = Eigen::Map<Matrix, 0, Strides>(
            static_cast<Scalar *>(info.ptr), info.shape[0], info.shape[1], strides);

        new (&m) Matrix(map);
    });
