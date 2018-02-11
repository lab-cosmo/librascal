#include "basic_types.h"
#include <iostream>

namespace py = pybind11;


Matrix cdist(Matrix &X,Matrix &Y){
    const int N = X.rows();
    const int M = X.cols();
    const int K = Y.rows();

    //std::cout << N<<M<<K << std::endl;
    // Allocate parts of the expression
    //Matrix D = Matrix::Zero(N,K);
    //std::cout << X << std::endl;
    //std::cout << Y << std::endl;
    //D = XX * Matrix::Ones(1,K) + Matrix::Ones(N,1) * YY - 2*X*Y.transpose();
    Matrix XX,YY,D;
    XX.resize(N,1);
    YY.resize(K,1);
    D.resize(N,K);
    XX = X.rowwise().squaredNorm();
    YY = Y.rowwise().squaredNorm();
    D = -2*X*Y.transpose();
    for( int ii= 0; ii < N; ii++ ){
        for( int jj = 0; jj < K; jj++ ){
            D(ii,jj) +=  XX(ii) + YY(jj);
        }}
    /*
    for( int ii= 0; ii < N; ii++ ){
        for( int jj = 0; jj < K; jj++ ){
            for( int kk = 0; kk < M; kk++ ){

                D(ii,jj) += (X(ii,kk)-Y(jj,kk))*(X(ii,kk)-Y(jj,kk));
    }}}
    */
    //D.array().Eigen::sqrt().matrix();
    //std::cout << D << std::endl;

    return Eigen::sqrt(D.array());
}

Matrix cdist_pf(Positions_f &X,Positions_f &Y){
    const int N = X.rows();
    const int M = Y.rows();

    Matrix D;
    D.resize(N,M);
    Vector x = Vector::Ones(N,1);Vector y = Vector::Ones(N,1);Vector z = Vector::Ones(N,1);

    #pragma omp parallel for shared(x,y,z)
    for ( int ii= 0; ii < N; ii++ ){
        D.row(ii) = (X(ii,0)*x-Y.col(0)).array().square() + (X(ii,1)*y-Y.col(1)).array().square() + (X(ii,2)*z-Y.col(2)).array().square();
           /* for( int jj= 0; jj < M; jj++ ){
            D(ii,jj) = (X(ii,0) - Y(jj,0))*(X(ii,0) - Y(jj,0)) + (X(ii,1) - Y(jj,1))*(X(ii,1) - Y(jj,1)) + (X(ii,2) - Y(jj,2))*(X(ii,2) - Y(jj,2));
        }*/
        }

    return Eigen::sqrt(D.array());
}


void export_proteus_basic_types(py::module& m){

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
        })
        .def("__getitem__", [](const Matrix &m, std::pair<ssize_t, ssize_t> i) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            return m(i.first, i.second);
        })
        .def("__setitem__", [](Matrix &m, std::pair<ssize_t, ssize_t> i, Scalar v) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            m(i.first, i.second) = v;
        })
        .def("get_shape", [](const Matrix &m) {
        std::pair<int, int> shape( m.rows(), m.cols() );
        return shape;
        })
        ;

    py::class_<Positions_c>(m, "Positions_c", py::buffer_protocol())
        .def_buffer([](Positions_c &m) -> py::buffer_info {
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
        .def("__init__", [](Positions_c &m, py::buffer b) {
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
        })
        .def("__getitem__", [](const Positions_c &m, std::pair<ssize_t, ssize_t> i) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            return m(i.first, i.second);
        })
        .def("__setitem__", [](Positions_c &m, std::pair<ssize_t, ssize_t> i, Scalar v) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            m(i.first, i.second) = v;
        })
        .def("get_shape", [](const Positions_c &m) {
        std::pair<int, int> shape( m.rows(), m.cols() );
        return shape;
        })
        ;
    py::class_<Positions_f>(m, "Positions_f", py::buffer_protocol())
        .def_buffer([](Positions_f &m) -> py::buffer_info {
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
        .def("__init__", [](Positions_f &m, py::buffer b) {
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
        })
        .def("__getitem__", [](const Positions_f &m, std::pair<ssize_t, ssize_t> i) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            return m(i.first, i.second);
        })
        .def("__setitem__", [](Positions_f &m, std::pair<ssize_t, ssize_t> i, Scalar v) {
            if (i.first >= m.rows() || i.second >= m.cols())
                throw py::index_error();
            m(i.first, i.second) = v;
        })
        .def("get_shape", [](const Positions_f &m) {
        std::pair<int, int> shape( m.rows(), m.cols() );
        return shape;
        })
        ;


    m.def("cdist",&cdist, R"pbdoc(
        Compute Euclidean distance between 2 set of vectors

        Some other explanation about the pdist function.
    )pbdoc");
    m.def("cdist_pf",&cdist_pf, R"pbdoc(
        Compute Euclidean distance between 2 set of vectors

        Some other explanation about the pdist function.
    )pbdoc");
}