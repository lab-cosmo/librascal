#include <pybind11/pybind11.h>
#include "basic_types.h"
#include <Eigen/Dense>

namespace py = pybind11;

Matrix cdist(const Eigen::Ref<const Eigen::MatrixX3d>& X,const Eigen::Ref<Eigen::MatrixX3d>& Y)
{  	
    const int N = X.rows();
    const int M = X.cols();
    const int K = Y.rows();
    
    Matrix  XX,YY,D;
    XX.resize(N,1);
    YY.resize(K,1);
    D.resize(N,K);
    XX = X.rowwise().squaredNorm();
    YY = Y.rowwise().squaredNorm();
    D = -2*X*Y.transpose();
    for( int ii= 0; ii < N; ii++ ){
    	for( int jj = 0; jj < K; jj++ ){
    		D(ii,jj) +=  XX(ii) + YY(jj);
    	}
    }
    return Eigen::sqrt(D.array());
}


void prot_dist_mat(py::module& m)
{
	m.def("cdist",&cdist);
}


