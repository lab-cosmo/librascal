
#include "basic_types.h"

namespace proteus {

  MatrixXdR cdist(Eigen::Ref<MatrixXdR> X,Eigen::Ref<MatrixXdR> Y) {
    const int N = X.rows();
    const int K = Y.rows();

    MatrixXdR  XX,YY,D;
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

}  // proteus
