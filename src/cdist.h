/**
 * @file   cdist.h
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  Implementation of a general distance matrix calculation
 *
 * Copyright © 2017 Till Junge
 *
 * Proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Proteus is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */




#include "basic_types.h"
#include <Eigen/Dense>

/**
 * cdist calculate the distance matrix between two different input vector containing an arbitrary number of points 
 * with an arbitrary number of dimensions. It is supposed to work with float and results for integer are not guaranteed.
 */
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


