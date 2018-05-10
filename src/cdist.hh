/**
 * @file   cdist.h
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  Implementation of a general distance matrix calculation
 *
 * Copyright © 2017 Felix Musil
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */




#include "basic_types.hh"
#include <Eigen/Dense>
#include <pybind11/eigen.h>


namespace py=pybind11;

namespace rascal {

/**
 * cdist calculate the distance matrix between two different input vector containing an arbitrary number of points 
 * with an arbitrary number of dimensions. It is supposed to work with float and results for integer are not guaranteed.
 */


  MatrixXdR cdist(Eigen::Ref<MatrixXdR> X,Eigen::Ref<MatrixXdR> Y);
//	MatrixXdR cdist(py::EigenDRef<Eigen::MatrixXd> X,py::EigenDRef<Eigen::MatrixXd> Y)
	
}

