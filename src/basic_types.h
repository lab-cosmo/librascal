/**
 * @file   basic_types.h
 *
 * @author Federico Giberti <federico.giberti@epfl.ch>
 *
 * @date   14 Mar 2018
 *
 * @brief  Implementation of base-type objects for Rascal
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


#include <Eigen/Dense>


namespace rascal {

/**
 * Dynamically allocated matrix of arbitrary dimension.
 */

using MatrixXdR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixXdC = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

using VecXb = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using VecXi = Eigen::Matrix<int, Eigen::Dynamic, 1>;



using Matrix3XdC = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>;
using Vecd_t = typename Eigen::Array<double, 3, 1>;
using Veci_t = typename Eigen::Array<int, 3, 1>;
using Mati_t = typename Eigen::Array<int, 3, 2>;
using Cell_t = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
using Vec3_t = Eigen::Matrix<double, 3,1>;
using Vec3i_t = Eigen::Matrix<int, 3, 1>;


using Dim_t = int;

}
