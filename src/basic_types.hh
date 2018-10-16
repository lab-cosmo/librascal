/**
 * @file   basic_types.hh
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

#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <Eigen/Dense>


namespace rascal {
  //! Dynamically allocated Column Major matrix used for storing positions.
  using Matrix3XdC = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>;
  //! Static double matrix for storing the lattice vectors.
  using Cell_t = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
  //! Static double vector to store a position, an angle etc..
  using Vec3_t = Eigen::Matrix<double, 3,1>;
  //! Reference vector defined using Eigen Map function
  using Vector_ref = Eigen::Map<Vec3_t>;
  /** Static integer vector to store the three indices
   * of either the atom in a triplet or the coordinates of a cell
   * in a supercell
   */
  using Vec3i_t = Eigen::Matrix<int, 3, 1>;
  /**
   * Static integer scalar to store the number of real space
   * dimensions of the system.
   */
  using Dim_t = int;

}

#endif /* BASIC_TYPES_H */