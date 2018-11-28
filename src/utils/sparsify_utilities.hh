/**
 * file   sparsify_utilities.hh
 *
 * @author  Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   15 August 2018
 *
 * @brief Header file (declarations) for sparsification utilities
 *
 * Copyright © 2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "basic_types.hh"

namespace rascal {
  namespace utils {
    using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic,
                                      Eigen::Dynamic, Eigen::RowMajor>;

    std::tuple<Eigen::ArrayXi, Eigen::ArrayXd>
    select_fps(const Eigen::Ref<const RowMatrixXd>& feature_matrix,
               int n_sparse=0, int i_first_point=0);

    std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXi, Eigen::ArrayXd>
    select_fps_voronoi(const Eigen::Ref<const RowMatrixXd>& feature_matrix,
                       int n_sparse=0, int i_first_point=0);
  } // namespace utils
} // namespace rascal
