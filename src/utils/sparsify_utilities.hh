/**
 * @file   sparsify_utilities.hh
 *
 * @author  Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   15 August 2018
 *
 * @brief Header file (declarations) for sparsification utilities
 *
 * Copyright  2018  Michele Ceriotti, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_UTILS_SPARSIFY_UTILITIES_HH_
#define SRC_UTILS_SPARSIFY_UTILITIES_HH_

#include "basic_types.hh"

#include <Eigen/Dense>

#include <tuple>

namespace rascal {
  namespace utils {

    using RowMatrixXd =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using FPSReturnTuple =
        std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd>;

    using FPSReturnTupleConst =
        std::tuple<const Eigen::ArrayXi, const Eigen::ArrayXd,
                   const Eigen::ArrayXd>;

    /**
     * Farthest Point Sampling selection of points given the feature matrix
     *
     * @param feature_matrix is a NxD matrix containing N inputs with D features
     *        each. Defaults to numpy row-major type
     * @param n_sparse specifies how many points should be selected. The default
     *        value of zero sorts the whole set of inputs
     * @param i_first_point indicates the index of the first FPS selection.
     *        Defaults to zero
     * @param restart is a tuple containing the outputs of a previous
     *        select_fps call (indices, distance decay, haussdorf dist),
     *        that can be used to restart the calculation.
     *        Defaults to empty arrays, signaling to start from scratch.
     *        If the distance arrays are empty, they are recalculated.
     * @return a tuple containing the list of n_sparse indices of the selected
     *        points, a list of the maximum minimum distance obtained at
     *        each stage, and a list of the Hausdorff distances between all the
     *        inputs and the selected set
     */
    FPSReturnTuple
    select_fps(const Eigen::Ref<const RowMatrixXd> & feature_matrix,
               int n_sparse = 0, int i_first_point = 0,
               const FPSReturnTupleConst & restart = std::make_tuple(
                   Eigen::ArrayXi(0), Eigen::ArrayXd(0), Eigen::ArrayXd(0)));

    /**
     * Farthest Point Sampling selection of points given the feature
     * matrix. Uses a Voronoi cell algorithm that can be faster when selecting
     * many points, or the dimensionality is relatively low.
     *
     * @param feature_matrix is a NxD matrix containing N inputs with D features
     *        each. Defaults to numpy row-major type
     * @param n_sparse specifies how many points should be selected. The default
     *        value of zero sorts the whole set of inputs
     * @param i_first_point indicates the index of the first FPS selection.
     *        Defaults to zero
     * @return a tuple containing the list of n_sparse indices of the selected
     *        points, a list of the maximum minimum distance obtained at each
     *        stage, the list of Hausdorff distances between all the inputs and
     *        the selected set, a list of the assignment of each of the input
     *        points to the FPS selected inputs, and the radius of each Voronoi
     *        cell
     */
    std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd, Eigen::ArrayXi,
               Eigen::ArrayXd>
    select_fps_voronoi(const Eigen::Ref<const RowMatrixXd> & feature_matrix,
                       int n_sparse = 0, int i_first_point = 0);
  }  // namespace utils
}  // namespace rascal
#endif  // SRC_UTILS_SPARSIFY_UTILITIES_HH_
