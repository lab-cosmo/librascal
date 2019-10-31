/**
 * @file   sparsify_fps.cc
 *
 * @author  Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   15 August 2018
 *
 * @brief Implementation of Farthest Point Sampling sparsification
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

#include "utils/sparsify_utilities.hh"

#define DO_TIMING 1
#ifdef DO_TIMING
#include <chrono>  // NOLINT
#include <iostream>
#include <tuple>
using hrclock = std::chrono::high_resolution_clock;
#endif

namespace rascal {
  namespace utils {
    FPSReturnTuple
    select_fps(const Eigen::Ref<const RowMatrixXd> & feature_matrix,
               int n_sparse, int i_first_point,
               const FPSReturnTupleConst & restart) {
      // number of inputs
      int n_inputs = feature_matrix.rows();

      // n. of sparse points. defaults to full sorting of the inputs
      if (n_sparse == 0) {
        n_sparse = n_inputs;
      }

      // TODO(ceriottm) <- use the exception mechanism
      // for librascal whatever it is
      if (n_sparse > n_inputs) {
        throw std::runtime_error("Cannot FPS more inputs than those provided");
      }

      /* return arrays */
      // FPS indices
      auto sparse_indices = Eigen::ArrayXi(n_sparse);
      // minmax distances (squared)
      auto sparse_minmax_d2 = Eigen::ArrayXd(n_sparse);

      // square moduli of inputs
      auto feature_x2 = Eigen::ArrayXd(n_inputs);
      // list of square distances to latest FPS point
      auto list_new_d2 = Eigen::ArrayXd(n_inputs);
      // list of minimum square distances to each input
      auto list_min_d2 = Eigen::ArrayXd(n_inputs);
      int i_new{};
      double d2max_new{};

      // computes the squared modulus of input points
      feature_x2 = feature_matrix.cwiseAbs2().rowwise().sum();

      // extracts (possibly empty) restart arrays
      auto i_restart = std::get<0>(restart);
      auto d_restart = std::get<1>(restart);
      auto ld_restart = std::get<2>(restart);
      ssize_t i_start_index = 1;
      if (i_restart.size() > 0) {
        if (i_restart.size() >= n_sparse) {
          throw std::runtime_error("Restart arrays larger than target ");
        }

        sparse_indices.head(i_restart.size()) = i_restart;
        i_start_index = i_restart.size();

        if (d_restart.size() > 0) {
          // restart the FPS calculation from previous run.
          // all information is available
          if (d_restart.size() != i_restart.size()) {
            throw std::runtime_error("Restart indices and distances mismatch");
          }

          // sets the part of the data we know already
          sparse_minmax_d2.head(i_restart.size()) = d_restart;
          list_min_d2 = ld_restart;
        } else {
          // distances are not available, so we recompute them.
          // note that this is as expensive as re-running a full
          // FPS, but it allows us to extend an existing FPS set
          list_new_d2 = feature_x2 + feature_x2(i_restart[0]) -
                        2 * (feature_matrix *
                             feature_matrix.row(i_restart[0]).transpose())
                                .array();
          list_min_d2 = list_new_d2;
          // this is basically the standard algorithm below, only that
          // it is run on the first i_start_index points. see below
          // for comments
          for (ssize_t i = 1; i < i_start_index; ++i) {
            // if the feature matrix has been expanded, the data will
            // not be selected in the same order, so we have to
            // override the selection
            i_new = i_restart[i];
            d2max_new = list_min_d2[i_new];
            /*d2max_new = list_min_d2.maxCoeff(&i_new);
             if (i_new != i_restart[i])
             throw std::runtime_error("Reconstructed distances \
              are inconsistent with restart array");*/
            sparse_indices(i) = i_new;
            sparse_minmax_d2(i - 1) = d2max_new;
            list_new_d2 =
                feature_x2 + feature_x2(i_new) -
                2 * (feature_matrix * feature_matrix.row(i_new).transpose())
                        .array();
            list_min_d2 = list_min_d2.min(list_new_d2);
          }
        }
      } else {
        // standard initialization
        // initializes arrays taking the first point provided in input
        sparse_indices(0) = i_first_point;
        //  distance square to the selected point
        list_new_d2 =
            feature_x2 + feature_x2(i_first_point) -
            2 * (feature_matrix * feature_matrix.row(i_first_point).transpose())
                    .array();
        list_min_d2 = list_new_d2;  // we only have this point....
      }

      for (ssize_t i = i_start_index; i < n_sparse; ++i) {
        // picks max dist and its index
        d2max_new = list_min_d2.maxCoeff(&i_new);
        sparse_indices(i) = i_new;
        sparse_minmax_d2(i - 1) = d2max_new;

        // compute distances^2 to the new point
        // TODO(ceriottm): encapsulate the distance calculation
        // into an interface function
        // dispatching to the proper distance/kernel needed
        list_new_d2 =
            feature_x2 + feature_x2(i_new) -
            2 * (feature_matrix * feature_matrix.row(i_new).transpose())
                    .array();

        // this actually returns a list with the element-wise minimum between
        // list_min_d2(i) and list_new_d2(i)
        list_min_d2 = list_min_d2.min(list_new_d2);
      }
      sparse_minmax_d2(n_sparse - 1) = 0;

      return std::make_tuple(sparse_indices, sparse_minmax_d2, list_min_d2);
    }

    /* ---------------------------------------------------------------------- */
    std::tuple<Eigen::ArrayXi, Eigen::ArrayXd, Eigen::ArrayXd, Eigen::ArrayXi,
               Eigen::ArrayXd>
    select_fps_voronoi(const Eigen::Ref<const RowMatrixXd> & feature_matrix,
                       int n_sparse, int i_first_point) {
      // number of inputs
      int n_inputs = feature_matrix.rows();
      // number of features
      int n_features = feature_matrix.cols();

      // defaults to full sorting of the inputs
      if (n_sparse == 0) {
        n_sparse = n_inputs;
      }

      // TODO(ceriottm) <- use the exception mechanism
      // for librascal whatever it is
      if (n_sparse > n_inputs) {
        throw std::runtime_error("Cannot FPS more inputs than those provided");
      }

      // return arrays
      // FPS indices
      auto sparse_indices = Eigen::ArrayXi(n_sparse);
      // minmax distances^2
      auto sparse_minmax_d2 = Eigen::ArrayXd(n_sparse);
      // size^2 of Voronoi cells
      auto voronoi_r2 = Eigen::ArrayXd(n_sparse);
      // assignment of points to Voronoi cells
      auto voronoi_indices = Eigen::ArrayXi(n_inputs);

      // work arrays
      // index of the maximum-d2 point in each cell
      auto voronoi_i_far = Eigen::ArrayXd(n_sparse);
      // square moduli of inputs
      auto feature_x2 = Eigen::ArrayXd(n_inputs);
      // list of distances^2 to latest FPS point
      auto list_new_d2 = Eigen::ArrayXd(n_inputs);
      // list of minimum distances^2 to each input
      auto list_min_d2 = Eigen::ArrayXd(n_inputs);
      // flags for "active" cells
      auto f_active = Eigen::ArrayXi(n_sparse);
      // list of dist^2/4 to previously selected points
      auto list_sel_d2q = Eigen::ArrayXd(n_sparse);
      // feaures of the latest FPS point
      auto feature_new = Eigen::VectorXd(n_features);
      // matrix of the features for the active point selection
      auto feature_sel = RowMatrixXd(n_sparse, n_features);

      int i_new{};
      double d2max_new{};

      // computes the squared modulus of input points
      feature_x2 = feature_matrix.cwiseAbs2().rowwise().sum();

      // initializes arrays taking the first point provided in input
      sparse_indices(0) = i_first_point;
      //  distance square to the selected point
      list_new_d2 =
          feature_x2 + feature_x2(i_first_point) -
          2 * (feature_matrix * feature_matrix.row(i_first_point).transpose())
                  .array();
      list_min_d2 = list_new_d2;  // we only have this point....

      voronoi_r2 = 0.0;
      voronoi_indices = 0;
      // picks the initial Voronoi radius and the farthest point index
      voronoi_r2(0) = list_min_d2.maxCoeff(&voronoi_i_far(0));

      feature_sel.row(0) = feature_matrix.row(i_first_point);

#ifdef DO_TIMING
      // timing code
      double tmax{0}, tactive{0}, tloop{0};
      int64_t ndist_eval{0}, npoint_skip{0}, ndist_active{0};
      auto gtstart = hrclock::now();
#endif
      for (int i = 1; i < n_sparse; ++i) {
#ifdef DO_TIMING
        auto tstart = hrclock::now();
#endif
        /*
         * find the maximum minimum distance and the corresponding point.  this
         * is our next FPS. The maxmin point must be one of the voronoi
         * radii. So we pick it from this smaller array. Note we only act on the
         * first i items as the array is filled incrementally picks max dist and
         * index of the cell
         */
        d2max_new = voronoi_r2.head(i).maxCoeff(&i_new);
        // the actual index of the fartest point
        i_new = voronoi_i_far(i_new);
#ifdef DO_TIMING
        auto tend = hrclock::now();
        tmax +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tstart)
                .count();
#endif
        // store properties of the new FPS selection
        sparse_indices(i) = i_new;
        sparse_minmax_d2(i - 1) = d2max_new;
        feature_new = feature_matrix.row(i_new);
        /*
         * we store indices of the selected features because we can then compute
         * some of the distances with contiguous array operations
         */
        feature_sel.row(i) = feature_new;

        /*
         * now we find the "active" Voronoi cells, i.e. those
         * that might change due to the new selection.
         */
        f_active = 0;

#ifdef DO_TIMING
        tstart = hrclock::now();
        ndist_active += i;
#endif
        /*
         * must compute distance of the new point to all the previous FPS.  some
         * of these might have been computed already, but bookkeeping could be
         * worse that recomputing (TODO: verify!)
         */
        list_sel_d2q.head(i) =
            feature_x2(i_new) -
            2 * (feature_sel.topRows(i) * feature_new).array();
        for (ssize_t j = 0; j < i; ++j) {
          list_sel_d2q(j) += feature_x2(sparse_indices(j));
        }
        list_sel_d2q.head(i) *= 0.25;  // triangle inequality: voronoi_r < d/2

        for (ssize_t j = 0; j < i; ++j) {
          /*
           * computes distances to previously selected points and uses triangle
           * inequality to find which voronoi sets might be affected by the
           * newly selected point divide by four so we don't have to do that
           * later to speed up the bound on distance to the new point
           */
          if (list_sel_d2q(j) < voronoi_r2(j)) {
            f_active(j) = 1;
            //! size of active cells will have to be recomputed
            voronoi_r2(j) = 0;
#ifdef DO_TIMING
          } else {
            ++npoint_skip;
#endif
          }
        }

#ifdef DO_TIMING
        tend = hrclock::now();
        tactive +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tstart)
                .count();

        tstart = hrclock::now();
#endif

        for (ssize_t j = 0; j < n_inputs; ++j) {
          int voronoi_idx_j = voronoi_indices(j);
          // only considers "active" points
          if (f_active(voronoi_idx_j) > 0) {
            /*
             * check if we can skip this check for point j. this is a tighter
             * bound on the distance, since |x_j-x_sel|<rvoronoi_sel
             */
            if (list_sel_d2q(voronoi_idx_j) < list_min_d2(j)) {
              double d2_j = feature_x2(i_new) + feature_x2(j) -
                            2 * feature_new.dot(feature_matrix.row(j));
#ifdef DO_TIMING
              ndist_eval++;
#endif
              /*
               * we have to reassign point j to the new selection. also, the
               * voronoi center is actually that of the new selection
               */
              if (d2_j < list_min_d2(j)) {
                list_min_d2(j) = d2_j;
                voronoi_indices(j) = voronoi_idx_j = i;
              }
            }
            // also must update the voronoi radius
            if (list_min_d2(j) > voronoi_r2(voronoi_idx_j)) {
              voronoi_r2(voronoi_idx_j) = list_min_d2(j);
              // stores the index of the FP of the cell
              voronoi_i_far(voronoi_idx_j) = j;
            }
          }
        }

#ifdef DO_TIMING
        tend = hrclock::now();
        tloop +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tstart)
                .count();
#endif
      }
      sparse_minmax_d2(n_sparse - 1) = 0;

#ifdef DO_TIMING
      auto gtend = hrclock::now();

      std::cout << "Skipped " << npoint_skip << " FPS centers of "
                << n_sparse * (n_sparse - 1) / 2 << " - "
                << npoint_skip * 100. / (n_sparse * (n_sparse - 1) / 2)
                << "%\n";
      std::cout << "Computed " << ndist_eval << " distances rather than "
                << n_inputs * n_sparse << " - "
                << ndist_eval * 100. / (n_inputs * n_sparse) << " %\n";

      std::cout << "Time total "
                << std::chrono::duration_cast<std::chrono::nanoseconds>(gtend -
                                                                        gtstart)
                           .count() *
                       1e-9
                << "\n";
      std::cout << "Time looking for max " << tmax * 1e-9 << "\n";
      std::cout << "Time looking for active " << tactive * 1e-9 << " with "
                << ndist_active << " distances\n";
      std::cout << "Time general loop " << tloop * 1e-9 << "\n";
#endif

      return std::make_tuple(sparse_indices, sparse_minmax_d2, list_min_d2,
                             voronoi_indices, voronoi_r2);
    }

  }  // namespace utils
}  // namespace rascal
