/**
 * @file   sparsify_utilities.hh
 *
 * @author  Michele Ceriotti <michele.ceriotti@gmail.com>
 *
 * @date   15 August 2018
 *
 * @brief Tests for sparsification routines
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

#include <boost/test/unit_test.hpp>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(UtilsSparsify);

  BOOST_AUTO_TEST_CASE(sparsify_fps_test) {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        feature_matrix(10, 4);
    feature_matrix << -1.75385227, -0.23250687, 0.57691232, -1.01651329,
        0.06937612, 0.87558439, 0.40366658, 0.53847493, 0.11625195, -0.40740825,
        0.02985093, 0.04758038, 1.13182798, 0.51740637, 1.93431393, 1.66922457,
        1.27045019, -0.43015544, -0.24192861, 0.3747168, -0.07650897,
        -0.2354746, 0.11482844, -0.07034375, 1.83259473, 0.04583692, 0.94401602,
        1.19349946, 0.29391326, 2.17371398, 0.42916834, 0.15842183, 0.73327368,
        0.86531859, -0.75413054, 0.84116485, 0.18576865, -0.00515566,
        -1.56720151, -0.57342492;

    std::vector<int> fps_order = {0, 6, 9, 7, 8};

    Eigen::ArrayXi ifps;
    Eigen::ArrayXd dfps;
    Eigen::ArrayXd ldmin2;
    std::tie(ifps, dfps, ldmin2) = utils::select_fps(feature_matrix, 5, 0);

    for (size_t i = 0; i < fps_order.size(); ++i) {
      BOOST_CHECK_EQUAL(ifps(i), fps_order[i]);
    }

    Eigen::ArrayXi v_ifps;
    Eigen::ArrayXd v_dfps;
    Eigen::ArrayXd v_ldmin2;
    Eigen::ArrayXi v_ivor;
    Eigen::ArrayXd v_rvor;
    std::tie(v_ifps, v_dfps, v_ldmin2, v_ivor, v_rvor) =
        utils::select_fps_voronoi(feature_matrix, 5, 0);

    for (size_t i = 0; i < fps_order.size(); ++i) {
      BOOST_CHECK_EQUAL(ifps(i), v_ifps[i]);
      BOOST_CHECK_EQUAL(dfps(i), v_dfps[i]);
    }

    for (size_t i = 0; i < (size_t)ldmin2.size(); ++i) {
      BOOST_CHECK_EQUAL(ldmin2(i), v_ldmin2[i]);
    }
  }

  BOOST_AUTO_TEST_SUITE_END();
}  // namespace rascal
