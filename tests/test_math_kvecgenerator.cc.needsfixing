/**
 * @file   test_kvecgenerator.cc
 *
 * @author Kevin Kazuki Huguenin-Dumittan <kevin.huguenin-dumittan@epfl.ch>
 *
 * @date   2021
 *
 * @brief Test the class used to generate the k-vectors, also called
 *        reciprocal space vectors, used in the k-space implementation
 *        of the spherical expansion coefficients of SOAP and LODE
 *
 * Copyright  2019  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "test_math.hh"

#include "rascal/math/kvec_generator.hh"

#include <boost/test/unit_test.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>

namespace rascal {

  BOOST_AUTO_TEST_SUITE(test_kvecgenerator);

  /* ---------------------------------------------------------------------- */
  /**
   * Test that the kvecgenerator produces the desired vectors for small
   * cells, where we can explicitly list the correct vectors
   */
  BOOST_AUTO_TEST_CASE(kvecgen_smallcells_test) {
    // Define reciprocal space cell to be used in the following tests
    // We use the identity matrix to avoid floating point errors
    Eigen::Matrix3d reciprocal_cell;
    reciprocal_cell << 1, 0, 0, 0, 1, 0, 0, 0, 1;

    // Store vectors that should be found, if the code runs correctly
    Eigen::RowVector3d v100, v010, v001, v110, v101, v011, v01m, v10m, v1m0,
        v111, v11m, v1m1, v1mm;
    v100 << 1, 0, 0;
    v010 << 0, 1, 0;
    v001 << 0, 0, 1;
    v110 << 1, 1, 0;
    v101 << 1, 0, 1;
    v011 << 0, 1, 1;
    v01m << 0, 1, -1;
    v10m << 1, 0, -1;
    v1m0 << 1, -1, 0;
    v111 << 1, 1, 1;
    v11m << 1, 1, -1;
    v1m1 << 1, -1, 1;
    v1mm << 1, -1, -1;

    // For small enough cutoff radii, namely for the intervals
    // [1,sqrt(2)), [sqrt(2), sqrt(3)), [sqrt(3),2) we can explicitly
    // list all vector that should be found by the algorithm.
    // Store these solutions in the correct order for each cutoff
    Eigen::MatrixXd kvec_expected_0(3, 3);
    kvec_expected_0 << v001, v010, v100;

    Eigen::MatrixXd kvec_expected_1(9, 3);
    kvec_expected_1 << v001, v01m, v010, v011, v1m0, v10m, v100, v101, v110;

    Eigen::MatrixXd kvec_expected_2(13, 3);
    kvec_expected_2 << v001, v01m, v010, v011, v1mm, v1m0, v1m1, v10m, v100,
        v101, v11m, v110, v111;

    Eigen::MatrixXd * kvecs_expected[3] = {&kvec_expected_0, &kvec_expected_1,
                                           &kvec_expected_2};

    // For various cutoff radii close to the thresholds after which
    // we would get further vectors, we compare the obtained vectors
    // to the correct solution. The variable i runs over the three intervals
    double cutoff_lims[4] = {1., math::SQRT_TWO, math::SQRT_THREE, 2.};
    double cutoff_eps = 1e-12;
    for (int i = 0; i < 3; i++) {
      math::Kvectors k_vectors_min(cutoff_lims[i] + cutoff_eps, reciprocal_cell,
                                   true);
      BOOST_CHECK_EQUAL(k_vectors_min.get_kvectors(), *kvecs_expected[i]);
      math::Kvectors k_vectors_max(cutoff_lims[i + 1] - cutoff_eps,
                                   reciprocal_cell, true);
      BOOST_CHECK_EQUAL(k_vectors_max.get_kvectors(), *kvecs_expected[i]);
    }

    // Check whether the code also works if we use a different, highly
    // distorted basis to describe the same lattice
    Eigen::Matrix3d cell_distorted_basis;
    cell_distorted_basis << 1, 0, 0, 5, 1, 0, 5, 5, 1;

    // Define the vectors that should be obtained from this procedure
    // Note that the obtained vectors should be the same as for the
    // undistorted lattice, but their order and representant out of
    // +k and -k will be different
    Eigen::MatrixXd kvec_expected_3(3, 3), kvec_expected_4(9, 3);
    kvec_expected_3 << v100, -v010, v001;
    kvec_expected_4 << v100, -v110, -v010, v1m0, v011, -v10m, v001, v101, -v01m;
    Eigen::MatrixXd * kvecs_distorted[2] = {&kvec_expected_3, &kvec_expected_4};

    // Run the algorithm and compare the obtained k-vectors to the
    // correct ones
    for (int i = 0; i < 2; i++) {
      math::Kvectors k_vectors_min(cutoff_lims[i] + cutoff_eps,
                                   cell_distorted_basis, true);
      BOOST_CHECK_EQUAL(k_vectors_min.get_kvectors(), *kvecs_distorted[i]);
      math::Kvectors k_vectors_max(cutoff_lims[i + 1] - cutoff_eps,
                                   cell_distorted_basis, true);
      BOOST_CHECK_EQUAL(k_vectors_max.get_kvectors(), *kvecs_distorted[i]);
    }
  }  // end test case

  /* ---------------------------------------------------------------------- */
  /**
   * @brief Check that all obtained k-vectors lie inside the ball
   * and that the norms computed by the kvecgenerator are correct.
   */
  BOOST_FIXTURE_TEST_CASE(kvecgen_allvecs_in_ball_test, KvecgenRefFixture) {
    // For each cell: get the required quantities for the tests
    double cutoff = 6.28;
    for (auto & cell : this->cells) {
      math::Kvectors kvectors(cutoff, cell);
      int nvec = kvectors.get_numvectors();
      auto kvecs = kvectors.get_kvectors();
      auto kvals = kvectors.get_kvector_norms();

      // Start checking the obtained norms
      double cutoffsq = cutoff * cutoff;
      for (int i = 0; i < nvec; i++) {
        auto kvec = kvecs.row(i);
        auto normsq = kvec.dot(kvec);
        BOOST_CHECK_LT(normsq, cutoffsq);
        BOOST_CHECK_CLOSE(normsq, kvals(i) * kvals(i), 1e-13);
      }  // for each k-vector
    }    // for each cell
  }      // end test case

  /* ---------------------------------------------------------------------- */
  /**
   * @brief Check that the number of obtained k-vectors is consistent with
   * the total number of searched vectors as well as the rigorous bound in
   * https://arxiv.org/pdf/math/0410522.pdf in section 3.2.
   */
  BOOST_FIXTURE_TEST_CASE(kvecgen_numvectors_reasonable_test,
                          KvecgenRefFixture) {
    // For each cell: get the required quantities for the tests
    const double cutoff = 6.28;
    for (auto & cell : this->cells) {
      // Get number of obtained k-vectors
      math::Kvectors kvectors(cutoff, cell);
      auto numvectors = kvectors.get_numvectors();

      // Preparation: Calculate various mathematical quantities needed
      // to generate the bounds. We start by recomputing the same auxilary
      // quantities as in the main file. Please check out the main file
      // src/rascal/math/kvec_generator.cc for more details.
      const auto basis_vecs = 2.0 * math::PI * cell.transpose().inverse();
      const auto M = basis_vecs * basis_vecs.transpose();
      auto detM = M.determinant();
      const auto n1max =
          floor(sqrt((M(1, 1) * M(2, 2) - M(1, 2) * M(1, 2)) / detM) * cutoff);
      const auto n2max =
          floor(sqrt((M(0, 0) * M(2, 2) - M(0, 2) * M(0, 2)) / detM) * cutoff);
      const auto n3max =
          floor(sqrt((M(0, 0) * M(1, 1) - M(0, 1) * M(0, 1)) / detM) * cutoff);
      const auto numvectors_searchbox =
          n3max + n2max * (2 * n3max + 1) +
          n1max * (2 * n2max + 1) * (2 * n3max + 1);

      // Compute the estimated number of k-vectors based on the
      // continuous approximation: the volume of the ball divided
      // by the volume of a unit cell. In practice, we divide by
      // a factor of 2 since we only return half of the vectors.
      const auto numvectors_estimate =
          2. * math::PI / 3. * math::pow(cutoff, 3) / sqrt(detM);

      // Compute a factor used for the bounds on the deviation between
      // the estimated computed above and the actual number of k-vectors.
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig_solver;
      const auto eigenvalues = eig_solver.compute(M).eigenvalues();
      const auto bound_factor =
          math::SQRT_THREE / 2. / cbrt(eigenvalues.maxCoeff());

      // RUN ACTUAL TESTS:
      // Number of found vectors has to be less than number in search box
      BOOST_CHECK_LT(numvectors, numvectors_searchbox);

      // Start generating upper and lower bounds from the paper
      // mentioned above.
      auto upperbound = math::pow((1. + bound_factor), 3) * numvectors_estimate;
      auto lowerbound = math::pow((1. - bound_factor), 3) * numvectors_estimate;
      BOOST_CHECK_LE(static_cast<double>(numvectors), upperbound);
      BOOST_CHECK_GE(static_cast<double>(numvectors), lowerbound);
    }  // for each cell
  }    // end test case

  /* ---------------------------------------------------------------------- */
  /**
   * @brief Test correct behavior under scaling
   *
   * If we double the cutoff as well as the basis vectors, the new
   * vectors should be the doubled versions of the original vectors.
   */
  BOOST_FIXTURE_TEST_CASE(kvecgen_scaling_test, KvecgenRefFixture) {
    // Get the required quantities for the tests
    double cutoff = math::PI * 2;
    double normsq;
    for (auto & cell : this->cells) {
      math::Kvectors kvectors(cutoff, cell, false);
      int nvec = kvectors.get_numvectors();
      Eigen::MatrixXd kvecs = kvectors.get_kvectors();
      Eigen::VectorXd kvals = kvectors.get_kvector_norms();

      // For multiple scaling factors, compare the k-vectors obtained
      // before and after scaling. We are only using powers of two to
      // eliminate unwanted behavior due to rounding errors.
      double scales[2] = {0.5, 2};
      for (auto & scale : scales) {
        math::Kvectors kvectors_scaled(cutoff * scale, cell / scale);
        int nvec_scaled = kvectors_scaled.get_numvectors();
        Eigen::MatrixXd kvecs_scaled = kvectors_scaled.get_kvectors() / scale;
        BOOST_CHECK_EQUAL(nvec, nvec_scaled);

        // Check that the obtained vectors agree using squared norm
        if (nvec == nvec_scaled) {
          for (int ik = 0; ik < nvec; ik++) {
            normsq = (kvecs.row(ik) - kvecs_scaled.row(ik)).squaredNorm();
            BOOST_CHECK_LT(normsq, cutoff * cutoff * 1e-14);
          }  // for all k-vectors (of this cell)
        }    // if number of vectors agree
      }      // for all scale factors
    }        // for each cell
  }          // end test case

  /* ---------------------------------------------------------------------- */
  /**
   * @brief Test correct behavior under rotations
   *
   * If we rotate a (real space) cell and compute its k-vectors, the
   * obtained k-vectors should be the same as the rotation applied to
   * the k-vectors of the original cell (stored in the same order)
   * The rotation matrices are generated randomly to avoid systematic
   * biases.
   */
  BOOST_AUTO_TEST_CASE(kvecgen_rotation_test) {
    // To make sure that rounding errors do not accidentally lead
    // to unwanted errors, use a cutoff and cells that are quaranteed
    // to be stable from Legendre's three squares theorem:
    // For an integer (reciprocal) lattice or a sublattice thereof,
    // the squared norm of k is an integer that is a sum of three
    // squares. The numbers in the sequence https://oeis.org/A004215
    // including the number 87 cannot be written as the sum of three
    // squares, thus ensuring that there are no k-vectors whose squared
    // norm satisfies 86 < k^2 < 88. This eliminates potential numerical
    // problems from vectors lying exactly on the boundary.
    // Note that the parameters are still chosen in a way that closely
    // represent typical cells and cutoffs in applications.
    double cutoff = sqrt(87);
    std::vector<Eigen::Matrix3d> cells{};
    Eigen::Matrix3d cell0;
    // Cubic cell
    cell0 << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    cells.push_back(cell0);
    // Tetragonal cell
    cell0 << 2., 0., 0., 0., 1., 0., 0., 0., 1.;
    cells.push_back(cell0);
    // Orthorhombic cell
    cell0 << 2., 0., 0., 0., 1., 0., 0., 0., 3.;
    cells.push_back(cell0);

    // Generate rotation matrices with fixed seed
    srand(4653056);
    size_t num_random = 100;
    std::vector<Eigen::Matrix3d> rotation_matrices{};
    for (size_t index = 0; index < num_random; index++) {
      auto random_matrix = Eigen::Matrix3d::Random();
      Eigen::Matrix3d rotation_matrix_0 =
          random_matrix.householderQr().householderQ();
      rotation_matrices.push_back(rotation_matrix_0);
    }

    // Identity to make sure the rotation matrices are orthogonal
    auto identity_3 = Eigen::MatrixXd::Identity(3, 3);

    for (auto & cell : cells) {
      math::Kvectors kvectors_obj(cutoff, cell, true);
      auto numvectors = kvectors_obj.get_numvectors();
      auto kvectors = kvectors_obj.get_kvectors();
      double total_squared_norm = kvectors.squaredNorm();

      // Start rotations
      for (auto & rotation : rotation_matrices) {
        // Make sure that the rotation matrix is orthogonal
        auto check_orthogonality = rotation * rotation.transpose();
        auto delta_orthogonality = identity_3 - check_orthogonality;
        BOOST_CHECK_LE(delta_orthogonality.array().abs().maxCoeff(), 1e-14);

        // Generate k-vectors of rotated frame
        math::Kvectors kvectors_rotated_obj(cutoff, cell * rotation, true);
        auto numvectors_rotated = kvectors_rotated_obj.get_numvectors();
        BOOST_CHECK_EQUAL(numvectors, numvectors_rotated);

        // Compare obtained vectors
        if (numvectors == numvectors_rotated) {
          auto kvectors_rotated = kvectors_rotated_obj.get_kvectors();
          kvectors_rotated = kvectors_rotated * rotation.transpose();
          auto diff_kvectors = kvectors - kvectors_rotated;
          // 2-norm
          BOOST_CHECK_LE(diff_kvectors.squaredNorm(),
                         total_squared_norm * 1e-15);
          // infinity-norm
          BOOST_CHECK_LE(diff_kvectors.array().abs().maxCoeff(), 1e-13);
        }  // if numvectors agree
      }    // for each rotation matrix
    }      // for each cell
  }        // end test case

  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
