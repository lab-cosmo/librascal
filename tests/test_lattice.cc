/**
 * @file   test_lattice.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief test implementation of lattice
 *
 * Copyright  2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#include "test_lattice.hh"

#include <boost/test/unit_test.hpp>

namespace rascal {
  constexpr static double lattice_tol{1e-10 * 100};

  BOOST_AUTO_TEST_SUITE(LatticeTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(lattice_constructor_test, ManagerFixtureLattice) {
    Vec3_t cell_lengths = lattice.get_cell_lengths();
    Vec3_t cell_angles = lattice.get_cell_angles();
    Vec3_t cell_lengths_true;
    Vec3_t cell_angles_true;
    cell_lengths_true << 6.190000000000000, 6.6053463194597155,
        7.383806606351496;
    // angle in radian
    cell_angles_true << 1.4313508192414794, 1.5423518764040736,
        1.1973182286702833;

    for (int ii{0}; ii < 3; ++ii) {
      auto error{std::abs(cell_lengths[ii] - cell_lengths_true[ii])};
      BOOST_CHECK_LE(error, lattice_tol);

      error = std::abs(cell_angles[ii] - cell_angles_true[ii]);
      BOOST_CHECK_LE(error, lattice_tol);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(cell_length_test, ManagerFixtureLattice) {
    Vec3_t cell_lengths = lattice.get_cell_lengths();
    Vec3_t cell_lengths_true;
    cell_lengths_true << 6.190000000000000, 6.6053463194597155,
        7.383806606351496;

    for (int ii{0}; ii < 3; ++ii) {
      BOOST_CHECK_EQUAL(cell_lengths[ii], cell_lengths_true[ii]);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(cell_angles_test, ManagerFixtureLattice) {
    Vec3_t cell_angles = lattice.get_cell_angles();
    Vec3_t cell_angles_true;
    // angle in radian
    cell_angles_true << 1.4313508192414794, 1.5423518764040736,
        1.1973182286702833;
    for (int ii{0}; ii < 3; ++ii) {
      auto error{std::abs(cell_angles[ii] - cell_angles_true[ii])};
      BOOST_CHECK_LE(error, lattice_tol);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(transformation_matrix_test, ManagerFixtureLattice) {
    Cell_t cartesian2scaled = lattice.get_cartesian2scaled_matrix();
    Cell_t cartesian2scaled_true;
    // clang-format off
    cartesian2scaled_true <<
       0.161550888529887,  0.00,                          0.000,
      -0.063306933553988,  0.162601626016260,             0.000,
       0.004192528814472, -0.022688598979013, 0.136798905608755;
    // clang-format on
    Cell_t scaled2cartesian = lattice.get_scaled2cartesian_matrix();
    Cell_t scaled2cartesian_true;
    // clang-format off
    scaled2cartesian_true <<
      6.190000000000000,             0.00 ,              0.00,
      2.409999999999999, 6.150000000000001,              0.00,
      0.210000000000000, 1.019999999999999, 7.310000000000000;
    // clang-format on
    for (int jj{0}; jj < 3; ++jj) {
      for (int ii{0}; ii < 3; ++ii) {
        BOOST_CHECK_CLOSE(cartesian2scaled_true(ii, jj),
                          cartesian2scaled(ii, jj), lattice_tol);
        BOOST_CHECK_CLOSE(scaled2cartesian_true(ii, jj),
                          scaled2cartesian(ii, jj), lattice_tol);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(get_cartesian2scaled_test, ManagerFixtureLattice) {
    Eigen::MatrixXd positions_test(3, 22);
    // clang-format off
    positions_test <<
      3.689540159937393, 5.123016813620886, 1.994119731169116,
      6.818437242389163, 2.630056617829216, 6.182500355729062,
      2.114977334498767, 6.697579639059512, 1.392155450018263,
      7.420401523540017, 2.432242071439904, 6.380314902118375,
      1.112656394115962, 7.699900579442317, 3.569715877854675,
      5.242841095703604, 3.122826344932127, 5.689730628626151,
      3.248684682453303, 5.563872291104976, 2.608353462112637,
      6.204203511445642, 5.035681855581504, 2.134827911489532,
      0.946910011088814, 6.223599755982222, 4.168634519120968,
      3.001875247950068, 1.980327734683430, 5.190182032387606,
      2.943861424421339, 4.226648342649697, 5.457161501166098,
      1.713348265904937, 1.501663178733906, 5.668846588337130,
      5.208365510425203, 1.962144256645833, 2.728127406527150,
      4.442382360543885, 2.839975217222644, 4.330534549848392,
      0.744216089807768, 6.426293677263268, 4.643695520786083,
      2.662204050783991, 1.250682335857938, 6.055217235712136,
      0.860905287815103, 6.444994283754972, 4.536108843695142,
      2.769790727874932, 5.609177455068640, 1.696722116501434,
      6.703053268421970, 0.602846303148105, 3.487609972580834,
      3.818289598989240, 1.436734374347541, 5.869165197222533,
      1.054504320562138, 6.251395251007936, 3.998423858825871,
      3.307475712744203, 5.323662899811682, 1.982236671758393;
    // clang-format on
    Vec3_t positions_sc;

    Eigen::MatrixXd positions_sc_true(3, 22);
    // clang-format off
    positions_sc_true <<
      0.296723741750796, 0.703639876645053, 0.267449366982580,
      0.732914251413269, 0.164593885207063, 0.835769733188786,
      0.235325758326898, 0.765037860068951, 0.062053748440047,
      0.938309869955802, 0.075557451185986, 0.924806167209864,
      0.099306843325001, 0.901056775070848, 0.252988672836491,
      0.747374945559359, 0.336207030045438, 0.664156588350411,
      0.361801281872678, 0.638562336523171, 0.396587390963031,
      0.603776227432818, 0.713451112366376, 0.286724809564553,
      0.125592877525700, 0.874583044405229, 0.658294016242431,
      0.341881905688497, 0.219086555224869, 0.781089366706060,
      0.351412276497280, 0.648763645433649, 0.735260445980754,
      0.264915475950175, 0.165043890527786, 0.835132031403143,
      0.814291210823212, 0.185884711107717, 0.419672726629966,
      0.580503195300962, 0.371065952685266, 0.629109969245663,
      0.000224293676929, 0.999951628253999, 0.635252465223814,
      0.364186600654445, 0.171091974809567, 0.828347091068692,
      0.117770901205896, 0.881668164672363, 0.620534725539691,
      0.378904340338568, 0.767329337218692, 0.232109728659567,
      0.916970351357315, 0.082468714520945, 0.477101227439239,
      0.522337838439020, 0.196543690061223, 0.802895375817036,
      0.144255037012604, 0.855184028865655, 0.546980008047315,
      0.452459057830944, 0.728271258524170, 0.271167807354089;
    // clang-format on
    for (int jj{0}; jj < 22; ++jj) {
      lattice.get_cartesian2scaled(positions_test.col(jj), positions_sc);
      for (int ii{0}; ii < 3; ++ii) {
        BOOST_CHECK_CLOSE(positions_sc_true(ii, jj), positions_sc[ii],
                          lattice_tol);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(get_scaled2cartesian_test, ManagerFixtureLattice) {
    Eigen::MatrixXd positions_sc_true(3, 22);
    // clang-format off
    positions_sc_true <<
      3.689540159937393, 5.123016813620886, 1.994119731169116,
      6.818437242389163, 2.630056617829216, 6.182500355729062,
      2.114977334498767, 6.697579639059512, 1.392155450018263,
      7.420401523540017, 2.432242071439904, 6.380314902118375,
      1.112656394115962, 7.699900579442317, 3.569715877854675,
      5.242841095703604, 3.122826344932127, 5.689730628626151,
      3.248684682453303, 5.563872291104976, 2.608353462112637,
      6.204203511445642, 5.035681855581504, 2.134827911489532,
      0.946910011088814, 6.223599755982222, 4.168634519120968,
      3.001875247950068, 1.980327734683430, 5.190182032387606,
      2.943861424421339, 4.226648342649697, 5.457161501166098,
      1.713348265904937, 1.501663178733906, 5.668846588337130,
      5.208365510425203, 1.962144256645833, 2.728127406527150,
      4.442382360543885, 2.839975217222644, 4.330534549848392,
      0.744216089807768, 6.426293677263268, 4.643695520786083,
      2.662204050783991, 1.250682335857938, 6.055217235712136,
      0.860905287815103, 6.444994283754972, 4.536108843695142,
      2.769790727874932, 5.609177455068640, 1.696722116501434,
      6.703053268421970, 0.602846303148105, 3.487609972580834,
      3.818289598989240, 1.436734374347541, 5.869165197222533,
      1.054504320562138, 6.251395251007936, 3.998423858825871,
      3.307475712744203, 5.323662899811682, 1.982236671758393;
    // clang-format o
    Vec3_t position;

    Eigen::MatrixXd positions_sc_test(3, 22);
    // clang-format off
    positions_sc_test <<
      0.296723741750796, 0.703639876645053, 0.267449366982580,
      0.732914251413269, 0.164593885207063, 0.835769733188786,
      0.235325758326898, 0.765037860068951, 0.062053748440047,
      0.938309869955802, 0.075557451185986, 0.924806167209864,
      0.099306843325001, 0.901056775070848, 0.252988672836491,
      0.747374945559359, 0.336207030045438, 0.664156588350411,
      0.361801281872678, 0.638562336523171, 0.396587390963031,
      0.603776227432818, 0.713451112366376, 0.286724809564553,
      0.125592877525700, 0.874583044405229, 0.658294016242431,
      0.341881905688497, 0.219086555224869, 0.781089366706060,
      0.351412276497280, 0.648763645433649, 0.735260445980754,
      0.264915475950175, 0.165043890527786, 0.835132031403143,
      0.814291210823212, 0.185884711107717, 0.419672726629966,
      0.580503195300962, 0.371065952685266, 0.629109969245663,
      0.000224293676929, 0.999951628253999, 0.635252465223814,
      0.364186600654445, 0.171091974809567, 0.828347091068692,
      0.117770901205896, 0.881668164672363, 0.620534725539691,
      0.378904340338568, 0.767329337218692, 0.232109728659567,
      0.916970351357315, 0.082468714520945, 0.477101227439239,
      0.522337838439020, 0.196543690061223, 0.802895375817036,
      0.144255037012604, 0.855184028865655, 0.546980008047315,
      0.452459057830944, 0.728271258524170, 0.271167807354089;
    // clang-format on
    for (int jj{0}; jj < 22; ++jj) {
      lattice.get_scaled2cartesian(positions_sc_test.col(jj), position);
      for (int ii{0}; ii < 3; ++ii) {
        BOOST_CHECK_CLOSE(positions_sc_true(ii, jj), position[ii], lattice_tol);
      }
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(crossproduct_test, ManagerFixtureLattice) {
    Vec3_t v_true;
    v_true << -0.23767568374958059, 0.24116413754264393, -0.063589340753780477;
    Vec3_t v1;
    v1 << 0.57942928, 0.82826958, 0.97551992;
    Vec3_t v2;
    v2 << 0.56475161, 0.69754369, 0.53459899;
    Vec3_t v3;

    lattice.crossproduct(v1, v2, v3);

    for (int ii{0}; ii < 3; ++ii) {
      BOOST_CHECK_CLOSE(v_true[ii], v3[ii], lattice_tol);
    }
  }

  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(set_reciprocal_vectors_test, ManagerFixtureLattice) {
    Cell_t reciprocical_vectors_true;  // 3,22
    // clang-format off
    reciprocical_vectors_true <<
       0.16155088852988697123, -0.0000000000000000,      0.0000000000000000,
      -0.06330693355398815669,  0.16260162601626021450, -0.0000000000000000,
       0.00419252881447217986, -0.02268859897901305545,  0.13679890560875518357;
    // clang-format on
    Cell_t reciprocical_vectors = lattice.get_reciprocal_vectors();

    Vec3_t reciprocal_lengths_true;
    reciprocal_lengths_true << 0.17356276881481585983, 0.16417692074942269453,
        0.13679890560875518357;
    Vec3_t reciprocal_lengths = lattice.get_reciprocal_lengths();

    for (int jj{0}; jj < 3; ++jj) {
      for (int ii{0}; ii < 3; ++ii) {
        BOOST_CHECK_CLOSE(reciprocical_vectors_true(ii, jj),
                          reciprocical_vectors(ii, jj), lattice_tol);
      }
      BOOST_CHECK_CLOSE(reciprocal_lengths_true[jj], reciprocal_lengths[jj],
                        lattice_tol);
    }
  }
  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // namespace rascal
