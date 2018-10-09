/**
 * file   test_lattice.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief test implementation of lattice
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "tests.hh"
#include "test_math.hh"

namespace rascal {
  constexpr static double math_tol{1e-14};

  BOOST_AUTO_TEST_SUITE(MathCEPHESTests);
  /* ---------------------------------------------------------------------- */
  BOOST_FIXTURE_TEST_CASE(math_hyp2f1_test, ManagerFixture_math) {
    bool vebose{true};
    for (int ii{0};ii<3;++ii) {
      double val{hyp2f1(numbers(0,ii),numbers(1,ii),numbers(2,ii),numbers(3,ii))};
      if (vebose){
        std::cout << std::setprecision(14) << val <<", ";
      }
      auto error{std::abs(val-results_hyp2f1[ii])};
      BOOST_CHECK_LE(error, math_tol);
    }

  }

  BOOST_FIXTURE_TEST_CASE(math_airy_test, ManagerFixture_math) {
    bool vebose{true};
    for (int ii{0};ii<3;++ii) {
      double Ai{0}, Aip{0}, Bi{0}, Bip{0};
      airy(numbers(0,ii),&Ai,&Aip,&Bi,&Bip);
      if (vebose){
        std::cout << std::setprecision(14) << Ai <<", "<<Aip <<", "<<Bi <<", "<<Bip <<", ";
      }
      auto error{std::abs(Ai-results_airy(ii,0))};
      BOOST_CHECK_LE(error, math_tol);
      error = std::abs(Aip-results_airy(ii,1));
      BOOST_CHECK_LE(error, math_tol);
      error = std::abs(Bi-results_airy(ii,2));
      BOOST_CHECK_LE(error, math_tol);
      error = std::abs(Bip-results_airy(ii,3));
      BOOST_CHECK_LE(error, math_tol);
    }

  }
  
  /* ---------------------------------------------------------------------- */
  BOOST_AUTO_TEST_SUITE_END();

}  // rascal
