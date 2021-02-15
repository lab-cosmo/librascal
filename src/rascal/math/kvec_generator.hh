/**
 * @file   kvec_generator.cc
 *
 * @author  Andrea Grisafi <andrea.grisafi@epfl.ch>
 * @author  Kevin Kazuki Huguenin-Dumittan <kevin.huguenin-dumittan@epfl.ch>
 *
 * @date   10 February 2021
 *
 * @brief implementation of k-vectors generation
 *
 * Copyright  2021  Andrea Grisafi, Kevin Kazuki Huguenin-Dumittan, COSMO(EPFL)
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

#ifndef SRC_RASCAL_MATH_KVEC_GENERATOR_HH_
#define SRC_RASCAL_MATH_KVEC_GENERATOR_HH_

#include "rascal/math/utils.hh"

namespace rascal {
  namespace math {
    class Kvectors {
     public:
      /*
       * Construct a Kvectors class 
       */
	Kvectors(int n) : nk(0), kval(Vector_t::Zero(n)), kvec(Matrix_t::Zero(n,3)){
        }
      void precompute(Matrix_t basisvectors, double kcut);
      int nk; //number of k-vectors
      Vector_t kval; //k-vectors modulii
      Matrix_t kvec; //k-vectors
    };
  }  // namespace math
}  // namespace rascal

#endif  
