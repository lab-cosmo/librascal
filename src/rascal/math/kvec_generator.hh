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
      explicit Kvectors(size_t n) : nk(0), kval(Vector_t::Zero(n)), kvec(Matrix_t::Zero(n,3)){
        }
      void precompute(size_t n1max, size_t n2max, size_t n3max, Matrix_t basisvectors, double kcut);
      size_t nk; //number of k-vectors
      Vector_t kval; //k-vectors modulii
      Matrix_t kvec; //k-vectors
      
	  // FUnctions for user
      size_t GetGridPointNumber() const; // Get number of grid points stored (e.g. to loop over them)
      size_t GetFullGridPointNumber() const; // Get full number of grid points = twice the above
      Matrix_t GetGrid() const;
      Matrix_t GetGridIndices() const;
      Vector_t GetGridNorms() const;


	  private:
        double rad=0;
		Matrix_t basisvecs = Matrix_t::Zero(3,3);

		size_t numstored=0;
		Matrix_t kvecs=Matrix_t::Zero(2,3);
		Matrix_t kvecindices = Matrix_t::Zero(2,3);
		Vector_t kvecnorms = Vector_t::Zero(1);

		// Auxiliary variables for development
		double detM=0;
		size_t n1_max=0;
		size_t n2_max=0;
		size_t n3_max=0;

    };
  }  // namespace math
}  // namespace rascal

#endif  
