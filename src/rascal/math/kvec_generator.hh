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
    private:
      // Quantities that will be provided by the user upon initialization
      double kcutoff{}; // cutoff radius
      Eigen::Matrix3d basisvecs{};// matrix containing the three basis vectors of the cell

      // Quantities that will be computed and be used in further steps (e.g. for LODE)
      size_t numvectors=0; //
      Matrix_t kvecs{};
      //Eigen::Matrix_Xi kvecindices{};
      Vector_t kvecnorms{};

      // Auxiliary variables for internal use
	  double detM{};
      size_t n1_max{};
      size_t n2_max{};
      size_t n3_max{};

	  

    public:
      /*
       * Construct a Kvectors class 
       */
      //explicit Kvectors(size_t n) : nk(0), kval(Vector_t::Zero(n)), kvec(Matrix_t::Zero(n,3)){}
      explicit Kvectors(size_t n){this->n1_max=n;}

      explicit Kvectors(double cutoff, Eigen::Matrix3d basisvectors) 
      {
        this->kcutoff = cutoff;
        this->basisvecs = basisvectors;
      } 
      void precompute(size_t n1max, size_t n2max, size_t n3max, Matrix_t basisvectors, double kcut);
      //size_t nk; //number of k-vectors
      //Vector_t kval; //k-vectors modulii
      //Matrix_t kvec; //k-vectors
      //void precompute2();
      //size_t nk; //number of k-vectors
     
    
	
      // FUNCTIONS FOR USER

      /** Get number of vectors found within cutoff radius
       * without double counting pair related by inversion
      */
      size_t get_numvectors() const {
        return this->numvectors;
      };

      /** Get number of vectors found within cutoff radius
       * counting pairs related by inversion separately
      */
      size_t get_numvectors_all() const {
        return 2*this->numvectors-1;
      };

      /** Get matrix containing all vectors within cutoff radius
       * without double counting, where row(i)= i-th vector 
      */
      Matrix_t get_kvectors() const {
        return this->kvecs;
      };

      /** Get vector containing the norm of all vectors
      */
      Vector_t get_norms() const {
        return this->kvecnorms;
      };



      // FUNCTIONS FOR DEVELOPERS

    };
  }  // namespace math
}  // namespace rascal

#endif  

