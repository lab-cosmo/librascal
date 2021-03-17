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

#include "rascal/math/kvec_generator.hh"
#include <iostream>

using namespace rascal::math;  // NOLINT


/*
  -----------------------------------------------------------
  Declare some simple functions that return stored properties
  -----------------------------------------------------------
*/

// Number of found wave vectors inside cutoff (without double counting)
size_t Kvectors::GetGridPointNumber() const
{
  return numstored;
}

// Number of found wave vectors including pairs delated by inversion (+k and -k)
// assuming that the origin was excluded to begin with
size_t Kvectors::GetFullGridPointNumber() const
{
  return 2*numstored;
}

// Return wave vectors 
Matrix_t Kvectors::GetGrid() const
{
  return kvecs;
}

// Return indices of wave vectors in terms of the reciprocal basis vectors
Matrix_t Kvectors::GetGridIndices() const
{
  return kvecindices;
}

// Return array containing the norms of the wave vectors
Vector_t Kvectors::GetGridNorms() const
{
  return kvecnorms;
}



// Compute the wave vectors
void Kvectors::precompute(size_t n1max, size_t n2max, size_t n3max,
  Matrix_t basisvecs, double kcut) {
 
  /*  INPUTS
    *   - basisvecs[3][3]: matrix containing the three basis vectors b1, b2, b3 of a 3D lattice, stored in the format b1 = basisvecs[0], b2=basisvecs[1] etc.
    *     a general vector in the lattice will then be of the form k = n1*b1 + n2*b2 + n3*b3, where n1-n3 are integers
    *   - kcut: length of the cutoff
    *
    *   OUTPUTS
    *   - Half the lattice points k excluding the origin lying in a sphere of radius kcut,
    *     where pairs related by k1 = -k2 are only included once.
    *     Mathematically, these are the points satisfying norm(k) <= kcut.
    *     More specifically, all points within the sphere are selected which are in the hemisphere
    *     lying in the b1-direction from the plane spanned by b2 and b3 that cuts the sphere in half.
  */

  // Auxiliary variables
  double kcutsq = kcut*kcut; // Computing the squared norm is faster than the norm
  double normsq; // Store the current squared norm

  /* EXECUTION:
    Begin loops to find the points within the search box contained in the sphere
    In order to avoid double counting pairs of points related by k2=-k1,
    the loops are chosen carefully, and separated into parts dealing with
    the cases
    - n1>0
    - n1=0, n2>0
    - n1=0, n2=0, n3>0
    in reverse order (order of increasing complexity of code)
  */


  // Define basis vectors to make them easier to use in the future
  Eigen::RowVector3d b1 = basisvecs.row(0);
  Eigen::RowVector3d b2 = basisvecs.row(1);
  Eigen::RowVector3d b3 = basisvecs.row(2);
  

  // Initialize current grid point
  // Include the origin for comparison with real space version of SOAP
  Eigen::RowVector3d kvec_new(0,0,0);
  kvec.row(0) = kvec_new;
  nk += 1;
 
  std::cout << "Current wave vector = " << kvec_new << "\n";

  // Step 1: Find all points of the form (0, 0, n3>0)
  for (size_t n3 = 1; n3 <= n3max; ++n3)
  {
    // Update the current grid point
    kvec_new += b3;
	normsq = kvec_new.dot(kvec_new);
	std::cout << "Current indices = (0, 0, " << n3 << ") \n";
	std::cout << "Current wave vector = " << kvec_new << "\n";
    if (normsq <= kcutsq) // Point inside sphere
    {
      kval(nk) = sqrt(normsq);
	  kvec.row(nk) = kvec_new;
      nk += 1;
    }
  }



  // Step 2: Find all points of the form (0, n2>0, n3)
  for (size_t n2 = 1; n2 <= n2max; ++n2)
  {
    // Update current vector for new n2 value
    // We subtract (n3max+1)*b3 s.t. we only have to add b3 at each iteration
    // to get the correct vector
	kvec_new = n2 * b2 - (n3max+1) * b3;

    for (size_t n3 = 0; n3 < 2*n3max+1; ++n3)
    {
      // Update the current grid point
      kvec_new += b3;
	  std::cout << "Current indices = (0, " << n2 << ", " << n3 << ") \n";
	  std::cout << "Current wave vector = " << kvec_new << "\n";
      normsq = kvec_new.dot(kvec_new);
      if (normsq <= kcutsq) // Point inside sphere
      {
        kval(nk) = sqrt(normsq);
        kvec.row(nk) = kvec_new;
        nk += 1;
      }

    } // end loop over n3
  } // end loop over n2



  // Step 3: Find all remaining points of the form (n1>0, n2, n3)
  for (size_t n1 = 1; n1 <= n1max; ++n1)
  {
    for (size_t n2 = 0; n2 < 2*n2max+1; ++n2)
    {
      // Update current vector for new n2 value
      // We subtract (n3max+1)*b3 s.t. we only have to add b3 at each iteration
      // to get the correct vector
      kvec_new = n1*b1 + n2*b2 - n2max*b2 - (n3max+1)*b3;
      for (size_t n3 = 0; n3 < 2*n3max+1; ++n3)
        {
        // Update the current grid point
        kvec_new += b3;
	    std::cout << "Current indices = (" << n1 << ", " << n2 << ", " << n3 << ") \n";
	    std::cout << "Current wave vector = " << kvec_new << "\n";
        normsq = kvec_new.dot(kvec_new);
        if (normsq <= kcutsq) // Point inside sphere
        {
          kval(nk) = sqrt(normsq);
          kvec.row(nk) = kvec_new;
		  nk += 1;
        }

      } // end loop over n3
    } // end loop over n2
  } // end loop over n1


  // Final adjustments
  std::cout << "Final adjustments" << "\n";
  //kvec.conservativeResize(nk+1,3);
  //kval.conservativeResize(nk+1);

}





