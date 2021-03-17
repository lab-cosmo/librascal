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



/** Class for finding and storing all reciprocal space vectors
 * of a lattice whose length is smaller than a cutoff.
 * INPUTS
 * - basisvecs[3][3]: 3x3 matrix containing the three basis
 *   vectors b1, b2, b3 of a 3D lattice, stored in the format
 *   b1 = basisvecs[0], b2=basisvecs[1] etc.
 *   a general vector in the lattice will then be of the form
 *   k = n1*b1 + n2*b2 + n3*b3, where n1,n2,n3 are integers
 * - kcut: length of the cutoff
 *
 * OUTPUTS
 * - The precompute method computes half the lattice points k
 *   excluding the origin lying in a sphere of radius kcut,
 *   where pairs related by k1 = -k2 are only included once.
 *   Mathematically, these are the vector s.t. norm(k) <= kcut.
 *   More specifically, all points within the sphere are
 *   selected which are in the hemisphere lying in the
 *   b1-direction from the plane spanned by b2 and b3 that
 *   cuts the sphere in half.
 * - If the cell is updated (e.g. during a NpT simultion)
 *   the wave vectors can be updated in two ways:
 *   1. Redo the same computations from scratch, finding all
 *      vectors s.t. norm(k) <= kcut again.
 *   2. Continuously adjust the vectors to the new cell.
 *      If the old cell had basis b1,b2,b3 and the new one
 *      the basis b1', b2', b3', then the vector
 *      k=n1b1+n2b2+n3b3 is updated to n1b1' + n2b2' + n3b3'
 *      Use this option during NpT simulations for small
 *      changes in the cell to ensure continuity of the
 *      ML descriptors.
*/



void Kvectors::precompute() {
  // PREPARATION
  // Generate optimal bounds for search space box 
  Eigen::Matrix3d bvecs = this->basisvecs;
  double kcut = this -> kcutoff;
  Matrix_t M = bvecs * bvecs.transpose(); // = inner product matrix M_ij = b_i * b_j
  double detM = M(0,0)*M(1,1)*M(2,2)+2*M(0,1)*M(1,2)*M(2,0);
  detM -= M(0,0)*M(1,2)*M(1,2);
  detM -= M(1,1)*M(0,2)*M(0,2);
  detM -= M(2,2)*M(0,1)*M(0,1);

  size_t n1max = floor(sqrt( (M(1,1)*M(2,2) - M(1,2)*M(1,2)) / detM) * kcut);
  size_t n2max = floor(sqrt( (M(0,0)*M(2,2) - M(0,2)*M(0,2)) / detM) * kcut);
  size_t n3max = floor(sqrt( (M(0,0)*M(1,1) - M(0,1)*M(0,1)) / detM) * kcut);
  
  this->n1_max = n1max;
  this->n2_max = n2max;
  this->n3_max = n3max;

// Total number of points that will be checked (used for initialization)
  size_t numtot = 1+ n3max + n2max * (2 * n3max + 1) + n1max * (2 * n2max + 1) * (2 * n3max + 1);
 
  // Auxiliary variables
  double kcutsq = kcut*kcut; // use squared norm for tests
  double normsq; // Store the current squared norm
  
  // Store basis vectors for increased readability
  // and ease of use
  Eigen::RowVector3d b1 = basisvecs.row(0);
  Eigen::RowVector3d b2 = basisvecs.row(1);
  Eigen::RowVector3d b3 = basisvecs.row(2);

  // Initialize current grid point
  Eigen::RowVector3d kvec_new(0,0,0);
  //Eigen::RowVector3i kvec_ind(0,0,0);

  // Initialize arrays in which to store results
  kvecs.resize(numtot,3);
  kvecnorms.resize(numtot);
  
  /* EXECUTION::w

    Begin loops to find the points within the search box
	contained in the sphere. In order to avoid double counting
	pairs of points related by k2=-k1, the loops are chosen
	carefully, and separated into parts dealing with the cases
    - n1>0
    - n1=0, n2>0
    - n1=0, n2=0, n3>0
    in reverse order (order of increasing complexity of code)
  */

  // Step 0: Include the origin for comparison with real space
  // version of SOAP (optional)
  //kvecindices.row(0) = kvec_ind;
  kvecs.row(0) = kvec_new;
  kvecnorms(0) = 0.;
  numvectors += 1;
 
  // Step 1: Find all points of the form (0, 0, n3>0)
  for (size_t n3 = 1; n3 <= n3max; ++n3)
  {
    // Update the current grid point
    kvec_new += b3;
	normsq = kvec_new.dot(kvec_new);
    
	if (normsq <= kcutsq) // Point inside sphere
    {
      kvecs.row(numvectors) = kvec_new;
      kvecnorms(numvectors) = sqrt(normsq);
      numvectors += 1;
    }
  } // end loop over n3


  // Step 2: Find all points of the form (0, n2>0, n3)
  for (size_t n2 = 1; n2 <= n2max; ++n2)
  {
    // Update current vector for new n2 value
    // We subtract (n3max+1)*b3 s.t. we only have to add b3
	// at each iteration to get the correct vector
	kvec_new = n2 * b2 - (n3max+1) * b3;

    for (size_t n3 = 0; n3 < 2*n3max+1; ++n3)
    {
      // Update the current grid point
      kvec_new += b3;
      normsq = kvec_new.dot(kvec_new);
      
	  if (normsq <= kcutsq) // Point inside sphere
      {
        kvecs.row(numvectors) = kvec_new;
        kvecnorms(numvectors) = sqrt(normsq);
        numvectors += 1;
      }

    } // end loop over n3
  } // end loop over n2



  // Step 3: Find all remaining points of the form (n1>0, n2, n3)
  for (size_t n1 = 1; n1 <= n1max; ++n1)
  {
    for (size_t n2 = 0; n2 < 2*n2max+1; ++n2)
    {
      // Update current vector for new n2 value
      // We subtract (n3max+1)*b3 s.t. we only have to add b3
	  // at each iteration to get the desired vector
      kvec_new = n1*b1 + n2*b2 - n2max*b2 - (n3max+1)*b3;
      for (size_t n3 = 0; n3 < 2*n3max+1; ++n3)
        {
        // Update the current grid point
        kvec_new += b3;
        normsq = kvec_new.dot(kvec_new);

        if (normsq <= kcutsq) // Point inside sphere
        {
          kvecs.row(numvectors) = kvec_new;
          kvecnorms(numvectors) = sqrt(normsq);
          numvectors += 1;
        }

      } // end loop over n3
    } // end loop over n2
  } // end loop over n1


  // Final adjustments: Get rid of the empty components
  kvecs.conservativeResize(numvectors,3);
  kvecnorms.conservativeResize(numvectors);
}

