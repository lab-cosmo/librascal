/**
 * @file   kvec_generator.cc
 *
 * @author  Andrea Grisafi <andrea.grisafi@epfl.ch>
 * @author  Kevin Kazuki Huguenin-Dumittan <kevin.huguenin-dumittan@epfl.ch>
 *
 * @date   10 February 2021
 *
 * @brief Generate the k-vectors (also called reciprocal or Fourier vectors)
 *        needed for the k-space implementation of LODE and SOAP.
 *        More specifically, these are all points of a (reciprocal space)
 *        lattice that lie within a ball of a specified cutoff radius.
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

#include <assert.h>

using namespace rascal::math;  // NOLINT

/** Function for finding and storing all reciprocal space vectors
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
 * Computes half of the lattice points k excluding the origin lying in
 * a ball of radius kcut, where pairs related by k1 = -k2 are only included
 * once. Mathematically, these are the vectors s.t. norm(k) <= kcut.
 * More specifically, all k-vectors within the ball are selected which are
 * lying in the b1-direction from the plane spanned by b2 and b3.
 *
 * DESCRIPTION OF ALGORITHM
 * Essentially, we find the vectors by performing a brute force search of
 * all vectors contained in a box. The first part of the code is used to
 * determine the optimal geometry of this box, while the second part loops
 * over all vectors within the box.
 *
 * The key idea is that a k-vector is parametrized by three integers via
 * k = n1*b1 + n2*b2 + n3*b3. Using these integers, and the Einstein summation
 * convention, we can write the squared norm of k as:
 * k*k = M_ij n_i n_j, where M_ij is the 3x3 matrix M_ij = b_i * b_j.
 * The rest of the algorithm is performed directly in this index-space, i.e.
 * the space of (n1, n2, n3).
 *
 * In the brute force part, we would then like to loop over all possible
 * values of (n1, n2, n3), where 0<n1<n1max, |n2|<n2max and |n3|<n3max
 * to ensure that only half of the points are found (+ some corrections
 * along the boundary n1=0, as described below before the loop starts) for
 * some n1max, n2max and n3max.
 *
 * The preparation part is used to determine the values of njmax that are
 * optimal, i.e. contains the smallest amount of redundancy. We can calculate
 * this analytically in the following way:
 * To get the bounds for n1, fix the value of n1 and consider the restriction
 * of the bilinear form M_ij n_i n_j onto the subspace spanned by (n2, n3).
 * The equation k*k = M_ij n_i n_j can then be recast into the form:
 * k*k + (constant part depending on n1) = quadratic form in (n2, n3).
 * For a solution to exist, the left hand side of this equation has to be
 * positive, which provides us with an inequality that results in a bound
 * for n1, namely n1max. Similarly, we obtain n2max and n3max.
 */
Kvectors::Kvectors(double cutoff, Eigen::Matrix3d basisvectors,
                   bool is_reciprocal_cell, bool need_origin) {
  // INITIALIZATION
  // Set internal variables apart from cell
  this->kcutoff = cutoff;
  this->include_origin = need_origin;

  /* Store provided basis vectors into the internal variable basisvecs.
   * We need to distinguish whether the provided vectors are already
   * reciprocal cell vectors or need to be converted first
   */
  if (not is_reciprocal_cell) {  // real space cell is given
    // Create reciprocal cell from real space cell
    Eigen::Matrix3d tcell = basisvectors.transpose();
    this->basisvecs = 2.0 * math::PI * tcell.inverse();
  } else {  // cell of reciprocal space is already given
    this->basisvecs = basisvectors;
  }

  /* PREPARATION: DEFINE SEARCH SPACE BOX
   * Determine the optimal bounds n1max, n2max, n3max that
   * define the search space box. Roughly speaking, our
   * goal will then be to find all vectors
   * k = n1*b1 + n2*b2 + n3*b3,
   * where |n1|<n1max, |n2|<n2max etc., whose norm is
   * smaller than kcut. In practice, only half of the
   * k-vectors are returned up to the identification
   * k2 = -k1, since such pairs can be grouped together in
   * the remaining part.
   */
  // Define some const shortcuts for quantities that shouldn't change
  // during the precomputation step
  const Eigen::Matrix3d bvecs = this->basisvecs;
  const double kcut = this->kcutoff;

  // Generate inner product matrix M_ij = b_i * b_j and the volume of
  // a cell in reciprocal space to compute the bounds of the optimal
  // search space box
  const Matrix_t M = bvecs * bvecs.transpose();
  const double kcell_vol = basisvecs.determinant();
  const size_t n1max =
      floor(sqrt(M(1, 1) * M(2, 2) - M(1, 2) * M(1, 2)) / kcell_vol * kcut);
  const size_t n2max =
      floor(sqrt(M(0, 0) * M(2, 2) - M(0, 2) * M(0, 2)) / kcell_vol * kcut);
  const size_t n3max =
      floor(sqrt(M(0, 0) * M(1, 1) - M(0, 1) * M(0, 1)) / kcell_vol * kcut);

  // Total number of points that will be checked (used for initialization)
  size_t numvectors_searchbox = n3max + n2max * (2 * n3max + 1) +
                                n1max * (2 * n2max + 1) * (2 * n3max + 1);
  // Add origin if needed
  if (this->include_origin) {
    numvectors_searchbox += 1;
  }

  /* Auxiliary variables: using the squared norm rather than norm itself is
   * computationally more efficient. Use one variable each for the squared
   * cutoff and the squared norm of the k-vectors (to be updated during search)
   */
  const double kcutsq = kcut * kcut;
  double normsq;

  // Store basis vectors for increased readability
  // and ease of use
  const Eigen::RowVector3d b1 = basisvecs.row(0);
  const Eigen::RowVector3d b2 = basisvecs.row(1);
  const Eigen::RowVector3d b3 = basisvecs.row(2);

  // Initialize current vector and number of found vectors
  size_t numvectors = 0;
  Eigen::RowVector3d kvec_new(0, 0, 0);

  // Initialize arrays in which to store results
  kvectors.resize(numvectors_searchbox, 3);
  kvector_norms.resize(numvectors_searchbox);

  /* START OF MAIN PART
   * Begin loops to find the points within the search box
   * contained in the ball. In order to avoid double counting
   * pairs of points related by k2=-k1, the loops are chosen
   * carefully, and separated into parts dealing with the cases
   * - n1>0
   * - n1=0, n2>0
   * - n1=0, n2=0, n3>0
   *  in reverse order (order of increasing complexity of code)
   */

  // Step 0: If desired (e.g. for SOAP), include origin:
  if (this->include_origin) {
    kvectors.row(0) = kvec_new;
    kvector_norms(0) = 0.;
    numvectors += 1;
  }

  // Step 1: Find all points of the form (0, 0, n3>0)
  for (size_t n3 = 1; n3 <= n3max; ++n3) {
    // Update the current vector
    kvec_new += b3;
    normsq = kvec_new.dot(kvec_new);

    // If vector is inside ball, add it
    if (normsq <= kcutsq) {
      kvectors.row(numvectors) = kvec_new;
      kvector_norms(numvectors) = sqrt(normsq);
      numvectors += 1;
    }  // end if
  }    // end loop over n3

  // Step 2: Find all points of the form (0, n2>0, n3)
  for (size_t n2 = 1; n2 <= n2max; ++n2) {
    // Update current vector for new n2 value
    // We subtract (n3max+1)*b3 s.t. we only have to add b3
    // at each iteration to get the correct vector
    kvec_new = n2 * b2 - (n3max + 1) * b3;

    for (size_t n3 = 0; n3 < 2 * n3max + 1; ++n3) {
      // Update the current vector
      kvec_new += b3;
      normsq = kvec_new.dot(kvec_new);

      // If vector is inside ball, add it
      if (normsq <= kcutsq) {
        kvectors.row(numvectors) = kvec_new;
        kvector_norms(numvectors) = sqrt(normsq);
        numvectors += 1;
      }  // end if
    }    // end loop over n3
  }      // end loop over n2

  // Step 3: Find all remaining points of the form (n1>0, n2, n3)
  for (size_t n1 = 1; n1 <= n1max; ++n1) {
    for (size_t n2 = 0; n2 < 2 * n2max + 1; ++n2) {
      // Update current vector for new n2 value
      // We subtract (n3max+1)*b3 s.t. we only have to add b3
      // at each iteration to get the desired vector
      kvec_new = n1 * b1 + n2 * b2 - n2max * b2 - (n3max + 1) * b3;

      for (size_t n3 = 0; n3 < 2 * n3max + 1; ++n3) {
        // Update the current vector
        kvec_new += b3;
        normsq = kvec_new.dot(kvec_new);

        // If vector is inside ball, add it
        if (normsq <= kcutsq) {
          kvectors.row(numvectors) = kvec_new;
          kvector_norms(numvectors) = sqrt(normsq);
          numvectors += 1;
        }  // end if
      }    // end loop over n3
    }      // end loop over n2
  }        // end loop over n1

  // Final adjustments: Get rid of the empty components
  assert(numvectors != 0);
  kvectors.conservativeResize(numvectors, 3);
  kvector_norms.conservativeResize(numvectors);
}
