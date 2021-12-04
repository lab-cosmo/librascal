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

#ifndef SRC_RASCAL_MATH_KVEC_GENERATOR_HH_
#define SRC_RASCAL_MATH_KVEC_GENERATOR_HH_

#include "rascal/math/utils.hh"

#include <Eigen/Dense>

namespace rascal {
  namespace math {
    /**
     * Produce the k-vectors (also called reciprocal space or Fourier
     * vectors) needed for the k-space implementation of LODE and SOAP.
     *
     * More specifically, if b1, b2 and b3 are basis vectors of a reciprocal
     * lattice, this class is used to get all vectors of the form:
     * \f{equation}{k = n_1b_1 + n_2b_2 + n_3b_3,\f}
     * where n1,n2,n3 are integers, which are inside a ball of a given
     * cutoff radius \f$k_\mathrm{cut}\f$.
     *
     * Note that for any radius > 0, the origin \f$k=(0,0,0)\f$ is always
     * in the ball and is the first entry of the returned vectors.
     * For applications where k=0 should be excluded, there is the option
     * to do so.
     * For any other k-vector k, symmetry typically allows us to group the
     * two vectors k and -k together. Thus, only one out of each such pair
     * is returned as a representant.
     *
     */
    class Kvectors {
     private:
      /* Store the quantities provided by the user that completely determine
       * the problem:
       * The first two quantities, the cutoff radius in reciprocal space as
       * well as the three basis vectors of the reciprocal cell, determine
       * the geometry, and thus are most important.
       * The optional parameter include_origin is stored as well which
       * determines whether the point (0,0,0) should be included or not.
       */
      double kcutoff{};
      Eigen::Matrix3d basisvecs{};
      bool include_origin{};

      /* Variables for quantities that can be returned by this class:
       * - the N x 3 Eigen matrix containing the k-vectors in the ball
       * - the Eigen vector of size N containing the norms of the k-vectors
       */
      Matrix_t kvectors{};
      Vector_t kvector_norms{};

     public:
      /** Constructor:
       * Upon initialization, directly compute all the k-vectors on the
       * lattice specified by three basisvectors lying within a ball.
       * @param cutoff Cutoff radius of reciprocal space ball.
       *     As a guideline, if the smearing parameter of the density
       *     is given by \f$sigma\f$, use \f$pi/sigma\f$ as cutoff.
       * @param basisvectors 3x3 matrix containing the three basis
       *     vectors a1, a2, a3 of the real space unit cell, stored as:
       *     a1 = basisvectors[0], a2=basisvectors[1], etc.
       *     These are then automatically converted to the corresponding
       *     reciprocal cell vectors b1, b2, b3 by the formulae:
       *     \f{equation}{b_1 = 2\pi \frac{a_1 \times a_2}
       *                 {a_1 \cdot (a_2 \times a_3)} \f}
       *     and similarly for the periodic permutations.
       *     It is also possible to provide the basis vectors of the
       *     reciprocal lattice directly by using the option below.
       * @param is_reciprocal_cell false by default.
       *     If set to true, it indicates that the cell provided
       *     in the variable basisvectors is already the reciprocal cell.
       *     Thus, the formula above for converting the a_i's to the b_i's
       *     is not used, and the provided cell is used directly.
       * @param need_origin false by default.
       *     If set to true, the list of found k-vectors will also
       *     contain the origin k=(0,0,0), which will then be stored
       *     as the very first element in the matrix of found vectors.
       */
      Kvectors(double cutoff, Eigen::Matrix3d basisvectors,
               bool is_reciprocal_cell = false, bool need_origin = false);

      /** Returns number of vectors found within cutoff radius
       * without double counting pairs related by inversion.
       * This is the actual number that can be used to iterate over all
       * vectors and its norms obtained from get_kvectors and get_norms.
       */
      size_t get_numvectors() const { return this->kvector_norms.size(); }

      /** Returns reference to the Eigen matrix containing all
       * k-vectors within cutoff radius without double counting, where
       * row(i) = i-th vector. If the origin k=(0,0,0) is included, it is
       * always stored in row(0).
       */
      Eigen::Ref<Matrix_t> get_kvectors() {
        return Eigen::Ref<Matrix_t>(this->kvectors);
      }

      /** Returns const reference to the Eigen vector containing the norm
       * of all k-vectors within the cutoff without double counting.
       */
      Eigen::Ref<const Vector_t> const get_kvector_norms() {
        return Eigen::Ref<const Vector_t>(this->kvector_norms);
      }
    };  // Kvectors class
  }     // namespace math
}  // namespace rascal

#endif  // SRC_RASCAL_MATH_KVEC_GENERATOR_HH_
