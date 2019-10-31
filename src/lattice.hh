/**
 * @file   lattice.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Lattice class used to compute face distances within the cell
 * and to scale and unscale positions
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

#ifndef SRC_LATTICE_HH_
#define SRC_LATTICE_HH_

#include "atomic_structure.hh"

#include <Eigen/Dense>

#include <cmath>

namespace rascal {
  /**
   * Class to store and change between different lattice representation (real
   * and reciprocal space in terms of lattice vectors, cell length and
   * angles). Also translate absolute to fractional and fractional to absolute
   * coordinates.
   */
  template <size_t Dim>
  class Lattice {
   public:
    // TODO(felix): implement the specialization for the 2D case
    static_assert(Dim == 3, "2D lattice is not implemented yet");
    using Cell_t = typename AtomicStructure<Dim>::Cell_t;
    using AtomTypes_t = typename AtomicStructure<Dim>::AtomTypes_t;
    using PBC_t = typename AtomicStructure<Dim>::PBC_t;
    using Positions_t = typename AtomicStructure<Dim>::Positions_t;
    using Vec_t = Eigen::Matrix<double, Dim, 1>;

    //! Default constructor
    Lattice() = default;

    //! Initializes the cell via set_cell function
    explicit Lattice(const Cell_t & cell) { this->set_cell(cell); }

    //! Copy constructor
    Lattice(const Lattice & other) = delete;

    //! Move constructor
    Lattice(Lattice && other) = default;

    //! Destructor
    virtual ~Lattice() = default;

    //! Copy assignment operator
    Lattice & operator=(const Lattice & other) = delete;

    //! Move assignment operator
    Lattice & operator=(Lattice && other) = default;

    /* ---------------------------------------------------------------------- */
    /**
     * Calculates the cell lengths and angles given the lattice vectors in a 3x3
     * matrix, as well as the reciprocal vectors and the transformation matrix.
     */
    void set_cell(const Cell_t & cell) {
      this->cell_vectors = cell;
      this->cell_lengths = cell.colwise().norm();
      this->cell_angles[0] =
          std::acos(this->cell_vectors.col(1).dot(this->cell_vectors.col(2)) /
                    this->cell_lengths[1] / this->cell_lengths[2]);
      this->cell_angles[1] =
          std::acos(this->cell_vectors.col(0).dot(this->cell_vectors.col(2)) /
                    this->cell_lengths[0] / this->cell_lengths[2]);
      this->cell_angles[2] =
          std::acos(this->cell_vectors.col(1).dot(this->cell_vectors.col(0)) /
                    this->cell_lengths[1] / this->cell_lengths[0]);
      this->set_transformation_matrix();
      this->set_reciprocal_vectors();
    }

    /* ---------------------------------------------------------------------- */
    //! Returns the cell lengths
    const Vec_t get_cell_lengths() { return this->cell_lengths; }

    //! Returns the cell angles
    const Vec_t get_cell_angles() { return this->cell_angles; }

    /**
     * Returns the transformation matrix to transform absolute cartesian
     * coordinates into fractional/scaled coordinates.
     */
    const Cell_t get_cartesian2scaled_matrix() {
      return this->cartesian2scaled;
    }

    /**
     * Returns the transformation matrix to transform fractional/scaled
     * coordinates into absolute cartesian coordinates.
     */
    const Cell_t get_scaled2cartesian_matrix() {
      return this->scaled2cartesian;
    }

    //! Returns the reciprocal space lattice vectors
    const Cell_t get_reciprocal_vectors() { return this->reciprocal_vectors; }

    //! Returns the reciprocal space cell lengths
    const Vec_t get_reciprocal_lengths() { return this->reciprocal_lengths; }

    /* ---------------------------------------------------------------------- */
    /**
     * Calculates the reciprocal space lattice vectors from the cartesian real
     * space lattice vectors and the reciprocal space cell lengths
     */
    void set_reciprocal_vectors() {
      Vec_t c_abg = cell_angles.array().cos();

      //! Cell volume
      double V{this->cell_lengths[0] * this->cell_lengths[1] *
               this->cell_lengths[2] *
               std::sqrt(1 - c_abg[0] * c_abg[0] - c_abg[1] * c_abg[1] -
                         c_abg[2] * c_abg[2] +
                         2 * c_abg[0] * c_abg[1] * c_abg[2])};

      double Vinv{1. / V};
      this->reciprocal_vectors *= Vinv;

      Vec_t recip1, recip2, recip3;

      this->crossproduct(this->cell_vectors.col(1), this->cell_vectors.col(2),
                         recip1);

      this->crossproduct(this->cell_vectors.col(2), this->cell_vectors.col(0),
                         recip2);

      this->crossproduct(this->cell_vectors.col(0), this->cell_vectors.col(1),
                         recip3);

      for (int ii{0}; ii < 3; ++ii) {
        this->reciprocal_vectors(ii, 0) *= recip1[ii];
        this->reciprocal_vectors(ii, 1) *= recip2[ii];
        this->reciprocal_vectors(ii, 2) *= recip3[ii];
      }

      this->reciprocal_lengths = this->reciprocal_vectors.colwise().norm();
    }
    /* ---------------------------------------------------------------------- */
    template <typename DerivedA, typename DerivedB>
    void crossproduct(const Eigen::MatrixBase<DerivedA> & v1,
                      const Eigen::MatrixBase<DerivedB> & v2, Vec_t & v3) {
      v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
      v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
      v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Calculates the transformation matrix to transform absolute cartesian
     * coordinates to fractional/scaled coordinates
     */
    void set_transformation_matrix() {
      Vec_t c_abg = cell_angles.array().cos();
      Vec_t s_abg = cell_angles.array().sin();

      //! Cell volume divided by a*b*c
      double V{std::sqrt(1 - c_abg[0] * c_abg[0] - c_abg[1] * c_abg[1] -
                         c_abg[2] * c_abg[2] +
                         2 * c_abg[0] * c_abg[1] * c_abg[2])};
      double Vinv{1 / V};

      //! compute transformation matrix from the cartesian system
      //! to the lattice coordinate system
      this->cartesian2scaled(0, 0) = 1.0 / this->cell_lengths[0];
      this->cartesian2scaled(1, 0) =
          -c_abg[2] / (this->cell_lengths[0] * s_abg[2]);
      this->cartesian2scaled(2, 0) = Vinv * (c_abg[0] * c_abg[2] - c_abg[1]) /
                                     (s_abg[2] * this->cell_lengths[0]);
      this->cartesian2scaled(1, 1) = 1.0 / (this->cell_lengths[1] * s_abg[2]);
      this->cartesian2scaled(2, 1) = Vinv * (c_abg[1] * c_abg[2] - c_abg[0]) /
                                     (s_abg[2] * this->cell_lengths[1]);
      this->cartesian2scaled(2, 2) = s_abg[2] * Vinv / this->cell_lengths[2];

      /**
       * compute transformation matrix from the lattice coordinate system to
       * cartesian
       */
      this->scaled2cartesian(0, 0) = this->cell_lengths[0];
      this->scaled2cartesian(1, 0) = this->cell_lengths[1] * c_abg[2];
      this->scaled2cartesian(2, 0) = this->cell_lengths[2] * c_abg[1];
      this->scaled2cartesian(1, 1) = this->cell_lengths[1] * s_abg[2];
      this->scaled2cartesian(2, 1) =
          this->cell_lengths[2] * (c_abg[0] - c_abg[1] * c_abg[2]) / s_abg[2];
      this->scaled2cartesian(2, 2) = this->cell_lengths[2] * V / s_abg[2];
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Calculates the fractional/scaled coordinates (position_sc) from the
     * absolute cartesian coordinates (position)
     */
    template <typename DerivedA, typename DerivedB>
    void get_cartesian2scaled(const Eigen::MatrixBase<DerivedA> & position,
                              Eigen::MatrixBase<DerivedB> & position_sc) {
      position_sc = this->cartesian2scaled.transpose() * position;
    }

    /** Calculates the absolute cartesian coordinates (position)
     * from the fractional/scaled coordinates (position_sc)
     */
    template <typename DerivedA, typename DerivedB>
    void get_scaled2cartesian(const Eigen::MatrixBase<DerivedA> & position_sc,
                              Eigen::MatrixBase<DerivedB> & position) {
      position = this->scaled2cartesian.transpose() * position_sc;
    }

   protected:
    //! lattice vectors
    Cell_t cell_vectors = Cell_t::Ones();
    //! Reciprocal lattice
    Cell_t reciprocal_vectors = Cell_t::Ones();
    //! Cell lengths
    Vec_t cell_lengths = Vec_t::Ones();
    //! Reciprocal Cell lengths
    Vec_t reciprocal_lengths = Vec_t::Ones();
    //! alpha(b,c) beta(a,c) gamma(a,b) in radian
    Vec_t cell_angles = Vec_t::Ones();
    //! transformation matrix from the lattice coordinate system to cartesian
    Cell_t scaled2cartesian = Cell_t::Zero();
    //! transformation matrix from the cartesian system to the lattice
    //! coordinate system
    Cell_t cartesian2scaled = Cell_t::Zero();
    double pi{EIGEN_PI};
  };
}  // namespace rascal

#endif  // SRC_LATTICE_HH_
