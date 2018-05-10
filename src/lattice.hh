/**
 * file   lattice.cc
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Implementation of the neighbourhood manager for lammps
 *        neighbourhood lists
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

#ifndef LATTICE_H
#define LATTICE_H

#include <Eigen/Dense>
#include <basic_types.hh>
#include <cmath>

namespace rascal {

//using Cell_t = typename rascal::Cell_t;
//using Vec3_t = typename rascal::Vec3_t;
using Vector_ref = Eigen::Map<Vec3_t>;

class Lattice {
  public:
    //! Default constructor
    Lattice() = default;

    Lattice(const Cell_t& cell) {
      this->set_cell(cell);
    };

    //! Copy constructor
    Lattice(const Lattice &other) = delete;

    //! Move constructor
    Lattice(Lattice &&other) = default;

    //! Destructor
    virtual ~Lattice() = default;

    //! Copy assignment operator
    Lattice& operator=(const Lattice &other) = delete;

    //! Move assignment operator
    Lattice& operator=(Lattice &&other) = default;

    void set_cell(const Cell_t& cell){
      this->cell_vectors = cell;
      this->cell_lenghts = cell.colwise().norm();
      this->cell_angles[0] = std::acos(this->cell_vectors.col(1).dot(this->cell_vectors.col(2))/this->cell_lenghts[1]/this->cell_lenghts[2]);
      this->cell_angles[1] = std::acos(this->cell_vectors.col(0).dot(this->cell_vectors.col(2))/this->cell_lenghts[0]/this->cell_lenghts[2]);
      this->cell_angles[2] = std::acos(this->cell_vectors.col(1).dot(this->cell_vectors.col(0))/this->cell_lenghts[1]/this->cell_lenghts[0]);
      this->set_transformation_matrix();
      this->set_reciprocal_vectors();
    }

    const Vec3_t get_cell_lengths() {
      return this->cell_lenghts;
    }

    const Vec3_t get_cell_angles() {
      return this->cell_angles;
    }
    
    const Cell_t get_cartesian2scaled_matrix() {
      return this->cartesian2scaled;
    }

    const Cell_t get_scaled2cartesian_matrix() {
      return this->scaled2cartesian;
    }

    const Cell_t get_reciprocal_vectors() {
      return this->reciprocal_vectors;
    }
    
    const Vec3_t get_reciprocal_lenghts() {
      return this->reciprocal_lenghts;
    }

    inline void set_reciprocal_vectors(){
      Vec3_t c_abg = cell_angles.array().cos();
      double V{ this->cell_lenghts[0] *this->cell_lenghts[1] *this->cell_lenghts[2] *  std::sqrt(1 - c_abg[0]*c_abg[0] - c_abg[1]*c_abg[1] - c_abg[2]*c_abg[2] + 2 * c_abg[0] * c_abg[1] * c_abg[2] )}; //! Cell volume
      double Vinv{1./V};
      this->reciprocal_vectors *= Vinv;
      Vec3_t recip1,recip2,recip3;
      this->crossproduct(this->cell_vectors.col(1),this->cell_vectors.col(2),recip1);
      this->crossproduct(this->cell_vectors.col(2),this->cell_vectors.col(0),recip2);
      this->crossproduct(this->cell_vectors.col(0),this->cell_vectors.col(1),recip3);
      for (int ii{0};ii<3;++ii) {
        this->reciprocal_vectors(ii,0) *= recip1[ii];
        this->reciprocal_vectors(ii,1) *= recip2[ii];
        this->reciprocal_vectors(ii,2) *= recip3[ii];
      }
      this->reciprocal_lenghts = this->reciprocal_vectors.colwise().norm();
    }

    template <typename DerivedA,typename DerivedB>
    inline void crossproduct(const Eigen::MatrixBase<DerivedA>& v1, const Eigen::MatrixBase<DerivedB>& v2, Vec3_t& v3){
      v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
      v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
      v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
    }

    inline void set_transformation_matrix(){
      Vec3_t c_abg = cell_angles.array().cos();
      Vec3_t s_abg = cell_angles.array().sin();
      double V{std::sqrt(1 - c_abg[0]*c_abg[0] - c_abg[1]*c_abg[1] - c_abg[2]*c_abg[2] + 2 * c_abg[0] * c_abg[1] * c_abg[2] )}; //! Cell volume divided by a*b*c
      double Vinv{1/V};
      //! compute transformation matrix from the cartesian system to the lattice coordinate system
      this->cartesian2scaled(0,0) = 1.0/this->cell_lenghts[0];
      this->cartesian2scaled(1,0) = -c_abg[2]/(this->cell_lenghts[0]*s_abg[2]);
      this->cartesian2scaled(2,0) = Vinv* ((c_abg[0]*c_abg[2]-c_abg[1]))/(s_abg[2]*this->cell_lenghts[0]);
      this->cartesian2scaled(1,1) = 1.0/(this->cell_lenghts[1]*s_abg[2]);
      this->cartesian2scaled(2,1) = Vinv*( c_abg[1]*c_abg[2]-c_abg[0] )/(s_abg[2]*this->cell_lenghts[1]);
      this->cartesian2scaled(2,2) = s_abg[2]*Vinv/this->cell_lenghts[2];
      //! compute transformation matrix from the lattice coordinate system to cartesian
      this->scaled2cartesian(0,0) = this->cell_lenghts[0];
      this->scaled2cartesian(1,0) = this->cell_lenghts[1] * c_abg[2];
      this->scaled2cartesian(2,0) = this->cell_lenghts[2] * c_abg[1];
      this->scaled2cartesian(1,1) = this->cell_lenghts[1] * s_abg[2];
      this->scaled2cartesian(2,1) = this->cell_lenghts[2] * (c_abg[0] - c_abg[1] * c_abg[2]) / s_abg[2];
      this->scaled2cartesian(2,2) = this->cell_lenghts[2] * V / s_abg[2];
    }

    template <typename DerivedA,typename DerivedB>
    inline void get_cartesian2scaled(const Eigen::MatrixBase<DerivedA>& position, Eigen::MatrixBase<DerivedB>& position_sc){
      position_sc = this->cartesian2scaled.transpose() * position;
    }

    template <typename DerivedA,typename DerivedB>
    inline void get_scaled2cartesian(const Eigen::MatrixBase<DerivedA>& position_sc, Eigen::MatrixBase<DerivedB>& position){
      position = this->scaled2cartesian.transpose() * position_sc;
    }

  protected:
    Cell_t cell_vectors = Cell_t::Ones();
    Cell_t reciprocal_vectors = Cell_t::Ones();
    Vec3_t cell_lenghts = Vec3_t::Ones();
    Vec3_t reciprocal_lenghts = Vec3_t::Ones();
    Vec3_t cell_angles = Vec3_t::Ones();// alpha(b,c) beta(a,c) gamma(a,b) in radian
    Cell_t scaled2cartesian = Cell_t::Zero(); //! transformation matrix from the lattice coordinate system to cartesian
    Cell_t cartesian2scaled = Cell_t::Zero(); //! transformation matrix from the cartesian system to the lattice coordinate system
    double pi{M_PI};


    constexpr static int Dim{3};
};


} // rascal

#endif /* LATTICE_H */
