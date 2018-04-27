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
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
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
#include <basic_types.h>
#include <cmath>

namespace proteus {

//using Cell_t = typename proteus::Cell_t;
//using Vec3_t = typename proteus::Vec3_t;
using Vector_ref = Eigen::Map<Vec3_t>;

class Lattice {
  public:
    //! Default constructor
    Lattice() = default;

    Lattice(const Cell_t& cell) {
      this->set_cell(cell);
      //this->cell_angles[0] = (cell[1].dot(cell[2])/this->cell_lenght[1]/this->cell_lenght[2]).array().acos();
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
      this->compute_transformation_matrix();
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

    inline void compute_transformation_matrix(){
      Vec3_t c_abg = cell_angles.array().cos();
      Vec3_t s_abg = cell_angles.array().sin();
      double V{std::sqrt(1 - c_abg[0]*c_abg[0] - c_abg[1]*c_abg[1] - c_abg[2]*c_abg[2] + 2 * c_abg[0] * c_abg[1] * c_abg[2] )}; // Cell volume divided by a*b*c
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
  /*
    inline Vector_ref get_cartesian2scaled(const Vector_ref& position){
      return this->cartesian2scaled.transpose().dot(position);
    }

    inline Vector_ref get_cartesian2scaled(const Vec3_t& position){
      return this->cartesian2scaled.transpose().dot(position);
    }
*/
    inline Eigen::MatrixXd get_cartesian2scaled(const Eigen::Ref<const Eigen::MatrixXd> positions,
                                                Eigen::Ref<Eigen::MatrixXd> positions_sc){
      positions_sc = this->cartesian2scaled.transpose() * positions;
    }

    inline void get_scaled2cartesian(const Eigen::Ref<const Eigen::MatrixXd> positions_sc,Eigen::Ref<Eigen::MatrixXd> positions){
      positions = this->scaled2cartesian.transpose() * positions_sc;
    }

  protected:
    Cell_t cell_vectors;
    Vec3_t cell_lenghts;
    Vec3_t cell_angles;// alpha(b,c) beta(a,c) gamma(a,b) in radian
    Cell_t scaled2cartesian = Cell_t::Zero(); //! transformation matrix from the lattice coordinate system to cartesian
    Cell_t cartesian2scaled = Cell_t::Zero(); //! transformation matrix from the cartesian system to the lattice coordinate system
    double pi{M_PI};


    constexpr static int Dim{3};
};


} // proteus

#endif /* LATTICE_H */
