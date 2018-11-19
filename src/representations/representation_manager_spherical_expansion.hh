/**
 * file   representation_manager_spherical_expansion.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 *
 * @date   19 October 2018
 *
 * @brief  Compute the spherical harmonics expansion of the local atom density
 *
 * Copyright © 2018 Max Veit, COSMO (EPFL), LAMMM (EPFL)
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

// Checklist (to get this working):
// * Implement associated Legendre polynomials
// * Orthogonalize the radial coefficients
// * Implement the accumulator
// * Write some basic tests

#ifndef BASIS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_H
#define BASIS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_H


#include "representations/representation_manager_base.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"
#include <algorithm>
#include <cmath>
#include <boost/math/special_functions/factorial.hpp>
#include <Eigen/Eigenvalues>
#include <vector>
#include <algorithm>

namespace rascal {

  namespace internal {

    const double SQRT_TWO = std::sqrt(2.0);
    const double PI = 3.1415926535897932384626433832795028841971; // OK I think that's enough

  }

  template<class StructureManager>
  class RepresentationManagerSphericalExpansion:
    public RepresentationManagerBase
  {
  public:
    // TODO make a traits mechanism
    using hypers_t = RepresentationManagerBase::hypers_t;
    using Property_t = Property<double, 1, 1, Eigen::Dynamic, Eigen::Dynamic>;
    //using Property_1_t = Property<double, 1, 1, Eigen::Dynamic, Eigen::Dynamic>;
    using Manager_t = StructureManager;

    //! Default constructor
    RepresentationManagerSphericalExpansion(Manager_t &sm, 
      double interaction_cutoff, double cutoff_smooth_width,
      size_t n_species, size_t max_radial, size_t max_angular)
      :structure_manager{sm}, interaction_cutoff{interaction_cutoff},
      cutoff_smooth_width{cutoff_smooth_width},
      max_radial{max_radial}, max_angular{max_angular},
      n_species{n_species}, soap_vectors{sm}
      {
        size_t radial_n;
        size_t radial_n1;
        size_t radial_n2;

        this->radial_sigmas.resize(this->max_radial, 1);
        for (radial_n = 0; radial_n < max_radial; radial_n += 1) {
          this->radial_sigmas[radial_n] = std::max(
                std::sqrt((float)radial_n), 1.0)
              * interaction_cutoff / (float)max_angular;
        }

        // Precompute common prefactors
        prefactors.resize(this->max_radial, 1);
        for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
            this->prefactors(radial_n) = std::sqrt(
              2.0 / std::tgamma(1.5 + radial_n)
              * std::pow(this->radial_sigmas[radial_n], 3.0 * 2.0*radial_n));
        }

        //radial_ortho_matrix.resize(this->max_radial, this->max_radial);
        Eigen::MatrixXd overlap(this->max_radial, this->max_radial);
        for (radial_n1 = 0; radial_n1 < this->max_radial; radial_n1++) {
          for (radial_n2 = 0; radial_n2 < this->max_radial; radial_n2++) {
            overlap(radial_n1, radial_n2) =
                std::pow(0.5/std::pow(this->radial_sigmas[radial_n1], 2)
                       + 0.5/std::pow(this->radial_sigmas[radial_n2], 2),
                         -0.5*(3.0 + radial_n1 + radial_n2))
                / (std::pow(this->radial_sigmas[radial_n1], radial_n1) +
                   std::pow(this->radial_sigmas[radial_n2], radial_n2))
                * std::tgamma(0.5 * (3.0 + radial_n1 + radial_n2))
                / std::pow(this->radial_sigmas[radial_n1] *
                           this->radial_sigmas[radial_n2], 1.5)
                * std::sqrt(std::tgamma(1.5 + radial_n1)
                          * std::tgamma(1.5 + radial_n2));
          }
        }
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(overlap);
        if (eigensolver.info() != Eigen::Success) {
          cout << "Warning: Could not diagonalize radial overlap matrix.";
          // TODO How are we doing error handling again?
          // And in a constructor, no less.  Maybe we should move this work
          // somewhere else, like in a "precompute" method?
        }
        // The orthogonalization matrix has eigenvalues 1/sqrt() of the
        // overlap matrix, and is otherwise its inverse
        eigenvalues = eigensolver.eigenvalues();
        Eigen::VectorXd eigs_invsqrt(this->max_radial);
        for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
          eigs_invsqrt(radial_n) = 1.0 / std::sqrt(eigenvalues(radial_n));
        }
        unitary = eigensolver.eigenvectors();
        //TODO check whether this is the right way to do this in Eigen
        radial_ortho_matrix = unitary.adjoint() * eigs_invsqrt.asDiagonal()
                              * unitary;
      }

    RepresentationManagerSphericalExpansion(Manager_t &sm,const hypers_t& hyper)
      :structure_manager{sm},central_decay{},
      interaction_cutoff{},
      interaction_decay{},coulomb_matrices{sm}
      {
        this->set_hyperparameters(hyper);
      }

    //! Copy constructor
    RepresentationManagerSphericalExpansion(
      const RepresentationManagerSphericalExpansion &other) = delete;

    //! Move constructor
    RepresentationManagerSphericalExpansion(
      RepresentationManagerSphericalExpansion &&other) = default;

    //! Destructor
    virtual ~RepresentationManagerSphericalExpansion()  = default;

    //! Copy assignment operator
    RepresentationManagerSphericalExpansion& operator=(
      const RepresentationManagerSphericalExpansion &other) = delete;

    //! Move assignment operator
    RepresentationManagerSphericalExpansion& operator=(
      RepresentationManagerSphericalExpansion && other) = default;

    //! compute representation
    void compute();

    //! getter for the representation
    Eigen::Map<Eigen::MatrixXd> get_representation_full() {
      auto Nb_centers{this->structure_manager.nb_clusters(1)};
      auto Nb_features{this->size};
      auto& raw_data{this->coulomb_matrices.get_raw_data()};
      Eigen::Map<Eigen::MatrixXd> representation(raw_data.data(),Nb_features,Nb_centers);
      std::cout <<"Sizes: "<< representation.size()<<", "<< Nb_features<<", "<<Nb_centers<<std::endl;
      return representation;
    }
    template <size_t Order, size_t Layer> 
    Eigen::Map<Eigen::MatrixXd> get_coulomb_matrix(const ClusterRefKey<Order, Layer>& center){
      return this->coulomb_matrices_full[center];
    }

    // TODO think of a generic input type for the hypers
    void set_hyperparameters(const hypers_t & );

    Manager_t& structure_manager;
    //hypers_t hyperparmeters;
    double interaction_cutoff;
    double cutoff_smooth_width;
    // TODO these are specific to the radial Gaussian basis
    Eigen::MatrixXd<Eigen::Dynamic, 1> radial_sigmas;
    Eigen::MatrixXd<Eigen::Dynamic, 1> prefactors;
    Eigen::MatrixXd<Eigen::Dynamic, Eigen::Dynamic> radial_ortho_matrix;
    size_t n_species;
    size_t max_radial;
    size_t max_angular;

    Property_t soap_vectors;

  protected:
  private:
  };


  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::compute(){

    this->soap_vectors.resize_to_zero();
    this->soap_vectors.set_nb_row(this->n_species * this->max_radial);
    this->soap_vectors.set_nb_col(std::pow(this->max_angular + 1, 2));

    size_t radial_n;
    size_t angular_l;
    size_t m_count;

    for (auto center : this->structure_manager) {

      Eigen::MatrixXd soap_vector = Eigen::MatrixXd::Zero(
          this->n_species * this->max_radial,
          std::pow(this->max_angular + 1, 2));

      // Start the accumulator with the central atom
      // All terms where l =/= 0 cancel
      double sigma2 = std::pow(this->get_gaussian_sigma(), 2);
      for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
        // TODO remember to update when we do multiple n_species!
        soap_vector(radial_n, 0) = prefactors(radial_n)
            * std::exp2(-0.5 * (1.0 - radial_n))
            * std::pow(1.0/sigma2 + 1.0/this->radial_sigmas[radial_n]**2,
                       -0.5 * (3.0 + radial_n))
            * std::tgamma(0.5 * (3.0 + radial_n))
            * prefactors(radial_n) / std::sqrt(4.0 * PI);
      }
      // TODO don't forget to orthogonalize!

      for (auto neigh : center) {
        auto dist{this->structure_manager.get_distance(neigh)};
        auto direction{this->structure_manager.get_direction_vector(neigh)};

        // Precompute exponential factor
        double exp_factor = std::exp(-0.5 * std::pow(
              dist / this->get_gaussian_sigma(), 2));

        // Precompute spherical harmonics TODO maths seems a little low-level
        // for here, move into its own function?
        // The cosine against the z-axis is just the z-component of the
        // direction vector
        double cos_theta = direction[2];
        // Only really useful with the multiple-angle formulae, whose cost
        // scales linearly with m (and includes combinatorial numbers)
        double phi = std::atan2(direction[1], direction[0]);
        //double inv_sqrt_xy = 1.0 / std::sqrt(std::pow(direction[0], 2) +
                                             //std::pow(direction[1], 2));
        //double cos_phi = direction[0] * inv_sqrt_xy;
        //double sin_phi = direction[1] * inv_sqrt_xy;
        MatrixXd harmonics(max_angular+1, 2*max_angular + 1);
        // Associated Legendre polynomials include the Condon-Shortley phase
        // TODO use the recurrence relations to compute P_l^m(cos theta)
        MatrixXd assoc_legendre_polynom(max_angular+1, 2*max_angular + 1);
        for(angular_l = 0; angular_l < this->max_angular + 1; angular_l++) {
          size_t m_array_idx = 0;
          for (m_count = -1.0*angular_l; m_count < 1; m_count++) {
            if (m_count == 0) {
              harmonics(angular_l, m_array_idx) = assoc_legendre_polynom(
                  angular_l, m_count);
            } else {
              harmonics(angular_l, m_array_idx) = assoc_legendre_polynom(
                  angular_l, m_count) * -1.0 * std::sin(m_count * phi)
                  * std::sqrt(2.0 * (2.0*angular_l + 1) / (4.0 * PI));
              harmonics(angular_l, m_array_idx - 2*m_count) = assoc_legendre_polynom(
                  angular_l, m_count) * std::cos(m_count * phi)
                  * boost::math::factorial<int>(angular_l + m_count)
                  / boost::math::factorial<int>(angular_l - m_count)
                  * std::sqrt(2.0 * (2.0*angular_l + 1.0) / (4.0 * PI));
              if ((m_count % 2) == 1) {
                // Factor of (-1)^m that goes with negative-m entries
                harmonics(angular_l, m_array_idx - 2*m_count) *= -1.0;
              }
            }
          }
        }

        //TODO how does this work with SpeciesFilter?
        for (size_t radial_n = 0; radial_n < this->max_radial; radial_n++) {
          for(size_t angular_l = 0; angular_l < this->max_angular + 1; angular_l++) {
            for (size_t idx_m = 0; idx_m < this->angular_l + 1; idx_m++) {
              // Compute SOAP coefficients (use your imagination for now)
            }
          }
        }
      } // for (neigh : center)
      this->soap_vectors.push_back(soap_vector)
    } // for (center : structure_manager)
  }

  //TODO what's the type of the argument?
  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::get_gaussian_sigma(){

  }

}

#endif /* BASIS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_H */

