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
#include "math/math_utils.hh"
#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace rascal {

  namespace internal {

  }

  /**
   * Handles the expansion of an environment in a spherical and radial basis.
   *
   * The local environment of each atom is represented by Gaussians of a
   * certain width (user-defined; can be constant, species-dependent, or
   * radially dependent).  This density field is expanded in an angular basis
   * of spherical harmonics (à la SOAP) and a radial basis of either Gaussians
   * (again, as in SOAP) or one of the more recent bases currently under
   * development.
   */
  template<class StructureManager>
  class RepresentationManagerSphericalExpansion:
    public RepresentationManagerBase
  {
  public:
    // TODO make a traits mechanism
    using hypers_t = RepresentationManagerBase::hypers_t;
    using Property_t = Property<double, 1, 1, Eigen::Dynamic, Eigen::Dynamic>;
    using Manager_t = StructureManager;

    /** Describes the radial basis options.  Analytical integration is
     *  available for the Gaussian-type basis; the polynomial basis is
     *  numerically evaluated but better conditioned.
     */
    enum RadialBasisType {
      GaussianType = 0;
      PolynomialBasis = 1;
    };

    /** Describes how the Gaussian sigma (width) is specified.  Options
     *  are constant, species-dependent, or radially-dependent.
     */
    enum GaussianSigmaType {
      Constant = 0;
      PerSpecies = 1;
      Radial = 2;
    };

    //! Default constructor.  Assumes constant Gaussian sigma (width).
    RepresentationManagerSphericalExpansion(Manager_t &sm, 
        double interaction_cutoff, double cutoff_smooth_width,
        size_t n_species, size_t max_radial, size_t max_angular,
        double gaussian_sigma)
        :structure_manager{sm}, interaction_cutoff{interaction_cutoff},
        cutoff_smooth_width{cutoff_smooth_width},
        max_radial{max_radial}, max_angular{max_angular},
        n_species{n_species}, soap_vectors{sm},
        constant_gaussian_sigma{gaussian_sigma}, gaussian_sigma_type{Constant},
        is_precomputed{false}
    {
      this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
      this->radial_sigmas.resize(this->max_radial, 1);
      this->radial_norm_factors.resize(this->max_radial, 1);
      this->radial_nl_factors.resize(this->max_radial, this->max_angular + 1);
      // TODO(max-veit) the type of Gaussian sigma should not be determined
      // using constructor overloading.  Better convert
      // RepManSphExpn::get_gaussian_sigma() to its own class & let that switch
      // on the Gaussian sigma type
    }

    /**
     * Set the hyperparameters of this descriptor from a json object.
     *
     * @param hypers structure (usually parsed from json) containing the
     *               options and hyperparameters
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the structure
     */
    void set_hyperparameters(const hypers_t& hypers) {
      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
      if (hypers.find("n_species") != hypers.end()){
        this->n_species = hypers.at("n_species");
      }
      this->interaction_cutoff = hypers.at("interaction_cutoff");
      this->cutoff_smooth_width = hypers.at("cutoff_smooth_width");
      if (hypers.at("gaussian_sigma_type").compare("Constant") == 0) {
        this->gaussian_sigma_type = GaussianSigmaType.Constant;
        this->constant_gaussian_sigma = hypers.at("gaussian_sigma");
      } else {
        throw std::logic_error("Requested Gaussian sigma type \'"
                               + hypers.at("gaussian_sigma_type")
                               + "\' has not been implemented.");
      }
      this->radial_ortho_matrix.resize(this->max_radial, this->max_radial);
      this->radial_sigmas.resize(this->max_radial, 1);
      this->radial_norm_factors.resize(this->max_radial, 1);
      this->radial_nl_factors.resize(this->max_radial, this->max_angular + 1);
    }

    /**
     * Set the hyperparameters of this descriptor from a json string
     *
     * @param hypers json string containing the options and hyperparameters
     *               The string will be parsed into a json container.
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the string
     *
     * @throw json.exception.parse_error if the string is not valid json
     */
    void set_hyperparameters(const std::string& hypers) {
      hypers_t hypers_structure = json::parse(hypers);
      this->set_hyperparameters(hypers_structure);
    }

    /**
     * Construct a new RepresentationManager using a hyperparameters container
     *
     * @param hypers container (usually parsed from json) for the options and
     *               hyperparameters
     *
     * @throw logic_error if an invalid option or combination of options is
     *                    specified in the container
     */
    RepresentationManagerSphericalExpansion(Manager_t &sm,
                                            const hypers_t& hyper)
        :structure_manager{sm}, interaction_cutoff{}, cutoff_smooth_width{},
        max_radial{}, max_angular{}, n_species{}, soap_vectors{sm},
        is_precomputed{false}
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

    //! Precompute radial Gaussian widths (NB specific to basis!)
    void precompute_radial_sigmas();

    //! Precompute radial orthogonalization matrix (NB specific to basis!)
    void precompute_radial_overlap();

    //! getter for the representation
    template <size_t Order, size_t Layer> 
    Eigen::Map<Eigen::MatrixXd> get_soap_vector(
        const ClusterRefKey<Order, Layer>& center) {
      return this->soap_vectors[center];
    }

    std::vector<precision_t&> get_representation_raw_data() {
      return this->soap_vector.get_raw_data();
    }

    /**
     * Return whether the radial sigmas and overlap matrix have been precomputed
     *
     * This only needs to be done once on initialization, but it is not done in
     * the constructor to avoid spending time and possibly throwing exceptions
     * there.  Instead, by default, it is done once on the first call of the
     * compute() method and the precomputed results are saved for subsequent
     * calls.
     */
    bool get_is_precomputed() {
      return this->is_precomputed;
    }

  protected:
  private:
    Manager_t& structure_manager;
    //hypers_t hyperparmeters;
    double interaction_cutoff;
    double cutoff_smooth_width;
    double constant_gaussian_sigma;
    GaussianSigmaType gaussian_sigma_type;
    // TODO(max-veit) these are specific to the radial Gaussian basis
    Eigen::MatrixXd<Eigen::Dynamic, 1> radial_sigmas;
    Eigen::MatrixXd<Eigen::Dynamic, 1> radial_norm_factors;
    Eigen::MatrixXd<Eigen::Dynamic, 1> radial_sigma_factors;
    Eigen::MatrixXd<Eigen::Dynamic, Eigen::Dynamic> radial_nl_factors;
    Eigen::MatrixXd<Eigen::Dynamic, Eigen::Dynamic> radial_ortho_matrix;
    size_t n_species;
    size_t max_radial;
    size_t max_angular;
    bool is_precomputed;

    Property_t soap_vectors;
  };


  /** Compute common prefactors for the radial Gaussian basis functions */
  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::
      precompute_radial_sigmas(){

    using math::PI;
    using math::SQRT_TWO;
    using std::pow;

    size_t radial_n;
    size_t angular_l;

    for (radial_n = 0; radial_n < this->max_radial; radial_n += 1) {
      this->radial_sigmas[radial_n] = std::max(
            std::sqrt((double)radial_n), 1.0)
          * this->interaction_cutoff / (double)this->max_radial;
    }

    // Precompute common prefactors
    for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
      this->radial_norm_factors(radial_n) = std::sqrt(
          2.0 / std::tgamma(1.5 + radial_n)
          * pow(this->radial_sigmas[radial_n], 3.0 + 2.0*radial_n));
      for (angular_l = 0; angular_l < this->max_angular + 1; angular_l++) {
        this->radial_nl_fators(radial_n, angular_l) =
            std::exp2(-0.5 * (1.0 + angular_l - radial_n))
            * std::tgamma(0.5 * (3.0 + angular_l + radial_n))
            / std::tgamma(1.5 + angular_l);
    }
  }

  /** Compute the radial overlap matrix for later orthogonalization.
   *
   *  @throw runtime_error if the overlap matrix cannot be diagonalized
   */
  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::
      precompute_radial_overlap(){

    using std::pow;
    using std::sqrt;
    using gamma = std::tgamma;

    size_t radial_n1;
    size_t radial_n2;

    //TODO(max-veit) see if we can replace the gammas with their natural logs,
    //since it'll overflow for relatively small n (n1 + n2 >~ 300)
    Eigen::MatrixXd overlap(this->max_radial, this->max_radial);
    for (radial_n1 = 0; radial_n1 < this->max_radial; radial_n1++) {
      for (radial_n2 = 0; radial_n2 < this->max_radial; radial_n2++) {
        overlap(radial_n1, radial_n2) =
            pow(0.5/pow(this->radial_sigmas[radial_n1], 2)
                   + 0.5/pow(this->radial_sigmas[radial_n2], 2),
                -0.5*(3.0 + radial_n1 + radial_n2))
            / (pow(this->radial_sigmas[radial_n1], radial_n1) +
               pow(this->radial_sigmas[radial_n2], radial_n2))
            * gamma(0.5 * (3.0 + radial_n1 + radial_n2))
            / pow(this->radial_sigmas[radial_n1] *
                       this->radial_sigmas[radial_n2], 1.5)
            * sqrt(gamma(1.5 + radial_n1)
                 * gamma(1.5 + radial_n2));
      }
    }

    // Compute the inverse square root of the overlap matrix
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(overlap);
    if (eigensolver.info() != Eigen::Success) {
      throw std::runtime_error("Warning: Could not diagonalize "
                               "radial overlap matrix.");
    }
    //TODO(max-veit) check whether this is the right way to do this in Eigen
    Eigen::MatrixXd eigenvalues = eigensolver.eigenvalues();
    Eigen::ArrayXd eigs_invsqrt(this->max_radial) =
        eigenvalues.array().sqrt().inverse();
    Eigen::MatrixXd unitary = eigensolver.eigenvectors();
    this->radial_ortho_matrix = unitary.adjoint() *
                                eigs_invsqrt.matrix().asDiagonal() * unitary;
    this->is_precomputed = true;
  }

  /**
   * Precompute everything that doesn't depend on the atomic structure
   * (only on the hyperparameters)
   *
   * Calls the methods precompute_radial_sigmas() and
   * precompute_radial_overlap(); if either of those methods throws exceptions,
   * they will be passed on from here.
   *
   * @throw runtime_error if the overlap matrix cannot be diagonalized
   */
  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::precompute() {
    this->precompute_radial_sigmas();
    this->precompute_radial_overlap();
    // Only if none of the above failed (threw exceptions)
    this->is_precomputed = true;
  }


  template<class Mngr>
  void RepresentationManagerSphericalExpansion<Mngr>::compute(){

    using math::PI;
    using std::pow;

    this->soap_vectors.resize_to_zero();
    this->soap_vectors.set_nb_row(this->n_species * this->max_radial);
    this->soap_vectors.set_nb_col(pow(this->max_angular + 1, 2));

    size_t radial_n;
    size_t angular_l;
    size_t m_count;
    size_t lm_collective_idx;
    double sigma2;
    double radial_sigma_factor;

    if (not this->is_precomputed) {
      this->precompute();
    }

    for (auto center : this->structure_manager) {

      Eigen::MatrixXd soap_vector = Eigen::MatrixXd::Zero(
          this->n_species * this->max_radial,
          pow(this->max_angular + 1, 2));
      Eigen::MatrixXd radial_integral(this->max_radial, this->max_angular + 1);

      // Start the accumulator with the central atom
      // All terms where l =/= 0 cancel
      sigma2 = pow(this->get_gaussian_sigma(center), 2);
      // TODO(max-veit) this is specific to the Gaussian radial basis
      // (along with the matching computation below)
      // And ditto on the gamma functions (potential overflow)
      for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
        // TODO(max-veit) remember to update when we do multiple n_species!
        radial_integral(radial_n, 0) = radial_norm_factors(radial_n)
            * radial_nl_fators(radial_n, 0)
            * pow(1.0/sigma2 + 1.0/this->radial_sigmas[radial_n]**2,
                       -0.5 * (3.0 + radial_n));
      }
      soap_vector.col(0) = this->radial_ortho_matrix *
                           radial_integral.col(0) / sqrt(4.0 * PI);

      for (auto neigh : center) {
        auto dist{this->structure_manager.get_distance(neigh)};
        auto direction{this->structure_manager.get_direction_vector(neigh)};
        double exp_factor = std::exp(-0.5 * pow(dist, 2) / sigma2);
        sigma2 = pow(this->get_gaussian_sigma(neigh), 2);

        // Note: the copy _should_ be optimized out (RVO)
        MatrixXd harmonics = math::compute_spherical_harmonics(
            direction, this->max_angular);

        // Precompute radial factors that also depend on the Gaussian sigma
        VectorXd radial_sigma_factors(this->max_radial);
        for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
          radial_sigma_factors(radial_n) = (pow(sigma2, 2) +
                                 sigma2*pow(this->radial_sigmas(radial_n), 2))
                                / pow(this->radial_sigmas(radial_n), 2);
        }

        for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
          // TODO(max-veit) this is all specific to the Gaussian radial basis
          // (doing just the angular integration would instead spit out
          // spherical Bessel functions below)
          //TODO(max-veit) how does this work with SpeciesFilter?
          for(angular_l = 0; angular_l < this->max_angular + 1; angular_l++) {
            radial_integral(radial_n, angular_l) =
                exp_factor * radial_nl_fators(radial_n, angular_l)
                * pow(1.0/sigma2 + 1.0/pow(this->radial_sigmas(radial_n), 2),
                      -0.5 * (3.0 + angular_l + radial_n))
                * pow(dist / sigma2, angular_l)
                * math::hyp1f1(0.5 * (3.0 + angular_l * radial_n),
                               1.5 + angular_l,
                               0.5 * pow(dist, 2)
                                   / radial_sigma_factors(radial_n));
          }
        }
        radial_integral = this->radial_ortho_matrix * radial_integral;

        for (radial_n = 0; radial_n < this->max_radial; radial_n++) {
          lm_collective_idx = 0;
          for(angular_l = 0; angular_l < this->max_angular + 1; angular_l++) {
            for (size_t m_array_idx = 0; m_array_idx < 2*this->angular_l + 1;
                m_array_idx++) {
              soap_vector(n, lm_collective_idx) +=
                  radial_integral(radial_n, angular_l)
                  * harmonics(angular_l, m_array_idx);
              ++lm_collective_idx;
            }
          }
        }
      } // for (neigh : center)
      this->soap_vectors.push_back(soap_vector)
    } // for (center : structure_manager)
  } // compute()

  /**
   * Return the width of the atomic Gaussian for each neighbour
   *
   * This is `sigma' in the definition `f(r) = A exp(r / (2 sigma^2))'.
   * The width may depend both on the atomic species of the neighbour as well
   * as the distance.
   *
   * @param pair Atom pair defining the neighbour, as e.g. returned by
   *             iteration over neighbours of a centre
   *
   * @throw logic_error if the requested sigma type has not been implemented
   *
   */
  template<class Mngr, size_t Order, size_t Layer>
  void RepresentationManagerSphericalExpansion<Mngr>::get_gaussian_sigma(
        ClusterRefKey<Order, Layer> & pair) {
    switch (this->gaussian_sigma_type) {
    case (Constant) : {
      return this->constant_gaussian_sigma;
      break;
    }
    case (PerSpecies) :
    case (Radial) : {
      throw std::logic_error("Requested a sigma type that has not yet "
                             "been implemented");
      break;
    }
    }
  }

}

#endif /* BASIS_REPRESENTATION_MANAGER_SPHERICAL_EXPANSION_H */

