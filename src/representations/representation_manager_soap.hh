/**
 * file   representation_manager_soap.hh
 *
 * @author Max Veit <max.veit@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @author Michael Willatt <michael.willatt@epfl.ch>
 * @author Andrea Grisafi <andrea.grisafi@epfl.ch>
 *
 * @date   12 March 2019
 *
 * @brief  Compute the spherical harmonics expansion of the local atom density
 *
 * Copyright © 2019 Max Veit, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_

#include "representations/representation_manager_base.hh"
#include "representations/representation_manager_spherical_expansion.hh"
#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "structure_managers/property_block_sparse.hh"
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
  }  // namespace internal

  template <class StructureManager>
  class RepresentationManagerSoap
      : public RepresentationManagerBase {
   public:
    using Manager_t = StructureManager;
    using hypers_t = RepresentationManagerBase::hypers_t;
    using key_t = std::vector<int>;
    using SparseProperty_t = BlockSparseProperty<double, 1, 0>;
    using data_t = typename SparseProperty_t::data_t;
    using input_data_t = typename SparseProperty_t::input_data_t;

    RepresentationManagerSoap(Manager_t & sm, const hypers_t & hyper)
        : structure_manager{sm}, soap_vectors{sm}, rep_expansion{sm, hyper} {
      this->set_hyperparameters(hyper);
    }

    void set_hyperparameters(const hypers_t & hypers) {
      this->max_radial = hypers.at("max_radial");
      this->max_angular = hypers.at("max_angular");
    }

    std::vector<precision_t> & get_representation_raw_data() {
      std::vector<precision_t> aa{};
      return aa;
    }

    data_t& get_representation_sparse_raw_data() {
      return this->soap_vectors.get_raw_data();
    }

    void compute() {
      // SparseProperty_t expansion_coefficients;
      rep_expansion.compute();
      auto expansion_coefficients{rep_expansion.soap_vectors};

      for (auto center : this->structure_manager) {
	key_t center_type{center.get_atom_type()};
        input_data_t coefficients{expansion_coefficients[center]};
        input_data2_t soap_vector;
	key2_t element_pair{};
	size_t n_row{pow(this->max_radial, 2)};
	size_t n_col{this->max_angular + 1};

	for (key1_t species1: coefficients) {
	  for (key1_t species2: coefficients) {
            soap_vector.emplace(
	      std::make_pair(element_pair, dense_t::Zero(n_row, n_col)));
	    size_t nn{0};
	    for (size_t n1 = 0; n1 < rep_expansion.max_radial; n1++) {
	      for (size_t n2 = 0; n2 < rep_expansion.max_radial; n2++) {
		size_t lm{0};
	        for (size_t l = 0; l < rep_expansion.max_angular+1; l++) {
	          for (size_t m = 0; m < 2*rep_expansion.max_angular+1; m++) {
	            soap_vector[element_pair][nn,l] += coefficients[species1][n1,lm] *
						       coefficients[species2][n2,lm];
		    lm++;
	          }
	        }
		nn++;
	      }
	    }
	  }
	}
      this->soap_vectors.push_back(soap_vector);
      }
    }

   protected:
   private:
    size_t max_radial{};
    size_t max_angular{};
    Manager_t & structure_manager;
    SparseProperty_t soap_vectors;
    RepresentationManagerSphericalExpansion<Manager_t> rep_expansion;



  }


}  // namespace rascal

#endif // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_SOAP_HH_
