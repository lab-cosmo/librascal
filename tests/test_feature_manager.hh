/**
 * file   .hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief  test  managers
 *
 * Copyright © 2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TEST_FEATURE_MANAGER_H
#define TEST_FEATURE_MANAGER_H


#include "test_structure.hh"
#include "test_representation_manager.hh"
#include "representations/feature_manager_base.hh"
#include "representations/feature_manager_dense.hh"


namespace rascal {

  template<typename StructureManager>
  struct MultipleStrictStructureManager
  {
    using Manager1_t = StructureManager; 
    using Manager2_t = AdaptorNeighbourList<Manager1_t>;
    using Manager_t = AdaptorStrict<Manager2_t>;

    MultipleStrictStructureManager(
          std::vector<std::string> filenames,
          std::vector<double> cutoffs) {
      for (auto filename : filenames){
        for (auto cutoff : cutoffs){
          Manager1_t m1{};
          m1.update(filename);
          Manager2_t m2{m1,cutoff};
          m2.update();
          Manager_t m3{m2,cutoff};
          m3.update();

          managers1.push_back(std::move(m1));
          managers2.push_back(std::move(m2));
          managers.push_back(std::move(m3));
        }
      }

    }

    ~MultipleStrictStructureManager() {}
    
    std::vector<Manager1_t> managers1{};
    std::vector<Manager2_t> managers2{};
    std::vector<Manager_t> managers{};
        
  };



template<template<typename> typename RepresentationManager>
  struct FeatureFixture: public MultipleStrictStructureManager<StructureManagerCenters>
  { 
    using Manager_t = typename MultipleStrictStructureManager<StructureManagerCenters>::Manager_t;
    using Representation_t = RepresentationManager<Manager_t>;
    FeatureFixture():
      MultipleStrictStructureManager{this->filenames,this->cutoffs}
    {
      for (Manager_t& manager : this->managers){
        Representation_t rp{manager,central_decay,
                interaction_cutoff,interaction_decay,size};
        rp.compute();
        this->representations.push_back(std::move(rp));
      }
      
    }

    ~FeatureFixture() {}
  
    std::vector<Representation_t> representations{};
    
    std::vector<std::string> filenames{
      {"./reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
       "./reference_data/simple_cubic_8.json"}
                                        };
    std::vector<double> cutoffs{{1,2,3,4,5,6}};

    double central_decay{0.5};
    double interaction_cutoff{3};
    double interaction_decay{0.5};
    size_t size{16};

  };

/* ---------------------------------------------------------------------- */
  
} // RASCAL

#endif /* TEST_FEATURE_MANAGER_H */