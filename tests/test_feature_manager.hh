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

#include <list>

namespace rascal {

  template<class StructureManager>
  struct MultipleStrictStructureManager
  {
    using Manager1_t = StructureManager; 
    using Manager2_t = AdaptorNeighbourList<Manager1_t>;
    using Manager_t = AdaptorStrict<Manager2_t>;

    MultipleStrictStructureManager() {
      std::vector<std::string> filenames{
      // "reference_data/CaCrP2O7_mvc-11955_symmetrized_.json",
      "reference_data/simple_cubic_8_.json"};
      std::vector<double> cutoffs{{3}};
      for (auto filename : filenames){
        for (auto cutoff : cutoffs){
          this->managers1.emplace_back();
          this->managers1.back().update(filename);
          this->managers2.emplace_back(managers1.back(),cutoff);
          this->managers2.back().update();
          this->managers.emplace_back(managers2.back(),cutoff);
          this->managers.back().update();
        }
      }
      
    }

    ~MultipleStrictStructureManager() {}

    std::list<Manager1_t> managers1{};
    std::list<Manager2_t> managers2{};
    std::list<Manager_t> managers{};
    
  };



  template<template<typename> typename RepresentationManager>
  struct FeatureFixture: public MultipleStrictStructureManager<StructureManagerCenters>
  { 
    using Parent = MultipleStrictStructureManager<StructureManagerCenters>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = RepresentationManager<Manager_t>;
   
    FeatureFixture() 
    :Parent{}
    { 

      for (Manager_t& manager : this->managers){
        
        Representation_t rep{manager,central_decay,
                 interaction_cutoff,interaction_decay,size};
        rep.compute();
        // this->representations.emplace_back(manager,central_decay,
        //         interaction_cutoff,interaction_decay,size);
        // std::cout << "############2222222222222222"<<std::endl;
        // std::cout << manager.nb_clusters(1)<< ", "<<
        //           manager.nb_clusters(2)  <<", "<<
        //          manager.size() <<std::endl;
        
        // // for (auto center : manager){
        // //   std::cout << "############ "<<std::endl;
        // //   for (auto neigh : center){
        // //     std::cout << manager.get_distance(neigh)<< ", ";
        // //   }
        // //   std::cout <<std::endl;
        // // }


        // this->representations.back().compute();
      }
      
    }

    ~FeatureFixture() {}
  
    std::list<Representation_t> representations{};

    double central_decay{0.5};
    double interaction_cutoff{3};
    double interaction_decay{0.5};
    size_t size{20};

  };

/* ---------------------------------------------------------------------- */
  
} // RASCAL

#endif /* TEST_FEATURE_MANAGER_H */