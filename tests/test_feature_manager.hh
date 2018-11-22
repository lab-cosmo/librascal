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

  struct TestFeatureData
  {
    TestFeatureData() = default;
    ~TestFeatureData() = default;

    std::vector<std::string> filenames{
      "reference_data/CaCrP2O7_mvc-11955_symmetrized_.json",
      "reference_data/simple_cubic_8_.json"
      };
    std::vector<double> cutoffs{{1,2,3}};

    std::list<json> hypers{
      { 
        {"central_decay",0.5},
        {"interaction_cutoff",10},
        {"interaction_decay",0.5},
        {"size",200}
      }
      };
  };


  template< typename T,
            template<typename,typename> class FeatureManager,
            class StructureManager,
            template<typename> class RepresentationManager,
            class BaseFixture>
  struct FeatureFixture
  :RepresentationFixture<StructureManager,RepresentationManager, BaseFixture>
  {
    using Parent = RepresentationFixture<StructureManager,
                                  RepresentationManager, BaseFixture>;
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = typename Parent::Representation_t;
    using Feature_t = FeatureManager<T,Representation_t>;
    using hypers_t = typename Representation_t::hypers_t;

    FeatureFixture() :Parent{}
    { 
      std::vector<size_t> Nfeatures{};
      auto& representations = this->representations;
      for (auto& hyper : this->hypers){
        for (auto& manager : this->managers_strict){
          representations.emplace_back(manager,hyper);
          representations.back().compute();
          Nfeatures.push_back(representations.back().get_n_feature());
          std::cout << representations.back().get_n_feature()<< ", ";
        }
      }
      
      Nfeature = *std::max_element(std::begin(Nfeatures), std::end(Nfeatures));
      
    }

    ~FeatureFixture() = default;
    
    size_t Nfeature{};
    std::list<Feature_t> features{};
    
  };

/* ---------------------------------------------------------------------- */
  
} // RASCAL

#endif /* TEST_FEATURE_MANAGER_H */