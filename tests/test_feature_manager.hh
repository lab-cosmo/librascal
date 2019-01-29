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
<<<<<<< HEAD
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
=======
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
>>>>>>> origin/master
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

<<<<<<< HEAD
=======
  /* ---------------------------------------------------------------------- */
  /**
   * A fixture providing access to different structures, read in from
   * json. These are used in a feature fixture to have data at hand.
   */
>>>>>>> origin/master
  struct TestFeatureData {
    TestFeatureData() = default;
    ~TestFeatureData() = default;

    std::vector<std::string> filenames{
      "reference_data/CaCrP2O7_mvc-11955_symmetrized.json",
      "reference_data/simple_cubic_8.json",
      "reference_data/small_molecule.json"
      };
<<<<<<< HEAD
    std::vector<double> cutoffs{{1, 2, 3}};

    std::list<json> hypers{
      {{"central_decay", 0.5},
      {"interaction_cutoff", 10},
      {"interaction_decay", 0.5},
      {"size", 120}}
      };
  };


  template< typename T,
            template<typename> class FeatureManager,
            class StructureManager,
            template<typename, Option ... opts > class RepresentationManager,
            class BaseFixture,
            Option ...options>
  struct FeatureFixture
  :RepresentationFixture<StructureManager, RepresentationManager,
                          BaseFixture, options...> {
    using Parent = RepresentationFixture<StructureManager,
                                    RepresentationManager,
                                    BaseFixture, options...>;
=======
    std::vector<double> cutoffs{{1., 2., 3.}};

    std::list<json> hypers{
      {{"central_decay", 0.5},
       {"interaction_cutoff", 10.},
       {"interaction_decay", 0.5},
       {"size", 120},
       {"sorting_algorithm", "distance"}},
      {{"central_decay", 0.5},
       {"interaction_cutoff", 10.},
       {"interaction_decay", 0.5},
       {"size", 120},
      {"sorting_algorithm", "row_norm"}}
    };
  };

  /* ---------------------------------------------------------------------- */
  /**
   * A templated Fixture, inherits from the ReperesentationFixture. It provides
   * access to different data structures. They are used to check the aggregation
   * of calculated feature data from multiple structures.
   */
  template<typename T,
           template<typename> class FeatureManager,
           class StructureManager,
           template<typename> class RepresentationManager,
           class BaseFixture>
  struct FeatureFixture
  :RepresentationFixture<StructureManager, RepresentationManager,
                         BaseFixture> {
    using Parent = RepresentationFixture<StructureManager,
                                         RepresentationManager,
                                         BaseFixture>;
>>>>>>> origin/master
    using Manager_t = typename Parent::Manager_t;
    using Representation_t = typename Parent::Representation_t;
    using Feature_t = FeatureManager<T>;
    using hypers_t = typename Representation_t::hypers_t;
<<<<<<< HEAD

    FeatureFixture() :Parent{} {
      std::vector<size_t> Nfeatures{};

      auto& representations = this->representations;

      for (auto& manager : this->managers_strict) {
        for (auto& hyper : this->hypers) {
          Ncenter += manager.size();

          representations.emplace_back(manager, hyper);
          representations.back().compute();
          Nfeatures.push_back(representations.back().get_feature_size());
        }
      }

      this->Nfeature =
           *std::max_element(std::begin(Nfeatures), std::end(Nfeatures));
    }

    ~FeatureFixture() = default;
    size_t Ncenter{0};
    size_t Nfeature{};
=======
    using precision_t = T;

    FeatureFixture() :Parent{} {
      std::vector<size_t> n_features{};

      auto & representations = this->representations;

      for (auto & manager : this->managers_strict) {
        for (auto & hyper : this->hypers) {
          n_center += manager.size();

          representations.emplace_back(manager, hyper);
          representations.back().compute();
          n_features.push_back(representations.back().get_feature_size());
        }
      }

      this->n_feature = *std::max_element(std::begin(n_features),
                                          std::end(n_features));
    }

    ~FeatureFixture() = default;
    size_t n_center{0};
    size_t n_feature{};
>>>>>>> origin/master
    std::list<Feature_t> features{};
  };

/* ---------------------------------------------------------------------- */

} // RASCAL

#endif /* TEST_FEATURE_MANAGER_H */
