// /**
//  * file   .hh
//  *
//  * @author Musil Felix <musil.felix@epfl.ch>
//  *
//  * @date   14 November 2018
//  *
//  * @brief  test  managers
//  *
//  * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
//  *
//  * Rascal is free software; you can redistribute it and/or
//  * modify it under the terms of the GNU Lesser General Public License as
//  * published by the Free Software Foundation, either version 3, or (at
//  * your option) any later version.
//  *
//  * Rascal is distributed in the hope that it will be useful, but
//  * WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  * Lesser General Public License for more details.
//  *
//  * You should have received a copy of the GNU Lesser General Public License
//  * along with this software; see the file LICENSE. If not, write to the
//  * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
//  * Boston, MA 02111-1307, USA.
//  */

// #ifndef TESTS_TEST_FEATURE_MANAGER_HH_
// #define TESTS_TEST_FEATURE_MANAGER_HH_

// #include "test_structure.hh"
// #include "test_representation_manager.hh"
// #include "representations/feature_manager_base.hh"
// #include "representations/feature_manager_dense.hh"
// #include "representations/feature_manager_block_sparse.hh"

// #include <list>

// namespace rascal {

//   /* ---------------------------------------------------------------------- */
//   /**
//    * A templated Fixture, inherits from the ReperesentationFixture. It provides
//    * access to different data structures. They are used to check the aggregation
//    * of calculated feature data from multiple structures.
//    */
//   template <typename T, template <typename> class FeatureManager,
//             template <typename> class RepresentationManager, class BaseFixture>
//   struct FeatureFixture
//       : RepresentationFixture<BaseFixture, RepresentationManager> {
//     using Parent = RepresentationFixture<BaseFixture, RepresentationManager>;
//     using Manager_t = typename Parent::Manager_t;
//     using Representation_t = typename Parent::Representation_t;
//     using Feature_t = FeatureManager<T>;
//     using Hypers_t = typename Representation_t::Hypers_t;
//     using Precision_t = T;

//     FeatureFixture() : Parent{} {
//       std::vector<size_t> n_features{};

//       auto & representations = this->representations;

//       for (auto & manager : this->managers) {
//         for (auto & hyper : this->hypers) {
//           n_center += manager->size();

//           representations.emplace_back(manager, hyper);
//           representations.back().compute();
//           // TODELETE
//           n_features.push_back(representations.back().get_feature_size());
//         }
//       }

//       this->n_feature =
//           *std::max_element(std::begin(n_features), std::end(n_features));
//     }

//     ~FeatureFixture() = default;
//     size_t n_center{0};
//     size_t n_feature{};
//     std::list<Feature_t> features{};
//   };

//   /* ---------------------------------------------------------------------- */
//   /**
//    * A templated Fixture, inherits from the ReperesentationFixture. It provides
//    * access to different data structures. They are used to check the aggregation
//    * of calculated feature data from multiple structures.
//    */
//   template <typename T, template <typename> class FeatureManager,
//             template <typename> class RepresentationManager, class BaseFixture>
//   struct SparseFeatureFixture
//       : RepresentationFixture<BaseFixture, RepresentationManager> {
//     using Parent = RepresentationFixture<BaseFixture, RepresentationManager>;
//     using Manager_t = typename Parent::Manager_t;
//     using Representation_t = RepresentationManager<Manager_t>;
//     using Key_t = typename Representation_t::Key_t;
//     using Feature_t = FeatureManager<T>;
//     using Hypers_t = typename Representation_t::Hypers_t;
//     using Precision_t = T;

//     SparseFeatureFixture() : Parent{} {
//       for (auto & hyper : this->hypers) {
//         this->representations.emplace_back();
//         for (auto & manager : this->managers) {
//           this->representations.back().emplace_back(manager, hyper);
//           this->representations.back().back().compute();
//         }
//         // TODELETE
//         this->inner_sizes.emplace_back(
//             this->representations.back().back().get_feature_size());
//       }
//     }

//     ~SparseFeatureFixture() = default;

//     std::vector<std::vector<Representation_t>> representations{};
//     std::vector<Feature_t> features{};
//     std::vector<size_t> inner_sizes{};
//   };

//   /* ---------------------------------------------------------------------- */

// }  // namespace rascal

// #endif  // TESTS_TEST_FEATURE_MANAGER_HH_
