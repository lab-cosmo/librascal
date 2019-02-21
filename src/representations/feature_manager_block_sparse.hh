/**
 * file feature_manager_dense.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief Generic manager aimed to aggregate the features computed
 *  with a representation on one or more atomic structures
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

#ifndef SRC_REPRESENTATIONS_FEATURE_MANAGER_DENSE_HH_
#define SRC_REPRESENTATIONS_FEATURE_MANAGER_DENSE_HH_

#include "representations/feature_manager_base.hh"
#include "representations/representation_manager_base.hh"

#include <unordered_map>
#include <set>

namespace rascal {
  /**
   * Handles the aggragation of features from compatible representation
   * managers.
   * The storage layout can be represented in 2 ways:
   * + [centers, sparse_keys, row_inner, col_inner]
   * + [centers, n_features]
   * The first view is what is actually stored in a vector with some helper
   * maps to mimic this layout.
   * centers corresponds to the number of central atoms accross all the
   * aggregated structures.
   * sparse_keys is the sparse dimension and can have a varying number of
   * actually stored keys.
   * row_inner and col_inner are the inner dense dimensions that are fixed for
   * all centers and sparse_keys.
   * The 2nd view is a consolidated (dense) copy of the 1st where missing
   * sparse keys have been filled with zeros.
   */
  template <typename T, typename key_t>
  class FeatureManagerBlockSparse : public FeatureManagerBase<T> {
   public:
    using Parent = FeatureManagerBase<T>;
    using RepresentationManager_t = typename Parent::RepresentationManager_t;
    using hypers_t = typename Parent::hypers_t;
    using Feature_Matrix_t = typename Parent::Feature_Matrix_t;
    using Feature_Matrix_ref = typename Parent::Feature_Matrix_ref;
    using precision_t = typename Parent::precision_t;

    using dense_t = Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic>;
    using dense_ref_t = Eigen::Map<dense_t>;
    using map_center_t = std::vector<std::pair<size_t, size_t>>;
    using map_sparse_t = std::vector<std::unordered_map<key_t, std::pair<size_t, std::pair<size_t, size_t>>>>;
    using keys_t = std::vector<std::list<key_t>>;
    using data_t = std::vector<precision_t>;

    /**
     * Constructor where hypers contains all relevant informations
     * to setup a new RepresentationManager.
     */
    FeatureManagerBlockSparse(size_t row_inner, size_t col_inner, hypers_t hypers)
        :feature_matrix{}, row_inner{row_inner}, col_inner{col_inner},  n_feature{row_inner*col_inner}, n_center{0}, hypers{hypers} {}

    //! Copy constructor
    FeatureManagerBlockSparse(const FeatureManagerBlockSparse & other) = delete;

    //! Move constructor
    FeatureManagerBlockSparse(FeatureManagerBlockSparse && other) = default;

    //! Destructor
    ~FeatureManagerBlockSparse() = default;

    //! Copy assignment operator
    FeatureManagerBlockSparse & operator=(const FeatureManagerBlockSparse & other) = delete;

    //! Move assignment operator
    FeatureManagerBlockSparse & operator=(FeatureManagerBlockSparse && other) = default;

    /**
     * resize the underlying data structure
     *
     * @param n_centers list (for each structures) of the number of centers
     * @param nb_unique_keys list (for each structures) of number of unique
     * keys.
     */
    void resize(const std::vector<size_t>& n_centers, const std::vector<size_t>& nb_unique_keys) {
      size_t n_elements{this->values.size()};
      auto n_inner_comp{this->row_inner*this->col_inner};
      size_t i_structure{0};
      // there are as many set of keys as structures
      for (const auto& nb_keys : nb_unique_keys) {
        n_elements += nb_keys * n_inner_comp * n_centers[i_structure];
        this->n_center += n_centers[i_structure];
        i_structure++;
      }
      this->resize(n_elements);
    }

    /**
     * resize the underlying data structure
     *
     * @param nb_unique_keys list (for each centers) of number of unique
     * keys.
     */
    void resize(const std::vector<size_t>& nb_unique_keys) {
      size_t n_elements{this->values.size()};
      auto n_inner_comp{this->row_inner*this->col_inner};
      // there are as many set of keys as centers
      for (const auto& nb_keys : nb_unique_keys) {
        n_elements += nb_keys * n_inner_comp;
        this->n_center++;
      }
      this->resize(n_elements);
    }

    void resize(const size_t& n_elements) {
      this->feature_matrix.resize(n_elements, 0.);
      this->map2centers.resize(this->n_center);
      this->map2sparse.resize(this->n_center);
      this->keys_list.resize(this->n_center);
    }

    //! move data from the representation manager property
    void push_back(RepresentationManager_t & rm) {
      auto & raw_data{rm.get_representation_raw_data()};
      auto n_center{rm.get_center_size()};
      auto unique_keys{rm.get_all_unique_keys()};

      this->n_center += n_center;
      this->feature_matrix.insert(this->feature_matrix.end(),
                                  std::make_move_iterator(raw_data.begin()),
                                  std::make_move_iterator(raw_data.end()));
    }

    //! return number of elements of the flattened array
    inline int size() { return this->feature_matrix.size(); }

    //! return the number of samples in the feature matrix
    inline int sample_size() { return this->size() / this->n_feature; }

    //! return the number of feature in the feature matrix
    inline int feature_size() { return this->n_feature; }

    //! get the shape of the feature matrix (Nrow,Ncol)
    inline std::tuple<int, int> shape() {
      return std::make_tuple(this->sample_size(), this->feature_size());
    }

    //! return the feature matrix as an Map over Eigen MatrixXd
    inline Feature_Matrix_ref get_feature_matrix() {
      return Feature_Matrix_ref(this->feature_matrix.data(),
                                this->feature_size(), this->sample_size());
    }

   protected:
    //! underlying data container for the feature matrix
    std::vector<T> feature_matrix;

    //! number of rows in the inner dimension
    size_t row_inner;
    //! number of columns in the inner dimension
    size_t col_inner;
    //! Number of samples in the feature matrix
    int n_center;

    /**
     * Contain all relevant information to initialize a compatible
     * RepresentationManager
     */
    hypers_t hypers;

    //! map center idx to position and size in values
    map_center_t map2centers{};
    //! map center idx + key to relative position and size in values
    map_sparse_t map2sparse{};
    //! list of the registered keys for each centers
    keys_t keys_list{};
  };

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_FEATURE_MANAGER_DENSE_HH_
