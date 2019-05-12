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

#ifndef SRC_REPRESENTATIONS_FEATURE_MANAGER_BLOCK_SPARSE_HH_
#define SRC_REPRESENTATIONS_FEATURE_MANAGER_BLOCK_SPARSE_HH_

#include "representations/feature_manager_base.hh"
#include "representations/representation_manager_base.hh"
#include "structure_managers/property_block_sparse.hh"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>

namespace rascal {
  /**
   * Handles the aggregation of features from compatible representation
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
  template <typename T>
  class FeatureManagerBlockSparse : public FeatureManagerBase<T> {
   public:
    using Parent = FeatureManagerBase<T>;
    using RepresentationManager_t = typename Parent::RepresentationManager_t;
    using Hypers_t = typename Parent::Hypers_t;
    using Feature_Matrix_t = typename Parent::Feature_Matrix_t;
    using Feature_Matrix_ref = typename Parent::Feature_Matrix_ref;
    using Precision_t = typename Parent::Precision_t;

    using Key_t = std::vector<int>;
    using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic>;
    using dense_ref_t = Eigen::Map<Dense_t>;
    using MapCenter_t = std::vector<std::pair<size_t, size_t>>;
    using map_sparse_t = std::vector<
        std::map<Key_t, std::pair<size_t, std::pair<size_t, size_t>>>>;
    using Keys_t = std::vector<std::list<Key_t>>;
    using Data_t = std::vector<Precision_t>;

    using Vector_ref =
        Eigen::Map<const Eigen::Matrix<Precision_t, Eigen::Dynamic, 1>>;

    /**
     * Constructor where hypers contains all relevant informations
     * to setup a new RepresentationManager.
     */
    FeatureManagerBlockSparse(size_t inner_size, Hypers_t hypers)
        : feature_matrix{}, inner_size{inner_size}, n_center{0}, hypers{
                                                                     hypers} {}

    //! Copy constructor
    FeatureManagerBlockSparse(const FeatureManagerBlockSparse & other) = delete;

    //! Move constructor
    FeatureManagerBlockSparse(FeatureManagerBlockSparse && other) = default;

    //! Destructor
    ~FeatureManagerBlockSparse() = default;

    //! Copy assignment operator
    FeatureManagerBlockSparse &
    operator=(const FeatureManagerBlockSparse & other) = delete;

    //! Move assignment operator
    FeatureManagerBlockSparse &
    operator=(FeatureManagerBlockSparse && other) = default;

    /**
     * resize the underlying data structure
     *
     * @param n_centers list (for each structures) of the number of centers
     * @param nb_unique_keys list (for each structures) of number of unique
     * keys.
     */
    void resize(const std::vector<size_t> & n_centers,
                const std::vector<size_t> & nb_unique_keys) {
      size_t n_elements{this->values.size()};
      auto n_inner_comp{this->get_inner_size()};
      size_t i_structure{0};
      // there are as many set of keys as structures
      for (const auto & nb_keys : nb_unique_keys) {
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
    void resize(const std::vector<size_t> & nb_unique_keys) {
      size_t n_elements{this->values.size()};
      auto n_inner_comp{this->get_inner_size()};
      // there are as many set of keys as centers
      for (const auto & nb_keys : nb_unique_keys) {
        n_elements += nb_keys * n_inner_comp;
        this->n_center++;
      }
      this->resize(n_elements);
    }

    void resize(const size_t & n_elements) {
      this->feature_matrix.resize(n_elements, 0.);
      this->map2centers.resize(this->n_center);
      this->map2sparse.resize(this->n_center);
      this->keys_list.resize(this->n_center);
    }

    void reserve(size_t & n_elements) {
      this->feature_matrix.reserve(n_elements);
      this->map2centers.reserve(this->n_center);
      this->map2sparse.reserve(this->n_center);
      this->keys_list.reserve(this->n_center);
    }
    //! move data from the representation manager property
    void push_back(RepresentationManager_t & rm) {
      const auto & raw_data{rm.get_representation_sparse_raw_data()};
      auto n_center{rm.get_center_size()};
      this->structure_stride.push_back(n_center + structure_stride.back());

      auto && new_center_start_id{this->feature_matrix.size()};

      for (size_t i_center{0}; i_center < n_center; i_center++) {
        this->map2sparse.emplace_back();
        this->keys_list.emplace_back();
        size_t key_start{0};
        for (const auto & element : raw_data[i_center]) {
          const auto & key{element.first};
          const auto & value{element.second};

          for (int i_col{0}; i_col < value.cols(); i_col++) {
            for (int i_row{0}; i_row < value.rows(); i_row++) {
              this->feature_matrix.push_back(value(i_row, i_col));
            }
          }
          this->unique_keys.emplace(key);
          this->keys_list.back().emplace_back(key);
          this->map2sparse.back().emplace(std::make_pair(
              key, std::make_pair(key_start,
                                  std::make_pair(value.rows(), value.cols()))));

          key_start += value.size();
        }
        size_t center_length{key_start};
        this->map2centers.emplace_back(
            std::make_pair(new_center_start_id, center_length));
        new_center_start_id += center_length;
      }
      this->n_center += n_center;
    }

    //! move data from the representation manager property
    void insert(size_t &, RepresentationManager_t &) {}

    //! return number of elements of the flattened array
    inline int size() { return this->feature_matrix.size(); }

    //! return the number of samples in the feature matrix
    inline int sample_size() { return this->n_center; }

    //! return the number of feature in the feature matrix
    inline int feature_size() {
      return this->unique_keys.size() * this->get_inner_size();
    }

    inline auto get_inner_size() { return this->inner_size; }
    //! get the shape of the feature matrix (Nrow,Ncol)
    inline std::tuple<int, int> shape() {
      return std::make_tuple(this->sample_size(), this->feature_size());
    }

    //! return the feature matrix as an Map over Eigen MatrixXd
    inline Feature_Matrix_t get_feature_matrix_dense() {
      Feature_Matrix_t mat =
          Feature_Matrix_t::Zero(this->feature_size(), this->sample_size());
      auto inner_size{this->get_inner_size()};
      // loop center
      for (int i_center{0}; i_center < this->sample_size(); i_center++) {
        auto & center_start{this->map2centers[i_center].first};
        size_t i_feature{0};
        // loop sparse key
        for (auto & key : this->unique_keys) {
          // if the key exist
          if (this->map2sparse[i_center].count(key) == 1) {
            auto key_start{center_start +
                           this->map2sparse[i_center][key].first};
            auto & key_shape{this->map2sparse[i_center][key].second};
            auto key_size{key_shape.first * key_shape.second};
            for (size_t i_val{0}; i_val < key_size; i_val++) {
              mat(i_feature, i_center) =
                  this->feature_matrix[key_start + i_val];
              i_feature++;
            }
          } else {
            i_feature += inner_size;
          }
        }  // keys
      }    // centers
      return mat;
    }

    inline Feature_Matrix_ref get_feature_matrix() {
      return Feature_Matrix_ref(this->feature_matrix.data(),
                                this->feature_matrix.size(), 1);
    }

    inline auto get_structure_stride() { return this->structure_stride; }

    inline const std::list<Key_t> & get_keys(const int & i_center) {
      return this->keys_list[i_center];
    }

    inline Vector_ref get_features(const int & i_center, const Key_t & key) {
      auto && center_start{this->map2centers[i_center].first};
      auto key_start{center_start + this->map2sparse[i_center][key].first};
      auto && key_shape{this->map2sparse[i_center][key].second};
      auto key_size{key_shape.first * key_shape.second};
      return Vector_ref(&this->feature_matrix[key_start], key_size);
    }

   protected:
    //! underlying data container for the feature matrix
    std::vector<T> feature_matrix;

    //! number of elements in the inner dimension
    size_t inner_size;
    //! Number of samples in the feature matrix
    int n_center;

    /**
     * Contain all relevant information to initialize a compatible
     * RepresentationManager
     */
    Hypers_t hypers;

    //! stride for the structures
    std::vector<int> structure_stride{0};
    //! map center idx to position and size in values
    MapCenter_t map2centers{};
    //! map center idx + key to relative position and size in values
    map_sparse_t map2sparse{};
    //! list of the registered keys for each centers
    Keys_t keys_list{};

    std::set<Key_t> unique_keys{};
  };

  template <typename Precision_t>
  using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic>;

  template <typename T>
  std::list<T> intersection_of(const std::list<T> & a, const std::list<T> & b) {
    std::list<T> rtn;
    std::unordered_multiset<T, internal::Hash<T>> st;
    std::for_each(a.begin(), a.end(), [&st](const T & k) { st.insert(k); });
    std::for_each(b.begin(), b.end(), [&st, &rtn](const T & k) {
      auto iter = st.find(k);
      if (iter != st.end()) {
        rtn.push_back(k);
        st.erase(iter);
      }
    });
    return rtn;
  }

  /**
   * Matrix multiplication between two block sparse feature matrices
   */
  template <typename Precision_t>
  decltype(auto) dot(FeatureManagerBlockSparse<Precision_t> & X1,
                     FeatureManagerBlockSparse<Precision_t> & X2) {
    auto n1_centers{X1.sample_size()};
    auto n2_centers{X2.sample_size()};
    Dense_t<Precision_t> mat =
        Dense_t<Precision_t>::Zero(n1_centers, n2_centers);
    for (int icenter1{0}; icenter1 < n1_centers; icenter1++) {
      auto && keys1{X1.get_keys(icenter1)};
      for (int icenter2{0}; icenter2 < n2_centers; icenter2++) {
        auto && keys2{X2.get_keys(icenter2)};
        auto keys_intersection{intersection_of(keys1, keys2)};
        for (const auto & key : keys_intersection) {
          auto vec1{X1.get_features(icenter1, key)};
          auto vec2{X2.get_features(icenter2, key)};
          mat(icenter1, icenter2) += vec1.dot(vec2);
        }
      }
    }
    return mat;
  }

  /**
   * Matrix multiplication between a block sparse feature matrix with itself
   */
  template <typename Precision_t>
  decltype(auto) dot(FeatureManagerBlockSparse<Precision_t> & X1) {
    auto n1_centers{X1.sample_size()};
    Dense_t<Precision_t> mat =
        Dense_t<Precision_t>::Zero(n1_centers, n1_centers);
    for (int icenter1{0}; icenter1 < n1_centers; icenter1++) {
      auto && keys1{X1.get_keys(icenter1)};
      mat(icenter1, icenter1) = 1.;

      for (int icenter2{icenter1 + 1}; icenter2 < n1_centers; icenter2++) {
        auto && keys2{X1.get_keys(icenter2)};
        auto keys_intersection{intersection_of(keys1, keys2)};
        for (const auto & key : keys_intersection) {
          auto vec1{X1.get_features(icenter1, key)};
          auto vec2{X1.get_features(icenter2, key)};
          mat(icenter1, icenter2) += vec1.dot(vec2);
        }
        mat(icenter2, icenter1) = mat(icenter1, icenter2);
      }
    }
    return mat;
  }

  /**
   * global cosine kernel with a power of zeta
   */
  template <typename Precision_t>
  decltype(auto)
  cosine_kernel_global(const int & zeta,
                       FeatureManagerBlockSparse<Precision_t> & X1,
                       FeatureManagerBlockSparse<Precision_t> & X2) {
    auto stride1{X1.get_structure_stride()};
    auto stride2{X2.get_structure_stride()};
    auto n1_frames{stride1.size() - 1};
    auto n2_frames{stride2.size() - 1};

    Dense_t<Precision_t> mat(n1_frames, n2_frames);
    // loop over the structures of the 1st set
    for (size_t iframe1{0}; iframe1 < n1_frames; iframe1++) {
      auto & iframe1_st{stride1[iframe1]};
      auto & iframe1_nd{stride1[iframe1 + 1]};
      // loop over the structures of the 2nd set
      for (size_t iframe2{0}; iframe2 < n2_frames; iframe2++) {
        auto & iframe2_st{stride2[iframe2]};
        auto & iframe2_nd{stride2[iframe2 + 1]};
        double kernel_val{0.};
        // loop over the centers of the the ith structure
        for (int icenter1{iframe1_st}; icenter1 < iframe1_nd; icenter1++) {
          auto && keys1{X1.get_keys(icenter1)};
          // loop over the centers of the the jth structure
          for (int icenter2{iframe2_st}; icenter2 < iframe2_nd; icenter2++) {
            auto && keys2{X2.get_keys(icenter2)};
            auto keys_intersection{intersection_of(keys1, keys2)};
            // compute the dot product (cosine) between 2 normalized sparse
            // vector
            double inner_val{0.};
            for (const auto & key : keys_intersection) {
              auto vec1{X1.get_features(icenter1, key)};
              auto vec2{X2.get_features(icenter2, key)};
              inner_val += vec1.dot(vec2);
            }
            // raise the environmental kernel value to a power
            kernel_val += std::pow(inner_val, zeta);
          }
        }
        double n_centers2{static_cast<double>((iframe1_nd - iframe1_st) *
                                              (iframe2_nd - iframe2_st))};
        mat(iframe1, iframe2) = kernel_val / n_centers2;
      }
    }
    return mat;
  }

  /**
   * global cosine kernel with a power of zeta with itself
   */
  template <typename Precision_t>
  decltype(auto)
  cosine_kernel_global(const int & zeta,
                       FeatureManagerBlockSparse<Precision_t> & X1) {
    auto stride1{X1.get_structure_stride()};
    auto stride2{X1.get_structure_stride()};
    auto n1_frames{stride1.size() - 1};
    auto n2_frames{stride2.size() - 1};

    Dense_t<Precision_t> mat(n1_frames, n2_frames);
    // loop over the structures of the 1st set
    for (size_t iframe1{0}; iframe1 < n1_frames; iframe1++) {
      auto & iframe1_st{stride1[iframe1]};
      auto & iframe1_nd{stride1[iframe1 + 1]};

      double kernel_diag{0.};
      // compute the similarity between structure i and itself
      for (int icenter1{iframe1_st}; icenter1 < iframe1_nd; icenter1++) {
        auto && keys_intersection{X1.get_keys(icenter1)};
        // add the diagonal term which is 1 because normalized feature vector
        kernel_diag += 1.;

        // loop over the centers again
        for (int icenter2{icenter1 + 1}; icenter2 < iframe1_nd; icenter2++) {
          // compute the dot product (cosine) between 2 normalized sparse vector
          double inner_val{0.};
          for (const auto & key : keys_intersection) {
            auto vec1{X1.get_features(icenter1, key)};
            auto vec2{X1.get_features(icenter2, key)};
            inner_val += vec1.dot(vec2);
          }
          // raise the environmental kernel value to a power and take into
          // account the upper and lower contribution
          kernel_diag += 2 * std::pow(inner_val, zeta);
        }

        double n_centers2{static_cast<double>((iframe1_nd - iframe1_st) *
                                              (iframe1_nd - iframe1_st))};
        mat(iframe1, iframe1) = kernel_diag / n_centers2;
      }

      // loop over the structures of the 2nd set
      for (size_t iframe2{iframe1 + 1}; iframe2 < n2_frames; iframe2++) {
        auto & iframe2_st{stride2[iframe2]};
        auto & iframe2_nd{stride2[iframe2 + 1]};
        double kernel_val{0.};
        // loop over the centers of the the ith structure
        for (int icenter1{iframe1_st}; icenter1 < iframe1_nd; icenter1++) {
          auto && keys1{X1.get_keys(icenter1)};
          // loop over the centers of the the jth structure
          for (int icenter2{iframe2_st}; icenter2 < iframe2_nd; icenter2++) {
            auto && keys2{X1.get_keys(icenter2)};
            auto keys_intersection{intersection_of(keys1, keys2)};
            // compute the dot product (cosine) between 2 normalized sparse
            // vector
            double inner_val{0.};
            for (const auto & key : keys_intersection) {
              auto vec1{X1.get_features(icenter1, key)};
              auto vec2{X1.get_features(icenter2, key)};
              inner_val += vec1.dot(vec2);
            }
            // raise the environmental kernel value to a power
            kernel_val += std::pow(inner_val, zeta);
          }
        }
        double n_centers2{static_cast<double>((iframe1_nd - iframe1_st) *
                                              (iframe2_nd - iframe2_st))};
        mat(iframe1, iframe2) = kernel_val / n_centers2;
        mat(iframe2, iframe1) = mat(iframe1, iframe2);
      }
    }
    return mat;
  }

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_FEATURE_MANAGER_BLOCK_SPARSE_HH_
