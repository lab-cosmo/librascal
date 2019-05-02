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
 * Copyright  2018 Musil Felix, COSMO (EPFL), LAMMM (EPFL)
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

namespace rascal {
  /**
   * Handles the aggragation of features from compatible representation managers
   * using a dense underlying data storage.
   */
  template <typename T>
  class FeatureManagerDense : public FeatureManagerBase<T> {
   public:
    using Parent = FeatureManagerBase<T>;
    using RepresentationManager_t = typename Parent::RepresentationManager_t;
    using Hypers_t = typename RepresentationManager_t::Hypers_t;
    using Feature_Matrix_t = typename Parent::Feature_Matrix_t;
    using Feature_Matrix_ref = typename Parent::Feature_Matrix_ref;
    using Precision_t = typename Parent::Precision_t;
    /**
     * Constructor where hypers contains all relevant informations
     * to setup a new RepresentationManager.
     */
    FeatureManagerDense(int n_feature, Hypers_t hypers)
        : feature_matrix{}, n_feature{n_feature}, n_center{0}, hypers{hypers} {}

    //! Copy constructor
    FeatureManagerDense(const FeatureManagerDense & other) = delete;

    //! Move constructor
    FeatureManagerDense(FeatureManagerDense && other) = default;

    //! Destructor
    ~FeatureManagerDense() = default;

    //! Copy assignment operator
    FeatureManagerDense & operator=(const FeatureManagerDense & other) = delete;

    //! Move assignment operator
    FeatureManagerDense & operator=(FeatureManagerDense && other) = default;

    void check_feature_compatibility(int & a, int & b) {
      if (a != b) {
        auto error = std::string("Incompatible number of features ") +
                     std::to_string(a) + std::string(" != ") +
                     std::to_string(b);
        throw std::length_error(error);
      }
    }

    //! pre-allocate memory
    void reserve(size_t & n_center) {
      this->feature_matrix.reserve(n_center * this->n_feature);
    }

    //! pre-allocate memory
    void resize(size_t & n_center) {
      this->feature_matrix.resize(n_center * this->n_feature, 0.);
      this->n_center += n_center;
    }

    void insert(size_t & i_center, RepresentationManager_t & rm) {
      auto & raw_data{rm.get_representation_raw_data()};
      auto n_feature{static_cast<int>(rm.get_feature_size())};
      auto n_centers{rm.get_center_size()};
      this->n_center += n_centers;
      this->check_feature_compatibility(n_feature, this->n_feature);

      auto iter{this->feature_matrix.begin()};
      iter += i_center * this->n_feature;
      this->feature_matrix.insert(iter,
                                  std::make_move_iterator(raw_data.begin()),
                                  std::make_move_iterator(raw_data.end()));
    }

    void assign(size_t & i_center, RepresentationManager_t & rm) {
      auto & raw_data{rm.get_representation_raw_data()};
      auto n_feature{static_cast<int>(rm.get_feature_size())};
      auto n_center{rm.get_center_size()};
      this->check_feature_compatibility(n_feature, this->n_feature);

      size_t i_end{n_center * this->n_feature};
      size_t i_start{i_center * this->n_feature};
      for (size_t i_el{0}; i_el < i_end; i_el++) {
        this->feature_matrix[i_start + i_el] = raw_data[i_el];
      }
    }

    //! move data from the representation manager property
    void push_back(RepresentationManager_t & rm) {
      auto & raw_data{rm.get_representation_raw_data()};
      auto n_center{rm.get_center_size()};
      int n_feature{static_cast<int>(rm.get_feature_size())};
      this->check_feature_compatibility(n_feature, this->n_feature);

      this->n_center += n_center;
      this->feature_matrix.insert(this->feature_matrix.end(),
                                  std::make_move_iterator(raw_data.begin()),
                                  std::make_move_iterator(raw_data.end()));
    }

    //! move data from a feature vector
    void push_back(std::vector<T> feature_vector) {
      int n_feature{static_cast<int>(feature_vector.size())};
      this->check_feature_compatibility(n_feature, this->n_feature);
      this->n_center += 1;
      this->feature_matrix.insert(
          this->feature_matrix.end(),
          std::make_move_iterator(feature_vector.begin()),
          std::make_move_iterator(feature_vector.end()));
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

    inline Feature_Matrix_t get_feature_matrix_dense() {
      return Feature_Matrix_t::Zero(2, 2);
    }

   protected:
    //! underlying data container for the feature matrix
    std::vector<T> feature_matrix;
    //! Number of feature.
    // TODO(felix) make it possible to change it after construction
    int n_feature;
    //! Number of samples in the feature matrix
    int n_center;
    /**
     * Contain all relevant information to initialize a compatible
     * RepresentationManager
     */
    Hypers_t hypers;
  };

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_FEATURE_MANAGER_DENSE_HH_
