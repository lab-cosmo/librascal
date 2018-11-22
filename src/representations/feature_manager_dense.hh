/**
 * file   
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief  
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


#ifndef FEATURE_MANAGER_DENSE_H
#define FEATURE_MANAGER_DENSE_H


#include "representations/feature_manager_base.hh"


namespace rascal {

template<typename T, typename RepresentationManager>
class FeatureManagerDense: public FeatureManagerBase {
  public:

    using RepresentationManager_t = RepresentationManager;
    using hypers_t = typename RepresentationManager::hypers_t;
    using Feature_Matrix_t = Eigen::MatrixXd;
    using Feature_Matrix_ref = Eigen::Map<Eigen::MatrixXd>;

    FeatureManagerDense(int n_feature, hypers_t hypers)
    :feature_matrix{},n_feature{n_feature},n_center{0},hypers{hypers}
    {}

    //! Copy constructor
    FeatureManagerDense(const FeatureManagerDense &other) = delete;

    //! Move constructor
    FeatureManagerDense(FeatureManagerDense &&other) = default;

    //! Destructor
    ~FeatureManagerDense() = default;

    //! Copy assignment operator
    FeatureManagerDense& operator=(const FeatureManagerDense &other) = delete;

    //! Move assignment operator
    FeatureManagerDense& operator=(FeatureManagerDense && other) = default;

    //! pre-allocate memory
    void reserve(int n_center){
      this->feature_matrix.reserve(n_center*this->n_feature);
    }

    //! move data from the representation manager property 
    void push_back(RepresentationManager_t& rm){
      auto& property{rm.get_property()};
      auto& raw_data{property.get_raw_data()};
      auto n_elem{property.get_nb_item()};
    
      int n_feature{property.get_nb_comp()};
      if (n_feature != this->n_feature){
        throw std::length_error("Incompatible number of features");
      }
      this->n_center += n_elem;
      this->feature_matrix.insert(this->feature_matrix.end(),
                        std::make_move_iterator(raw_data.begin()), 
                        std::make_move_iterator(raw_data.end())
                        );
    }

    //! move data from a feature vector
    void push_back(std::vector<T> feature_vector){
      int n_feature{feature_vector.size()};
      if (n_feature != this->n_feature){
        throw std::length_error("Incompatible number of features");
      }
      this->n_center += 1;
      this->feature_matrix.insert(this->feature_matrix.end(),
                    std::make_move_iterator(feature_vector.begin()), 
                    std::make_move_iterator(feature_vector.end()));
    }

    inline int get_nb_comp(){
      return this->feature_matrix.size();
    }

    inline int get_nb_center(){
      return this->feature_matrix.size()/this->n_feature;
    }

    inline Feature_Matrix_ref get_feature_matrix(){
      return Feature_Matrix_ref(this->feature_matrix.data(),
                                this->n_feature,this->n_center);
    }

  protected:
    std::vector<T> feature_matrix;
    int n_feature;
    int n_center;
    hypers_t hypers;
  
};


} // rascal

#endif /* FEATURE_MANAGER_DENSE_H */
