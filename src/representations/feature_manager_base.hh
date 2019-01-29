/**
 * file  feature_manager_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief base class for storage of features from multiple managers
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


#ifndef FEATURE_MANAGER_BASE_H
#define FEATURE_MANAGER_BASE_H

#include "structure_managers/structure_manager_base.hh"
#include "representations/representation_manager_base.hh"
#include "json_io.hh"

#include <vector>
#include <stdexcept>
#include <iterator>


namespace rascal {
  /**
   * Base class of the Feature Managers. Defines the basic interface and some
   * common short hand types.
   */
  template<typename T>
  class FeatureManagerBase {
   public:
    using RepresentationManager_t = RepresentationManagerBase;
    using hypers_t = typename RepresentationManagerBase::hypers_t;
    using precision_t = T;
    using Feature_Matrix_t = Eigen::Matrix<precision_t, Eigen::Dynamic,
                                    Eigen::Dynamic, Eigen::ColMajor>;
    using Feature_Matrix_ref = Eigen::Map<const Feature_Matrix_t>;

    //! Default constructor
    FeatureManagerBase() = default;

    //! Copy constructor
    FeatureManagerBase(const FeatureManagerBase &other) = delete;

    //! Move constructor
    FeatureManagerBase(FeatureManagerBase &&other) = delete;

    //! Destructor
    virtual ~FeatureManagerBase() = default;

    //! Copy assignment operator
    FeatureManagerBase& operator=(const FeatureManagerBase & other) = delete;

    //! Move assignment operator
    FeatureManagerBase& operator=(FeatureManagerBase && other) = delete;

    //! pre-allocate memory
    virtual void reserve(size_t&) = 0;

    //! move data from the representation manager property
    virtual void push_back(RepresentationManager_t&) = 0;

    //! return number of elements of the flattened array
    virtual inline int size() = 0;

    //! return the number of samples in the feature matrix
    virtual inline int sample_size() = 0;

    //! return the number of feature in the feature matrix
    virtual inline int feature_size() = 0;

    //! get the shape of the feature matrix (Nrow,Ncol)
    virtual inline std::tuple<int, int> shape() = 0;

    //! expose the feature matrix
    virtual inline Feature_Matrix_ref get_feature_matrix() = 0;
  };


} // rascal

#endif /* FEATURE_MANAGER_BASE_H */
