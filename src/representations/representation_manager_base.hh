/**
 * file   representation_manager_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  base class for representation managers
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

#ifndef SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_BASE_HH_
#define SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_BASE_HH_

#include "structure_managers/structure_manager_base.hh"
#include "json_io.hh"

#include <string>
#include <vector>
#include <iostream>

namespace rascal {

  class RepresentationManagerBase {
   public:
    //! type for the hyper parameter class
    using hypers_t = json;
    //! type for representation
    using precision_t = double;
    //! type used to register the valid key and values of hypers_t
    using reference_hypers_t = std::map<std::string, std::vector<std::string>>;

    RepresentationManagerBase() = default;

    //! Copy constructor
    RepresentationManagerBase(const RepresentationManagerBase & other) = delete;

    //! Move constructor
    RepresentationManagerBase(RepresentationManagerBase && other) = default;

    //! Destructor
    virtual ~RepresentationManagerBase() = default;

    //! Copy assignment operator
    RepresentationManagerBase &
    operator=(const RepresentationManagerBase & other) = delete;

    //! Move assignment operator
    RepresentationManagerBase &
    operator=(RepresentationManagerBase && other) = default;

    //! Pure Virtual Function to set hyperparameters of the representation
    virtual void set_hyperparameters(const hypers_t &) = 0;

    //! Pure Virtual Function to set hyperparameters of the representation
    void check_hyperparameters(const reference_hypers_t &, const hypers_t &);

    //! Compute the representation using a StructureManager
    virtual void compute() = 0;

    //! get the raw data of the representation
    virtual std::vector<precision_t> & get_representation_raw_data() = 0;

    //! get the size of a feature vector
    virtual size_t get_feature_size() = 0;

    //! get the number of centers for the representation
    virtual size_t get_center_size() = 0;

    //! returns a string representation of the current options values
    //! in alphabetical order
    std::string get_options_string();

    //! returns a string representation of the current hypers dict
    std::string get_hypers_string();

    //! stores all the hyper parameters of the representation
    hypers_t hypers{};
    //! stores the hyperparameters that change
    //! the behaviour of the representation
    std::map<std::string, std::string> options{};
  };

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_REPRESENTATION_MANAGER_BASE_HH_
