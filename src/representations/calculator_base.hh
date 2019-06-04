/**
 * file   calculator_base.hh
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 September 2018
 *
 * @brief  base class for representation managers
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

#ifndef SRC_REPRESENTATIONS_CALCULATOR_BASE_HH_
#define SRC_REPRESENTATIONS_CALCULATOR_BASE_HH_

#include "structure_managers/structure_manager_base.hh"
#include "structure_managers/property_block_sparse.hh"
#include "json_io.hh"

#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <unordered_map>
#include <Eigen/Dense>

namespace rascal {

  class CalculatorBase {
   public:
    //! type for the hyper parameter class
    using Hypers_t = json;
    //! type for representation
    using Precision_t = double;
    //! type used to register the valid key and values of Hypers_t
    using ReferenceHypers_t = std::map<std::string, std::vector<std::string>>;

    // TODO(felix) make sure it is not need anymore
    using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic,
                                  Eigen::RowMajor>;
    using InputData_t =
        internal::InternallySortedKeyMap<std::vector<int>, Dense_t>;
    using Data_t = std::vector<InputData_t>;

    CalculatorBase() = default;

    //! Copy constructor
    CalculatorBase(const CalculatorBase & other) = delete;

    //! Move constructor
    CalculatorBase(CalculatorBase && other) = default;

    //! Destructor
    virtual ~CalculatorBase() = default;

    //! Copy assignment operator
    CalculatorBase & operator=(const CalculatorBase & other) = delete;

    //! Move assignment operator
    CalculatorBase & operator=(CalculatorBase && other) = default;

    //! Pure Virtual Function to set hyperparameters of the representation
    virtual void set_hyperparameters(const Hypers_t &) = 0;

    //! Pure Virtual Function to set hyperparameters of the representation
    void check_hyperparameters(const ReferenceHypers_t &, const Hypers_t &);

    //! Compute the representation using a StructureManager
    virtual void compute() = 0;

    //! get the raw data of the representation
    virtual std::vector<Precision_t> & get_representation_raw_data() = 0;

    virtual Data_t & get_representation_sparse_raw_data() = 0;

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
    Hypers_t hypers{};
    //! stores the hyperparameters that change
    //! the behaviour of the representation
    std::map<std::string, std::string> options{};
  };

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_BASE_HH_
