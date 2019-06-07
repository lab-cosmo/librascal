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

    using Key_t = std::vector<int>;

    // TODO(felix) make sure it is not need anymore
    using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic,
                                  Eigen::RowMajor>;
    using InputData_t =
        internal::InternallySortedKeyMap<Key_t, Dense_t>;
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

    //! return the name of the calculator
    virtual std::string get_name() = 0;

    //! Compute the representation using a StructureManager
    // virtual void compute() = 0;

    /**
     * use the property interface to get a property from the manager to
     * fill it with a new representation.
     * it return a reference to the typed property of property_name
     */
    template<class StructureManager, template <class> class Property>
    inline decltype(auto) get_property(
          std::shared_ptr<StructureManager>& manager,
          const std::string& property_name) {
      // check if the property already exist and create it if it does not
      if (not manager->has_property()) {
        manager->template create_property<
            Property<StructureManager>>(property_name);
      }
      return manager->template get_validated_property_ref<Property<StructureManager>>(property_name);
    }


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
