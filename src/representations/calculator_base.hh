/**
 * @file   calculator_base.hh
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

#include "json_io.hh"
#include "structure_managers/property_block_sparse.hh"
#include "structure_managers/structure_manager_base.hh"

#include <Eigen/Dense>

#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

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

    //!
    void check_hyperparameters(const ReferenceHypers_t &, const Hypers_t &);

    //! return the name of the calculator
    const std::string & get_name() const { return this->name; }

    //! return the name of the calculator's gradients
    std::string get_gradient_name() const {
      return this->name + std::string("_gradients");
    }

    const std::string & get_prefix() const { return this->default_prefix; }

    //! set the name of the calculator
    void set_name(const std::string & name) { this->name = name; }
    //! set the prefix for the default naming of the representation
    void set_default_prefix(const std::string & default_prefix) {
      this->default_prefix = default_prefix;
    }
    /**
     * identifier is a user defined name for the representation used to
     * register the computed representation
     * if not provided use the hypers to generate a unique identifier
     */
    void set_name(const Hypers_t & hyper) {
      if (hyper.count("identifier") == 1) {
        this->set_name(this->default_prefix +
                       hyper["identifier"].get<std::string>());
      } else {
        this->set_name(this->default_prefix + hyper.dump());
      }
    }

    /**
     * Computes the representation associated to the input structure
     * manager. It is templated so it can't be virtual but it is still
     * expected.
     */
    // template<class StructureManager>
    // virtual void compute(StructureManager& ) = 0;

    //! returns a string representation of the current options values
    //! in alphabetical order
    std::string get_options_string();

    //! returns a string representation of the current hypers dict
    std::string get_hypers_string();

    //! name of the calculator
    std::string name{""};
    //! default prefix of the calculator
    std::string default_prefix{""};

    //! stores all the hyper parameters of the representation
    Hypers_t hypers{};
    //! stores the hyperparameters that change
    //! the behaviour of the representation
    std::map<std::string, std::string> options{};
  };

}  // namespace rascal

#endif  // SRC_REPRESENTATIONS_CALCULATOR_BASE_HH_
