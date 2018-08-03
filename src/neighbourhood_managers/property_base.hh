/**
 * file   property_base.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   03 Aug 2018
 *
 * @brief implementation of non-templated base class for Properties,
 *        Properties are atom-, pair-, triplet-, etc-related values
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef PROPERTY_BASE_H
#define PROPERTY_BASE_H

#include <typeinfo>

#include "basic_types.hh"

namespace rascal {

  class PropertyBase
  {
  public:
    //! Default constructor
    PropertyBase() = delete;

    //! Copy constructor
    PropertyBase(const PropertyBase &other) = delete;

    //! Move constructor
    PropertyBase(PropertyBase &&other) = default;

    //! Destructor
    virtual ~PropertyBase() = default;

    //! Copy assignment operator
    PropertyBase& operator=(const PropertyBase &other) = delete;

    //! Move assignment operator
    PropertyBase& operator=(PropertyBase &&other) = default;

    //! return runtime info about the stored (e.g., numerical) type
    virtual std::type_info & get_type_info() = 0;

    //! returns the number of degrees of freedom stored per cluster
    inline Dim_t get_nb_comp() const {return this->nb_comp;}

    //! returns the number of rows stored per cluster
    inline Dim_t get_nb_row() const {return this->nb_row;}

    //! returns the number of columns stored per cluster
    inline Dim_t get_nb_col() const {return this->nb_col;}

    //! returns the cluster order
    inline Dim_t get_order() const {return this->order;}


  protected:
    const Dim_t nb_row;  //!< number of rows stored
    const Dim_t nb_col;  //!< number of columns stored
    const Dim_t nb_comp; //!< number of dofs stored
    const size_t order;  //!< order of the clusters
  private:
  };

}  // rascal

#endif /* PROPERTY_BASE_H */
