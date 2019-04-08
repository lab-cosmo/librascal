/**
 * file   filter_base.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   25 Mar 2019
 *
 * @brief untemplated base class for filter type for use with
 * TupleStandardisation for easier looping over heterogeneous data
 *
 * Copyright © 2019 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_STRUCTURE_MANAGERS_FILTER_BASE_HH_
#define SRC_STRUCTURE_MANAGERS_FILTER_BASE_HH_

namespace rascal {
  class FilterBase {
   public:
    //! Default constructor
    FilterBase() = default;

    //! Copy constructor
    FilterBase(const FilterBase & other) = delete;

    //! Move constructor
    FilterBase(FilterBase && other) = default;

    //! Destructor
    virtual ~FilterBase() = default;

    //! Copy assignment operator
    FilterBase & operator=(const FilterBase & other) = delete;

    //! Move assignment operator
    FilterBase & operator=(FilterBase && other) = default;

    /**
     * Virtual base class to be called e.g. in species map to clear existing
     * filters when updating
     */
    virtual void reset_initial_state() = 0;

   protected:
    // intentionally left blank
   private:
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_FILTER_BASE_HH_
