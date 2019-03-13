/**
 * file   updatable_base.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   13 Mar 2019
 *
 * @brief base class to provide update functionality across StructureManager,
 * SpeciesManager, AdaptorFilter to be used in a tree-update along the stacking.
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
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_STRUCTURE_MANAGERS_UPDATABLE_BASE_HH_
#define SRC_STRUCTURE_MANAGERS_UPDATABLE_BASE_HH_

namespace rascal {

  //! base class for updatable functionality
  class Updatable {
   public:
   protected:
    /**
     * Trigger to update the tree of stacked StructureManager as well as
     * Adaptors and possible SpeciesManager and AdaptorFilters.
     */
    void update_children() const = 0;
    /**
     * When the underlying structure changes, all computations are potentially
     * invalid. This function triggers the setting of the statue variable to
     * `false` along the tree. Should only be called from the root of the tree
     * (usually a StructureManager).
     */
    void send_changed_structure_signal() = 0;

    //! Setter function for update statue variable
    inline void set_is_up_to_date(const bool sig) {
      this->set_is_up_to_date = sig;
    }

    //! Getter function for update status variable.
    inline bool get_is_up_to_date() const { return this->is_up_to_date; }

    //! status variable for update
    bool is_up_to_date;

   private:
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_UPDATABLE_BASE_HH
