/**
 * file   updateable_base.hh
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

#ifndef SRC_STRUCTURE_MANAGERS_UPDATEABLE_BASE_HH_
#define SRC_STRUCTURE_MANAGERS_UPDATEABLE_BASE_HH_

#include <memory>
#include <vector>

namespace rascal {

  //! base class for updatable functionality
  class Updateable {
   public:
    using Children_t = std::weak_ptr<Updateable>;

    //! Default constructor sets the status variable for update to false.
    Updateable() : updated{false} {};

    //! Copy constructor
    Updateable(const Updateable & other) = delete;

    //! Move constructor
    Updateable(Updateable && other) = default;

    //! Destructor
    virtual ~Updateable() = default;

    //! Copy assignment operator
    Updateable & operator=(const Updateable & other) = delete;

    //! Move assignment operator
    Updateable & operator=(Updateable && other) = default;

    /**
     * Trigger to update the tree of stacked StructureManager as well as
     * Adaptors and possible SpeciesManager and AdaptorFilters.
     */
    virtual void update_children() = 0;

    /**
     * Add a stacked adaptor as a child node of the current StructureManager
     * tree of children. `Children_t` is a weak pointer to a manager or adaptor.
     */
    void add_child(Children_t child) { this->children.emplace_back(child); }

    // //! Create a new shared pointer to `Updatable`
    // std::shared_ptr<Updateable> get_updateable_shared_ptr() {
    //   return this->shared_from_this();
    // }

    // //! Create a new weak pointer to `Updatable`
    // std::weak_ptr<Updateable> get_updateable_weak_ptr() {
    //   return std::weak_ptr<Updateable>(this->shared_from_this());
    // }

    // virtual void update_adaptor() = 0;
    /**
     * When the underlying structure changes, all computations are potentially
     * invalid. This function triggers the setting of the statue variable to
     * `false` along the tree. Should only be called from the root of the tree
     * (usually a StructureManager).
     */
    // virtual void send_changed_structure_signal() = 0;

    //! Setter function for update statue variable
    inline void set_update_status(const bool sig) { this->updated = sig; }

    //! Getter function for update status variable.
    inline bool get_update_status() const { return this->updated; }

   protected:
    //! List of children which are stacked on top of current object.
    std::vector<Children_t> children{};

    /**
     * Status variable for the state of update. Avoids updating adaptors, when
     * the underlying structure did not change, if .update() is called.
     */
    bool updated;

   private:
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_UPDATEABLE_BASE_HH_
