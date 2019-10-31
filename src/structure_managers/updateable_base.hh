/**
 * @file   updateable_base.hh
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

  /**
   * Base class providing/defining the interface of an updatable object.
   * Adaptors, StructureManagers and ManagerSpecies are,
   * for instance updatable, objects since the underlying atomic structure
   * can change, e.g. along a MD run, with its 'neighbour list'.
   *
   * Updatable objects are linked together to form a tree with a root, some
   * branches and some number of leaves. Here is an illustration of such
   * structure:
   *                        ROOT
   *               /                  \
   *          brancheA1              brancheB1
   *            /    \                      \
   *    brancheA11 brancheA12              brancheB11
   *      /   \        /   \                /
   *  leave1 leave2 leave3 leave4        leave5
   *
   * The update mechanism allows to triger an update from
   * anywhere in the tree that will send the update signal to the root which
   * will in turn update itself and the whole tree.
   * The sending of the upward signal to the root is handled with the update
   * function implemented in each class inheriting Udatable. Since this
   * function involves the templated parent structure manager and that it is
   * itself templated it can't be defined explicitly in Updatable.
   * So this class mainly handles the parent to child relationship while
   * the child to parent relashionship is implementated in the classes
   * inheriting from Udatable.
   */
  class Updateable {
   public:
    using Children_t = std::weak_ptr<Updateable>;

    //! Default constructor sets the status variable for update to false.
    Updateable() : updated{false} {}

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

    /**
     * Implements how the updateable updates itself.
     */
    virtual void update_self() = 0;

    /**
     * When the underlying structure changes, all computations are potentially
     * invalid. This function triggers the setting of the status variable to
     * `false` along the tree to the managers and the properties it holds.
     */
    virtual void send_changed_structure_signal() = 0;

    //! Setter function for update statue variable
    void set_update_status(const bool sig) { this->updated = sig; }

    //! Getter function for update status variable.
    bool get_update_status() const { return this->updated; }

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
