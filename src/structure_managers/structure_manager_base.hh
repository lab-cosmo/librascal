/**
 * file   structure_manager_base.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   06 Aug 2018
 *
 * @brief  Polymorphic base class for structure managers
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef STRUCTURE_MANAGER_BASE_H
#define STRUCTURE_MANAGER_BASE_H

#include <string>
#include <memory>

namespace rascal {

  //! polymorphic base class type for StructureManagers
  class StructureManagerBase {
   public:
    //! access to specific properties with a string
    inline decltype(auto) get_property(std::string name);
    //! only? virtual function needed.
    virtual size_t nb_clusters(size_t cluster_size) const = 0;
    virtual ~StructureManagerBase() = default;
    virtual void update_children() = 0;
    // virtual void update_adaptor() = 0;

   protected:

    };
}  // rascal

#endif /* STRUCTURE_MANAGER_BASE_H */
