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

#ifndef BASIS_REPRESENTATION_MANAGER_BASE_H
#define BASIS_REPRESENTATION_MANAGER_BASE_H

#include "structure_managers/structure_manager_base.hh"
#include "json_io.hh"
#include <string>
#include <vector>


namespace rascal {

  template <class RepresentationImplementation>
  struct RepresentationManager_traits
  {};

  // template <class RepresentationImplementation>
  class RepresentationManagerBase
  {
  public:


    using hypers_t = json;
    
    RepresentationManagerBase() = default;

    //! Copy constructor
    RepresentationManagerBase(const RepresentationManagerBase &other) = delete;

    //! Move constructor
    RepresentationManagerBase(RepresentationManagerBase &&other) = default;

    //! Destructor
    virtual ~RepresentationManagerBase()  = default;

    //! Copy assignment operator
    RepresentationManagerBase& operator=(const RepresentationManagerBase &other) = delete;

    //! Move assignment operator
    RepresentationManagerBase& operator=(RepresentationManagerBase && other) = default;

    // Resolves the mismatch between the expected traits 
    // and the effective traits of the Structure Manager
    template <class Mngr>
    void check_traits_compatibility(Mngr &structure_manager);

    // Pure Virtual Function to set hyperparameters of the representation
    // TODO think of a generic input type for the hypers
    // make class similar to atomic_structure 
    // to handle the hyper io (json format from python and files)
    // virtual void set_hyperparameters(const hypers_t & ) = 0;
        
  

  protected:
  private:
  };

}

#endif /* REPRESENTATION_MANAGER_BASE_H */

