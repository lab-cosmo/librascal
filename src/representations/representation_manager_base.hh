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



namespace rascal {

template <class ManagerImplementation>
  struct RepresentationManager_traits
  {};

  class RepresentationManagerBase
  {
  public:
    

    //! Default constructor
    RepresentationManagerBase() = delete;

    //! Construct from file
    RepresentationManagerBase(FILE *);

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

    // Initialize hyperparameters
    void read_hyperparameters(const FILE *);
    void set_random_hyperparameters(const int);

    

  protected:
  private:
  }

}

#endif /* BASIS_REPRESENTATION_MANAGER_BASE_H */

