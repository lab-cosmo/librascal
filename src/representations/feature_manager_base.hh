/**
 * file   .cc
 *
 * @author Musil Felix <musil.felix@epfl.ch>
 *
 * @date   14 November 2018
 *
 * @brief  
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


#ifndef FEATURE_MANAGER_BASE_H
#define FEATURE_MANAGER_BASE_H

#include "structure_managers/structure_manager_base.hh"
#include "json_io.hh"

#include <vector>
#include <stdexcept>
#include <iterator>
		

namespace rascal {

class FeatureManagerBase {
  public:
    FeatureManagerBase() = default;

    //! Copy constructor
    FeatureManagerBase(const FeatureManagerBase &other) = delete;

    //! Move constructor
    FeatureManagerBase(FeatureManagerBase &&other) = delete;

    //! Destructor
    ~FeatureManagerBase() = default;

    //! Copy assignment operator
    FeatureManagerBase& operator=(const FeatureManagerBase &other) = delete;

    //! Move assignment operator
    FeatureManagerBase& operator=(FeatureManagerBase && other) = delete;
  
};


} // rascal

#endif /* FEATURE_MANAGER_BASE_H */
