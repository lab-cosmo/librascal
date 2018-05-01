/**
 * file   neighbourhood_manager_cell.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef NEIGHBOURHOOD_BOX_H
#define NEIGHBOURHOOD_BOX_H

#include "neighbourhood_managers/neighbourhood_manager_cell.hh"
#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <lattice.hh> //! this is a header enabeling a nice for i,j in zip(a1,a2) kind of loops. see for more details https://github.com/cshelton/zipfor
#include <basic_types.h>
#include <neighbourhood_managers/field.hh>



namespace proteus {

  

} //proteus

#endif /* NEIGHBOURHOOD_BOX_H */