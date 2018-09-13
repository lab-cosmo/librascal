/**
 * file   test_lattice.hh
 *
 * @author  Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Test implementation of lattice.cc
 *
 * Copyright © 2018  Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TEST_LATTICE_H
#define TEST_LATTICE_H

#include "lattice.hh"
#include "basic_types.hh"


namespace rascal {

  struct ManagerFixture_lattice
  {
    ManagerFixture_lattice(){
      Cell_t cell;
      cell << 6.19,2.41,0.21,
              0.00,6.15,1.02,
              0.00,0.00,7.31;
      lattice.set_cell(cell);

    }



    ~ManagerFixture_lattice() {

    }

    Lattice lattice{};

  };
}  // rascal

#endif /* TEST_LATTICE_H */
