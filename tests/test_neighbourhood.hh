/**
 * file   test_neighbourhood.hh
 *
 * @author Till Junge <till.junge@altermail.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Common headers for tests related to `NeighbourhoodManager`
 *
 * @section LICENSE
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef TEST_NEIGHBOURHOOD_H
#define TEST_NEIGHBOURHOOD_H

#include "neighbourhood_managers/neighbourhood_manager_lammps.hh"

namespace proteus {

  struct ManagerFixture
  {
    using Manager_t = NeighbourhoodManagerLammps;
    constexpr static int nb{3};
    constexpr static int dim{3};
    using ptr_t = double**;
    ManagerFixture()
      :firstneigh{new int*[nb]},
       x{new double*[nb]},
       f{new double*[nb]},
       vatom{new double*[nb]},
       manager{} {
         manager.reset_impl(inum, tot_num, ilist, numneigh,
                            static_cast<int**>(firstneigh),
                            ptr_t(x),
                            ptr_t(f), type, eatom,
                            static_cast<double**>(vatom));
      firstneigh[0] = new int[2];
      firstneigh[0][0] = 1;
      firstneigh[0][1] = 2;
      firstneigh[1] = new int;
      firstneigh[1][0] = 0;
      firstneigh[2] = new int;
      firstneigh[2][0] = 0;

      for (int i{0} ; i < nb ; ++i) {
        x[i] = new double[dim];
        f[i] = new double[dim];
        for (int j{0} ; j < dim ; ++j) {
          x[i][j] = tx[i][j];
          f[i][j] = tf[i][j];
        }
      }
    }
    ~ManagerFixture() {
      delete[] firstneigh[0];
      delete firstneigh[1];
      delete firstneigh[2];
      delete[] firstneigh;
      delete[] vatom;
      for (int i{0} ; i < nb ; ++i) {
        delete[] x[i];
        delete[] f[i];
      }
      delete[] x;
      delete[] f;
    }
    double tx[nb][nb] = {{0,0,0},{1,0,0},{0,1,0}};
    double tf[nb][nb] = {{1,1,0},{-1,0,0},{0,-1,0}};

    int inum{nb};
    int tot_num{nb}; //includes ghosts
    int ilist[nb]{0,1,2};
    int numneigh[nb]{2,1,1};
    int ** firstneigh;
    double **x;
    double **f;
    int type[nb]{1,1,1};
    double  eatom[3]{2,1,1};
    double ** vatom;
    int nb_pairs;
    Manager_t manager;

  };


}  // proteus

#endif /* TEST_NEIGHBOURHOOD_H */
