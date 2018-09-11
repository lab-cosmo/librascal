/**
 * file   test_structure.hh
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

#ifndef TEST_NEIGHBOURHOOD_H
#define TEST_NEIGHBOURHOOD_H

#include "structure_managers/structure_manager_base.hh"
#include "structure_managers/structure_manager_lammps.hh"
#include "structure_managers/structure_manager_cell.hh"
#include "structure_managers/structure_manager_chain.hh"
#include "structure_managers/structure_manager_json.hh"
#include "structure_managers/structure_manager_centers.hh"

namespace rascal {

  // TODO: this is not a general case of a manager fixture. Should not be so
  // complicated
  template<class ManagerImplementation>
  struct ManagerFixture
  {
    ManagerFixture():
      pbc{{true,true,true}}, cutoff_max{3}, center_ids(22),
      cell(3, 3), positions(3, 22), numbers(22)
    {
      cell <<
	6.19, 2.41, 0.21,
        0.00, 6.15, 1.02,
        0.00, 0.00, 7.31;
      positions <<
        3.689540159937393, 5.123016813620886, 1.994119731169116,
        6.818437242389163, 2.630056617829216, 6.182500355729062,
        2.114977334498767, 6.697579639059512, 1.392155450018263,
        7.420401523540017, 2.432242071439904, 6.380314902118375,
        1.112656394115962, 7.699900579442317, 3.569715877854675,
        5.242841095703604, 3.122826344932127, 5.689730628626151,
        3.248684682453303, 5.563872291104976, 2.608353462112637,
        6.204203511445642, 5.035681855581504, 2.134827911489532,
        0.946910011088814, 6.223599755982222, 4.168634519120968,
        3.001875247950068, 1.980327734683430, 5.190182032387606,
        2.943861424421339, 4.226648342649697, 5.457161501166098,
        1.713348265904937, 1.501663178733906, 5.668846588337130,
        5.208365510425203, 1.962144256645833, 2.728127406527150,
        4.442382360543885, 2.839975217222644, 4.330534549848392,
        0.744216089807768, 6.426293677263268, 4.643695520786083,
        2.662204050783991, 1.250682335857938, 6.055217235712136,
        0.860905287815103, 6.444994283754972, 4.536108843695142,
        2.769790727874932, 5.609177455068640, 1.696722116501434,
        6.703053268421970, 0.602846303148105, 3.487609972580834,
        3.818289598989240, 1.436734374347541, 5.869165197222533,
        1.054504320562138, 6.251395251007936, 3.998423858825871,
        3.307475712744203, 5.323662899811682, 1.982236671758393;
      numbers << 20, 20, 24, 24, 15, 15, 15, 15,  8,  8,  8,
        8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8;
      center_ids << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15, 16, 17, 18, 19, 20, 21;
      manager.update(positions,numbers,center_ids,cell,pbc,cutoff_max);

    }

    ~ManagerFixture() {

    }

    ManagerImplementation manager{};
    std::array<bool, 3> pbc;
    double cutoff_max;
    VecXi center_ids;
    Eigen::MatrixXd cell;
    Eigen::MatrixXd positions; // 3, 22
    VecXi numbers;
  };

  /* ---------------------------------------------------------------------- */
  template<class ManagerImplementation>
  struct ManagerFixtureSimple
  {
    ManagerFixtureSimple():
      pbc{{true,true,true}}, cutoff{1.}, center_ids(natoms),
      cell(dim, dim), positions(dim, natoms), atom_types(natoms)
    {}

    ~ManagerFixtureSimple() {}

    ManagerImplementation manager{};
    std::array<bool, 3> pbc;
    double cutoff;
    VecXi center_ids;
    Eigen::MatrixXd cell;
    Eigen::MatrixXd positions;
    VecXi atom_types;
    int natoms;
    int dim;
  };

  /* ---------------------------------------------------------------------- */
  template <>
  struct ManagerFixture<StructureManagerLammps>
  {
    using Manager_t = StructureManagerLammps;
    constexpr static int nb{3};
    constexpr static int dim{3};
    using ptr_t = double**;

    ManagerFixture()
      :firstneigh{new int*[nb]},
       x{new double*[nb]},
       f{new double*[nb]},
       vatom{new double*[nb]},
       manager{} {
         manager.update(inum, tot_num, ilist, numneigh,
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

    ManagerFixture( ManagerFixture &) = delete;
    ManagerFixture & operator=(const ManagerFixture&) = delete;
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
    double tx[nb][nb] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    double tf[nb][nb] = {{1, 1, 0}, {-1, 0, 0}, {0, -1, 0}};

    int inum{nb};
    int tot_num{nb}; //includes ghosts
    int ilist[nb]{0, 1, 2};
    int numneigh[nb]{2, 1, 1};
    int ** firstneigh;
    double **x;
    double **f;
    int type[nb]{1, 1, 1};
    double  eatom[3]{2, 1, 1};
    double ** vatom;
    Manager_t manager;

  };

  /* ---------------------------------------------------------------------- */
  template <>
  struct ManagerFixture<StructureManagerChain>
  {
    using Manager_t = StructureManagerChain;

    ManagerFixture()
      : manager_chain{}, cutoff{1.0} {
        manager_chain.read_structure_from_json("simple_cubic_8.json");
        manager_chain.update(cutoff);
      }

    ~ManagerFixture() {BOOST_TEST_MESSAGE("teardown ManagerChain fixture");}

    Manager_t manager_chain;
    double cutoff;
  };

  /* ---------------------------------------------------------------------- */
  template <>
  struct ManagerFixture<StructureManagerJson>
  {
    using Manager_t = StructureManagerJson;
    ManagerFixture()
      : manager_json{}, cutoff{1.1} {
      manager_json.read_structure_from_json("simple_cubic_8.json");
      manager_json.update();
    }

    ~ManagerFixture () {BOOST_TEST_MESSAGE("teardown ManagerJson fixture");}

    Manager_t manager_json;
    double cutoff;
  };

  /* ---------------------------------------------------------------------- */
  template<>
  struct ManagerFixture<StructureManagerCenters>
  {

    using Manager_t = StructureManagerCenters;

    ManagerFixture():
      pbc{{true,true,true}}, cell(3, 3), positions(3, 22), numbers(22),
      cutoff{3.}
    {
      cell <<
        6.19, 2.41, 0.21,
        0.00, 6.15, 1.02,
        0.00, 0.00, 7.31;

      positions <<
        3.689540159937393, 5.123016813620886, 1.994119731169116,
        6.818437242389163, 2.630056617829216, 6.182500355729062,
        2.114977334498767, 6.697579639059512, 1.392155450018263,
        7.420401523540017, 2.432242071439904, 6.380314902118375,
        1.112656394115962, 7.699900579442317, 3.569715877854675,
        5.242841095703604, 3.122826344932127, 5.689730628626151,
        3.248684682453303, 5.563872291104976, 2.608353462112637,
        6.204203511445642, 5.035681855581504, 2.134827911489532,
        0.946910011088814, 6.223599755982222, 4.168634519120968,
        3.001875247950068, 1.980327734683430, 5.190182032387606,
        2.943861424421339, 4.226648342649697, 5.457161501166098,
        1.713348265904937, 1.501663178733906, 5.668846588337130,
        5.208365510425203, 1.962144256645833, 2.728127406527150,
        4.442382360543885, 2.839975217222644, 4.330534549848392,
        0.744216089807768, 6.426293677263268, 4.643695520786083,
        2.662204050783991, 1.250682335857938, 6.055217235712136,
        0.860905287815103, 6.444994283754972, 4.536108843695142,
        2.769790727874932, 5.609177455068640, 1.696722116501434,
        6.703053268421970, 0.602846303148105, 3.487609972580834,
        3.818289598989240, 1.436734374347541, 5.869165197222533,
        1.054504320562138, 6.251395251007936, 3.998423858825871,
        3.307475712744203, 5.323662899811682, 1.982236671758393;

      numbers << 20, 20, 24, 24, 15, 15, 15, 15, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8;

      manager.update(positions, numbers, cell, pbc);
    }

    ~ManagerFixture() {

    }

    Manager_t manager{};
    std::array<bool, 3> pbc;
    Eigen::MatrixXd cell;
    Eigen::MatrixXd positions; // 3, 22
    VecXi numbers;

    double cutoff;

    int natoms{22};
  };

  /* ---------------------------------------------------------------------- */
  template<>
  struct ManagerFixtureSimple<StructureManagerCenters>
  {

    using Manager_t = StructureManagerCenters;

    ManagerFixtureSimple():
      pbc{{true,false,false}}, cell(3, 3), positions(3, 9), numbers(9),
      cutoff{1.5}
    {
      cell <<
        2., 0., 0.,
        0., 2., 0.,
        0., 0., 2.;

      positions <<
        0.5, 0.5, 1.5, 0.5, 1.5, 0.5, 1.5, 0.5, 1.5,
        0.4, 0.5, 0.5, 1.5, 1.5, 0.5, 0.5, 1.5, 1.5,
        0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.8;

      numbers << 1, 1, 1, 1, 1, 1, 1, 1, 1;

      manager.update(positions, numbers, cell, pbc);
    }

    ~ManagerFixtureSimple() {

    }

    Manager_t manager{};
    std::array<bool, 3> pbc;
    Eigen::MatrixXd cell;
    Eigen::MatrixXd positions;
    VecXi numbers;

    double cutoff;

    const int natoms{9};
  };

}  // rascal

#endif /* TEST_NEIGHBOURHOOD_H */
