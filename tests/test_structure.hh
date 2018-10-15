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
#include "structure_managers/adaptor_strict.hh"
#include "structure_managers/adaptor_increase_maxorder.hh"
#include "structure_managers/adaptor_neighbour_list.hh"

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
  template<class ManagerImplementation>
  struct ManagerFixtureNeighbourComparison
  {
    ManagerFixtureNeighbourComparison():
      pbc{{true,true,true}}, cutoff{1.}, center_ids(natoms),
      cell_1(dim, dim), cell_2(dim, dim),
      positions_1(dim, natoms), positions_2(dim, natoms), atom_types(natoms)
    {}

    ~ManagerFixtureNeighbourComparison() {}

    ManagerImplementation manager_1{};
    ManagerImplementation manager_2{};
    std::array<bool, 3> pbc;
    double cutoff;
    VecXi center_ids;
    Eigen::MatrixXd cell_1;
    Eigen::MatrixXd cell_2;
    Eigen::MatrixXd positions_1;
    Eigen::MatrixXd positions_2;
    VecXi atom_types;
    int natoms;
    int dim;
  };

  /* ---------------------------------------------------------------------- */
  template<class ManagerImplementation>
  struct ManagerFixtureNeighbourCheckFcc
  {
    ManagerFixtureNeighbourCheckFcc():
      pbc{{true,true,true}}, cutoff{1.},
      center_ids_1(natoms_1), center_ids_2(natoms_2),
      cell_1(dim, dim), cell_2(dim, dim),
      positions_1(dim, natoms_1), positions_2(dim, natoms_1),
      atom_types_1(natoms_2), atom_types_2(natoms_2)
    {}

    ~ManagerFixtureNeighbourCheckFcc() {}

    ManagerImplementation manager_1{};
    ManagerImplementation manager_2{};
    std::array<bool, 3> pbc;
    double cutoff;
    VecXi center_ids_1;
    VecXi center_ids_2;
    Eigen::MatrixXd cell_1;
    Eigen::MatrixXd cell_2;
    Eigen::MatrixXd positions_1;
    Eigen::MatrixXd positions_2;
    VecXi atom_types_1;
    VecXi atom_types_2;
    int natoms_1;
    int natoms_2;
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
      positions(22, 3), numbers(22), cell(3, 3), pbc{{true,true,true}},
      // positions(3,8), numbers(8), cell(3, 3), pbc{{true,true,true}},
      cutoff{2.}
    {
      cell <<
        6.19, 2.41, 0.21,
        0.00, 6.15, 1.02,
        0.00, 0.00, 7.31;
      // cell.transposeInPlace();
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

      positions.transposeInPlace();
      numbers << 20, 20, 24, 24, 15, 15, 15, 15, 8, 8, 8,
        8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8;
      // cell <<
      //   2., 0., 0.,
      //   0., 2., 0.,
      //   0., 0., 2.;

      // positions <<
      //   0.4, 1.4, 0.4, 1.4, 0.4, 1.4, 0.4, 1.4,
      //   0.4, 0.4, 1.4, 1.4, 0.4, 0.4, 1.4, 1.4,
      //   0.4, 0.4, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4;

      // numbers << 1, 2, 3, 4, 5, 6, 7, 8;

      manager.update(positions, numbers, cell,
                     Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});
    }

    ~ManagerFixture() {

    }

    Manager_t manager{};
    Eigen::MatrixXd positions;
    VecXi numbers;
    Eigen::MatrixXd cell;
    std::array<int, 3> pbc;
    double cutoff;

    //int natoms{22};
  };

  /* ---------------------------------------------------------------------- */
  /**
   * Comparison of two cells to check if the zeroth level the AdaptorMaxOrder
   * (building the neighbourlist) works properly
   */
  template<>
  struct ManagerFixtureNeighbourComparison<StructureManagerCenters>
  {
    using Manager_t = StructureManagerCenters;

    ManagerFixtureNeighbourComparison():
      pbc{{true, true, true}}, cell_1(3, 3), cell_2(3, 3),
      positions_1(3, 2), positions_2(3, 2), numbers(2), cutoff{0.7}
    {
      /**
       * hcp crystal with lattice parameter a = 1, c = sqrt(8/3), defined in two
       * unit cells: basal and prismatic 1. The neighbourlist is built with the
       * same cutoff. The test checks, if all atoms have the same number of
       * neighbours.
       */
      auto a{1.};
      auto c{std::sqrt(8./3.)};

      cell_1 <<
        a,  -0.5*a ,            0.,
        0., std::sqrt(3.)/2.*a, 0.,
        0.,  0.,                c;

      cell_2 <<
        a,   0.,         0.5*a,
        0.,  c,             0.,
        0.,  0.,  std::sqrt(3.)/2.*a;

      auto p_1 = 2./3. * cell_1.col(0)
        + 1./3. * cell_1.col(1)
        + 1./2. * cell_1.col(2);

      positions_1 <<
        0.0, p_1[0],
        0.0, p_1[1],
        0.0, p_1[2];

      auto p_2 = -1./3. * cell_2.col(0)
        + 1./2. * cell_2.col(1)
        + 2./3. * cell_2.col(2);

      positions_2 <<
        0.0, p_2[0],
        0.0, p_2[1],
        0.0, p_2[2];

      numbers << 1, 1;

      manager_1.update(positions_1, numbers, cell_1,
                       Eigen::Map<Eigen::Matrix<int, 3, 1>>
                       {pbc.data()});

      manager_2.update(positions_2, numbers, cell_2,
                       Eigen::Map<Eigen::Matrix<int, 3, 1>>
                       {pbc.data()});
    }
    Manager_t manager_1{};
    Manager_t manager_2{};
    std::array<int, 3> pbc;
    Eigen::MatrixXd cell_1;
    Eigen::MatrixXd cell_2;
    Eigen::MatrixXd positions_1;
    Eigen::MatrixXd positions_2;
    VecXi numbers;

    double cutoff;

    const int natoms{2};

  };

  /* ---------------------------------------------------------------------- */
  /**
   * Comparison of two fcc cells to check if the zeroth level the AdaptorMaxOrder
   * (building the neighbourlist) works properly
   */
  template<>
  struct ManagerFixtureNeighbourCheckFcc<StructureManagerCenters>
  {
    using Manager_t = StructureManagerCenters;

    ManagerFixtureNeighbourCheckFcc():
      pbc{{true, true, true}},
      cell_1(3, 3), cell_2(3, 3),
      positions_1(3, 1), positions_2(3, 4),
      numbers_1(1), numbers_2(4),
      cutoff{0.7}, // start with zero neighbours
      natoms_1{1}, natoms_2{4}
    {
      /**
       * fcc unit cells: first cell consists of only one atom, which is
       * repeated, second cell is the conventional 4 atoms. This test checks, if
       * the found number of neighbours with increasing cutoff is the same for
       * the atom at position (0, 0, 0).
       */
      auto a{1.};

      cell_1 <<
        a,  0.5*a, 0.5*a,
        0., 0.5*a, 0.   ,
        0., 0.,    0.5*a;

      cell_2 <<
        a,   0., 0.,
        0.,  a,  0.,
        0.,  0., a ;

      positions_1 <<
        0.,
        0.,
        0.;

      auto p_2 = 0.5 * cell_2.col(0) + 0.5 * cell_2.col(1);
      auto p_3 = 0.5 * cell_2.col(0) + 0.5 * cell_2.col(2);
      auto p_4 = 0.5 * cell_2.col(1) + 0.5 * cell_2.col(2);

      positions_2 <<
        0.0, p_2[0], p_3[0], p_4[0],
        0.0, p_2[1], p_3[1], p_4[1],
        0.0, p_2[2], p_3[2], p_4[2];

      numbers_1 << 1;
      numbers_2 << 1, 1, 1, 1;

      manager_1.update(positions_1, numbers_1, cell_1,
                       Eigen::Map<Eigen::Matrix<int, 3, 1>>
                       {pbc.data()});

      manager_2.update(positions_2, numbers_1, cell_2,
                       Eigen::Map<Eigen::Matrix<int, 3, 1>>
                       {pbc.data()});
    }
    Manager_t manager_1{};
    Manager_t manager_2{};
    std::array<int, 3> pbc;
    Eigen::MatrixXd cell_1;
    Eigen::MatrixXd cell_2;
    Eigen::MatrixXd positions_1;
    Eigen::MatrixXd positions_2;
    VecXi numbers_1;
    VecXi numbers_2;

    double cutoff;

    const int natoms_1{1};
    const int natoms_2{4};

  };

  /* ---------------------------------------------------------------------- */
  /**
   * A simple manager using ManagerCenters to check the neighbourlist algorithm
   * with simple positions and a periodicity only in x-direction.
   *
   */
  template<>
  struct ManagerFixtureSimple<StructureManagerCenters>
  {

    using Manager_t = StructureManagerCenters;

    ManagerFixtureSimple():
      pbc{{true, false, false}}, cell(3, 3), positions(3, 8), numbers(8),
      cutoff{2.1}
    {
      cell <<
        2., 0., 0.,
        0., 2., 0.,
        0., 0., 2.;

      positions <<
        0.4, 1.4, 0.4, 1.4, 0.4, 1.4, 0.4, 1.4,
        0.4, 0.4, 1.4, 1.4, 0.4, 0.4, 1.4, 1.4,
        0.4, 0.4, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4;

      numbers << 1, 1, 1, 1, 1, 1, 1, 1;

      manager.update(positions, numbers, cell,
                     Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});
    }

    ~ManagerFixtureSimple() {}

    Manager_t manager{};
    std::array<int, 3> pbc;
    Eigen::MatrixXd cell;
    Eigen::MatrixXd positions;
    VecXi numbers;

    double cutoff;

    const int natoms{8};
  };

  // template<>
  // struct ManagerFixtureSimple<StructureManagerCenters>
  // {

  //   using Manager_t = StructureManagerCenters;

  //   ManagerFixtureSimple():
  //     pbc{{true,false,false}}, cell(3, 3), positions(3, 2), numbers(2),
  //     cutoff{5.0}
  //   {
  //     cell <<
  //       1.99, 0.99, 0.,
  //       0.99, 1.99, 0.,
  //       0., 0., 1.99;

  //     positions <<
  //       0.5, 1.5,
  //       0.5, 1.5,
  //       0.5, 0.5;

  //     numbers << 1, 1;//, 1, 1, 1, 1, 1, 1;

  //     manager.update(positions, numbers, cell,
  //                    Eigen::Map<Eigen::Matrix<int, 3, 1>>{pbc.data()});
  //   }

  //   ~ManagerFixtureSimple() {}

  //   Manager_t manager{};
  //   std::array<int, 3> pbc;
  //   Eigen::MatrixXd cell;
  //   Eigen::MatrixXd positions;
  //   VecXi numbers;

  //   double cutoff;

  //   const int natoms{2};
  // };

}  // rascal

#endif /* TEST_NEIGHBOURHOOD_H */
