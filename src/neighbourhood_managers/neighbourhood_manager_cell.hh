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


#ifndef NEIGHBOURHOOD_MANAGER_CELL_H
#define NEIGHBOURHOOD_MANAGER_CELL_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <lattice.hh> //! this is a header enabeling a nice for i,j in zip(a1,a2) kind of loops. see for more details https://github.com/cshelton/zipfor
#include <basic_types.h>
#include <neighbourhood_managers/field.hh>

using namespace std;

namespace proteus {
  //! forward declaration for traits
  class NeighbourhoodManagerCell;

  //! traits specialisation for Lammps manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerCell> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{2};
  };
  class NeighbourhoodManagerCell: public NeighbourhoodManagerBase<NeighbourhoodManagerCell>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerCell>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerCell>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;
    using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;
    using AtomVectorField_t = Field<NeighbourhoodManagerCell, double, 1, 3>;
    //! Default constructor
    NeighbourhoodManagerCell() = default;

    //! Copy constructor
    NeighbourhoodManagerCell(const NeighbourhoodManagerCell &other) = delete;

    //! Move constructor
    NeighbourhoodManagerCell(NeighbourhoodManagerCell &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerCell() = default;

    //! Copy assignment operator
    NeighbourhoodManagerCell& operator=(const NeighbourhoodManagerCell &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerCell& operator=(NeighbourhoodManagerCell &&other) = default;

    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}


    // return position vector
    inline Vector_ref get_position(const AtomRef_t& atom) {
      //cout << "get_position" << endl;
      int index{atom.get_index()};
      double * xval{this->positions.col(index).data()};//&this-
      return Vector_ref(xval);
      //return this->positions[index];
    }

    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->centers.size();
    }
    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const {
      return this->centers[i_atom_id].get_index();
    }

    /*
    void set_positions(const Eigen::MatrixXd& pos){
      int dim{this->dim()};
      int Natom{pos.cols()};

      this->positions = new double*[Natom];
      for (int i{0}; i < Natom; ++i){
        this->positions[i] = new double[dim];
        for (int j{0}; j < dim; ++j){
          this->positions[i][j] = pos(j,i);
        }
      }
    }
    */

    void set_positions(const Eigen::Ref<const Eigen::MatrixXd> pos){
      this->positions = pos;
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      auto && ij_atom_id{this->neighlist[i_atom_id][j_atom_id].get_index()};
      return ij_atom_id;
    }


    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level, MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      return this->neighlist[cluster.get_atoms().back().get_index()].size();
    }

    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

    //void build(const Eigen::MatrixXd& positions, const Eigen::MatrixXd& cell,const std::array<bool,3>& pbc, const double& cutoff_max);
    void build(const Eigen::Ref<const Eigen::MatrixXd>, const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);


    void update(const Eigen::Ref<const Eigen::MatrixXd>, const Eigen::Ref<const Eigen::MatrixXd> cell,const std::array<bool,3>& pbc, const double& cutoff_max);


  protected:
    std::vector<AtomRef_t> centers; //! list of center indices. after build it points to neighpos
    std::vector<std::vector<AtomRef_t>> neighlist; //! list of list of neighbour indices. fisrt dimension is in the same order as centers. it points to neighpos
    std::vector<std::array<double,3>> neighpos; //! list of positions for the neighboor list. contains positions of the centers and then the positions of the neighboors.
    Matrix3XdC positions; //! list of pointers. after build it points to neighpos
    Lattice lattice;
    std::array<bool,3> pbc;
  private:
  };

  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerCell::build(const Eigen::Ref<const Eigen::MatrixXd>  positions,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
    {
      Eigen::Index Natom{positions.cols()};
      const int Ncenter{Natom};
      //this->centers.resize(Ncenter);
      this->neighlist.resize(Ncenter);

      cout << "Natom " << Natom << endl;
      const int dim{traits::Dim};
      using Vecd = typename Eigen::Array<double, dim, 1>;
      using Veci = typename Eigen::Array<int, dim, 1>;
      using Mati = typename Eigen::Array<int, dim, 2>;

      cout << "init centers " << endl;
      int count{0};
      this->set_positions(positions);

      for (int id{0} ; id<Ncenter ; ++id) {
          this->centers.push_back(NeighbourhoodManagerCell::AtomRef_t(this->get_manager(),id));
          count += 1;
      }

      Cell_t lat = cell;
      this->lattice.set_cell(lat);
      

      
      //this->set_positions(positions);
      //this->positions = positions.data();

      
    }

  /* ---------------------------------------------------------------------- */

  void NeighbourhoodManagerCell::update(const Eigen::Ref<const Eigen::MatrixXd>,
                                        const Eigen::Ref<const Eigen::MatrixXd> cell,
                                        const std::array<bool,3>& pbc, const double& cutoff_max)
  {
    bool some_condition{false};
    if (some_condition){
      NeighbourhoodManagerCell::build(positions,cell,pbc,cutoff_max);
    }
  }
  /* ---------------------------------------------------------------------- */

  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerCell::
  get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 2, "This class cas only handle single atoms and pairs");
    static_assert(MaxLevel == traits::MaxLevel, "Wrong maxlevel");

    auto atoms{cluster.get_atoms()};
    auto i{atoms.front().get_index()};
    auto j{cluster.get_index()};
    auto main_offset{this->neighlist[i][j]};
    return main_offset;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms

  template <>
  inline int NeighbourhoodManagerCell:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }


  //----------------------------------------------------------------------------//
  // check in muSpectre src/common/ccoord_operations.hh (class Pixels<dim>)
    //! get the i-th pixel in a grid of size sizes
    /*
    template <size_t dim>
    constexpr std::array<int, dim> get_ccoord(const std::array<int, dim> & resolutions,
                                              Dim_t index) {
      std::array<int, dim> retval{{0}};
      int factor{1};
      for (int i = dim-1; i >=0; --i) {
        retval[i] = index/factor%resolutions[i];
        if (i != 0 ) {
          factor *= resolutions[i];
        }
      }
      return retval;
    }
    */

}  // proteus

#endif /* NEIGHBOURHOOD_MANAGER_CELL_H */
