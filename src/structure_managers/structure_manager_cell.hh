/**
 * file   structure_manager_cell.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager linked cell
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


#ifndef STRUCTURE_MANAGER_CELL_H
#define STRUCTURE_MANAGER_CELL_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "lattice.hh"
#include "basic_types.hh"

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <stdexcept>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

namespace rascal {
  namespace internal {

    // conversion of a linear id to a Dim_t id
    template<int Dim>
    inline void lin2mult(const Dim_t& index,
                         const Eigen::Ref<const Vec3i_t> shape,
                         Eigen::Ref<Vec3i_t> retval) {
      // TODO: what does the factor do?
      Dim_t factor{1};

      for (Dim_t i{0}; i < Dim; ++i) {
        retval[i] = index / factor%shape[i];
        if (i != Dim-1) {
          factor *= shape[i];
        }
      }
    }

    template<int Dim>
    inline Dim_t mult2lin(const Eigen::Ref<const Vec3i_t> coord,
                          const Eigen::Ref<const Vec3i_t> shape) {
      Dim_t index{0};
      Dim_t factor{1};
      for (Dim_t i = 0; i < Dim; ++i) {
        index += coord[i]*factor;
        if (i != Dim-1 ) {
          factor *= shape[i];
        }
      }
      return index;
    }

    // https://stackoverflow.com/questions/828092/python-style-integer-division-modulus-in-c
    // TODO more efficient implementation without if (would be less general) !
    // div_mod function returning python like div_mod, i.e. signed integer
    // division truncates towards negative infinity, and signed integer modulus
    // has the same sign the second operand.
    // Integral is just a name, it does not mean an actual integration
    template<class Integral>
    void div_mod(const Integral& x, const Integral & y,
                 std::array<int, 2> & out) {
      const Integral quot = x/y;
      const Integral rem  = x%y;

      if (rem != 0 && (x<0) != (y<0)) {
        out[0] = quot-1;
        out[1] = rem+y;
      } else {
        out[0] = quot;
        out[1] = rem;
      }
    }

    //! warning: works only with negative x if |x| < y
    // TODO Does not pass tests
    template<class Integral>
    void branchless_div_mod(const Integral & x, const Integral & y,
                            std::array<int, 2> & out) {
      const Integral quot = (x-y) / y;
      const Integral rem  = (x+y) % y;
      out[0] = quot;
      out[1] = rem;
    }

    // maybe get_box_index from _chain.cc

  }  // internal

  //! forward declaration for traits
  class StructureManagerCell;

  //! traits specialisation for Cell manager
  template <>
  struct StructureManager_traits<StructureManagerCell> {
    constexpr static int Dim{3};
    constexpr static size_t MaxOrder{2};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDirectionVectors{false};
    constexpr static bool HasDistances{false};
    using LayerByOrder = std::index_sequence<0, 0>;
  };
  class StructureManagerCell:
    public StructureManager<StructureManagerCell>
  {
  public:
    using traits = StructureManager_traits<StructureManagerCell>;
    using Parent = StructureManager<StructureManagerCell>;
    using Vector_ref = typename Parent::Vector_ref;
    using Vector_t = typename Parent::Vector_t;
    using AtomRef_t = typename Parent::AtomRef;
    template <size_t Order>
    using ClusterRef_t = typename Parent::template ClusterRef<Order>;
    using AtomVectorField_t =
      typename StructureManagerCell::Property_t<double, 1, 3>;

    //! access to periodic boundary conditions
    using PBC_ref = Eigen::Map<Eigen::Matrix<int, 1, traits::Dim>>;

    //! Default constructor
    StructureManagerCell()
      : particles{}, centers{}, positions{}, shifted_position{}, lattice{},
        cell{}, pbc{}, part2bin{}, boxes{}, number_of_neighbours{0},
        neighbour_bin_id{}, number_of_neighbours_stride{},
        neighbour_atom_index{}, particle_types{} {}

    //! Copy constructor
    StructureManagerCell(const StructureManagerCell &other) = delete;

    //! Move constructor
    StructureManagerCell(StructureManagerCell &&other) = default;

    //! Destructor
    virtual ~StructureManagerCell() = default;

    //! Copy assignment operator
    StructureManagerCell
    & operator=(const StructureManagerCell &other) = delete;

    //! Move assignment operator
    StructureManagerCell
    & operator=(StructureManagerCell &&other) = default;

    class Box;

    // return position vector for atom
    inline Vector_ref get_position(const AtomRef_t & atom) {
      auto index{atom.get_index()};
      auto * xval{this->positions.col(index).data()};
      return Vector_ref(xval);
    }

    // return position vector for atom_index
    inline Vector_ref get_position(const size_t & atom_index) {
      auto * xval{this->positions.col(atom_index).data()};
      return Vector_ref(xval);
    }

    //! Returns a map of size traits::Dim with 0/1 for periodicity
    inline decltype(auto) get_periodic_boundary_conditions() {
      return this->pbc;
    }

    inline Vector_ref get_shift(const int& i_bin_id, const int& shift_index);

    // return position vector atom is the neighbour atom. center_atom is the
    // current center. j_linear_id is the index of the current neighbour
    // iterator.
    template<size_t Order, size_t Layer>
    inline Vector_t get_neighbour_position(const ClusterRefKey<Order, Layer> &
                                           cluster) {
      static_assert(Order > 1,
                    "Only possible for Order > 1.");
      static_assert(Order <= traits::MaxOrder,
                    "this implementation should only work up to MaxOrder.");


      // TODO: why is there a j_linear_id and a j_atom_id, which is the same?
      // TODO: this is the wrong index -> Repair by Felix?
      auto && j_linear_id = cluster.back();
      auto & i_atom_id{cluster.front()}; // center_atom index
      auto & i_bin_id{this->part2bin[i_atom_id]};
      auto & shift_index{this->neighbour_bin_id[i_bin_id][j_linear_id].get_index()};
      auto & j_atom_id{cluster.back()}; // neighbour atom index
      auto && shifted_position{this->positions.col(j_atom_id)};
      auto && tmp{this->cell * this->get_shift(i_bin_id,shift_index)};
      return shifted_position + tmp;
    }

    // return number of center in the list
    inline size_t get_size() const {
      return this->centers.size();
    }

    // return the index-th neighbour of cluster
    template<size_t Order, size_t Layer>
    inline int get_cluster_neighbour(const ClusterRefKey<Order, Layer>
                                     & cluster,
                                     size_t index) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.back()};
      auto && i_bin_id{this->part2bin[i_atom_id]};
      auto && ij_atom_id{this->neighbour_atom_index[i_bin_id][index].get_index()};
      return ij_atom_id;
    }


    // return the atom_index of the index-th atom in manager
    inline int get_cluster_neighbour(const Parent & /*cluster */ ,
                                     size_t index) const {
      return this->centers[index].get_index();
    }

    //! return atom type
    inline size_t get_atom_type(const AtomRef_t& atom) {
      auto && index{atom.get_index()};
      return this->particle_types[index];
    }

    //! return atom type from atom index
    inline size_t get_atom_type(const size_t& index) {
      return this->particle_types[index];
    }

    // return the number of neighbours of a given atom
    template<size_t Order, size_t Layer>
    inline size_t get_cluster_size(const ClusterRefKey<Order, Layer>
                                   & cluster) const {
      static_assert(Order <= traits::MaxOrder,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.back()};
      auto && box_id{this->part2bin[i_atom_id]};
      auto && size{this->neighbour_atom_index[box_id].size()};
      return size;
    }

    template<size_t Order>
    inline size_t get_offset_impl(const std::array<size_t, Order>
                                  & counters) const;

    size_t get_nb_clusters(size_t cluster_size) const;

    void update(const Eigen::Ref<const Eigen::MatrixXd> positions,
                const Eigen::Ref<const Eigen::VectorXi>  particle_types,
                const Eigen::Ref<const Eigen::VectorXi> center_ids,
                const Eigen::Ref<const Eigen::MatrixXd> cell,
                const std::array<bool,3>& pbc,
                const double & cutoff_max);

    //Box get_box(const int& bin_id);

    size_t get_box_nb();

  protected:

    void build(const Eigen::Ref<const Eigen::MatrixXd> positions,
               const Eigen::Ref<const Eigen::VectorXi>  particle_types,
               const Eigen::Ref<const Eigen::VectorXi> center_ids,
               const Eigen::Ref<const Eigen::MatrixXd> cell,
               const std::array<bool,3>& pbc,
               const double& cutoff_max);

    void set_positions(const Eigen::Ref<const Eigen::MatrixXd> pos){
      this->positions = pos;
    }

    std::vector<AtomRef_t> particles;
    std::vector<AtomRef_t> centers; //!
    Matrix3XdC positions; //!
    Vector_t shifted_position;
    Lattice lattice;
    Cell_t cell; // to simplify get_neighbour_position()
    std::array<int, 3> pbc;
    std::vector<int> part2bin; //!
    std::vector<Box> boxes;
    size_t number_of_neighbours;
    std::vector<std::vector<AtomRef_t>> neighbour_bin_id;
    std::vector<size_t> number_of_neighbours_stride;
    std::vector<std::vector<AtomRef_t>> neighbour_atom_index;
    std::vector<int> particle_types;

  private:
  };

  /* ---------------------------------------------------------------------- */

  class StructureManagerCell::Box {
  public:
    using Manager_t = StructureManager<StructureManagerCell>;
    using AtomRef_t = typename StructureManagerCell::AtomRef_t;
    using Vector_t = typename StructureManagerCell::Vector_t;
    using Vector_ref = typename StructureManagerCell::Vector_ref;
    //! Default constructor
    Box() = default;

    //! constructor
    Box(Manager_t& manager,
        const Vec3i_t& coord,
        const std::array<bool, 3>& pbc,
        const Vec3i_t& neigh_search,
        const Vec3i_t& nbins_c);
    //const std::array<std::array<Dim_t, 3>,2>& neigh_bounds,


    //! copy constructor
    Box(const Box & other) = default;
    //! assignment operator
    Box & operator=(const Box & other) = default;

    virtual ~Box() = default;

    constexpr static int dim() {return StructureManagerCell::dim();}

    inline void push_particle_back(const int& part_index);

    inline size_t get_number_of_particles();

    inline size_t get_number_of_neighbours();

    inline size_t get_number_of_neighbour_box();

    inline void set_number_of_neighbours(const size_t& neigh_nb);

    inline  int get_neighbour_bin_index(const int& j_index);

    inline size_t get_particle_index(const int& index);

    inline Vector_ref get_neighbour_bin_shift(const int& neigh_bin_index){
      auto * xval{this->neighbour_bin_shift[neigh_bin_index].col(0).data()};
      return Vector_ref(xval);
    }

  protected:
    Manager_t & manager;
    std::vector<AtomRef_t> particles;
    //stores double for the dot product with the Cell vectors
    std::vector<Vector_t,Eigen::aligned_allocator<Vector_t>> neighbour_bin_shift;
    std::vector<int> neighbour_bin_index;
    size_t number_of_neighbours;
    Vec3i_t coordinates;
  };

  /* ---------------------------------------------------------------------- */
  // buildup
  template<size_t Order>
  inline size_t StructureManagerCell::
  get_offset_impl(const std::array<size_t, Order> & counters) const {
    static_assert (Order == 1, "this manager can only give the offset "
                   "(= starting index) for a pair iterator, given the i atom "
                   "of the pair");
    return this->number_of_neighbours_stride[counters.front()];
  }

  /* ---------------------------------------------------------------------- */
  inline Vector_ref StructureManagerCell::
  get_shift(const int& i_bin_id, const int& neigh_bin_index){
    return this->boxes[i_bin_id].get_neighbour_bin_shift(neigh_bin_index);
  }

}  // rascal

#endif /* STRUCTURE_MANAGER_CELL_H */
