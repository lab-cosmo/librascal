/**
 * file   neighbourhood_manager_lammps.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief Neighbourhood manager for lammps neighbourhood lists
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


#ifndef NEIGHBOURHOOD_MANAGER_LAMMPS_H
#define NEIGHBOURHOOD_MANAGER_LAMMPS_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"

#include <stdexcept>
#include <vector>

namespace rascal {
  //! forward declaration for traits
  class NeighbourhoodManagerLammps;

  //! traits specialisation for Lammps manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerLammps> {
    constexpr static int Dim{3};
    constexpr static size_t MaxLevel{2};
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    using DepthByDimension = std::index_sequence<0, 0>;
  };


  //----------------------------------------------------------------------------//
  class NeighbourhoodManagerLammps: public NeighbourhoodManagerBase<NeighbourhoodManagerLammps>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerLammps>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerLammps>;
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;
    template <size_t Level>
    using ClusterRef_t = typename Parent::template ClusterRef<Level>;

    //! Default constructor
    NeighbourhoodManagerLammps() = default;

    //! Copy constructor
    NeighbourhoodManagerLammps(const NeighbourhoodManagerLammps & other) = delete;

    //! Move constructor
    NeighbourhoodManagerLammps(NeighbourhoodManagerLammps && other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerLammps() = default;

    //! Copy assignment operator
    NeighbourhoodManagerLammps & operator=(const NeighbourhoodManagerLammps & other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerLammps & operator=(NeighbourhoodManagerLammps && other) = default;

    /**
     * resetting is required every time the list changes. Here, this
     * is implemented without explicit dependency to lammmps. The
     * signature could be simplified by including lammps as a
     * dependency, but it is unclear that the convenience would
     * outweigh the hassle of maintaining the dependency.
     *
     * @param inum Property `inum` in the lammps `NeighList` structure
     *
     * @param tot_num sum of the properties `nlocal` and `nghost` in the
     *                lammps `Atom` structure
     *
     * @param ilist Property `ilist` in the lammps `NeighList` structure
     *
     * @param numneigh Property `numneigh` in the lammps `NeighList` structure
     *
     * @param firstneigh Property `firstneigh` in the lammps `NeighList` structure
     *
     * @param x Property `x` in the lammps `Atom` structure
     *
     * @param f Property `f` in the lammps `Atom` structure
     *
     * @param type Property `type` in the lammps `Atom` structure
     *
     * @param eatom per-atom energy
     *
     * @param vatom per-atom virial
     */
    void update(const int & inum, const int & tot_num,
                int * ilist, int * numneigh, int ** firstneigh,
                double ** x, double ** f, int * type,
                double * eatom, double ** vatom);



    // return position vector
    inline Vector_ref get_position(const int & atom_index) {
      auto * xval{this->x[atom_index]};
      return Vector_ref(xval);
    }

    // return position vector
    template<size_t Level, size_t Depth>
    inline Vector_ref get_neighbour_position(const ClusterRefBase<Level, Depth> & cluster) {
      static_assert(Level > 1,
                    "Only possible for Level > 1.");
      static_assert(Level <= traits::MaxLevel,
                    "Level to large, not available.");

      return this->get_position(cluster.back());
    }

    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->inum;
    }

    // return the number of neighbours of a given atom
    template<size_t Level>
    inline size_t get_cluster_size(const ClusterRef<Level> & cluster)
      const {
      static_assert(Level < traits::MaxLevel,
                    "this implementation only handles atoms and pairs");
      return this->numneigh[cluster.back()];
    }

    // return the global id of an atom
    template<size_t Level>
    inline size_t get_atom_id(const ClusterRef_t<Level> & cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.back()};
      return this->firstneigh[std::move(i_atom_id)][j_atom_id];
    }

    // return the global id of an atom
    inline size_t get_atom_id(const Parent & /*cluster*/,
                              size_t i_atom_id) const {
      return (size_t) this->ilist[i_atom_id];
    }

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<size_t Level>
    inline size_t get_offset_impl(const ClusterRefBase<Level,
                               cluster_depth<Level>(traits::DepthByDimension{})>
                               & cluster) const;

    size_t get_nb_clusters(int cluster_size);

  protected:
    int inum{};
    int tot_num{}; //includes ghosts
    int * ilist{};
    int * numneigh{};
    int ** firstneigh{};
    double **x{}; //! pointer to pointer
    double **f{};
    int * type{};
    double * eatom{};
    double ** vatom{};
    int nb_pairs{};
    std::vector<int> offsets{};

  private:
  };


  /* ---------------------------------------------------------------------- */
  // template<size_t Level, size_t Depth>
  // inline int NeighbourhoodManagerLammps::
  // get_offset_impl(const ClusterRefBase<Level, Depth>& cluster) const {
  //   static_assert(Level == 2, "This class cas only handle single atoms and pairs");

  //   auto i{cluster.front()};
  //   auto j{cluster.get_index()};
  //   auto main_offset{this->offsets[i]};
  //   return main_offset + j;
  // }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms // TODO: is not correct anymore
  template <size_t Level>
  inline size_t NeighbourhoodManagerLammps:: template
  get_offset_impl(const ClusterRefBase<Level,
                  cluster_depth<Level>(traits::DepthByDimension{})> & cluster) const {
    return cluster.get_cluster_index(cluster_depth<Level>(traits::DepthByDimension{}));
  }

}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_LAMMPS_H */
