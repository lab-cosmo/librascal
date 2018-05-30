/**
 * file   neighbourhood_manager_chain.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   30 May 2018
 *
 * @brief Neighbourhood manager for polyalanine chain, reading
 *        structure from json file
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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


#ifndef NEIGHBOURHOOD_MANAGER_CHAIN_H
#define NEIGHBOURHOOD_MANAGER_CHAIN_H

#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"
#include <stdexcept>
#include <vector>

namespace rascal {
  //! forward declaration for traits
  class NeighbourhoodManagerChain;

  //! traits specialisation for Chain manager
  template <>
  struct NeighbourhoodManager_traits<NeighbourhoodManagerChain> {
    constexpr static int Dim{3};
    constexpr static int MaxLevel{2}; // only triplets needed for angle
  };
  class NeighbourhoodManagerChain:
    public NeighbourhoodManagerBase<NeighbourhoodManagerChain>
  {
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManagerChain>;
    using Parent = NeighbourhoodManagerBase<NeighbourhoodManagerChain>;
    using Vector_ref = typename Parent::Vector_ref;
    using AtomRef_t = typename Parent::AtomRef;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template ClusterRef<Level, MaxLevel>;

    //! Default constructor
    NeighbourhoodManagerChain() = default;

    //! Copy constructor
    NeighbourhoodManagerChain(const NeighbourhoodManagerChain &other) = delete;

    //! Move constructor
    NeighbourhoodManagerChain(NeighbourhoodManagerChain &&other) = default;

    //! Destructor
    virtual ~NeighbourhoodManagerChain() = default;

    //! Copy assignment operator
    NeighbourhoodManagerChain& operator=(const NeighbourhoodManagerChain &other) = delete;

    //! Move assignment operator
    NeighbourhoodManagerChain& operator=(NeighbourhoodManagerChain &&other) = default;

    //! something like this!
    void reset_impl()
    // void reset_impl(const int & inum, const int & tot_num,
    //                 int * ilist, int * numneigh, int ** firstneigh,
    //                 double ** x, double ** f, int * type,
    //                 double * eatom, double ** vatom);
    // required for the construction of vectors, etc
    constexpr static int dim() {return traits::Dim;}


    // return position vector
    inline Vector_ref get_position(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      auto * xval{this->x[index]};
      return Vector_ref(xval);
    }

    // return force vector
    inline Vector_ref get_force(const AtomRef_t& atom) {
      return Vector_ref(this->f[atom.get_index()]);
    }

    // return position vector
    inline int get_atom_type(const AtomRef_t& atom) {
      auto index{atom.get_index()};
      return this->type[index];
    }


    // return number of I atoms in the list
    inline size_t get_size() const {
      return this->inum;
    }

    // return the number of neighbours of a given atom
    template<int Level, int MaxLevel>
    inline size_t get_cluster_size(const ClusterRef_t<Level,
				   MaxLevel>& cluster) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      return this->numneigh[cluster.get_atoms().back().get_index()];
    }

    // return the number of atoms forming the next higher cluster with this one
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const {
      static_assert(Level == traits::MaxLevel-1,
                    "this implementation only handles atoms and pairs");
      auto && i_atom_id{cluster.get_atoms().back().get_index()};
      return this->firstneigh[std::move(i_atom_id)][j_atom_id];
    }

    // return the number of neighbours of a given atom
    inline size_t get_atom_id(const Parent& /*cluster*/,
                              int i_atom_id) const {
      return this->ilist[i_atom_id];
    }

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

    size_t get_nb_clusters(int cluster_size);

  protected:
    int inum{}; // total number of atoms in structure
    //int tot_num{}; //includes ghosts
    //int * ilist{};
    //int * numneigh{};
    //int ** firstneigh{};
    // positions double **x{};
    //atom_type //int * type{};

  private:
  };


  /* ---------------------------------------------------------------------- */
  // adjust for triplets
  template<int Level, int MaxLevel>
  inline int NeighbourhoodManagerChain::
  get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const {
    static_assert(Level == 2, "This class cas only handle single atoms and pairs");
    static_assert(MaxLevel == traits::MaxLevel, "Wrong maxlevel");

    auto atoms{cluster.get_atoms()};
    auto i{atoms.front().get_index()};
    auto j{cluster.get_index()};
    auto main_offset{this->offsets[i]};
    return main_offset + j;
  }

  /* ---------------------------------------------------------------------- */
  // specialisation for just atoms
  template <>
  inline int NeighbourhoodManagerChain:: template
  get_offset_impl<1, 2>(const ClusterRef_t<1, 2>& cluster) const {
    return cluster.get_atoms().back().get_index();
  }

}  // rascal

#endif /* NEIGHBOURHOOD_MANAGER_CHAIN_H */
