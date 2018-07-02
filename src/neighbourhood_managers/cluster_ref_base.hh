/**
 * file   ClusterRefBase.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   21 Jun 2018
 *
 * @brief  a base class for getting access to clusters
 *
 * Copyright © 2018 Till Junge, Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef CLUSTERREFBASE_H
#define CLUSTERREFBASE_H

#include <array>

namespace rascal {
  namespace internal {
    template <int Depth, int HiDepth,typename T, int... Ints>
    std::array<T, Depth>
    head_helper(const std::array<T, HiDepth> & arr,
                std::index_sequence<Ints...>) {
      return std::array<T, Depth> {arr[Ints]...};
    }

    template <int Depth, int HiDepth,typename T>
    std::array<T, Depth> head(const std::array<T, HiDepth> & arr) {
      return head_helper(arr, std::make_index_sequence<Depth>{});
    }

  }  // internal

  template<int Level, int Depth=2>
  class ClusterRefBase
  {
  public:
    //! Default constructor
    ClusterRefBase() = delete;

    //! direct constructor
    ClusterRefBase(std::array<int, Level> atom_indices,
                   std::array<int, Depth> cluster_indices={}):
      atom_indices{atom_indices}, cluster_indices{cluster_indices} {}

    // //! constructor from higher depth
    // template<int HiDepth>
    // ClusterRefBase(const ClusterRefBase<Level, HiDepth> & other):
    //   atom_indices{other.atom_indices},
    //   cluster_indices{internal::head<Depth>(other.cluster_indices)} {
    //     static_assert(HiDepth >= Depth,
    //                   "You are trying to access a property that "
    //                   "does not exist at this low a level in the "
    //                   "adaptor stack.");
    //   }

    //! Copy constructor
    ClusterRefBase(const ClusterRefBase &other) = default; //

    //! Move constructor
    ClusterRefBase(ClusterRefBase &&other) = default;

    //! Destructor
    virtual ~ClusterRefBase() = default;

    //! Copy assignment operator
    ClusterRefBase& operator=(const ClusterRefBase &other) = delete;

    //! Move assignment operator
    ClusterRefBase& operator=(ClusterRefBase &&other) = default;

    const std::array<int, Level> & get_atom_indices() const {return this->indices;}

    const int & front() const{return this->atom_indices.front();}
    const int & back() const{return this->atom_indices.back();}

    inline int get_cluster_index(int depth) const {return this->cluster_indices[depth];}

  protected:
    std::array<int, Level> atom_indices;
    /**
     * cluster indices by depth level, highest depth, means last
     * adaptor, and mean last entry (.back())
     */
    std::array<int, Depth> cluster_indices;
  private:
  };

} // rascal

#endif /* CLUSTERREFBASE_H */
