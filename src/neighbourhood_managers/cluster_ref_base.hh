/**
 * file   cluster_ref_base.hh
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

#include <tuple>
#include <array>

namespace rascal {

  /**
   * Depth calculations and manipulations
   */
  /**
   * Computes depth by cluster dimension for new adaptor layer
   */
  template <size_t MaxLevel, class T>
  struct DepthIncreaser{};

  template <size_t MaxLevel, int... Ints>
  struct DepthIncreaser<MaxLevel,
                        std::integer_sequence<int, Ints...>>{
    using type = std::integer_sequence<int, (Ints+1)...>;
  };

  template <size_t MaxLevel, int... Ints>
  using DepthIncreaser_t = typename DepthIncreaser<MaxLevel,
                                                   std::integer_sequence<int, Ints...>>::type;

  /**
   * Extends depth by cluster for an additional cluster dimension
   */
  template <size_t MaxLevel, class T>
  struct DepthExtender{};

  template <size_t MaxLevel, int... Ints>
  struct DepthExtender<MaxLevel,
                       std::integer_sequence<int, Ints...>>{
    using type = std::integer_sequence<int, Ints..., 0>;
  };

  template <size_t MaxLevel, int... Ints>
  using DepthExtender_t = typename DepthExtender<MaxLevel,
                                                 std::integer_sequence<int, Ints...>>::type;

  /**
   * Dynamic access to all depths by cluster dimension (probably not
   * necessary)
   */
  template <size_t MaxLevel, int... Ints>
  constexpr std::array<int, MaxLevel>
  get_depths(std::integer_sequence<int, Ints...>) {
    return std::array<int, MaxLevel>{Ints...};
  }

  namespace internal {
    template <int head, int... tail>
    struct Min {
      constexpr static int value{ head < Min<tail...>::value ? head : Min<tail...>::value};
    };

    template <int head>
    struct Min<head> {
      constexpr static int value{head};
    };

    template <class Sequence>
    struct MinExtractor {};

    template <int... Ints>
    struct MinExtractor<std::integer_sequence<int, Ints...>> {
      constexpr static int value {Min<Ints...>::value};
    };

    template <int Level, class Sequence, int... Ints>
    struct HeadExtractor {};

    template <int... seq>
    struct HeadExtractorTail {
      using type = std::integer_sequence<int, seq...>;
    };
      

    template <int Level, int head, int... tail, int... seq>
    struct HeadExtractor<Level, std::integer_sequence<int, seq...>, head, tail...> {
      using Extractor_t = std::conditional_t
        <(Level > 1),
        HeadExtractor<Level-1,
                      std::integer_sequence<int, seq..., head>,
                      tail...>,
        HeadExtractorTail<seq..., head>>;
      using type = typename Extractor_t::type; 
    };


  }  // internal

  template <int Level, int... Ints>
  constexpr int compute_cluster_depth(const std::integer_sequence<int, Ints...> &) {
    using ActiveDimensions = typename internal::HeadExtractor
      <Level, std::integer_sequence<int>, Ints...>::type;
    return internal::MinExtractor<ActiveDimensions>::value;
  }

  /**
   * Dynamic access to depth by cluster dimension (possibly not
   * necessary)
   */
  template <size_t MaxLevel, int... Ints>
  constexpr int
  get_depth(size_t index, std::integer_sequence<int, Ints...>) {
    constexpr int arr[] {Ints...};
    return arr [index];
  }


  /**
   * Static access to depth by cluster dimension (e.g., for defining
   * template parameter `NbRow` for a property
   */
  template <size_t index, int... Ints>
  constexpr int get(std::integer_sequence<int, Ints...>) {
    return get<index>(std::make_tuple(Ints...));
  }

  
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

  template<int Level, int Depth>
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
