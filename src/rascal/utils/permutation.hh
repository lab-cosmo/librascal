/**
 * @file   permutation.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   08 Jan 2020
 *
 * @brief  Compile-time permutation (for evaluating symmetry functions in
 * arbitrary cluster order)
 *
 * Copyright © 2020 Till Junge
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */
#include <type_traits>
#ifndef SRC_RASCAL_UTILS_PERMUTATION_HH_
#define SRC_RASCAL_UTILS_PERMUTATION_HH_

namespace rascal {

  template <size_t Size, size_t First, size_t Second, size_t Third = Size - 1>
  struct Permutation {
    static_assert((First != Second) && (Size > First) && (Size > Second) &&
                      (Size > Third) &&
                      ((Size == 2) ||
                       ((Size == 3) && (Second != Third) && (First != Third))),
                  "Not a valid  pair or triplet permutation");
    static constexpr int leading() { return First; }
    static constexpr int second() { return Second; }
    template <bool IsTriplet = (Size == 3)>
    static constexpr std::enable_if_t<IsTriplet, int> third() {
      return Third;
    }

    template <typename StructureManager, typename Cluster>
    static int leading(StructureManager & manager, const Cluster & cluster) {
      return manager.get_atom_index(cluster.get_atom_tag_list(leading()));
    }

    template <typename StructureManager, typename Cluster>
    static int second(StructureManager & manager, const Cluster & cluster) {
      return manager.get_atom_index(cluster.get_atom_tag_list(second()));
    }

    template <class Derived, bool IsPair = (Size == PairOrder),
              std::enable_if_t<IsPair, int> = 0>
    static auto
    flip_direction(const Eigen::MatrixBase<Derived> & direction_vector)
        -> decltype(1 * direction_vector) {
      static_assert(Size == PairOrder, "IsPair is a SFINAE, don't touch");
      if (First > Second) {
        return -1 * direction_vector;
      } else {
        return 1 * direction_vector;
      };
    }

    template <typename StructureManager, typename Cluster,
              bool IsTriplet = (Size == 3)>
    static std::enable_if_t<IsTriplet, int> third(StructureManager & manager,
                                                  const Cluster & cluster) {
      return manager.get_atom_index(cluster.get_atom_tag_list(third()));
    }
  };

}  // namespace rascal

#endif  // SRC_RASCAL_UTILS_PERMUTATION_HH_
