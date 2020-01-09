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
namespace rascal {

  template <size_t Size, size_t First, size_t Second, size_t Third>
  struct Permutation {
    static_assert((First != Second) && (First != Third) && (Second != Third) &&
                      (Size > First) && (Size > Second) && (Size > Third) &&
                      ((Size == 2) || (Size == 3)),
                  "Not a valid  pair or triplet permutation");
    static constexpr int leading() { return First; }
    static constexpr int second() { return Second; }
    template <bool IsTriplet = (size == 3)>
    static constexpr std::enable_if_t<IsTriplet, int> third() {
      return Third;
    }

    template <typename StructureManager, typename Cluster>
    static int leading(StructureManager & manager, const Cluster & cluster) {
      return manager.get_atom_index(pair.get_atom_tag_list(leading()));
    }
    template <typename StructureManager, typename Cluster>
    static int second(StructureManager & manager, const Cluster & cluster) {
      return manager.get_atom_index(pair.get_atom_tag_list(second()));
    }
    template <typename StructureManager, typename Cluster>
    template <bool IsTriplet = (size == 3)>
    static std::enable_if_t<IsTriplet, int>
    third(StructureManager & manager, const Cluster & cluster) {
      return manager.get_atom_index(pair.get_atom_tag_list(third()));
    }
  };

}  // namespace rascal

Permutation<3,0,1,2>;
Permutation<3,1,0,2>;Permutation<3,1,2,0>;Permutation<3,2,1,0>;Permutation<3,2,0,1>;