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

  /**
   * Refers to the pattern of repeated species in the permuted cluster, as it is
   * evaluated by the symmetry function
   */
  enum class RepeatedSpecies {
    Unknown,    // has not been evaluated yet
    Not,        // all species in this cluster are unique
    All,        // all atoms in this cluster are same species
    FirstTwo,   // the first two atoms of this cluster are of same species
    SecondTwo,  // the second two atoms of this cluster are of same species
    OuterTwo    // the first and last atom in this cluster are of same species
  };

  // // Proposed for usage in test for iterating RepeatedSpecies
  // static constexpr RepeatedSpecies AllRepSpecies[] = {
  //     RepeatedSpecies::Unknown,   RepeatedSpecies::Not,
  //     RepeatedSpecies::All,       RepeatedSpecies::FirstTwo,
  //     RepeatedSpecies::SecondTwo, RepeatedSpecies::OuterTwo};

  constexpr size_t nb_triplet_orderings(const RepeatedSpecies rep) {
    switch (rep) {
    case RepeatedSpecies::Not: {
      return 1;
      break;
    }
    case RepeatedSpecies::FirstTwo: {
      // fall-through
    }
    case RepeatedSpecies::OuterTwo: {
      return 2;
      break;
    }
    case RepeatedSpecies::SecondTwo: {
      return 0;  // yes, this case has nothing to compute
      break;
    }
    case RepeatedSpecies::All: {
      return 3;
      break;
    }
    default: {
      throw std::runtime_error("Unknown species repetition");
    }
    }
  }

  template <size_t Size_, size_t First, size_t Second, size_t Third = Size_ - 1>
  struct Permutation {
    static_assert((First != Second) && (Size_ > First) && (Size_ > Second) &&
                      (Size_ > Third) &&
                      ((Size_ == 2) ||
                       ((Size_ == 3) && (Second != Third) && (First != Third))),
                  "Not a valid  pair or triplet permutation");
    static constexpr int leading() { return First; }
    static constexpr int second() { return Second; }
    static constexpr size_t Size{Size_};
    template <bool IsTriplet = (Size_ == 3)>
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

    template <class Derived, bool IsPair = (Size_ == PairOrder),
              std::enable_if_t<IsPair, int> = 0>
    static auto
    flip_direction(const Eigen::MatrixBase<Derived> & direction_vector)
        -> decltype(1 * direction_vector) {
      static_assert(Size_ == PairOrder, "IsPair is a SFINAE, don't touch");
      if (First > Second) {
        return -1 * direction_vector;
      } else {
        return 1 * direction_vector;
      }
    }

    template <typename StructureManager, typename Cluster,
              bool IsTriplet = (Size_ == 3)>
    static std::enable_if_t<IsTriplet, int> third(StructureManager & manager,
                                                  const Cluster & cluster) {
      return manager.get_atom_index(cluster.get_atom_tag_list(third()));
    }

    template <RepeatedSpecies RepSpecies>
    constexpr static std::array<std::array<size_t, 3>,
                                nb_triplet_orderings(RepSpecies)>
    get_triplet_orderings() {
      std::array<std::array<size_t, 3>, nb_triplet_orderings(RepSpecies)>
          retval{};

      switch (RepSpecies) {
      case RepeatedSpecies::Not: {
        retval[0][0] = leading();
        retval[0][1] = second();
        retval[0][2] = third();
        break;
      }
      case RepeatedSpecies::FirstTwo: {
        retval[0][0] = leading();
        retval[0][1] = second();
        retval[0][2] = third();

        retval[1][0] = second();
        retval[1][1] = leading();
        retval[1][2] = third();
        break;
      }
      case RepeatedSpecies::OuterTwo: {
        retval[0][0] = leading();
        retval[0][1] = second();
        retval[0][2] = third();

        retval[1][0] = third();
        retval[1][1] = second();
        retval[1][2] = leading();
        break;
      }
      case RepeatedSpecies::All: {
        retval[0][0] = leading();
        retval[0][1] = second();
        retval[0][2] = third();

        retval[1][0] = second();
        retval[1][1] = third();
        retval[1][2] = leading();

        retval[1][0] = third();
        retval[1][1] = leading();
        retval[1][2] = second();
        break;
      }
      default:
        throw std::runtime_error("Unknown species repetition");
        break;
      }

      return retval;
    }

    constexpr static std::array<double, 3>
    apply_ordering(std::array<double, 3> values,
                   std::array<size_t, 3> ordering) {
      std::array<double, 3> ret_val;
      ret_val[0] = values[ordering[0]];
      ret_val[1] = values[ordering[1]];
      ret_val[2] = values[ordering[2]];
      return ret_val;
      //      std::array<double, 3> ret_val{} : ret_val[0] =
      //      values[ordering[0]];
    }
  };  // Permutation
  template <size_t Size_, size_t First, size_t Second, size_t Third>
  constexpr size_t Permutation<Size_, First, Second, Third>::Size;

}  // namespace rascal

#endif  // SRC_RASCAL_UTILS_PERMUTATION_HH_
