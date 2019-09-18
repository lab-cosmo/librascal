/**
 * @file   for_each_at_order.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief defines a generic tool to apply a function to each cluster of a given
 * order in a manager
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with librascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_UTILS_FOR_EACH_AT_ORDER_HH_
#define SRC_UTILS_FOR_EACH_AT_ORDER_HH_

#include <cstddef>
#include <utility>

namespace rascal {

  namespace utils {

    namespace detail {

      /**
       * General case
       */
      template <size_t CurrentOrder, size_t ApplicationOrder, class Cluster,
                class Function, class... Args>
      struct ForEachAtOrderHelper {
        using NextLoopHelper =
            ForEachAtOrderHelper<CurrentOrder + 1, ApplicationOrder, Cluster,
                                 Function, Args...>;

        static void loop(Cluster & cluster, Function && function,
                         Args &&... args) {
          for (auto nextcluster : cluster) {
            NextLoopHelper::loop(nextcluster, std::forward<Function>(function),
                                 std::forward<Args>(args)...);
          }
        }
      };

      /**
       * recursion end
       */
      template <size_t Order, class Cluster, class Function, class... Args>
      struct ForEachAtOrderHelper<Order, Order, Cluster, Function, Args...> {
        static void loop(Cluster & cluster, Function && function,
                         Args &&... args) {
          function(cluster, std::forward<Args>(args)...);
        }
    };

    }  // namespace detail

    template <size_t ApplicationOrder, class Manager, class Function,
              class... Args>
    void for_each_at_order(Manager & manager, Function && function,
                           Args &&... args) {
      using LoopHelper =
          detail::ForEachAtOrderHelper<1, ApplicationOrder, Manager, Function,
                                       Args...>;
      for (auto && atom : manager) {
        LoopHelper::loop(atom, function, args...);
      }
    }

  }  // namespace utils

}  // namespace rascal

#endif  // SRC_UTILS_FOR_EACH_AT_ORDER_HH_
