/**
 * file   utility.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 *
 * @date   16 Jul 2018
 *
 * @brief  utilities for rascal
 *
 * Copyright © 2018 Markus Stricker, COSMO (EPFL), LAMMM (EPFL)
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

#include <utility>

namespace rascal {
  /* ---------------------------------------------------------------------- */
  /**
   * Utility for applying a function to individual tuple elements. The function
   * `f` which is applied should be defined close to the usage.
   */
  template<typename Func, typename Last>
  void for_each_impl(Func&& f, Last&& last) {
    f(last);
  }

  template<typename Func, typename Head, typename ... Tail>
  void for_each_impl(Func&& f, Head&& head, Tail&&...tail) {
    f(head);
    for_each_impl( std::forward<Func>(f), tail...);
  }

  template<typename Func, size_t ... Indices, typename ... Args>
  void for_each_helper(Func&& f,
                       std::index_sequence<Indices...>,
                       std::tuple<Args...>&& tup) {
    for_each_impl(std::forward<Func>(f),
                  std::forward<Args>(std::get<Indices>(tup))...);
  }

  template<typename Func, typename ... Args>
  void for_each(std::tuple<Args...>& tup, Func&& f) {
    for_each_helper(std::forward<Func>(f),
                    std::index_sequence_for<Args...>{},
                    std::forward<std::tuple<Args...>>(tup));
  }
  /* ---------------------------------------------------------------------- */

}  // rascal
