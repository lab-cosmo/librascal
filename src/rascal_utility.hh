/**
 * file   rascal_utility.hh
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

#ifndef RASCAL_UTILITY_H
#define RASCAL_UTILITY_H

#include <utility>


namespace rascal {
  namespace internal {
    /* ---------------------------------------------------------------------- */
    /** Helper functions to apply a functor to all items in a tuple.
     * The actual function is for_each, other functions construct the
     * template loop over the items in the tuple.
     */
    template<typename Func, typename Last>
    inline void for_each_impl(Func&& f, Last&& last) {
      f(last);
    }

    template<typename Func, typename Head, typename ... Tail>
    inline void for_each_impl(Func&& f, Head&& head, Tail&&...tail) {
      f(head);
      for_each_impl( std::forward<Func>(f), tail...);
    }

    template<typename Func, size_t ... Indices, typename ... Args>
    inline void for_each_helper(Func&& f,
                         std::index_sequence<Indices...>,
                         std::tuple<Args...>&& tup) {
      for_each_impl(std::forward<Func>(f),
                    std::forward<Args>(std::get<Indices>(tup))...);
    }

    /**
     * Utility for applying a function to individual tuple elements.
     * `tup` is a tuple that can be templated with an arbitrary number
     * of arguments. `f` is the function that should be applied to each
     * element of the tuple.
     */
    template<typename Func, typename ... Args>
    inline void for_each(std::tuple<Args...>& tup, Func&& f) {
      for_each_helper(std::forward<Func>(f),
                      std::index_sequence_for<Args...>{},
                      std::forward<std::tuple<Args...>>(tup));
    }
    /* ---------------------------------------------------------------------- */
    /* A collection of useful functors to be applied to tuples somewhere
     * else in the code. */
     /**
     * Functor for resetting properties to zero size
     */
    struct ResizePropertyToZero {
      template<typename T>
      void operator() (T& t) {
        t.resize_to_zero();
      }
    };
  }  // internal
}  // rascal

#endif /* RASCAL_UTILITY_H */