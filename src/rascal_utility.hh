/**
 * file   rascal_utility.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @date   16 Jul 2018
 *
 * @brief  utilities for rascal
 *
 * Copyright © 2018 Markus Stricker, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
 *
 * Rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * Rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef RASCAL_UTILITY_H
#define RASCAL_UTILITY_H

// Detects which compiler is used
#if defined(__clang__)
#define CLANG_COMPILER
#elif defined(__GNUC__) || defined(__GNUG__)
#define GCC_COMPILER
#endif

#include "representations/representation_manager_base.hh"


#include <utility>
#include <string>
#include <regex> // NOLINT
#include <tuple>
#include <map>

namespace rascal {
  namespace internal {

    /* ---------------------------------------------------------------------- */
    /**
     * Helper functions to apply a functor to all items in a tuple. The actual
     * function is ``for_each``, other functions construct the template loop over
     * the items in the tuple.
     */
    template<typename Func, typename Last>
    inline void for_each_impl(Func&& f, Last&& last) {
      f(last);
    }

    template<typename Func, typename Head, typename ... Tail>
    inline void for_each_impl(Func&& f, Head&& head, Tail&&...tail) {
      f(head);
      for_each_impl(std::forward<Func>(f), tail...);
    }

    template<typename Func, size_t ... Indices, typename ... Args>
    inline void for_each_helper(Func&& f,
                                std::index_sequence<Indices...>,
                                std::tuple<Args...>&& tup) {
      for_each_impl(std::forward<Func>(f),
                    std::forward<Args>(std::get<Indices>(tup))...);
    }

    /**
     * Utility for applying a function to individual tuple elements. `tup` is a
     * tuple that can be templated with an arbitrary number of arguments. `f` is
     * the function that should be applied to each element of the tuple.
     */
    template<typename Func, typename ... Args>
    inline void for_each(std::tuple<Args...>& tup, Func&& f) {
      for_each_helper(std::forward<Func>(f),
                      std::index_sequence_for<Args...>{},
                      std::forward<std::tuple<Args...>>(tup));
    }

    /* ---------------------------------------------------------------------- */
    // useful functors to be applied to tuples somewhere else in the code

    //! Functor for resetting properties to zero size
    struct ResizePropertyToZero {
      template<typename T> void operator() (T& t) { t.resize_to_zero();}
    };

    /* ---------------------------------------------------------------------- */
    // inspiered from
    // https://blog.molecular-matters.com/2015/12/11/getting-
    // the-type-of-a-template-argument-as-string-without-rtti/
    /**
     * Transforms the template typename to a string.
     * This functionality is compiler dependant so for the moment
     * clang and gcc are compatible.
     * @template T type that should be stringifyied
     * @returns std::string name of the type
     */
    template <typename T>
    struct GetTypeNameHelper {
      static const std::string GetTypeName() {
        // The output of of Pretty Function depends on the compiler
        // the #define strings is a pain to split
#if defined(GCC_COMPILER)
        #define FUNCTION_MACRO __PRETTY_FUNCTION__
        #define PREFIX "static const string rascal::internal::GetTypeNameHelper<T>::GetTypeName() [with T = " // NOLINT
        #define SUFFIX_1 "; std::__cxx11::string = std::__cxx11::basic_string<char>]" // NOLINT
        #define SUFFIX_2 ""
        #define NUM_TYPE_REPEATS 1
#elif defined(CLANG_COMPILER)
        #define FUNCTION_MACRO __PRETTY_FUNCTION__
        #define PREFIX "static const std::string rascal::internal::GetTypeNameHelper<" // NOLINT
        #define SUFFIX_1 ">::GetTypeName() [T ="
        #define SUFFIX_2 "]"
        #define NUM_TYPE_REPEATS 2
#else
        #error "No implementation for current compiler"
#endif

        const size_t funcNameLength{sizeof(FUNCTION_MACRO) - 1u};
        const size_t prefixLength{sizeof(PREFIX) - 1u};
        const size_t suffixLength{sizeof(SUFFIX_1) - 1u +
                                            sizeof(SUFFIX_2) - 1u};
        const size_t typeLength{(funcNameLength -
                        (prefixLength + suffixLength)) / NUM_TYPE_REPEATS};
        std::string typeName{FUNCTION_MACRO + prefixLength, typeLength};
        return typeName;
        #undef FUNCTION_MACRO
        #undef PREFIX
        #undef SUFFIX_1
        #undef SUFFIX_2
        #undef NUM_TYPE_REPEATS
      }
    };
    //! return a pretty form of the template typename
    template <typename T>
    std::string GetTypeName() {
      std::string full_typeName = GetTypeNameHelper<T>::GetTypeName();

      std::string tn2{std::regex_replace(full_typeName,
                                    std::regex("rascal::"), "")};
      std::string tn3{std::regex_replace(tn2, std::regex("<"), "_")};
      std::string tn4{std::regex_replace(tn3, std::regex(">"), "")};
      std::string tn5{std::regex_replace(tn4, std::regex(" "), "")};
      std::string tn6{std::regex_replace(tn5, std::regex(","), "_")};
      return tn6;
    }

  }  // internal
}  // rascal

#endif /* RASCAL_UTILITY_H */
