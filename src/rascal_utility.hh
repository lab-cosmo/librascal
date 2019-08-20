/**
 * file   rascal_utility.hh
 *
 * @author Markus Stricker <markus.stricker@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 * @date   16 Jul 2018
 *
 * @brief  utilities for rascal
 *
 * Copyright  2018 Markus Stricker, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_RASCAL_UTILITY_HH_
#define SRC_RASCAL_UTILITY_HH_

// Detects which compiler is used
#if defined(__clang__)
#define CLANG_COMPILER
#elif defined(__GNUC__) || defined(__GNUG__)
#define GCC_COMPILER
#endif

#include "representations/representation_manager_base.hh"

#include <utility>
#include <string>
#include <regex>  // NOLINT
#include <tuple>
#include <map>
#include <fstream>

namespace rascal {
  namespace internal {

    /* ---------------------------------------------------------------------- */
    /**
     * Utility to deduce the type of a Manager with a list of Adaptors
     */
    template <typename ManagerImplementation,
              template <class> class AdaptorImplementation,
              template <class> class... AdaptorImplementationPack>
    struct AdaptorTypeStacker {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using type =
          typename AdaptorTypeStacker<Manager_t,
                                      AdaptorImplementationPack...>::type;
    };

    template <typename ManagerImplementation,
              template <class> class AdaptorImplementation>
    struct AdaptorTypeStacker<ManagerImplementation, AdaptorImplementation> {
      using type = AdaptorImplementation<ManagerImplementation>;
    };
    /* ---------------------------------------------------------------------- */
    /**
     * Utilities to combine Enum to flatten the nested switch cases
     * The caveat is that the enum list need to be finished with the
     * End_
     * taken from:
     * https://www.fluentcpp.com/2017/06/27/how-to-collapse-nested-switch-statements/
     * // NOLINT
     */

    //! get the underlying value of the enum
    template <typename Enum>
    constexpr size_t enumValue(Enum e) {
      return static_cast<size_t>(e);
    }

    //! compute the lenght of the enum assuming the last element is End_
    template <typename Enum>
    constexpr size_t enumSize() {
      return enumValue(Enum::End_);
    }

    //! combine 2 enum into a new integer making sure types are properly
    //! provided
    // template <typename Enum1, typename Enum2>
    // struct CombineEnums {
    //   constexpr size_t operator()(Enum1 e1, Enum2 e2) {
    //     return enumValue(e1) * enumSize<Enum2>() + enumValue(e2);
    //   }
    // };
    // the above code does not compile with gcc 5 and 6 (4 and 7 works though)
    // and clang has no problem. it has to do with the implementation of the
    // standard see for more details
    // https://stackoverflow.com/questions/16493652/constexpr-not-working-if-the-function-is-declared-inside-class-scope
    // // NOLINT
    template <typename Enum1, typename Enum2>
    constexpr size_t combineEnums(Enum1 e1, Enum2 e2) {
      return enumValue(e1) + enumSize<Enum1>() * enumValue(e2);
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Implementation of the generation of an index sequence from Min to Max
     */

    template <size_t N, size_t... Seq>
    constexpr std::index_sequence<N + Seq...>
    add_to_sequence(std::index_sequence<Seq...>) {
      return {};
    }

    template <size_t Min, size_t Max>
    using make_index_range =
        decltype(add_to_sequence<Min>(std::make_index_sequence<Max - Min>()));

    /* ---------------------------------------------------------------------- */
    //! Implementation of index_apply
    template <class F, size_t... Is>
    constexpr auto index_apply_impl(F func, std::index_sequence<Is...>) {
      return func(std::integral_constant<size_t, Is>{}...);
    }

    /**
     * index_apply and its implementation take a callable of type F and return
     * it evaluated with an index sequence starting from Min to Max
     */
    template <size_t Min, size_t Max, class F>
    constexpr auto index_apply(F func) {
      return index_apply_impl(func, make_index_range<Min, Max>{});
    }

    /**
     * Extract the elements of a tuple from Min to Max
     */
    template <size_t Min, size_t Max, class Tuple>
    constexpr auto take_range(Tuple t) {
      static_assert(Min <= Max, "Min should be smaller or equal than Max");
      return index_apply<Min, Max>(
          [&](auto... Is) { return make_tuple(std::get<Is>(t)...); });
    }

    /**
     * Evaluate a callable of type F with the elements of the Tuple
     */
    template <class Tuple, class F>
    constexpr auto apply(F func, Tuple t) {
      return index_apply<0, std::tuple_size<Tuple>{}>(
          [&](auto... Is) { return func(std::get<Is>(t)...); });
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Helper functions to apply a functor to all items in a tuple. The actual
     * function is ``for_each``, other functions construct the template loop
     * over the items in the tuple.
     */
    template <typename Func, typename Last>
    inline void for_each_impl(Func && f, Last && last) {
      f(last);
    }

    template <typename Func, typename Head, typename... Tail>
    inline void for_each_impl(Func && f, Head && head, Tail &&... tail) {
      f(head);
      for_each_impl(std::forward<Func>(f), tail...);
    }

    template <typename Func, size_t... Indices, typename... Args>
    inline void for_each_helper(Func && f, std::index_sequence<Indices...>,
                                std::tuple<Args...> && tup) {
      for_each_impl(std::forward<Func>(f),
                    std::forward<Args>(std::get<Indices>(tup))...);
    }

    /**
     * Utility for applying a function to individual tuple elements. `tup` is a
     * tuple that can be templated with an arbitrary number of arguments. `f` is
     * the function that should be applied to each element of the tuple.
     */
    template <typename Func, typename... Args>
    inline void for_each(std::tuple<Args...> & tup, Func && f) {
      for_each_helper(std::forward<Func>(f), std::index_sequence_for<Args...>{},
                      std::forward<std::tuple<Args...>>(tup));
    }

    /* ---------------------------------------------------------------------- */
    // useful functors to be applied to tuples somewhere else in the code

    //! Functor for resetting properties to zero size
    struct ResizePropertyToZero {
      template <typename T>
      void operator()(T & t) {
        t.resize_to_zero();
      }
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
#define PREFIX                                                                 \
  "static const string rascal::internal::GetTypeNameHelper<T>::GetTypeName() " \
  "[with T = "  // NOLINT
#define SUFFIX_1                                                               \
  "; std::__cxx11::string = std::__cxx11::basic_string<char>]"  // NOLINT
#define SUFFIX_2 ""
#define NUM_TYPE_REPEATS 1
#elif defined(CLANG_COMPILER)
#define FUNCTION_MACRO __PRETTY_FUNCTION__
#define PREFIX "static const std::string rascal::internal::GetTypeNameHelper<"
#define SUFFIX_1 ">::GetTypeName() [T ="
#define SUFFIX_2 "]"
#define NUM_TYPE_REPEATS 2
#else
#error "No implementation for current compiler"
#endif

        const size_t funcNameLength{sizeof(FUNCTION_MACRO) - 1u};
        const size_t prefixLength{sizeof(PREFIX) - 1u};
        const size_t suffixLength{sizeof(SUFFIX_1) - 1u + sizeof(SUFFIX_2) -
                                  1u};
        const size_t typeLength{
            (funcNameLength - (prefixLength + suffixLength)) /
            NUM_TYPE_REPEATS};
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

      std::string tn2{
          std::regex_replace(full_typeName, std::regex("rascal::"), "")};
      std::string tn3{std::regex_replace(tn2, std::regex("<"), "_")};
      std::string tn4{std::regex_replace(tn3, std::regex(">"), "")};
      std::string tn5{std::regex_replace(tn4, std::regex(" "), "")};
      std::string tn6{std::regex_replace(tn5, std::regex(","), "_")};
      return tn6;
    }

    /**
     * Reads a binary file and puts it into a vector
     * Taken from
     * https://stackoverflow.com/questions/15138353/how-to-read-a-binary-file-into-a-vector-of-unsigned-chars
     * // NOLINT
     */
    template <typename BINARY>
    void read_binary_file(const std::string & filename,
                          std::vector<BINARY> & vec) {
      // open the file:
      std::ifstream file(filename, std::ios::binary);

      // Stop eating new lines in binary mode!!!
      file.unsetf(std::ios::skipws);

      // get its size:
      std::streampos fileSize;

      file.seekg(0, std::ios::end);
      fileSize = file.tellg();
      file.seekg(0, std::ios::beg);

      // reserve capacity
      vec.reserve(fileSize);

      // read the data:
      vec.insert(vec.begin(), std::istream_iterator<BINARY>(file),
                 std::istream_iterator<BINARY>());
    }

  }  // namespace internal
}  // namespace rascal

#endif  // SRC_RASCAL_UTILITY_HH_
