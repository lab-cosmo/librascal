/**
 * file   named_tuple.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   18 March 2019
 *
 * @brief implement a named tuple object used in the structure manager factory
 *
 * Copyright © 2018 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_NAMED_TUPLE_HH_
#define SRC_NAMED_TUPLE_HH_


// Named tuple for C++
// Example code from http://vitiy.info/
// Written by Victor Laskin (victor.laskin@gmail.com)

// Parts of code were taken from: https://gist.github.com/Manu343726/081512c43814d098fe4b


namespace rascal {
  namespace internal {

    template<int N, typename... Ts>
    using NthTypeOf =
        typename std::tuple_element<N, std::tuple<Ts...>>::type;

    template<char x, char... xs>
    struct hash_calc {
        static constexpr std::uint64_t apply () {
          return  (hash_calc<xs...>::apply() ^ x) * 16777619u;
        };
    };

    template<char x>
    struct hash_calc<x> {
        static constexpr std::uint64_t apply () {
          return  2166136261u;
        };
    };

    template<char... xs>
    constexpr std::uint64_t hash () {
        return hash_calc<xs...>::apply();
    }

    namespace detail {
      using hash_type = std::uint64_t;

      constexpr hash_type fnv_basis = 14695981039346656037ull;
      constexpr hash_type fnv_prime = 109951162821ull;

      // FNV-1a 64 bit hash
      constexpr hash_type sid_hash(const char *str, hash_type hash = fnv_basis) noexcept
      {
          return *str ? sid_hash(str + 1, (hash ^ *str) * fnv_prime) : hash;
      }
    }

    template<std::uint64_t HashValue>
    struct Hash {
      using value = std::integral_constant<std::uint64_t, HashValue>;
    };

    template <typename Hash, class T>
    struct named_value {

      const std::decay_t<T> value;

      using hash = typename Hash::value; ///< name

      named_value(T&& t)
        : value{std::forward<T>(t)}
        { };

      const std::decay_t<T>& get() const {
        return this->value;
      }
    };

    /// Named tuple is just tuple of named params
    template <typename... Params>
    struct named_tuple {
      const std::tuple<Params...> data;
      static constexpr const std::size_t size{sizeof...(Params)};

      template <typename... Args>
      named_tuple(Args&&... args)
        : data{std::make_tuple(args...)} {}

      static const std::size_t error = -1;

      template<std::size_t I = 0, typename HashValue>
      constexpr typename std::enable_if<I == size, const std::size_t>::type
      static get_element_index() {
        return error;
      }

      template<std::size_t I = 0, typename HashValue>
      constexpr typename std::enable_if<I < size, const std::size_t>::type
      static get_element_index() {
        using Value_t = NthTypeOf<I, std::decay_t<Params>...>;
        return (std::is_same<typename Value_t::hash, HashValue>::value) ? I : get_element_index<I + 1, HashValue>();
      }

      template<typename HashValue>
      const auto& get() const {
        constexpr std::size_t index = get_element_index<0, HashValue>();
        static_assert((index != error), "Wrong named tuple key");
        const auto& param = std::get< index >(this->data);
        return param.get();
      }

      template<size_t index>
      const auto& get_by_index() const {
        static_assert((index <= size), "Wrong named tuple index");
        const auto& param = std::get< index >(this->data);
        return param.get();
      }

      template<size_t index>
      const auto& get_item_by_index() const {
        static_assert((index <= size), "Wrong named tuple index");

        const auto& param = std::get< index >(this->data);
        return param;
      }

      template<typename Hash>
      const auto& operator[](Hash&& /*param*/) {
        return get<typename Hash::value>();
      }

    };


    template<size_t index, typename... Params>
    const auto& get(const named_tuple<Params...>& tup) noexcept {
      return tup.template get_item_by_index<index>();
    }

  } // namespace internal

  // only clang/gcc compatible
  template <typename Char, Char... Cs>
  constexpr auto operator""_hash() {
    return internal::Hash<internal::hash<Cs...>()>{};
  };

  template <char... Cs>
  constexpr auto const_hash() {
    return internal::Hash<internal::hash<Cs...>()>{};
  }

  template <typename Hash, typename Arg>
  decltype(auto) make_named_value(Hash&& /*hash*/, Arg&& arg) {
      return internal::named_value<Hash,Arg>(std::forward<Arg>(arg));
  }

  template <typename... Args>
  decltype(auto) make_named_tuple(Args&&... args) {
      return internal::named_tuple<Args...>(std::forward<Args>(args)...);
  }

  namespace internal {
    template <typename Tuple1, size_t... Indices1, typename Tuple2, size_t... Indices2>
    decltype(auto) named_tuple_cat_impl(Tuple1&& tup1, Tuple2&& tup2,
                    std::index_sequence<Indices1...>, std::index_sequence<Indices2...>) {
      return make_named_tuple(
          get<Indices1>(std::forward<Tuple1>(tup1))...,
          get<Indices2>(std::forward<Tuple2>(tup2))...
        );
    }
  } // namespace internal

  template< class T >
  class named_tuple_size;

  template< class... Types >
  class named_tuple_size< internal::named_tuple<Types...> >
    : public std::integral_constant<std::size_t, sizeof...(Types)> { };

  template <typename Tuple1, typename Tuple2>
  decltype(auto) named_tuple_cat(Tuple1&& tup1, Tuple2&& tup2)
  {
    return internal::named_tuple_cat_impl(
    std::forward<Tuple1>(tup1),
    std::forward<Tuple2>(tup2),
    std::make_index_sequence<named_tuple_size<std::decay_t<Tuple1>>::value>{},
    std::make_index_sequence<named_tuple_size<std::decay_t<Tuple2>>::value>{}
    );
  }

#define param(x) internal::Hash<internal::detail::sid_hash(x)>{}
} // namespace rascal

#endif  // SRC_NAMED_TUPLE_HH_
