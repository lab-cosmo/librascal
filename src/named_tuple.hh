// Named tuple for C++
// Example code from http://vitiy.info/
// Written by Victor Laskin (victor.laskin@gmail.com)

// Parts of code were taken from: https://gist.github.com/Manu343726/081512c43814d098fe4b

namespace foonathan {
    namespace string_id {
        namespace detail
        {
            using hash_type = std::uint64_t;

            constexpr hash_type fnv_basis = 14695981039346656037ull;
            constexpr hash_type fnv_prime = 109951162821ull;

            // FNV-1a 64 bit hash
            constexpr hash_type sid_hash(const char *str, hash_type hash = fnv_basis) noexcept
            {
                return *str ? sid_hash(str + 1, (hash ^ *str) * fnv_prime) : hash;
            }
        }
    }
} // foonathan::string_id::detail


namespace fn_detail {

      /// Named parameter (could be empty!)
    template <typename Hash, class... Ts>
    struct named_param {

      const std::tuple<std::decay_t<Ts>...> data;

      using hash = Hash;                                                              ///< key

      template<int N>
      using NthTypeOf =
        typename std::tuple_element<N, std::tuple<std::decay_t<Ts>...>>::type;

      named_param(Ts&&... ts)
        : data{std::make_tuple(std::forward<Ts>(ts)...)}
        { };

      template<size_t index>
      const NthTypeOf<index>& get() const {
        return std::get<index>(this->data);
      }

      template <class P>
      named_param<Hash,P> operator=(P&& p){
        return named_param<Hash,P>(std::forward<P>(p));
      };
    };

    template <typename Hash>
    using make_named_param = named_param<Hash>;




    /// Named tuple is just tuple of named params
    template <typename... Params>
    struct named_tuple
    {
        const std::tuple<Params...> data;

        template <typename... Args>
        named_tuple(Args&&... args)
          : data{std::make_tuple(args...)} {}

        static const std::size_t error = -1;

        template<std::size_t I = 0, typename Hash>
        constexpr typename std::enable_if<I == sizeof...(Params), const std::size_t>::type
        static get_element_index() {
            return error;
        }

        template<std::size_t I = 0, typename Hash>
        constexpr typename std::enable_if<I < sizeof...(Params), const std::size_t>::type
        static get_element_index() {
            using elementType = typename std::tuple_element<I, std::tuple<Params...>>::type;
            //return (typeid(typename elementType::hash) == typeid(Hash)) ? I : get_element_index<I + 1, Hash>();
            return (std::is_same<typename elementType::hash, Hash>::value) ? I : get_element_index<I + 1, Hash>();
        }

        template<typename Hash>
        const auto& get() const {
          constexpr std::size_t index = get_element_index<0, Hash>();
          static_assert((index != error), "Wrong named tuple key");
          const auto& param = std::get< index >(this->data);
          return param.template get<0>();
        }

        template<size_t index>
        const auto& get_by_index() const {
          // std::cout << sizeof...(Params) << std::endl;
          static_assert((index <= sizeof...(Params)), "Wrong named tuple index");

          const auto& param = std::get< index >(this->data);
          return param.template get<0>();
        }

        template<size_t index>
        const auto& get_item_by_index() const {
          // std::cout << sizeof...(Params) << std::endl;
          static_assert((index <= sizeof...(Params)), "Wrong named tuple index");

          const auto& param = std::get< index >(this->data);
          return param;
        }

        template<typename NP>
        const auto& operator[](NP&& /*param*/) {
          return get<typename NP::hash>();
        }

    };

    template<int N, typename... Ts>
    using NthTypeOf =
        typename std::tuple_element<N, std::tuple<Ts...>>::type;

    // template<size_t index, typename... Params>
    // constexpr const NthTypeOf<index,Params...>& get(const named_tuple<Params...>& tup) noexcept
    // {
    //   std::cout << tup.data;
    //   return std::get<index>(tup.data);
    // }

    template<size_t index, typename... Params>
    const auto& get(const named_tuple<Params...>& tup) noexcept
    {
      // return std::get<NthTypeOf<index,std::decay_t<Params>...>>(tup.data);
      // return tup.template get_item_by_index<index>();
      using hash = typename NthTypeOf<index,std::decay_t<Params>...>::hash;
      return tup.template get<hash>();
    }
}

template <typename... Args>
decltype(auto) make_named_tuple(Args&&... args)
{
    return fn_detail::named_tuple<Args...>(std::forward<Args>(args)...);
}

