/**
 * file   property_base.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   03 April 2019
 *
 * @brief implementation of a property container that has sparse keys
 *
 * Copyright © 2019 Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_PROPERTY_BLOCK_SPARSE_HH_
#define SRC_STRUCTURE_MANAGERS_PROPERTY_BLOCK_SPARSE_HH_

#include "structure_managers/property_base.hh"
#include "structure_managers/cluster_ref_key.hh"

#include <unordered_map>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <iterator>
#include <type_traits>

namespace rascal {

  namespace internal {

    /**
     * custom hash function for vector, list...
     * https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
     * // NOLINT
     */
    template <class KeyType>
    struct Hash {
      using result_type = size_t;
      using argument_type = KeyType;
      result_type operator()(argument_type const & vec) const {
        result_type seed{vec.size()};
        for (const auto & i : vec) {
          seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
      }
    };

    template <bool IsSorted>
    struct Sorted {};

    /**
     * Special key container that ensures its content is sorted.
     */
    template <class KeyType>
    struct SortedKey {
      using Key_t = KeyType;
      Key_t data;
      using Value_t = typename Key_t::value_type;

      explicit SortedKey(const Key_t & key) : data{std::move(key)} {
        if (data.size() > 1) {
          std::sort(data.begin(), data.end());
        }
      }

      SortedKey(const Sorted<false> &, const Key_t & key) : SortedKey{key} {}

      SortedKey(const Sorted<true> &, const Key_t & key)
          : data{std::move(key)} {}

      Key_t copy_sort(const Key_t & key) {
        Key_t skey{key};
        if (key.size() > 1) {
          std::sort(skey.begin(), skey.end());
        }
        return skey;
      }

      //! access or insert specified element. use with caution !
      inline Value_t & operator[](const size_t & id) { return this->data[id]; }

      inline const Key_t & get_key() const { return data; }
    };

    template <class K, class V>
    class InternallySortedKeyMap {
     public:
      using MyMap_t = std::map<K, V>;
      // using Map_t = std::unordered_map<K, std::array<int, 3>, Hash<K>>;
      using Map_t = std::map<K, std::tuple<int, int, int>>;
      using Precision_t = typename V::value_type;
      using Data_t = Eigen::Array<Precision_t, Eigen::Dynamic, 1>;
      //! the data holder.
      Data_t data{};
      Map_t map{};

      // some member types
      using key_type = typename MyMap_t::key_type;
      using mapped_type = V;
      using value_type = std::pair<const K, mapped_type>;
      using size_type = typename MyMap_t::size_type;
      using reference = typename Eigen::Map<V>;
      using const_reference = typename Eigen::Map<const V>;

      using SortedKey_t = SortedKey<key_type>;

      // member typedefs provided through inheriting from std::iterator
      template <typename Value>
      class Iterator
          : public std::iterator<
                std::bidirectional_iterator_tag,
                std::pair<K, typename std::remove_const<Value>::type>> {
       public:
        using Self_t = Iterator<Value>;

        // to handle properly const and normal cases
        using It_t = typename std::conditional<std::is_const<Value>::value,
                                               typename Map_t::const_iterator,
                                               typename Map_t::iterator>::type;

        using MyData_t = typename std::conditional<std::is_const<Value>::value,
                                                   const Data_t, Data_t>::type;
        // Map<const Matrix> is already write-only so remove the const
        // which is used to determine the cv of the iterator
        using Value_t = typename std::remove_const<Value>::type;

        Iterator(MyData_t & data, It_t map_iterator)
            : data{data}, map_iterator{map_iterator} {}

        Self_t & operator++() {
          map_iterator++;
          // this->update_current_data();
          return *this;
        }
        Self_t operator++(int) {
          Self_t retval = *this;
          ++(*this);
          return retval;
        }
        std::pair<K, Value_t> operator*() const {
          auto && el{*this->map_iterator};
          auto && key{el.first};
          auto && pos{el.second};
          return std::make_pair(key,
                                Value_t(&data[std::get<0>(pos)],
                                        std::get<1>(pos), std::get<2>(pos)));
        }
        Self_t & operator--() {
          map_iterator--;
          // this->update_current_data();
          return *this;
        }
        Self_t operator--(int) {
          Self_t retval = *this;
          --(*this);
          return retval;
        }
        bool operator==(const Self_t & rhs) const {
          return map_iterator == rhs.map_iterator;
        }
        bool operator!=(const Self_t & rhs) const {
          return map_iterator != rhs.map_iterator;
        }

       protected:
        MyData_t & data;
        It_t map_iterator;
      };

      using iterator = Iterator<reference>;
      using const_iterator = Iterator<const const_reference>;

      iterator begin() noexcept { return iterator(data, map.begin()); }
      iterator end() noexcept { return iterator(data, map.end()); }

      const_iterator begin() const noexcept {
        return const_iterator(data, map.begin());
      }
      const_iterator end() const noexcept {
        return const_iterator(data, map.end());
      }

      //! Default constructor
      InternallySortedKeyMap() = default;

      //! Copy constructor
      InternallySortedKeyMap(const InternallySortedKeyMap & other) = default;

      //! Move constructor
      InternallySortedKeyMap(InternallySortedKeyMap && other) = default;

      //! Destructor
      ~InternallySortedKeyMap() = default;

      //! Copy assignment operator
      InternallySortedKeyMap &
      operator=(const InternallySortedKeyMap & other) = default;

      //! Move assignment operator
      InternallySortedKeyMap &
      operator=(InternallySortedKeyMap && other) = default;

      /**
       * Returns a reference to the mapped value of the element with key
       * equivalent to key. If no such element exists, an exception of type
       * std::out_of_range is thrown.
       * The elements of the key are sorted in ascending order.
       *
       */
      reference at(const key_type & key) {
        SortedKey_t skey{key};
        return this->at(skey);
      }
      const_reference at(const key_type & key) const {
        SortedKey_t skey{key};
        return this->at(skey);
      }
      //! access or insert specified element
      reference operator[](const key_type & key) {
        SortedKey_t skey{key};
        return this->operator[](skey);
      }
      const_reference operator[](const key_type & key) const {
        SortedKey_t skey{key};
        return this->operator[](skey);
      }

      /**
       * Same as above but does not try to sort since we know it already is.
       */
      reference at(const SortedKey_t & skey) {
        auto & pos{this->map.at(skey.get_key())};
        return reference(&this->data[std::get<0>(pos)], std::get<1>(pos),
                         std::get<2>(pos));
      }

      const_reference at(const SortedKey_t & skey) const {
        auto & pos{this->map.at(skey.get_key())};
        return const_reference(&this->data[std::get<0>(pos)], std::get<1>(pos),
                               std::get<2>(pos));
      }
      //! access or insert specified element
      reference operator[](const SortedKey_t & skey) {
        auto & pos{this->map[skey.get_key()]};
        return reference(&this->data[std::get<0>(pos)], std::get<1>(pos),
                         std::get<2>(pos));
      }
      const_reference operator[](const SortedKey_t & skey) const {
        auto & pos{this->map[skey.get_key()]};
        return const_reference(&this->data[std::get<0>(pos)], std::get<1>(pos),
                               std::get<2>(pos));
      }

      /**
       * resize the underlying data to the proper size and can initialize
       * the elements
       */
      template <typename Key_List>
      void resize(const Key_List & keys, const int & n_row, const int & n_col,
                  const Precision_t & val) {
        this->resize(keys, n_row, n_col);
        this->data = val;
      }

      template <typename Key_List>
      void resize(const Key_List & keys, const int & n_row, const int & n_col) {
        std::vector<SortedKey_t> skeys{};
        for (auto && key : keys) {
          SortedKey_t skey{key};
          skeys.push_back(skey);
        }
        this->resize(skeys, n_row, n_col);
      }

      void resize(const std::vector<SortedKey_t> & skeys, const int & n_row,
                  const int & n_col) {
        int new_size{0};
        for (auto && skey : skeys) {
          if (this->count(skey) == 0) {
            auto && key{skey.get_key()};
            this->map[key] = std::make_tuple(new_size, n_row, n_col);
            new_size += static_cast<int>(n_row * n_col);
          }
        }
        this->data.resize(new_size);
      }

      //! Returns the number of elements with key that compares equivalent to
      //! the specified argument, which is either 1 or 0 since this container
      //! does not allow duplicates.
      decltype(auto) count(const key_type & key) {
        SortedKey_t skey{key};
        return this->count(skey);
      }

      decltype(auto) count(const SortedKey_t & skey) {
        return this->map.count(skey.get_key());
      }

      //! Erases all elements from the container. After this call, size()
      //! returns zero.
      void clear() noexcept {
        this->data.clear();
        this->map.clear();
      }

      /**
       * returns a vector of the valid keys of the map
       */
      std::vector<key_type> get_keys() const {
        std::vector<key_type> keys{};
        std::transform(this->map.begin(), this->map.end(),
                       std::back_inserter(keys), RetrieveKey());
        return keys;
      }

      /**
       * Normalize the whole vector
       */
      inline void normalize() {
        using ref = typename Eigen::Map<Eigen::VectorXd>;
        auto data_ref{ref(&data[0], data.size())};
        data_ref /= data_ref.norm();
      }

      inline void multiply_offdiagonal_elements_by(const double & fac) {
        for (const auto & el : this->map) {
          auto && pair_type{el.first};
          auto && pos{el.second};
          if (pair_type[0] != pair_type[1]) {
            auto block{reference(&this->data[std::get<0>(pos)],
                                 std::get<1>(pos), std::get<2>(pos))};
            block *= fac;
          }
        }
      }

     private:
      /**
       * Functor to get a key from a map
       */
      struct RetrieveKey {
        template <typename T>
        typename T::first_type operator()(T keyValuePair) const {
          return keyValuePair.first;
        }
      };
    };

  }  // namespace internal
  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename Precision_t, size_t Order, size_t PropertyLayer,
            class Manager, typename Key>
  class BlockSparseProperty : public PropertyBase {
   public:
    using Parent = PropertyBase;
    using Manager_t = Manager;
    using traits = typename Manager::traits;

    using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic,
                                  Eigen::RowMajor>;
    using DenseRef_t = Eigen::Map<Dense_t>;
    using sizes_t = std::vector<size_t>;
    using Key_t = Key;
    using Keys_t = std::set<Key_t>;
    using keys_list_t = std::vector<std::set<Key_t>>;
    using InputData_t = internal::InternallySortedKeyMap<Key_t, Dense_t>;
    using Data_t = std::vector<InputData_t>;
    using DataOrders_t = std::array<std::vector<InputData_t>, Order>;

   protected:
    DataOrders_t values{};  //!< storage for properties

   public:
    //! constructor
    BlockSparseProperty(Manager_t & manager,
                        std::string metadata = "no metadata")
        : Parent{static_cast<StructureManagerBase &>(manager),
                 0,
                 0,
                 Order,
                 PropertyLayer,
                 metadata} {}

    //! Default constructor
    BlockSparseProperty() = delete;

    //! Copy constructor
    BlockSparseProperty(const BlockSparseProperty & other) = delete;

    //! Move constructor
    BlockSparseProperty(BlockSparseProperty && other) = default;

    //! Destructor
    virtual ~BlockSparseProperty() = default;

    //! Copy assignment operator
    BlockSparseProperty & operator=(const BlockSparseProperty & other) = delete;

    //! Move assignment operator
    BlockSparseProperty & operator=(BlockSparseProperty && other) = default;

    /* ---------------------------------------------------------------------- */
    //! return runtime info about the stored (e.g., numerical) type
    const std::type_info & get_type_info() const final {
      return typeid(Precision_t);
    };

    /**
     * the case consider_ghost_atoms == true is limited to cluster_index
     * objects.
     */
    template <size_t Order_ = Order, std::enable_if_t<(Order_ == 1), int> = 0>
    size_t get_validated_property_length(bool consider_ghost_atoms) {
      if (consider_ghost_atoms) {
        if (traits::MaxOrder < 2) {
          throw std::runtime_error(
              "consider_ghost_atoms is true,"
              " but can only be use for underlying manager with"
              " MaxOrder at least 2.");
        }
        if (not(this->get_manager().get_consider_ghost_neighbours())) {
          throw std::runtime_error(
              "consider_ghost_atoms is true,"
              " but underlying manager does not have ghost atoms in"
              " cluster_indices_container. Turn consider_ghost_neighbours"
              " on, to consider ghost atoms with independent property values"
              " from their corresponding central atoms.");
        }
        return this->get_manager().size_with_ghosts();
      }
      return this->get_manager().size();
    }
    template <size_t Order_ = Order,
              std::enable_if_t<not(Order_ == 1), int> = 0>
    size_t get_validated_property_length(bool /*consider_ghost_atoms*/ = false) {
      return this->base_manager.nb_clusters(Order);
    }

    //! Adjust size of values (only increases, never frees)
    template <size_t Order_ = Order, std::enable_if_t<(Order_ > 1), int> = 0>
    inline void resize(bool consider_ghost_atoms = false) {
      this->values.resize<Order_ - 1>(consider_ghost_atoms);
      size_t new_size{this->get_validated_property_length(consider_ghost_atoms)};
      this->values[Order_ - 1].resize(new_size);
    }

    template <size_t Order_ = Order, std::enable_if_t<(Order_ == 1), int> = 0>
    inline void resize(bool consider_ghost_atoms = false) {
      size_t new_size{this->get_validated_property_length(consider_ghost_atoms)};
      this->values[Order_ - 1].resize(new_size);
    }

    inline size_t size() const {
      size_t size{0};
      for (size_t i_order{0}; i_order < Order; ++i_order) {
        size += this->values[i_order].size();
      }
      return size;
    }

    //! clear all the content of the property
    inline void clear() {
      for (size_t i_order{0}; i_order < Order; ++i_order) {
        this->values[i_order].clear();
      }
    }

    inline Manager_t & get_manager() {
      return static_cast<Manager_t &>(this->base_manager);
    }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerOrder, size_t CallerLayer, size_t Order_ = Order,
              std::enable_if_t<(CallerOrder <= Order_), int> = 0>
    inline decltype(auto)
    operator[](const ClusterRefKey<CallerOrder, CallerLayer> & id) {
      static_assert(CallerOrder <= Order, "should be CallerOrder <= Order");
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator()(id.get_cluster_index(CallerLayer), CallerOrder-1);
    }

    /**
     * Access a property of order 1 with a clusterRef of order 2
     */
    template <size_t CallerOrder, size_t CallerLayer, size_t Order_ = Order,
        std::enable_if_t<(Order_ == 1) and (CallerOrder > 1), int> = // NOLINT
                  0>  // NOLINT
    inline decltype(auto)
    operator[](const ClusterRefKey<CallerOrder, CallerLayer> & id) {
      return this->operator()(this->get_manager().get_atom_index(
          id.get_internal_neighbour_atom_tag()), 0);
    }


    inline InputData_t & operator()(const size_t & index,
                                    const size_t & i_order) {
      return this->values[i_order][index];
    }

     //! Accessor for property by index for dynamically sized properties
    inline decltype(auto) operator()(const size_t & index) {
      return this->operator[](index);
    }

    template <size_t CallerOrder, size_t CallerLayer>
    inline decltype(auto)
    operator()(const ClusterRefKey<CallerOrder, CallerLayer> & id,
               const Key_t & key) {
      static_assert(CallerOrder <= Order, "should be CallerOrder <= Order");
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator()(id.get_cluster_index(CallerLayer), CallerOrder-1, key);
    }

    //! Accessor for property by index for dynamically sized properties
    inline DenseRef_t operator()(const size_t & index,
                                  const size_t & i_order,
                                  const Key_t & key) {
      return DenseRef_t(&this->values[i_order][index].at(key)(0, 0),
                         this->values[i_order][index].at(key).rows(),
                         this->values[i_order][index].at(key).cols());
    }

    // inline size_t get_dense_feature_size(const size_t & index) {
    //   auto keys{this->values[0][index].get_keys()};
    //   return this->get_nb_comp() * keys.size();
    // }

    //! Accessor for property by cluster index and return a dense
    //! representation of the property associated to this cluster
    template <size_t CallerOrder, size_t CallerLayer>
    inline Dense_t get_dense_row(const ClusterRefKey<CallerOrder, CallerLayer> & id) {
      static_assert(CallerOrder <= Order, "should be CallerOrder <= Order");
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->get_dense_row(id.get_cluster_index(CallerLayer), CallerOrder - 1);
    }

    inline Dense_t get_dense_row(const size_t & index, const size_t & i_order) {
      auto keys = this->values[i_order][index].get_keys();
      Dense_t feauture_row = Dense_t::Zero(this->get_nb_comp(), keys.size());
      size_t i_col{0};
      for (const auto & key : keys) {
        size_t i_row{0};
        for (int i_pos{0}; i_pos < this->values[i_order][index][key].size(); i_pos++) {
          feauture_row(i_row, i_col) = this->values[i_order][index][key](i_pos);
          i_row++;
        }
        i_col++;
      }
      return feauture_row;
    }

    inline Dense_t get_dense_rep() {
      size_t n_elements{this->size()};
      int inner_size{this->get_nb_comp()};
      Keys_t all_keys{};
      for (size_t i_order{0}; i_order < Order; ++i_order) {
        size_t n_center{this->values[i_order].size()};
        for (size_t i_center{0}; i_center < n_center; i_center++) {
          auto keys = this->values[i_order][i_center].get_keys();
          for (auto & key : keys) {
            all_keys.insert(key);
          }
        }
      }
      Dense_t features =
          Dense_t::Zero(inner_size * all_keys.size(), n_elements);

      for (size_t i_order{0}; i_order < Order; ++i_order) {
        size_t n_center{this->values[i_order].size()};
        for (size_t i_center{0}; i_center < n_center; i_center++) {
          int i_feat{0};
          for (const auto & key : all_keys) {
            if (this->values[i_order][i_center].count(key) == 1) {
              for (int i_pos{0}; i_pos < this->values[i_order][i_center][key].size();
                  i_pos++) {
                features(i_feat, i_center) = this->values[i_order][i_center][key](i_pos);
                i_feat++;
              }
            } else {
              i_feat += inner_size;
            }
          }
        }
      }
      return features;
    }

    //! getter to the underlying data storage
    // A hack here is fine because the function will disapear soon
    inline std::vector<InputData_t> & get_raw_data() { return this->values[0]; }
    //! get number of different distinct element in the property
    //! (typically the number of center)
    inline size_t get_nb_item() const {
      return this->size();
    }

    template <size_t CallerOrder, size_t CallerLayer>
    inline decltype(auto)
    get_keys(const ClusterRefKey<CallerOrder, CallerLayer> & id) const {
      static_assert(CallerOrder <= Order, "should be CallerOrder <= Order");
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");
      return this->values[CallerOrder - 1][id.get_cluster_index(CallerLayer)].get_keys();
    }

  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_BLOCK_SPARSE_HH_
