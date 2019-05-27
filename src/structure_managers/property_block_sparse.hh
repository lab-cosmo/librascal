/**
 * file   property_base.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   03 Aug 2018
 *
 * @brief implementation of non-templated base class for Properties, Properties
 *        are atom-, pair-, triplet-, etc-related values
 *
 * Copyright © 2018 Till Junge, Felix Musil, COSMO (EPFL), LAMMM (EPFL)
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

    template <class K, class V>
    class InternallySortedKeyMap {
     public:
      using MyMap_t = std::map<K, V>;
      // using Map_t = std::unordered_map<K, std::array<int, 3>, Hash<K>>;
      using Map_t = std::map<K, std::array<int, 3>>;
      using Data_t = std::vector<typename V::value_type>;
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
          return std::make_pair(key, Value_t(&data[pos[0]], pos[1], pos[2]));
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
        key_type skey{this->copy_sort(key)};
        auto & pos{this->map.at(skey)};
        return reference(&this->data[pos[0]], pos[1], pos[2]);
      }
      const_reference at(const key_type & key) const {
        key_type skey{this->copy_sort(key)};
        auto & pos{this->map.at(skey)};
        return const_reference(&this->data[pos[0]], pos[1], pos[2]);
      }
      //! access or insert specified element
      reference operator[](const key_type & key) {
        key_type skey{this->copy_sort(key)};
        auto & pos{this->map[skey]};
        return reference(&this->data[pos[0]], pos[1], pos[2]);
      }
      const_reference operator[](const key_type & key) const {
        key_type skey{this->copy_sort(key)};
        auto & pos{this->map[skey]};
        return const_reference(&this->data[pos[0]], pos[1], pos[2]);
      }
      template <typename Key_List>
      void resize(const Key_List & keys, const int & n_row, const int & n_col) {
        using Array_t = typename Map_t::mapped_type;
        int new_size{0};
        for (auto & key : keys) {
          key_type skey{this->copy_sort(key)};
          if (this->map.count(skey) == 0) {
            Array_t val{new_size, n_row, n_col};
            this->map[skey] = val;
            new_size += static_cast<int>(n_row * n_col);
          }
        }
        this->data.resize(new_size, 0.);
      }

      //! Returns the number of elements with key that compares equivalent to
      //! the specified argument, which is either 1 or 0 since this container
      //! does not allow duplicates.
      template <class Key>
      decltype(auto) count(const Key & key) {
        key_type skey{this->copy_sort(key)};
        return this->map.count(skey);
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

      key_type copy_sort(const key_type & key) {
        key_type skey{key};
        if (key.size() > 1) {
          std::sort(skey.begin(), skey.end());
        }
        return skey;
      }

      key_type copy_sort(key_type && key) {
        key_type skey{key};
        if (key.size() > 1) {
          std::sort(skey.begin(), skey.end());
        }
        return skey;
      }
    };

  }  // namespace internal
  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename Precision_t, size_t Order, size_t PropertyLayer,
            typename Key = std::vector<int>>
  class BlockSparseProperty : public PropertyBase {
   public:
    using Parent = PropertyBase;
    using Dense_t = Eigen::Matrix<Precision_t, Eigen::Dynamic, Eigen::Dynamic>;
    using dense_ref_t = Eigen::Map<Dense_t>;
    using sizes_t = std::vector<size_t>;
    using Key_t = Key;
    using Keys_t = std::set<Key_t>;
    using keys_list_t = std::vector<std::set<Key_t>>;
    using InputData_t = internal::InternallySortedKeyMap<Key_t, Dense_t>;
    using Data_t = std::vector<InputData_t>;

    //! constructor
    BlockSparseProperty(StructureManagerBase & manager,
                        std::string metadata = "no metadata")
        : Parent{manager, 0, 0, Order, PropertyLayer, metadata} {}

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

    //! Adjust size so that each center are accessible
    void resize() {
      auto order = this->get_order();
      auto new_size = this->base_manager.nb_clusters(order);
      this->values.resize(new_size);
      this->center_sizes.resize(new_size);
    }

    template <size_t CallerLayer>
    void initialize_to_zeros(const ClusterRefKey<Order, CallerLayer> & id,
                             Key_t & key) {
      auto && index{id.get_cluster_index(CallerLayer)};
      Dense_t mat = Dense_t::Zero(this->get_nb_row(), this->get_nb_col());
      this->values[index][key] = mat;
    }

    size_t size() const { return this->values.size(); }

    /**
     * clear all the content of the property
     */
    void clear() {
      this->values.clear();
      this->center_sizes.clear();
    }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    inline decltype(auto)
    operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    //! Accessor for property by cluster index and return a sparse
    //! representation of the property associated to this cluster
    inline InputData_t & operator[](const size_t & index) {
      return this->values[index];
    }

    inline size_t get_dense_feature_size(const size_t & index) {
      auto keys = this->values[index].get_keys();
      return this->get_nb_comp() * keys.size();
    }

    //! Accessor for property by cluster index and return a dense
    //! representation of the property associated to this cluster
    template <size_t CallerLayer>
    inline Dense_t get_dense_row(const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->get_dense_row(id.get_cluster_index(CallerLayer));
    }

    inline Dense_t get_dense_row(const size_t & index) {
      auto keys = this->values[index].get_keys();
      Dense_t feauture_row = Dense_t::Zero(this->get_nb_comp(), keys.size());
      size_t i_col{0};
      for (const auto & key : keys) {
        size_t i_row{0};
        for (int i_pos{0}; i_pos < this->values[index][key].size(); i_pos++) {
          feauture_row(i_row, i_col) = this->values[index][key](i_pos);
          i_row++;
        }
        i_col++;
      }
      return feauture_row;
    }

    inline Dense_t get_dense_rep() {
      auto n_center{this->get_nb_item()};
      Keys_t all_keys{};
      for (size_t i_center{0}; i_center < n_center; i_center++) {
        auto keys = this->values[i_center].get_keys();
        for (auto & key : keys) {
          all_keys.insert(key);
        }
      }
      Dense_t features =
          Dense_t::Zero(this->get_nb_comp() * all_keys.size(), n_center);

      for (size_t i_center{0}; i_center < n_center; i_center++) {
        int i_feat{0};
        for (const auto & key : all_keys) {
          if (this->values[i_center].count(key) == 1) {
            for (int i_pos{0}; i_pos < this->values[i_center][key].size();
                 i_pos++) {
              features(i_feat, i_center) = this->values[i_center][key](i_pos);
              i_feat++;
            }
          }
        }
      }
      return features;
    }

    //! Accessor for property by index for dynamically sized properties
    inline decltype(auto) operator()(const size_t & index) {
      return this->operator[](index);
    }

    template <size_t CallerLayer>
    inline decltype(auto)
    operator()(const ClusterRefKey<Order, CallerLayer> & id,
               const Key_t & key) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator()(id.get_cluster_index(CallerLayer), key);
    }

    //! Accessor for property by index for dynamically sized properties
    inline dense_ref_t operator()(const size_t & index, const Key_t & key) {
      return dense_ref_t(&this->values[index].at(key)(0, 0),
                         this->values[index].at(key).rows(),
                         this->values[index].at(key).cols());
    }

    //! getter to the underlying data storage
    inline Data_t & get_raw_data() { return this->values; }
    //! get number of different distinct element in the property
    //! (typically the number of center)
    inline size_t get_nb_item() const { return this->values.size(); }

    /**
     * Accessor for last pushed entry for dynamically sized properties
     */
    inline decltype(auto) back() { return this->values.back(); }

    //! push back data associated to a new center atom
    inline void push_back(const InputData_t & ref) {
      for (const auto & element : ref) {
        const auto & value{element.second};
        if (value.size() != this->get_nb_comp()) {
          auto error{std::string("Size should match: ") +
                     std::to_string(value.size()) + std::string(" != ") +
                     std::to_string(this->get_nb_comp())};
          throw std::length_error(error);
        }
      }

      this->values.push_back(ref);

      size_t n_keys{0};
      for (const auto & element : ref) {
        n_keys++;
      }
      this->center_sizes.push_back(n_keys * this->get_nb_comp());
    }

    template <size_t CallerLayer>
    inline decltype(auto)
    get_keys(const ClusterRefKey<Order, CallerLayer> & id) const {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");
      return this->values[id.get_cluster_index(CallerLayer)].get_keys();
    }

   protected:
    Data_t values{};  //!< storage for properties
    sizes_t center_sizes{};
  };

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_BLOCK_SPARSE_HH_
