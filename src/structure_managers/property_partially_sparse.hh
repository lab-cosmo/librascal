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

#ifndef SRC_STRUCTURE_MANAGERS_PROPERTY_PARTIALLY_SPARSE_HH_
#define SRC_STRUCTURE_MANAGERS_PROPERTY_PARTIALLY_SPARSE_HH_

#include "structure_managers/property_base.hh"
#include "structure_managers/cluster_ref_key.hh"

#include <unordered_map>
#include <set>
#include <list>

namespace rascal {



  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename precision_t, typename key_t, size_t Order, size_t PropertyLayer>
  class PartiallySparseProperty : public PropertyBase {
   public:
    using Parent = PropertyBase;
    using dense_t = Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic>;
    using dense_ref_t = Eigen::Map<dense_t>;
    using map_center_t = std::vector<std::pair<size_t, size_t>>;
    using map_sparse_t = std::vector<std::unordered_map<key_t, std::pair<size_t, std::pair<size_t, size_t>>>>;
    using keys_t = std::vector<std::list<key_t>>;
    using data_t = std::vector<precision_t>;
    using input_data_t = std::unordered_map<key_t, dense_t>;

    //! constructor
    PartiallySparseProperty(StructureManagerBase & manager, std::string metadata = "no metadata")
        : Parent{manager, 0, 0, Order, PropertyLayer, metadata} {}

    //! Default constructor
    PartiallySparseProperty() = delete;

    //! Copy constructor
    PartiallySparseProperty(const PartiallySparseProperty & other) = delete;

    //! Move constructor
    PartiallySparseProperty(PartiallySparseProperty && other) = default;

    //! Destructor
    virtual ~PartiallySparseProperty() = default;

    //! Copy assignment operator
    PartiallySparseProperty & operator=(const PartiallySparseProperty & other) = delete;

    //! Move assignment operator
    PartiallySparseProperty & operator=(PartiallySparseProperty && other) = default;

    /* ---------------------------------------------------------------------- */
    //! return runtime info about the stored (e.g., numerical) type
    const std::type_info & get_type_info() const final { return typeid(precision_t); };


    //! Adjust size of values (only increases, never frees)
    // void resize() {
    //   auto order = this->get_order();
    //   auto new_size = this->base_manager.nb_clusters(order);
    //   this->values.resize(new_size);
    // }

    size_t size() const { return this->map2center.size(); }

    /**
     * clear all the content of the property
     */
    void clear() {
      this->values.clear();
      this->map2centers.clear();
      this->map2sparse.clear();
      this->keys_list.clear();

    }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    inline decltype(auto) operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    //! Accessor for property by index for dynamically sized properties
    inline dense_ref_t operator[](const size_t & index) {
      auto && start_id{this->map2centers[index].first};
      auto && length{this->map2centers[index].second};
      return dense_ref_t(&this->values[start_id], length, 1);
    }

    //! Accessor for property by index for dynamically sized properties
    inline decltype(auto) operator()(const size_t & index) {
      return this->operator[](index);
    }

    template <size_t CallerLayer>
    inline decltype(auto) operator()(const ClusterRefKey<Order, CallerLayer> & id, const key_t & key) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator()(id.get_cluster_index(CallerLayer), key);
    }

    //! Accessor for property by index for dynamically sized properties
    inline dense_ref_t operator()(const size_t & index, const key_t & key) {
      auto && center_start_id{this->map2centers[index].first};
      auto && sparse_start_id{center_start_id + this->map2sparse[index][key].first};
      auto && sparse_shape{this->map2sparse[index][key].second};
      return dense_ref_t(&this->values[sparse_start_id], sparse_shape.first, sparse_shape.second);
    }

    //! getter to the underlying data storage
    inline data_t& get_raw_data() { return this->values; }

    //! get number of different distinct element in the property
    //! (typically the number of center)
    inline size_t get_nb_item() const {
      return this->map2centers.size();
    }

    /**
     * Accessor for last pushed entry for dynamically sized properties
     */
    inline decltype(auto) back() {
      auto && index{this->map2centers.size()-1};
      return this->operator[](index);
    }

    //! push back data associated to a new center atom
    inline void push_back(input_data_t& ref) {
      // get where the new center starts
      auto && new_center_start_id{this->values.size()};
      this->map2sparse.emplace_back();
      this->keys_list.emplace_back();
      size_t key_start{0};
      for(const auto& element : ref) {
        auto&& key{element.first};
        auto&& value{element.second};

        this->keys_list.back().emplace_back(key);
        this->map2sparse.back().emplace(
          std::make_pair(key, std::make_pair(key_start, std::make_pair(value.rows(), value.cols())))
        );

        key_start += value.size();

        for (int jj{0}; jj < value.cols(); ++jj) {
          for (int ii{0}; ii < value.rows(); ++ii) {
            this->values.push_back(value(ii, jj));
          }
        }
      }

      auto center_length{key_start};
      map2centers.push_back(
        std::make_pair(new_center_start_id, center_length)
      );
    }

    template <size_t CallerLayer>
    inline decltype(auto) get_keys(const ClusterRefKey<Order, CallerLayer> & id) const {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");
      return this->keys_list[id.get_cluster_index(CallerLayer)];
    }


   protected:
    data_t values{};  //!< storage for properties
    map_center_t map2centers{};
    map_sparse_t map2sparse{};
    keys_t keys_list{};
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_PARTIALLY_SPARSE_HH_
