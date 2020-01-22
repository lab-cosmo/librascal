/**
 * @file   property_lookup_keys.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   22 Jan 2020
 *
 * @brief  Property for storing LookupkeysRefKey entries
 *
 * Copyright © 2020 Till Junge
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
 * along with Rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_RASCAL_STRUCTURE_MANAGERS_PROPERTY_LOOKUP_KEYS_HH_
#define SRC_RASCAL_STRUCTURE_MANAGERS_PROPERTY_LOOKUP_KEYS_HH_

#include "rascal/structure_managers/cluster_ref_key.hh"
#include "rascal/structure_managers/property_base.hh"

#include <array>
#include <vector>

namespace rascal {

  template <class StoredClusterRefKey, size_t Order, class Manager,
            size_t NbKeys>
  class PropertyLookupKeys : public PropertyBase {
   public:
    static_assert(Order > 1, "Current implementation does not handle atom "
                             "properties in any sensible way");
    using Parent = PropertyBase;
    using Manager_t = Manager;
    using value = std::array<StoredClusterRefKey, NbKeys>;
    using reference = value &;
    using const_reference = const reference;
    constexpr static size_t StoredOrder{value::order()};

    //! Default constructor
    PropertyLookupKeys() = delete;

    //! constructor
    PropertyLookupKeys(Manager & manager,
                       const std::string & metadata = "no metadata",
                       const bool  /*dummy_exclude_ghosts*/ = false)
        : Parent{manager,
                 NbKeys,
                 1,
                 Order,
                 manager.template cluster_layer_from_order<Order>(),
                 metadata} {}

    //! Copy constructor
    PropertyLookupKeys(const PropertyLookupKeys & other) = delete;

    //! Move constructor
    PropertyLookupKeys(PropertyLookupKeys && other) = default;

    //! Destructor
    virtual ~PropertyLookupKeys() = default;

    //! Copy assignment operator
    PropertyLookupKeys & operator=(const PropertyLookupKeys & other) = delete;

    //! Move assignment operator
    PropertyLookupKeys & operator=(PropertyLookupKeys && other) = default;

    const std::string & get_type_info() const {
      return this->type_name;
    }

    Manager_t & get_manager() {
      return static_cast<Manager_t &>(this->base_manager);
    }

    /**
     * Adjust size of values (only increases, never frees).
     *
     * Uses SFINAE to differenciate behavior between values of Order.
     * Order > 1 then the size of the property is directly taken
     * from the manager.
     * Order == 0 then the size of the property is 1.
     * Order == 1 then the size of the property depends on the bolean
     * exclude_ghosts. By default it is set to false so that property size is
     * always larger than what could be needed.
     */

    //! Adjust size of values (only increases, never frees).
    void resize() {
      auto && new_size{this->base_manager.nb_clusters(Order)};
      this->values.resize(new_size * this->get_nb_comp());
    }

    //! Returns the size of one component
    size_t size() const { return this->values.size(); }

    /**
     * shortens the vector so that the manager can push_back into it (capacity
     * not reduced)
     */
    void clear() { this->values.clear(); }

    template<class StoredClusterRefKeyArray>
    void push_back(StoredClusterRefKeyArray && val) {
      this->values.push_back(std::forward<StoredClusterRefKeyArray>(val));
    }

    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    reference operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      // You are trying to access a property that does not exist at this depth
      // in the adaptor stack.
      assert(static_cast<int>(CallerLayer) >= this->get_property_layer());
      return this->operator[](id.get_cluster_index(this->get_property_layer()));
    }

    //! Accessor for property by index for dynamically sized properties
    reference operator[](const size_t & index) { return this->values[index]; }

    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    const_reference
    operator[](const ClusterRefKey<Order, CallerLayer> & id) const {
      // You are trying to access a property that does not exist at this depth
      // in the adaptor stack.
      assert(static_cast<int>(CallerLayer) >= this->get_property_layer());
      return this->operator[](id.get_cluster_index(this->get_property_layer()));
    }

    //! Accessor for property by index for dynamically sized properties
    const_reference operator[](const size_t & index) const  {
      return this->values[index];
    }

   protected:
    std::vector<std::array<StoredClusterRefKey, NbKeys>> values{};
    std::string type_name{typeid(PropertyLookupKeys).name()};
  };
}  // namespace rascal

#endif /* SRC_RASCAL_STRUCTURE_MANAGERS_PROPERTY_LOOKUP_KEYS_HH_ */
