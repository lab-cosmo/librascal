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


namespace rascal {



  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename precision_t, typename key_t, size_t Order, size_t PropertyLayer>
  class PartiallySparseProperty : public PropertyBase {
    using Parent = PropertyBase;
    using dense_t = Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic>;
    using dense_ref_t = Eigen::Map<Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic>>;
    using sparse_t = std::unordered_map<key_t, dense_t>;
    using type = std::vector<sparse_t>;

    // using Value = internal::Value<T, Eigen::Dynamic, Eigen::Dynamic>;

   public:
    // using value_type = typename Value::type;
    // using reference = typename Value::reference;

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
    void resize() {
      auto order = this->get_order();
      auto new_size = this->base_manager.nb_clusters(order);
      this->values.resize(new_size);
    }

    //! Adjust size of values (only increases, never frees)
    size_t size() const { return this->values.size(); }

    /**
     * shortens the vector so that the manager can push_back into it (capacity
     * not reduced)
     */
    void resize_to_zero() { this->values.resize(0); }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    inline sparse_t operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    //! Accessor for property by index for dynamically sized properties
    sparse_t& operator[](const size_t & index) {
      return this->values[index];
    }

    //! Accessor for property by index for dynamically sized properties
    sparse_t& operator()(const size_t & index) {
      return this->values[index];
    }

    //! Accessor for property by index for dynamically sized properties
    dense_ref_t& operator()(const size_t & index, const key_t & key) {
      return dense_ref_t(this->values[index][key], , );
    }

    //! getter to the underlying data storage
    inline type& get_raw_data() { return this->values; }

    //! get number of different distinct element in the property
    //! (typically the number of center)
    inline size_t get_nb_item() const {
      return values.size();
    }

    /**
     * Accessor for last pushed entry for dynamically sized properties
     */
    sparse_t& back() {
      return this->values.back();
    }

    inline void push_back(sparse_t& ref) {
      this->values.push_back(ref);
    }




   protected:
    type values{};  //!< storage for properties
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_PARTIALLY_SPARSE_HH_
