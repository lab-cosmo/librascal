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
    using dense_t = Eigen::Map<Eigen::Matrix<precision_t, Eigen::Dynamic, Eigen::Dynamic>>;
    using dense_ref_t = dense_t;
    using full_t = std::unordered_map<std::vector<precision_t>>;

    // using Value = internal::Value<T, Eigen::Dynamic, Eigen::Dynamic>;

   public:
    // using value_type = typename Value::type;
    // using reference = typename Value::reference;

    //! constructor
    PartiallySparseProperty(StructureManagerBase & manager, Dim_t nb_row,
                  Dim_t nb_col = 1, std::string metadata = "no metadata")
        : Parent{manager, nb_row, nb_col, Order, PropertyLayer, metadata} {}

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

    //! Fill sequence, used for *_cluster_indices initialization
    inline void fill_sequence() {
      this->resize();
      for (size_t i{0}; i < this->values.size(); ++i) {
        values[i] = i;
      }
    }

    //! Adjust size of values (only increases, never frees)
    void resize() {
      auto order = this->get_order();
      auto n_components = this->get_nb_comp();
      auto new_size = this->base_manager.nb_clusters(order) * n_components;
      this->values.resize(new_size);
    }

    //! Adjust size of values (only increases, never frees)
    size_t size() const { return this->values.size() / this->get_nb_comp(); }

    /**
     * shortens the vector so that the manager can push_back into it (capacity
     * not reduced)
     */
    void resize_to_zero() { this->values.resize(0); }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    inline reference operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    //! Accessor for property by index for dynamically sized properties
    reference operator[](const size_t & index) {
      return Value::get_ref(this->values[index * this->get_nb_comp()],
                            this->get_nb_row(), this->get_nb_col());
    }

    //! getter to the underlying data storage
    inline std::vector<precision_t> & get_raw_data() { return this->values; }

    //! get number of different distinct element in the property
    //! (typically the number of center)
    inline size_t get_nb_item() const {
      return values.size() / this->get_nb_comp();
    }

    /**
     * Accessor for last pushed entry for dynamically sized properties
     */
    reference back() {
      auto && index{this->values.size() - this->get_nb_comp()};
      return Value::get_ref(this->values[index * this->get_nb_comp()],
                            this->get_nb_row(), this->get_nb_col());
    }

   protected:
    std::vector<precision_t> values{};  //!< storage for properties
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_PARTIALLY_SPARSE_HH_
