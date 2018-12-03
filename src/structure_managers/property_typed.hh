/**
 * file   property_typed.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   06 Aug 2018
 *
 * @brief Implements intermediate property class for which the type of stored
 *          objects is known, but not the size
 *
 * Copyright © 2018 Federico Giberti, Till Junge, COSMO (EPFL), LAMMM (EPFL)
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

#ifndef PROPERTY_TYPED_H
#define PROPERTY_TYPED_H

#include "structure_managers/property_base.hh"
#include "structure_managers/cluster_ref_key.hh"

namespace rascal {

  /* ---------------------------------------------------------------------- */
  namespace internal {

    /**
     * Structure providing access to a ``property`` and also access data stored
     * in a property
     */
    template <typename T, Dim_t NbRow, Dim_t NbCol>
    struct Value {
      using type = Eigen::Map<Eigen::Matrix<T, NbRow, NbCol>>;
      using reference = type;

      //! get a reference to specific value at row and colum
      static reference get_ref(T & value, int nb_row, int nb_col) {
        return type(&value, nb_row, nb_col);
      }

      //! get a reference
      static reference get_ref(T & value) {
        return type(&value);
      }

      //! push back data into ``property``
      static void push_in_vector(std::vector<T> & vec, reference ref) {
        for (size_t j{0}; j < NbCol; ++j) {
          for (size_t i{0}; i < NbRow; ++i) {
            vec.push_back(ref(i, j));
          }
        }
      }

      //! Used for extending cluster_indices
      template<typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref) {
        static_assert(Derived::RowsAtCompileTime == NbRow,
                      "NbRow has incorrect size.");
        static_assert(Derived::ColsAtCompileTime == NbCol,
                      "NbCol has incorrect size.");
        for (size_t j{0}; j < NbCol; ++j) {
          for (size_t i{0}; i < NbRow; ++i) {
            vec.push_back(ref(i, j));
          }
        }
      }
    };

    /* ---------------------------------------------------------------------- */
    //! specialisation for scalar properties
    template <typename T>
    struct Value<T, 1, 1> {
      constexpr static Dim_t NbRow{1};
      constexpr static Dim_t NbCol{1};
      using type = T;
      using reference = T&;

      //! get a reference to a scalar value
      static reference get_ref(T & value) {return value;}

      //! push a scalar in a vector
      static void push_in_vector(std::vector<T> & vec, reference ref) {
        vec.push_back(ref);
      }

      //! Used for extending cluster_indices
      template<typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref) {
        static_assert(Derived::RowsAtCompileTime == NbRow,
                      "NbRow has incorrect size.");
        static_assert(Derived::ColsAtCompileTime == NbCol,
                      "NbCol has incorrect size.");
        vec.push_back(ref(0, 0));
      }
    };

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_t = typename Value<T, NbRow, NbCol>::type;

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_ref = typename Value<T, NbRow, NbCol>::reference;

  }  // internal

  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename T, size_t Order, size_t PropertyLayer>
  class TypedProperty: public PropertyBase {
    using Parent = PropertyBase;
    using Value = internal::Value<T, Eigen::Dynamic, Eigen::Dynamic>;

   public:
    using value_type = typename Value::type;
    using reference = typename Value::reference;

    //! constructor
    TypedProperty(StructureManagerBase & manager, Dim_t nb_row,
                  Dim_t nb_col = 1, std::string metadata = "no metadata"):
      Parent{manager, nb_row, nb_col, Order, PropertyLayer, metadata}
    {}

    //! Default constructor
    TypedProperty() = delete;

    //! Copy constructor
    TypedProperty(const TypedProperty & other) = delete;

    //! Move constructor
    TypedProperty(TypedProperty && other) = default;

    //! Destructor
    virtual ~TypedProperty() = default;

    //! Copy assignment operator
    TypedProperty & operator=(const TypedProperty & other) = delete;

    //! Move assignment operator
    TypedProperty & operator=(TypedProperty && other) = default;

    /* ---------------------------------------------------------------------- */
    //! return runtime info about the stored (e.g., numerical) type
    const std::type_info & get_type_info() const override final {
      return typeid(T);
    };

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

    /**
     * shortens the vector so that the manager can push_back into it (capacity
     * not reduced)
     */
    void resize_to_zero() {
      this->values.resize(0);
    }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template<size_t CallerLayer>
    inline reference operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    //! Accessor for property by index for dynamically sized properties
    reference operator[](const size_t & index) {
      return Value::get_ref(this->values[index * this->get_nb_comp()],
                            this->get_nb_row(),
                            this->get_nb_col());
    }

   protected:
    std::vector<T> values{}; //!< storage for properties
  };

}  // rascal

#endif /* PROPERTY_TYPED_H */
