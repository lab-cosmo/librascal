/**
 * file   property.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Properties of atom-, pair-, triplet-, etc-related values
 *
 * @section LICENSE
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

#ifndef NEIGHBOURHOOD_PROPERTY_H
#define NEIGHBOURHOOD_PROPERTY_H

#include "structure_managers/property_typed.hh"
#include <basic_types.hh>
#include <cassert>
#include <Eigen/Dense>

#include <type_traits>

namespace rascal {

  //! Forward declaration of traits to use `Property`.
  template <class Manager>
  struct StructureManager_traits;

  /**
   * Class definition of `property`. This is a container of data (specifiable),
   * which can be access with clusters directly, without the need for dealing
   * with indices.
   */
  template <typename T,
            size_t Order,
            size_t PropertyLayer,
            Dim_t NbRow = 1, Dim_t NbCol = 1>
  class Property: public TypedProperty<T, Order, PropertyLayer> {
    static_assert((std::is_arithmetic<T>::value ||
                   std::is_same<T, std::complex<double>>::value),
                  "can currently only handle arithmetic types");

   public:
    using Parent = TypedProperty<T, Order, PropertyLayer>;

    using Value = internal::Value<T, NbRow, NbCol>;
    static_assert(std::is_same<Value, internal::Value<T, NbRow, NbCol>>::value,
                  "type alias failed");

    using value_type = typename Value::type;
    using reference = typename Value::reference;

    static constexpr bool
    IsStaticallySized{ (NbCol != Eigen::Dynamic) and
                       (NbRow != Eigen::Dynamic) };

    //! Empty type for tag dispatching to differenciate between
    //! the Dynamic and Static size case
    struct DynamicSize {};
    struct StaticSize {};

    //! Default constructor
    Property() = delete;

    //! Constructor with Manager
    Property(StructureManagerBase & manager,
              std::string metadata = "no metadata")
      :Parent{manager, NbRow, NbCol, metadata}
    {}

    //! Copy constructor
    Property(const Property & other) = delete;

    //! Move constructor
    Property(Property && other) = default;

    //! Destructor
    virtual ~Property() = default;

    //! Copy assignment operator
    Property & operator=(const Property & other) = delete;

    //! Move assignment operator
    Property & operator=(Property && other) = delete;

    /* ---------------------------------------------------------------------- */
    /**
     * Cast operator: takes polymorphic base class reference, and returns
     * properly casted fully typed and sized reference, or throws a runttime
     * error
     */
    // TODO(felix) Need to make an equivalent for dynamic sized property
    static inline Property & check_compatibility(PropertyBase & other) {
      // check ``type`` compatibility
      if (not(other.get_type_info().hash_code() == typeid(T).hash_code())) {
        std::stringstream err_str{};
        err_str << "Incompatible types: '" << other.get_type_info().name()
                << "' != '" << typeid(T).name() << "'.";
        throw std::runtime_error(err_str.str());
      }

      // check ``order`` compatibility
      if (not(other.get_order() == Order)) {
        std::stringstream err_str{};
        err_str << "Incompatible property order: input is of order "
                << other.get_order() << ", this property is of order "
                << Order << ".";
        throw std::runtime_error(err_str.str());
      }

      // check property ``layer`` compatibility
      if (not(other.get_property_layer() == PropertyLayer)) {
        std::stringstream err_str{};
        err_str << "At wrong layer in stack: input is at layer "
                << other.get_property_layer() << ", this property is at layer "
                << PropertyLayer << ".";
        throw std::runtime_error(err_str.str());
      }

      // check size compatibility
      if (not((other.get_nb_row() == NbRow) and
               (other.get_nb_col() == NbCol))) {
        std::stringstream err_str{};
        err_str << "Incompatible sizes: input is " << other.get_nb_row() << "×"
                << other.get_nb_col() << ", but should be " << NbRow << "×"
                << NbCol  << ".";
        throw std::runtime_error(err_str.str());
      }
      return static_cast<Property& > (other);
    }

    /* ---------------------------------------------------------------------- */
    /**
     * allows to add a value to `Property` during construction of the
     * neighbourhood.
     */
    inline void push_back(reference ref) {
      // use tag dispatch to use the proper definition
      // of the push_in_vector function
      this->push_back(ref,
          std::conditional_t<(IsStaticallySized),
                              StaticSize, DynamicSize>{});
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Function for adding Eigen-based matrix data to `property`
     */
    template<typename Derived>
    inline void push_back(const Eigen::DenseBase<Derived> & ref) {
      // use tag dispatch to use the proper definition
      // of the push_in_vector function
      this->push_back(ref,
          std::conditional_t<(IsStaticallySized),
                              StaticSize, DynamicSize>{});
    }

    /* ---------------------------------------------------------------------- */
    /**
     * Property accessor by cluster ref
     */
    template<size_t CallerLayer>
    inline reference operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that "
                    "does not exist at this depth in the "
                    "adaptor stack.");

      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    /**
     * Accessor for property by index for properties
     */
    inline reference operator[](const size_t & index) {
      // use tag dispatch to use the proper definition
      // of the get function
      return this->get(index,
                        std::conditional_t<(IsStaticallySized),
                          StaticSize, DynamicSize>{});
    }

   protected:
    inline void push_back(reference ref, StaticSize) {
        Value::push_in_vector(this->values, ref);
    }

    inline void push_back(reference ref, DynamicSize) {
      Value::push_in_vector(this->values, ref,
                            this->get_nb_row(), this->get_nb_col());
    }

    template<typename Derived>
    inline void push_back(const Eigen::DenseBase<Derived> & ref, StaticSize) {
      static_assert(Derived::RowsAtCompileTime == NbRow,
                      "NbRow has incorrect size.");
      static_assert(Derived::ColsAtCompileTime == NbCol,
                      "NbCol has incorrect size.");
      Value::push_in_vector(this->values, ref);
    }
    template<typename Derived>
    inline void push_back(const Eigen::DenseBase<Derived> & ref, DynamicSize) {
      Value::push_in_vector(this->values, ref,
                              this->get_nb_row(), this->get_nb_col());
    }

    inline reference get(const size_t & index, StaticSize) {
      return Value::get_ref(this->values[index*NbRow*NbCol]);
    }

    inline reference get(const size_t & index, DynamicSize) {
      return get_ref(this->values[index*this->get_nb_comp()]);
    }

    //! get a reference
    reference get_ref(T & value) {
      return reference(&value, this->get_nb_row(), this->get_nb_col());
    }
  };

}  // rascal

#endif /* NEIGHBOURHOOD_PROPERTY_H */
