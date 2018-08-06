/**
 * file   property.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Properties of atom-, pair-, triplet-, etc-related values
 *
 * @section LICENSE
 *
 * Copyright © 2018 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef NEIGHBOURHOOD_PROPERTY_H
#define NEIGHBOURHOOD_PROPERTY_H

#include <Eigen/Dense>
#include <type_traits>
#include <vector>
#include <basic_types.hh>

#include "structure_managers/property_typed.hh"
#include "structure_managers/cluster_ref_base.hh"

namespace rascal {


  // Forward declaration of traits to use `Property`.
  template <class Manager>
  struct StructureManager_traits;


  template <class StructureManager, typename T,
            size_t Order,
            Dim_t NbRow = 1, Dim_t NbCol = 1,
            size_t ActiveLayer=compute_cluster_layer<Order>(typename StructureManager::traits::LayerByDimension{})>
  class Property: public TypedProperty<T>
  {
    static_assert((std::is_arithmetic<T>::value or
                   std::is_same<T, std::complex<double>>::value),
                  "can currently only handle arithmetic types");
  public:
    using traits = StructureManager_traits<StructureManager>;

    using Parent = TypedProperty<T>;
    constexpr static size_t NbComp{NbRow*NbCol};

    using Value = internal::Value<T, NbRow, NbCol>;
    static_assert(std::is_same<Value, internal::Value<T, NbRow, NbCol>>::value,
                  "type alias failed");

    using value_type = typename Value::type;
    using reference = typename Value::reference;

    static constexpr bool IsStaticallySized{(NbCol != Eigen::Dynamic)
        and (NbRow != Eigen::Dynamic)};

    //! Default constructor
    Property() = delete;

    //! Constructor with Manager
    Property(StructureManager & manager)
      :Parent{manager, NbRow, NbCol, Order}
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

    /**
     * allows to add a value to `Property` during construction of the
     * neighbourhood.
     */
    inline void push_back(reference ref) {
      Value::push_in_vector(this->values, ref);
    }

    /**
     * Function for adding Eigen-based matrix data to `property`
     */
    template<typename Derived>
    inline void
    push_back(const Eigen::DenseBase<Derived> & ref) {
      static_assert(Derived::RowsAtCompileTime==NbRow,
                    "NbRow has incorrect size.");
      static_assert(Derived::ColsAtCompileTime==NbCol,
                    "NbCol has incorrect size.");

      Value::push_in_vector(this->values, ref);
    }

    /**
     * Not sure about the actual implementation of this one.
     */
    template<size_t CallerLayer>
    reference operator[](const ClusterRefBase<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= ActiveLayer,
                    "You are trying to access a property that "
                    "does not exist at this depth in the "
                    "adaptor stack.");

      return
        Value::get_ref(this->values[id.get_cluster_index(CallerLayer)*NbComp]);
    }

    /**
     * Accessor for property by index for statically sized properties
     */
    reference operator[](const size_t & index) {
      return Value::get_ref(this->values[index*NbComp]);
    }

  protected:
  private:
  };

}  // rascal

#endif /* NEIGHBOURHOOD_PROPERTY_H */
