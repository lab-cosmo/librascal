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

namespace rascal {

  namespace internal {

    template <typename T, Dim_t NbRow, Dim_t NbCol>
    struct Value {
      using type = Eigen::Map<Eigen::Matrix<T, NbRow, NbCol>>;
      using reference = type;

      static reference get_ref(T & value, int nb_row, int nb_col) {
        return type(&value, nb_row, nb_col);
      }

      static reference get_ref(T & value) {
        return type(&value);
      }

      static void push_in_vector(std::vector<T> & vec, reference ref) {
        for (size_t j{0}; j < NbCol; ++j) {
          for (size_t i{0}; i < NbRow; ++i) {
            vec.push_back(ref(i,j));
          }
        }
      }

      // Used for extending cluster_indices
      template<typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref) {
        static_assert(Derived::RowsAtCompileTime==NbRow,
                      "NbRow has incorrect size.");
        static_assert(Derived::ColsAtCompileTime==NbCol,
                      "NbCol has incorrect size.");
        for (size_t j{0}; j < NbCol; ++j) {
          for (size_t i{0}; i < NbRow; ++i) {
            vec.push_back(ref(i,j));
          }
        }
      }
    };

    //! specialisation for scalar properties
    template <typename T>
    struct Value<T, 1, 1> {
      constexpr static Dim_t NbRow{1};
      constexpr static Dim_t NbCol{1};
      using type = T;
      using reference = T&;

      static reference get_ref(T & value) {return value;}

      static void push_in_vector(std::vector<T> & vec, reference ref) {
        vec.push_back(ref);
      }

      // Used for extending cluster_indices
      template<typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref) {
        static_assert(Derived::RowsAtCompileTime==NbRow,
                      "NbRow has incorrect size.");
        static_assert(Derived::ColsAtCompileTime==NbCol,
                      "NbCol has incorrect size.");
        vec.push_back(ref(0,0));
      }
    };

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_t = typename Value<T, NbRow, NbCol>::type;

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_ref = typename Value<T, NbRow, NbCol>::reference;

  }  // internal

  // Forward declaration of traits to use `Property`.
  template <class Manager>
  struct NeighbourhoodManager_traits;


  template <class NeighbourhoodManager, typename T,
            size_t Level,
            Dim_t NbRow = 1, Dim_t NbCol = 1>
  class Property
  {
    static_assert((std::is_arithmetic<T>::value or
                   std::is_same<T, std::complex<double>>::value),
                  "can currently only handle arithmetic types");
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManager>;
    constexpr static size_t NbComp{NbRow*NbCol};

    using Value = internal::Value<T, NbRow, NbCol>;

    using value_type = typename Value::type;
    using reference = typename Value::reference;

    static constexpr bool IsStaticallySized{(NbCol != Eigen::Dynamic)
        and (NbRow != Eigen::Dynamic)};

    //! Default constructor
    Property() = delete;

    //! Constructor with Manager
    Property(NeighbourhoodManager & manager)
      :manager{manager}, values{}
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

    //! Adjust size of values (only increases, never frees)
    void resize() {
      this->values.resize(this->manager.nb_clusters(Level) * NbComp);
    }

    /**
     * shortens the vector so that the manager can push_back into it
     * (capacity not reduced)
     */
    void resize_to_zero() {
      this->values.resize(0);
    }

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
     * Fill sequence for *_cluster_indices to initialize
     */
    inline void fill_sequence() {
      this->resize();
      for (size_t i{0}; i<this->values.size(); ++i) {
        values[i] = i;
      }
    }

    /**
     * Not sure about the actual implementation of this one.
     */
    template<size_t CallerDepth>
    reference operator[](const ClusterRefBase<Level, CallerDepth> & id) {
      constexpr auto ActiveDepth{
        compute_cluster_depth<Level>(typename traits::DepthByDimension{})};
      static_assert(CallerDepth >= ActiveDepth,
                    "You are trying to access a property that "
                    "does not exist at this low a level in the "
                    "adaptor stack.");

      return
        Value::get_ref(this->values[id.get_cluster_index(CallerDepth)*NbComp]);
    }

    /**
     * Accessor for property by index for statically sized properties
     */
    template <bool Static = IsStaticallySized>
    reference operator[](const std::enable_if_t<Static, size_t> & index) {
      static_assert(Static == IsStaticallySized,
                    "SFINAE, don't set 'Static'");
      return Value::get_ref(this->values[index*NbComp]);
    }

    /**
     * Accessor for property by index for dynamically sized properties
     */
    template <bool Dynamic = !IsStaticallySized>
    std::enable_if_t<Dynamic, reference> operator[](const size_t & index) {
      static_assert(Dynamic == !IsStaticallySized,
                    "SFINAE, don't set 'Dynamic'");
      return Value::get_ref(this->values[index*this->get_nb_comp()],
                            this->get_nb_row(),
                            this->get_nb_col());
    }

  protected:
    NeighbourhoodManager & manager;
    std::vector<T> values;
  private:
  };

}  // rascal

#endif /* NEIGHBOURHOOD_PROPERTY_H */
