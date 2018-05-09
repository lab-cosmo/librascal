/**
 * file   field.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   05 Apr 2018
 *
 * @brief  Fields of atom-, pair-, triplet-, etc-related values
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

#ifndef NEIGHBOURHOOD_FIELD_H
#define NEIGHBOURHOOD_FIELD_H

#include <Eigen/Dense>
#include <type_traits>
#include <vector>
#include <basic_types.h>

namespace rascal {

  namespace internal {

    template <typename T, size_t NbRow, size_t NbCol>
    struct Value {
      using type = Eigen::Map<Eigen::Matrix<T, NbRow, NbCol>>;
      using reference = type; //Map<MatrixXd>( resultC, resultEigen.rows(), resultEigen.cols() ) = resultEigen;

      static reference get_ref(T & value) {
        return type(&value);
      }

      static void push_in_vector(std::vector<T> & vec, reference ref) {
        for (int j = 0; j < NbCol; ++j) {
          for (int i = 0; i < NbRow; ++i) {
            vec.push_back(ref(i,j));
          }
        }
      }
    };

    //! specialisation for scalar fields
    template <typename T>
    struct Value<T, 1, 1> {
      using type = T;
      using reference = T&;

      static reference get_ref(T & value) {
        return value;
      }

      static void push_in_vector(std::vector<T> & vec, reference ref) {
        vec.push_back(ref);
      }
    };

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_t = typename Value<T, NbRow, NbCol>::type;

    template<typename T, size_t NbRow, size_t NbCol>
    using Value_ref = typename Value<T, NbRow, NbCol>::reference;

    inline void lin2mult(const Dim_t& index, const Eigen::Ref<const Vec3i_t> shape, Eigen::Ref< Vec3i_t> retval)  {
      int dim{3};
      //Vec3i_t retval;
      Dim_t factor{1};
      for (Dim_t i{0}; i < dim; ++i) {
        retval[i] = index/factor%shape[i];
        if (i != dim-1 ) {
          factor *= shape[i];
        }
      }
      //return retval;
    }

    inline Dim_t mult2lin( const Eigen::Ref<const Vec3i_t> coord, const Eigen::Ref<const Vec3i_t> shape)  {
      int dim{3};
      Dim_t index{0};
      Dim_t factor{1};
      for (Dim_t i = 0; i < dim; ++i) {
        index += coord[i]*factor;
        if (i != dim-1 ) {
          factor *= shape[i];
        }
      }
      return index;
    }

    // https://stackoverflow.com/questions/828092/python-style-integer-division-modulus-in-c
    // TODO more efficient implementation without if (would be less general)
    //! div_mod function returning python like div_mod, i.e. signed integer division truncates towards negative infinity, and signed integer modulus has the same sign the second operand.
    template<class Integral>
    void div_mod(const Integral& x, const Integral& y, std::array<int,2>& out)
    {
      const Integral quot = x / y;
      const Integral rem  = x % y;
      if (rem != 0 && (x < 0) != (y < 0)){
        out[0] = quot - 1;
        out[1] = rem + y;
      }
      else {
        out[0] = quot;
        out[1] = rem;
      }
    }
    
  }  // internal

  template <class NeighbourhoodManager, typename T,
            size_t ClusterSize,
            size_t NbRow = 1, size_t NbCol = 1>
  class Field
  {
    static_assert((std::is_arithmetic<T>::value or
                   std::is_same<T, std::complex<double>>::value),
                  "can currently only handle arithmetic types");
  public:
    using traits = NeighbourhoodManager_traits<NeighbourhoodManager>;
    constexpr static size_t NbDof{NbRow*NbCol};
    using Cluster_t =
      typename NeighbourhoodManager::template ClusterRef
      <ClusterSize, traits::MaxLevel>;

    using Value = internal::Value<T, NbRow, NbCol>;

    using value_type = typename Value::type;
    using reference = typename Value::reference;

    //! Default constructor
    Field() = delete;

    //! constructor with Manager
    Field(NeighbourhoodManager & manager)
      :manager{manager}
    {}

    //! Copy constructor
    Field(const Field &other) = delete;

    //! Move constructor
    Field(Field &&other) = delete;

    //! Destructor
    virtual ~Field() = default;

    //! Copy assignment operator
    Field& operator=(const Field &other) = delete; //! Disallow copying

    //! Move assignment operator
    Field& operator=(Field &&other) = delete; //! Disallow copying

    //! adjust size (only increases, never frees)
    void resize() {
      this->values.resize(this->manager.get_nb_clusters(ClusterSize) * NbDof);
    }

    /**
     * shortens the vector so that the manager can push_back into it
     * (capacity not reduced)
     */
    void resize_to_zero() {
      this->values.resize(0);
    }

    /**
     * allows to add a value to field during construction of the neighbourhood.
     */
    inline void push_back(reference ref) {
      Value::push_in_vector(this->values, ref);
    }

    reference operator[](const Cluster_t& id) {
      return Value::get_ref(this->values[id.get_global_index()*NbDof]);
    }

  protected:
    NeighbourhoodManager & manager;
    std::vector<T> values;
  private:
  };

}  // rascal

#endif /* NEIGHBOURHOOD_FIELD_H */