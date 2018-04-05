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
 * proteus is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * proteus is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Emacs; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <Eigen/Dense>

#include <type_traits>
#include <vector>

namespace proteus {

  namespace internal {

    template <typename T, size_t NbRow, size_t NbCol>
    struct Value {
      using type = Eigen::Map<Eigen::Matrix<T, NbRow, NbCol>>;
      using reference = type;

      static reference get_ref(T & value) {
        return type(&value);
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
    };

    template<typename t, size_t NbRow, size_t NbCol>
    using Value_t = typename Value<T, NbRow, NbCol>::type;

    template<typename t, size_t NbRow, size_t NbCol>
    using Value_ref = typename Value<T, NbRow, NbCol>::reference;

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
    using reference =typename Value::reference;

    //! Default constructor
    Field() = delete;

    //! constructor with Manager
    Field(NeighbourhoodManager & manager): manager{manager},
    {}

    //! Copy constructor
    Field(const Field &other) = delete;

    //! Move constructor
    Field(Field &&other) = delete;

    //! Destructor
    virtual ~Field() = default;

    //! Copy assignment operator
    Field& operator=(const Field &other) = delete;

    //! Move assignment operator
    Field& operator=(Field &&other) = delete;

    //! adjust size (only increases, never frees)
    void resize() {
      this->values.resize(manager.get_nb_clusters(ClusterSize) * NbDof)};
    }

    reference operator[](const Cluster_t& id) {
      return Value::ret_ref(this->values[id.get_global_index()*NbDof];)

  protected:
    NeighbourhoodManager & manager;
    std::vector<T> values;
  private:
  };


}  // proteus
