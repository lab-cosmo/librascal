/**
 * file   property_typed.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   06 Aug 2018
 *
 * @brief Implements intermediate property class for which the type of stored
 *          objects is known, but not the size
 *
 * Copyright  2018 Federico Giberti, Till Junge, Felix Musil, COSMO (EPFL),
 * LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_PROPERTY_TYPED_HH_
#define SRC_STRUCTURE_MANAGERS_PROPERTY_TYPED_HH_

#include "rascal_utility.hh"
#include "math/math_utils.hh"
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
      using value_type = Eigen::Matrix<T, NbRow, NbCol>;
      using reference = Eigen::Map<Eigen::Matrix<T, NbRow, NbCol>>;
      using const_reference =
          const Eigen::Map<const Eigen::Matrix<T, NbRow, NbCol>>;

      //! get a reference to specific value at row and colum
      static reference get_ref(T & value, int nb_row, int nb_col) {
        return reference(&value, nb_row, nb_col);
      }

      //! get a reference
      static reference get_ref(T & value) { return reference(&value); }

      //! push back data into ``property``
      static void push_in_vector(std::vector<T> & vec, reference ref) {
        for (size_t j{0}; j < NbCol; ++j) {
          for (size_t i{0}; i < NbRow; ++i) {
            vec.push_back(ref(i, j));
          }
        }
      }

      //! Dynamic size overloading of push back data into ``property``
      static void push_in_vector(std::vector<T> & vec, reference ref,
                                 Dim_t & nb_row, Dim_t & nb_col) {
        for (Dim_t j{0}; j < nb_col; ++j) {
          for (Dim_t i{0}; i < nb_row; ++i) {
            vec.push_back(ref(i, j));
          }
        }
      }

      //! Used for extending cluster_indices
      template <typename Derived>
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

      //! Dynamic size overloading of push back data into ``property``
      template <typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref,
                                 const Dim_t & nb_row, const Dim_t & nb_col) {
        for (Dim_t j{0}; j < nb_col; ++j) {
          for (Dim_t i{0}; i < nb_row; ++i) {
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
      using value_type = T;
      using reference = T &;
      using const_reference = const T &;

      //! get a reference to a scalar value
      static reference get_ref(T & value) { return value; }

      //! get a reference to a scalar value
      static const_reference get_ref(const T & value) { return value; }

      //! push a scalar in a vector
      static void push_in_vector(std::vector<T> & vec, reference ref) {
        vec.push_back(ref);
      }

      //! Used for extending cluster_indices
      template <typename Derived>
      static void push_in_vector(std::vector<T> & vec,
                                 const Eigen::DenseBase<Derived> & ref) {
        static_assert(Derived::RowsAtCompileTime == NbRow,
                      "NbRow has incorrect size.");
        static_assert(Derived::ColsAtCompileTime == NbCol,
                      "NbCol has incorrect size.");
        vec.push_back(ref(0, 0));
      }
    };

  }  // namespace internal

  /* ---------------------------------------------------------------------- */
  /**
   * Typed ``property`` class definition, inherits from the base property class
   */
  template <typename T, size_t Order, size_t PropertyLayer, class Manager>
  class TypedProperty : public PropertyBase {
   public:
    using Parent = PropertyBase;
    using Value_t = internal::Value<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Manager_t = Manager;
    using Self_t = TypedProperty<T, Order, PropertyLayer, Manager>;
    using traits = typename Manager::traits;
    using Matrix_t = math::Matrix_t;

    using value_type = typename Value_t::value_type;
    using reference = typename Value_t::reference;

    //! constructor
    TypedProperty(Manager_t & manager, Dim_t nb_row, Dim_t nb_col = 1,
                  std::string metadata = "no metadata")
        : Parent{static_cast<StructureManagerBase &>(manager),
                 nb_row,
                 nb_col,
                 Order,
                 PropertyLayer,
                 metadata},
          type_id{internal::GetTypeNameHelper<Self_t>::GetTypeName()} {}

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
    //! return info about the type
    const std::string & get_type_info() const { return this->type_id; }

    Manager_t & get_manager() {
      return static_cast<Manager_t &>(this->base_manager);
    }

    template <size_t Order_ = Order, std::enable_if_t<(Order_ == 1), int> = 0>
    size_t get_validated_property_length(bool consider_ghost_atoms) {
      if (consider_ghost_atoms) {
        if (traits::MaxOrder < 2) {
          throw std::runtime_error(
              "consider_ghost_atoms is true,"
              " but can only be use for underlying manager with"
              " MaxOrder at least 2.");
        }
        if (not(this->get_manager().get_consider_ghost_neighbours())) {
          throw std::runtime_error(
              "consider_ghost_atoms is true,"
              " but underlying manager does not have ghost atoms in"
              " cluster_indices_container. Turn consider_ghost_neighbours"
              " on, to consider ghost atoms with independent property values"
              " from their corresponding central atoms.");
        }
        return this->get_manager().size_with_ghosts();
      }
      return this->get_manager().size();
    }
    template <size_t Order_ = Order,
              std::enable_if_t<not(Order_ == 1), int> = 0>
    size_t get_validated_property_length(bool = false) {
      return this->base_manager.nb_clusters(Order);
    }

    /* Fill sequence, used for *_cluster_indices initialization
     * if consdier_ghost_atoms is true, ghost atoms also can have
     * their own propery value independent from its correpsonding central atom.
     * This function is used for all Order 1 ManagerImplementations
     */
    inline void fill_sequence(bool consider_ghost_atoms = false) {
      // adjust size of values (only increases, never frees)
      this->resize(consider_ghost_atoms);
      for (size_t i{0}; i < this->values.size(); ++i) {
        values[i] = i;
      }
      // fill_sequence happens in update_self so if it is called it means that
      // the property was not up to date (no need to check here if it should
      // be updated)
      this->set_updated_status(true);
    }

    //! Adjust size of values (only increases, never frees)
    inline void resize(bool consider_ghost_atoms = false) {
      auto n_components = this->get_nb_comp();
      size_t new_size =
          this->get_validated_property_length(consider_ghost_atoms) *
          n_components;
      this->values.resize(new_size);
    }

    //! Returns the size of one component
    size_t size() const { return this->values.size() / this->get_nb_comp(); }

    /**
     * shortens the vector so that the manager can push_back into it (capacity
     * not reduced)
     */
    void clear() { this->values.clear(); }

    /* ---------------------------------------------------------------------- */
    //! Property accessor by cluster ref
    template <size_t CallerLayer>
    inline reference operator[](const ClusterRefKey<Order, CallerLayer> & id) {
      static_assert(CallerLayer >= PropertyLayer,
                    "You are trying to access a property that does not exist at"
                    "this depth in the adaptor stack.");
      return this->operator[](id.get_cluster_index(CallerLayer));
    }

    template <size_t CallerOrder, size_t CallerLayer, size_t Order_ = Order>
    inline std::enable_if_t<(Order_ == 1) and (CallerOrder > 1),  // NOLINT
                            reference>                            // NOLINT
    operator[](const ClusterRefKey<CallerOrder, CallerLayer> & id) {
      return this->operator[](
          static_cast<Manager_t &>(this->base_manager)
              .get_cluster_index(id.get_internal_neighbour_atom_tag()));
    }

    //! Accessor for property by index for dynamically sized properties
    reference operator[](const size_t & index) {
      return Value_t::get_ref(this->values[index * this->get_nb_comp()],
                              this->get_nb_row(), this->get_nb_col());
    }

    // //! getter to the underlying data storage
    // inline std::vector<T> & get_raw_data() { return this->values; }

    inline void fill_dense_feature_matrix(Eigen::Ref<Matrix_t> features) {
      size_t n_center{this->get_nb_item()};
      auto n_cols{this->get_nb_comp()};
      auto mat = reference(this->values.data(), n_cols, n_center);
      for (size_t i_center{0}; i_center < n_center; i_center++) {
        for (int i_pos{0}; i_pos < n_cols; i_pos++) {
          // the storage order is swapped here because mat is ColMajor
          features(i_center, i_pos) = mat(i_pos, i_center);
        }
      }
    }

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
      return Value_t::get_ref(this->values[index * this->get_nb_comp()],
                              this->get_nb_row(), this->get_nb_col());
    }

    inline Matrix_t get_dense_feature_matrix() {
      auto nb_centers{this->get_nb_item()};
      auto nb_features{this->get_nb_comp()};
      Matrix_t features(nb_centers, nb_features);
      this->fill_dense_feature_matrix(features);
      return features;
    }

   protected:
    std::string type_id{};
    std::vector<T> values{};  //!< storage for properties
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_PROPERTY_TYPED_HH_
