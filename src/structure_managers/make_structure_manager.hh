/**
 * file   make_structure_manager.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   05 Feb 2019
 *
 * @brief implements the structure manager factory
 *
 * Copyright © 2018 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_MAKE_STRUCTURE_MANAGER_HH_
#define SRC_STRUCTURE_MANAGERS_MAKE_STRUCTURE_MANAGER_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"

#include "rascal_utility.hh"
#include "json_io.hh"
#include "atomic_structure.hh"

#include <type_traits>

namespace rascal {
  /**
   * Factory function to make a structure manager
   * @tparams Adaptor partial type of the adaptor
   * @params Manager input structure manager
   * @params args additional argument for the constructor
   */
  template <typename Manager>
  std::shared_ptr<Manager> make_structure_manager() {
    return std::make_shared<Manager>();
  }

  /**
   * Factory function to make an adapted structure manager
   * @tparams Adaptor partial type of the adaptor
   * @params Manager input structure manager
   * @params args additional argument for the constructor
   */
  template <template <class> class Adaptor, typename Manager, typename... Args>
  std::shared_ptr<Adaptor<Manager>>
  make_adapted_manager(std::shared_ptr<Manager> & arg, const Args &... args) {
    auto manager{std::make_shared<Adaptor<Manager>>(arg, args...)};
    arg->add_child(manager->get_weak_ptr());
    return manager;
  }

  /**
   * Factory function to make an adapted structure manager
   * @tparams Adaptor partial type of the adaptor
   * @params Manager input structure manager
   * @params adaptor_hypers additional argument for the constructor given
   *         in a dictionary like containner, e.g. json type.
   */
  template <template <class> class Adaptor, typename Manager, typename Hypers_t>
  std::shared_ptr<Adaptor<Manager>>
  make_adapted_manager_hypers(std::shared_ptr<Manager> & arg,
                              const Hypers_t & adaptor_hypers) {
    auto manager{std::make_shared<Adaptor<Manager>>(arg, adaptor_hypers)};
    arg->add_child(manager->get_weak_ptr());
    return manager;
  }

  /* ---------------------------------------------------------------------- */

  namespace internal {
    /**
     * helper to make an adapted manager taking a list of dictionary like
     * objects containing the initialization of some adaptors and which one
     * to pick.
     * @tparams Adaptor partial type of the adaptor
     * @tparams HyperPos position of the relevant initialization_arguments for
     *          the Adaptor
     * @params manager a structure manager
     * @params list of dictionary like objects
     *
     */
    template <template <class> class Adaptor, size_t HyperPos>
    struct make_adapted_manager_hypers_util {
      template <typename Manager_t, typename Hypers_t>
      static constexpr auto apply(Manager_t & manager,
                                  const Hypers_t & adaptors_hypers) {
        const Hypers_t & adaptor_hypers{
            adaptors_hypers.at(HyperPos).at("initialization_arguments")};
        return make_adapted_manager_hypers<Adaptor>(manager, adaptor_hypers);
      }
    };

    /**
     * Given an instanciated structure manager, recursively stack adaptors
     * on it
     * while instanciating them and create backward links (add_child).
     * It assumes the adaptors arguments needed for constructions
     * are gathered in a list of dictionary like objects, e.g. json type.
     */
    template <size_t CurrentPosition, typename ManagerImplementation,
              template <class> class AdaptorImplementation,
              template <class> class... AdaptorImplementationPack>
    struct AdaptorFactory_hypers {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      constexpr static size_t NextPosition{CurrentPosition + 1};

      using NextFactory_t = AdaptorFactory_hypers<NextPosition, Manager_t,
                                                  AdaptorImplementationPack...>;

      //! General case
      template <typename... Args, typename Hypers_t>
      AdaptorFactory_hypers(ImplementationPtr_t & m,
                            const Hypers_t & adaptors_hypers)
          : manager{make_adapted_manager_hypers_util<
                AdaptorImplementation,
                CurrentPosition>::apply(m, adaptors_hypers)},
            next_stack{manager, adaptors_hypers} {}

      ManagerPtr_t manager;
      NextFactory_t next_stack;

      decltype(auto) get_manager() { return this->next_stack.get_manager(); }
    };

    //! End of recursion
    template <size_t CurrentPosition, typename ManagerImplementation,
              template <class> class AdaptorImplementation>
    struct AdaptorFactory_hypers<CurrentPosition, ManagerImplementation,
                                 AdaptorImplementation> {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      using type = Manager_t;
      constexpr static size_t NextPosition{CurrentPosition + 1};

      template <typename... Args, typename Hypers_t>
      AdaptorFactory_hypers(ImplementationPtr_t & m,
                            const Hypers_t & adaptors_hypers)
          : manager{make_adapted_manager_hypers_util<
                AdaptorImplementation,
                CurrentPosition>::apply(m, adaptors_hypers)} {}

      ManagerPtr_t manager;

      ManagerPtr_t get_manager() { return this->manager; }
    };
  }  // namespace internal

  /**
   * Factory function to make an adapted structure manager from a list of
   * template parameters and constructor arguments.
   * @tparams AdaptorImplementationPack list of partial types of the adaptors
   * @tparams Manager type of the structure manager root
   * @tparams Hypers_t type of the dictionary like container that aggregate
   *          the parameters to construct the successive adapted managers
   * @params structure_inputs info to get an atomic structure using
   *         the AtomicStructure class
   * @params adaptor_hypers arguments for the constructructor of the adapted
   *         managers given in the same order as AdaptorImplementationPack
   *         in a dictionary like containner, e.g. json type.
   */
  template <typename Manager,
            template <class> class... AdaptorImplementationPack,
            typename Hypers_t>
  decltype(auto) make_structure_manager_stack(const Hypers_t & structure_inputs,
                                              const Hypers_t & adaptor_inputs) {
    // check input consistency
    if (not adaptor_inputs.is_array()) {
      throw std::runtime_error(
          "adaptor_inputs should be a list of dictionary like objects "
          "containing the parameters to construct the adaptors.");
    }

    if (adaptor_inputs.size() != sizeof...(AdaptorImplementationPack)) {
      throw std::runtime_error("adaptor_inputs should have as many dictionary "
                               "of parameters as there are adaptors to build.");
    }

    // instanciate the base manager
    auto manager_base = make_structure_manager<Manager>();
    // build the stack of adaptors
    auto factory = internal::AdaptorFactory_hypers<
        0, Manager, AdaptorImplementationPack...>(manager_base, adaptor_inputs);
    // get the manager with the full stack
    auto manager = factory.get_manager();
    // give a structure to the underlying base manager
    // and update all the stack of adaptors
    if (structure_inputs.size() > 0) {
      AtomicStructure<3> structure{};
      structure.set_structure(structure_inputs);
      manager->update(structure);
    }

    return manager;
  }

  /* ---------------------------------------------------------------------- */
  //! Utility to hold a list of Adaptors partial types
  template <template <class> class... AdaptorImplementation>
  struct AdaptorTypeHolder;
  /**
   * Utility to hold the fully typed structure manager with adaptors as stack
   * and as a list.
   *
   * StructureManagerTypeHolder<StructureManagerCenters, AdaptorNeighbourList,
   *                            AdaptorStrict>
   *
   * -> type_list
   *    std::tuple<StructureManagerCenters,
   *               AdaptorTypeHolder<AdaptorNeighbourList, AdaptorStrict>>
   *
   * -> type
   *    AdaptorStrict<AdaptorNeighbourList<StructureManagerCenters>>
   */
  template <typename Manager, template <class> class... AdaptorImplementation>
  struct StructureManagerTypeHolder {
    // handle the case without adaptors with conditional_t
    using type_list = std::conditional_t<
        sizeof...(AdaptorImplementation) == 0, std::tuple<Manager>,
        std::tuple<Manager, AdaptorTypeHolder<AdaptorImplementation...>>>;

    using type =
        std::conditional_t<sizeof...(AdaptorImplementation) == 0, Manager,
                           typename internal::AdaptorTypeStacker<
                               Manager, AdaptorImplementation...>::type>;
  };

  namespace detail {
    template <template <typename Manager,
                        template <class> class... AdaptorImplementation>
              class Collection,
              typename SM, typename AdaptorTypeHolder_>
    struct InjectTypeHolderHelper;

    template <template <typename Manager,
                        template <class> class... AdaptorImplementation>
              class Collection,
              typename SM, template <class> class... Ti>
    struct InjectTypeHolderHelper<Collection, SM, AdaptorTypeHolder<Ti...>> {
      using type = Collection<SM, Ti...>;
    };

    template <template <typename Manager,
                        template <class> class... AdaptorImplementation>
              class Collection,
              typename StructureManagerTypeHolder_>
    struct InjectTypeHolderUtil;

    template <template <typename Manager,
                        template <class> class... AdaptorImplementation>
              class Collection,
              typename... T>
    struct InjectTypeHolderUtil<Collection, std::tuple<T...>> {
      using type = typename InjectTypeHolderHelper<Collection, T...>::type;
    };
  }  // namespace detail
     /**
      * Utility class holding the fully typed Collection class in type member
      *
      * @tparam Collection a class templated by a structure manager and a list
      * of adaptors
      *
      * @tparam StructureManagerTypeHolder_ a
      *                  StructureManagerTypeHolder::type_list
      *
      * This utility does not help directly for templated function,
      * so to handle this case the function should inserted in a functor.
      * C++17 would allow to avoid the functor
      * see https://stackoverflow.com/a/49291186/11609484.
      */
  template <template <typename Manager,
                      template <class> class... AdaptorImplementation>
            class Collection,
            typename StructureManagerTypeHolder_>
  struct TypeHolderInjector {
    using type = typename detail::InjectTypeHolderUtil<
        Collection, StructureManagerTypeHolder_>::type;
  };

  namespace internal {
    /**
     * Allow to provide a StructureManagerTypeHolder instead of the list types
     * to make_structure_manager_stack.
     *
     * @tparams SM structure manager type
     * @tparams Ti list of adaptor partial types
     */
    template <typename SM, typename AdaptorTypeHolder_>
    struct make_structure_manager_stack_with_hypers_util;

    template <typename SM, template <class> class... Ti>
    struct make_structure_manager_stack_with_hypers_util<
        SM, AdaptorTypeHolder<Ti...>> {
      using Manager_t = typename AdaptorTypeStacker<SM, Ti...>::type;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;

      template <typename Hypers_t>
      static ManagerPtr_t apply(const Hypers_t & structure_inputs,
                                const Hypers_t & adaptor_inputs) {
        return make_structure_manager_stack<SM, Ti...>(structure_inputs,
                                                       adaptor_inputs);
      }
    };
  }  // namespace internal

  /**
   * TODO(felix) use TypeHolderInjector to simplify this thing
   *
   * Factory function to make a manager with its types provided with
   * a StructureManagerTypeHolder and arguments packaged in two json object.
   */
  template <typename StructureManagerTypeHolder_>
  struct make_structure_manager_stack_with_hypers_and_typeholder;

  template <typename... T>
  struct make_structure_manager_stack_with_hypers_and_typeholder<
      std::tuple<T...>> {
    template <typename Hypers_t>
    static decltype(auto) apply(const Hypers_t & structure_inputs,
                                const Hypers_t & adaptor_inputs) {
      return internal::make_structure_manager_stack_with_hypers_util<
          T...>::apply(structure_inputs, adaptor_inputs);
    }
  };

  /**
   * Factory function to stack adaptors on a structure managers with a valid
   *  structure already registered.
   *
   * @tparams Manager type of the base manager, e.g. StructureManagerCenters
   * @tparams AdaptorImplementationPack list of adaptors to stack on the base
   * manager type
   * @params manager_base a structure manager with a structure inside
   * @params args list of arguments to build the adaptor (packed in tuples and
   *  in the same order as in AdaptorImplementationPack)
   * @return shared pointer to the fully built structure manager
   */
  template <typename Manager,
            template <class> class... AdaptorImplementationPack,
            typename Hypers_t>
  decltype(auto) stack_adaptors(std::shared_ptr<Manager> & manager_base,
                                const Hypers_t & hypers) {
    // build the stack of adaptors
    auto factory = internal::AdaptorFactory_hypers<
        0, Manager, AdaptorImplementationPack...>(manager_base, hypers);
    // get the manager with the full stack
    auto manager = factory.get_manager();

    // update all the stack of adaptors
    manager->update();

    return manager;
  }

  namespace internal {
    template <typename StructureManager, int TargetLevel>
    struct UnderlyingManagerExtractor {
      using Manager_t = typename StructureManager::ManagerImplementation_t;
      using type = UnderlyingManagerExtractor<Manager_t, TargetLevel - 1>;

      explicit UnderlyingManagerExtractor(
          std::shared_ptr<StructureManager> & sm)
          : manager{sm->get_previous_manager()}, next_stack{manager} {}

      std::shared_ptr<Manager_t> manager;
      type next_stack;

      decltype(auto) get_manager() { return this->next_stack.get_manager(); }
    };

    template <typename StructureManager>
    struct UnderlyingManagerExtractor<StructureManager, 0> {
      explicit UnderlyingManagerExtractor(
          std::shared_ptr<StructureManager> & sm)
          : manager{sm} {}

      std::shared_ptr<StructureManager> manager;

      decltype(auto) get_manager() { return this->manager->get_shared_ptr(); }
    };

  }  // namespace internal

  template <int TargetLevel, typename StructureManager>
  decltype(auto)
  extract_underlying_manager(std::shared_ptr<StructureManager> manager) {
    static constexpr int n_step_below =
        StructureManager::traits::StackLevel - TargetLevel;
    static_assert(
        n_step_below >= 0,
        "TargetLevel is larger than the number of manager in the stack");
    auto test{
        internal::UnderlyingManagerExtractor<StructureManager, n_step_below>(
            manager)};

    return test.get_manager();
  }

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_MAKE_STRUCTURE_MANAGER_HH_
