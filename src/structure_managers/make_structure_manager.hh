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

namespace rascal {
  /**
   * Factory function to make a structure manager
   * @tparams Adaptor partial type of the adaptor
   * @params Manager input structure manager
   * @params args additional argument for the constructructor
   */
  template <typename Manager>
  std::shared_ptr<Manager> make_structure_manager() {
    return std::make_shared<Manager>();
  }

  /**
   * Factory function to make an adapted structure manager
   * @tparams Adaptor partial type of the adaptor
   * @params Manager input structure manager
   * @params args additional argument for the constructructor
   */
  template <template <class> class Adaptor, typename Manager, typename... Args>
  std::shared_ptr<Adaptor<Manager>>
  make_adapted_manager(std::shared_ptr<Manager> & arg, const Args &... args) {
    auto manager{std::make_shared<Adaptor<Manager>>(arg, args...)};
    arg->add_child(manager->get_weak_ptr());
    return manager;
  }

  template <template <class> class Adaptor, typename Manager>
  std::shared_ptr<Adaptor<Manager>>
  make_adapted_manager_json(std::shared_ptr<Manager> & arg, const json &adaptor_hypers) {
    auto manager{std::make_shared<Adaptor<Manager>>(arg, adaptor_hypers)};
    arg->add_child(manager->get_weak_ptr());
    return manager;
  }

  /* ---------------------------------------------------------------------- */

  namespace internal {
    /**
     * Make an adapted structure manager with the elements of the Tuple
     * from Min to Max and the provided Args
     */
    template <template <class> class Adaptor, size_t Min, size_t Max>
    struct make_adapted_manager_util {
      template <typename Manager, typename... Argst>
      static constexpr auto apply(Manager & manager,
                                  std::tuple<Argst...> & tuple) {
        static_assert(Min < Max, "Min should be smaller than Max");
        return index_apply<Min, Max>([&](auto... Is) {
          return make_adapted_manager<Adaptor>(manager, std::get<Is>(tuple)...);
        });
      }
    };
    // case Min == Max
    template <template <class> class Adaptor, size_t Val>
    struct make_adapted_manager_util<Adaptor, Val, Val> {
      template <typename Manager, typename... Argst>
      static constexpr auto apply(Manager & manager,
                                  std::tuple<Argst...> & tuple) {
        return make_adapted_manager<Adaptor>(manager);
      }
    };

    //! unpack the tuple elements and give them to the update function of the
    //! manager
    template <class ManagerPtr, class Tuple>
    constexpr void apply_update_util(ManagerPtr & manager, Tuple t) {
      return index_apply<0, std::tuple_size<Tuple>{}>(
          [&](auto... Is) { manager->update(std::get<Is>(t)...); });
    }

    /**
     * Given an instanciated structure manager, recursively stack adaptors
     * on it
     * while instanciating them and create backward links (add_child).
     * It assumes the adaptors arguments needed for constructions
     * are gathered in tuples.
     */
    template <size_t CurrentPosition, typename ManagerImplementation,
              template <class> class AdaptorImplementation,
              template <class> class... AdaptorImplementationPack>
    struct AdaptorFactory {
      // static_assert(CurrentPosition < TotalLenght, "Problem in the
      // recursion")
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      constexpr static size_t NextPosition{
          CurrentPosition + Manager_t::traits::AdaptorInitiParams};

      using type =
          AdaptorFactory<NextPosition, Manager_t, AdaptorImplementationPack...>;

      //! General case
      template <typename... Args>
      AdaptorFactory(ImplementationPtr_t & m, std::tuple<Args...> & tuple)
          : manager{make_adapted_manager_util<AdaptorImplementation,
                                              CurrentPosition,
                                              NextPosition>::apply(m, tuple)},
            next_stack{manager, tuple} {}

      ManagerPtr_t manager;
      type next_stack;

      decltype(auto) get_manager() { return this->next_stack.get_manager(); }
    };

    //! End of recursion
    template <size_t CurrentPosition, typename ManagerImplementation,
              template <class> class AdaptorImplementation>
    struct AdaptorFactory<CurrentPosition, ManagerImplementation,
                          AdaptorImplementation> {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      using type = Manager_t;
      constexpr static size_t NextPosition{
          CurrentPosition + Manager_t::traits::AdaptorInitiParams};

      template <typename... Args>
      AdaptorFactory(ImplementationPtr_t & m, std::tuple<Args...> & tuple)
          : manager{make_adapted_manager_util<AdaptorImplementation,
                                              CurrentPosition,
                                              NextPosition>::apply(m, tuple)} {}

      ManagerPtr_t manager;

      ManagerPtr_t get_manager() { return this->manager; }
    };

  }  // namespace internal

  /**
   * Factory function to build a structure managers.
   *
   * @tparams Manager type of the base manager, e.g. StructureManagerCenters
   * @tparams AdaptorImplementationPack list of adaptors to stack on the base
   * manager type
   * @params structure argument to update the structure of the
   * base structure manager. It should be a tuple if several arguments are
   * needed.
   * @params adaptor_inputs list of arguments to build the adaptor
   * in the same order as AdaptorImplementationPack (the order of the
   * constructors signature need to be followed too).
   * @return shared pointer to the fully built structure manager
   */
  template <typename Manager,
            template <class> class... AdaptorImplementationPack,
            typename Structure, typename... Args>
  decltype(auto) make_structure_manager_stack(Structure structure,
                                              Args... adaptor_inputs) {
    auto tuple{std::make_tuple(adaptor_inputs...)};
    // instanciate the base manager
    auto manager_base = make_structure_manager<Manager>();
    // build the stack of adaptors
    auto factory =
        internal::AdaptorFactory<0, Manager, AdaptorImplementationPack...>(
            manager_base, tuple);
    // get the manager with the full stack
    auto manager = factory.get_manager();
    // give a structure to the underlying base manager
    // and update all the stack of adaptors
    manager->update(structure);

    return manager;
  }


  namespace internal {
    /**
     */
    template <template <class> class Adaptor, size_t HyperPos>
    struct make_adapted_manager_hypers_util {
      template <typename Manager_t, typename Hypers_t>
      static constexpr auto apply(Manager_t & manager,
                                  const Hypers_t & adaptors_hypers) {
        const Hypers_t& adaptor_hypers{adaptors_hypers.at(HyperPos).at("initialization_arguments")};
        return make_adapted_manager_json<Adaptor>(manager, adaptor_hypers);
      }
    };

    /**
     * Given an instanciated structure manager, recursively stack adaptors
     * on it
     * while instanciating them and create backward links (add_child).
     * It assumes the adaptors arguments needed for constructions
     * are gathered in tuples.
     */
    template <size_t CurrentPosition, typename ManagerImplementation,
              template <class> class AdaptorImplementation,
              template <class> class... AdaptorImplementationPack>
    struct AdaptorFactory_hypers {
      // static_assert(CurrentPosition < TotalLenght, "Problem in the
      // recursion")
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;
      constexpr static size_t NextPosition{CurrentPosition + 1};

      using NextFactory_t =
          AdaptorFactory_hypers<NextPosition, Manager_t, AdaptorImplementationPack...>;

      //! General case
      template <typename... Args, typename Hypers_t>
      AdaptorFactory_hypers(ImplementationPtr_t & m, const Hypers_t & adaptors_hypers)
          : manager{make_adapted_manager_hypers_util<AdaptorImplementation,
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
      AdaptorFactory_hypers(ImplementationPtr_t & m, const Hypers_t & adaptors_hypers)
          : manager{make_adapted_manager_hypers_util<AdaptorImplementation,
                                              CurrentPosition>::apply(m, adaptors_hypers)} {}

      ManagerPtr_t manager;

      ManagerPtr_t get_manager() { return this->manager; }
    };
  }

  template <typename Manager,
            template <class> class... AdaptorImplementationPack, typename Hypers_t>
  decltype(auto) make_structure_manager_stack_hypers(const Hypers_t& structure_inputs,
                                              const Hypers_t& adaptor_inputs) {
    // check input consistency
    if (not adaptor_inputs.is_array()) {
      throw std::runtime_error("adaptor_inputs should be a list of dictionary like objects containing the parameters to construct the adaptors.");
    }

    if (adaptor_inputs.size() != sizeof...(AdaptorImplementationPack)) {
      throw std::runtime_error("adaptor_inputs should have as many dictionary of parameters as there are adaptors to build.");
    }

    AtomicStructure<3> structure{};
    structure.set_structure(structure_inputs);
    // instanciate the base manager
    auto manager_base = make_structure_manager<Manager>();
    // build the stack of adaptors
    auto factory =
        internal::AdaptorFactory_hypers<0, Manager, AdaptorImplementationPack...>(manager_base, adaptor_inputs);
    // get the manager with the full stack
    auto manager = factory.get_manager();
    // give a structure to the underlying base manager
    // and update all the stack of adaptors
    manager->update(structure);

    return manager;
  }


  //! same as make_structure_manager_stack but with a tuple input
  template <typename Manager,
            template <class> class... AdaptorImplementationPack,
            typename... Structure, typename... AdaptorInputs>
  decltype(auto) make_structure_manager_stack_with_tuple(
      std::tuple<std::tuple<Structure...>, std::tuple<AdaptorInputs...>> &
          tuple) {
    auto & structure{std::get<0>(tuple)};
    auto & adaptor_inputs{std::get<1>(tuple)};
    // instanciate the base manager
    auto manager_base = make_structure_manager<Manager>();
    // build the stack of adaptors
    auto factory =
        internal::AdaptorFactory<0, Manager, AdaptorImplementationPack...>(
            manager_base, adaptor_inputs);
    // get the manager with the full stack
    auto manager = factory.get_manager();
    // give a structure to the underlying base manager
    // and update all the stack of adaptors
    internal::apply_update_util(manager, structure);

    return manager;
  }

  /* ---------------------------------------------------------------------- */
  //! Utility to hold a list of Adaptors partial types
  template <template <class> class... AdaptorImplementation>
  struct AdaptorTypeHolder;
  //! Utility to hold a fully typed structure manager and a list of Adaptors
  template <typename Manager, template <class> class... AdaptorImplementation>
  struct StructureManagerTypeHolder {
    using type_list =
        std::tuple<Manager, AdaptorTypeHolder<AdaptorImplementation...>>;
    using type =
        typename internal::AdaptorTypeStacker<Manager,
                                              AdaptorImplementation...>::type;
  };

  namespace internal {
    /**
     * Allow to provide a StructureManagerTypeHolder instead of the list types
     * to make_structure_manager_stack.
     */
    template <typename MI, typename AdaptorTypeHolder_>
    struct make_structure_manager_stack_with_tuple_util;

    template <typename MI, template <class> class... Ti>
    struct make_structure_manager_stack_with_tuple_util<
        MI, AdaptorTypeHolder<Ti...>> {
      using Manager_t = typename AdaptorTypeStacker<MI, Ti...>::type;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;

      template <typename Tuple>
      static ManagerPtr_t apply(Tuple & tuple) {
        return make_structure_manager_stack_with_tuple<MI, Ti...>(tuple);
      }
    };

    /**
     * Allow to provide a StructureManagerTypeHolder instead of the list types
     * to make_structure_manager_stack.
     */
    template <typename MI, typename AdaptorTypeHolder_>
    struct make_structure_manager_stack_with_hypers_util;

    template <typename MI, template <class> class... Ti>
    struct make_structure_manager_stack_with_hypers_util<
        MI, AdaptorTypeHolder<Ti...>> {
      using Manager_t = typename AdaptorTypeStacker<MI, Ti...>::type;
      using ManagerPtr_t = std::shared_ptr<Manager_t>;

      template <typename Hypers_t>
      static ManagerPtr_t apply(const Hypers_t& structure_inputs,
                                const Hypers_t& adaptor_inputs) {
        return make_structure_manager_stack_hypers<MI, Ti...>(structure_inputs, adaptor_inputs);
      }
    };
  }  // namespace internal
  /**
   * Factory function to make a manager with its types provided with
   * a StructureManagerTypeHolder and arguments packaged in two json object.
   */
  template <typename StructureManagerTypeHolder_>
  struct make_structure_manager_stack_with_hypers_and_typeholder;

  template <typename... T>
  struct make_structure_manager_stack_with_hypers_and_typeholder<
      std::tuple<T...>> {
    template <typename Hypers_t>
    static decltype(auto) apply(const Hypers_t& structure_inputs,
                                const Hypers_t& adaptor_inputs) {
      return internal::make_structure_manager_stack_with_hypers_util<
          T...>::apply(structure_inputs, adaptor_inputs);
    }
  };

  /**
   * Factory function to make a manager with its types provided with
   * a StructureManagerTypeHolder and arguments packaged such as
   * tuple<tuple<StructureArgs>,tuple<AdaptorConstructorArgs>>.
   */
  template <typename StructureManagerTypeHolder_>
  struct make_structure_manager_stack_with_tuple_and_typeholder;

  template <typename... T>
  struct make_structure_manager_stack_with_tuple_and_typeholder<
      std::tuple<T...>> {
    template <typename Tuple>
    static decltype(auto) apply(Tuple & tuple) {
      return internal::make_structure_manager_stack_with_tuple_util<
          T...>::apply(tuple);
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
            typename... Args>
  decltype(auto) stack_adaptors(std::shared_ptr<Manager> & manager_base,
                                Args... args) {
    // build the stack of adaptors
    auto factory =
        internal::AdaptorFactory<0, Manager, AdaptorImplementationPack...>(
            manager_base, std::make_tuple(args...));
    // get the manager with the full stack
    auto manager = factory.get_manager();

    // TODO(felix) make sure that when updating adaptors without a new
    // structure that the whole tree is not regenerated again but only
    // the parts that needs it
    // update all the stack of adaptors
    manager->update();

    return manager;
  }

  template <typename StructureManagerPtr, int TargetLevel>
  struct UnderlyingManagerExtractor {
    using ManagerPtr_t =
        typename StructureManagerPtr::element_type::ImplementationPtr_t;
    using type = UnderlyingManagerExtractor<ManagerPtr_t, TargetLevel + 1>;

    explicit UnderlyingManagerExtractor(StructureManagerPtr & sm)
        : manager{sm->get_previous_manager()}, next_stack{manager} {}

    ManagerPtr_t manager;
    type next_stack;

    decltype(auto) get_manager() { return this->next_stack.get_manager(); }
  };

  template <typename StructureManagerPtr>
  struct UnderlyingManagerExtractor<StructureManagerPtr, 0> {
    using ManagerPtr_t = StructureManagerPtr;

    explicit UnderlyingManagerExtractor(StructureManagerPtr & sm)
        : manager{sm} {}

    ManagerPtr_t manager;

    decltype(auto) get_manager() { return this->manager->get_shared_ptr(); }
  };

  template <int TargetLevel, typename StructureManagerPtr>
  decltype(auto) extract_underlying_manager(StructureManagerPtr manager) {
    static_assert(TargetLevel < 0, "target_level should be negative");
    auto test{
        UnderlyingManagerExtractor<StructureManagerPtr, TargetLevel>(manager)};

    return test.get_manager();
  }

}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_MAKE_STRUCTURE_MANAGER_HH_
