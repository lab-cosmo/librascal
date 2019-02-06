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

#ifndef MAKE_STRUCTURE_MANAGER_H
#define MAKE_STRUCTURE_MANAGER_H

#include "structure_managers/structure_manager.hh"
#include "structure_managers/property.hh"
#include "rascal_utility.hh"

namespace rascal {

  /**
   * Given an instanciated structure manager, recurcivelly stack adaptors on it
   * while instanciating them and create backward links (add_child).
   * It assumes the adaptors arguments needed for constructions
   * are gathered in tuples.
   */
  template< typename ManagerImplementation,  template<class> class AdaptorImplementation,
                              template<class> class ... AdaptorImplementationPack>
  struct AdaptorFactory {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr = std::shared_ptr<Manager_t>;
      using type = AdaptorFactory<Manager_t,AdaptorImplementationPack...>;

      //! General case
      template<typename Arg, typename ...Args>
      AdaptorFactory(ImplementationPtr_t& m, Arg& arg, Args& ...args)
      :manager{std::make_shared<Manager_t>(m,arg)},
      next_stack{manager, args...} {
        m->add_child(this->manager->get_weak_ptr());
        std::cout<<this->manager->get_name()<<std::endl;
      }

      ManagerPtr manager;
      type next_stack;

      decltype(auto) get_manager() {
          return this->next_stack.get_manager();
      }
  };

  template< typename ManagerImplementation,  template<class> class AdaptorImplementation>
  struct AdaptorFactory<ManagerImplementation,AdaptorImplementation> {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr = std::shared_ptr<Manager_t>;
      using type = Manager_t;

      //! End of recursion
      template<typename Arg>
      AdaptorFactory(ImplementationPtr_t& m, Arg& arg)
      :manager{std::make_shared<Manager_t>(m,arg)} {
        m->add_child(this->manager->get_weak_ptr());
        std::cout<<this->manager->get_name()<<std::endl;
      }

      ManagerPtr manager;

      decltype(auto) get_manager() {
          return this->manager;
      }
  };


  /**
   * Factory function to build a structure managers.
   *
   * @tparams Manager type of the base manager, e.g. StructureManagerCenters
   * @tparams AdaptorImplementationPack list of adaptors to stack on the base
   * manager type
   * @params arg argument to update the structure of the base structure manager
   * @params args list of arguments to build the adaptor (packed in tuples and
   *  in the same order as in AdaptorImplementationPack)
   * @return shared pointer to the fully built structure manager
   */
  template <typename Manager, template<class> class ... AdaptorImplementationPack, typename Arg, typename ...Args >
  decltype(auto) make_structure_manager(Arg arg, Args ...args) {
    // instanciate the base manager
    auto manager_base = std::make_shared<Manager>();

    // build the stack of adaptors
    auto factory =
      AdaptorFactory<Manager, AdaptorImplementationPack... >(manager_base, args...);
    // get the manager with the full stack
    auto manager = factory.get_manager();

    // std::cout  <<std::endl;
    // std::cout << manager->get_name() <<std::endl;

    // give aa structure to the underlying base manager
    // and update all the stack of adaptors
    manager->update(arg);

    return manager;
  }

  template<typename StructureManagerPtr, int TargetLevel>
  struct UnderlyingManagerExtractor {
    using ManagerPtr_t = typename StructureManagerPtr::element_type::ImplementationPtr_t;
    using type = UnderlyingManagerExtractor<ManagerPtr_t, TargetLevel + 1>;

    UnderlyingManagerExtractor(StructureManagerPtr& sm)
    :manager{sm->get_previous_manager()}, next_stack{manager} {}

    ManagerPtr_t manager;
    type next_stack;

    decltype(auto) get_manager() {
        return this->next_stack.get_manager();
    }
  };

  template<typename StructureManagerPtr>
  struct UnderlyingManagerExtractor<StructureManagerPtr, 0> {
    using ManagerPtr_t = StructureManagerPtr;

    UnderlyingManagerExtractor(StructureManagerPtr& sm)
    :manager{sm} {}

    ManagerPtr_t manager;

    decltype(auto) get_manager() {
        return this->manager->get_shared_ptr();
    }
  };

  template<int TargetLevel, typename StructureManagerPtr>
  decltype(auto) extract_underlying_manager(StructureManagerPtr manager) {

    static_assert(TargetLevel < 0, "target_level should be negative");
    auto test{UnderlyingManagerExtractor<StructureManagerPtr, TargetLevel>(manager)};

    return test.get_manager();
  }


}  // namespace rascal

#endif /* MAKE_STRUCTURE_MANAGER_H */