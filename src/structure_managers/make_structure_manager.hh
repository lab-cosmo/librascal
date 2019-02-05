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

  template< typename ManagerImplementation,  template<class> class AdaptorImplementation,
                              template<class> class ... AdaptorImplementationPack>
  struct AdaptorFactory {
      using Manager_t = AdaptorImplementation<ManagerImplementation>;
      using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
      using ManagerPtr = std::shared_ptr<Manager_t>;
      using type = AdaptorFactory<Manager_t,AdaptorImplementationPack...>;

      template<typename Arg, typename ...Args>
      AdaptorFactory(ImplementationPtr_t& m, Arg& arg, Args& ...args)
      :manager{std::make_shared<Manager_t>(m,arg)},
      next_stack{manager, args...}
      {
        std::cout<<this->manager->get_name()<<std::endl;
      }

      template<typename Arg1, typename Arg2>
      AdaptorFactory(ImplementationPtr_t& m, Arg1& arg1, Arg2& arg2)
      :manager{std::make_shared<Manager_t>(m,arg1)},
      next_stack{manager, arg2}
      {
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

      template<typename Arg>
      AdaptorFactory(ImplementationPtr_t& m, Arg& arg)
      :manager{std::make_shared<Manager_t>(m,arg)}
      {
        std::cout<<this->manager->get_name()<<std::endl;
      }

      ManagerPtr manager;

      decltype(auto) get_manager() {
          return this->manager;
      }
  };

  template <typename ManagerImplementation, template<class> class AdaptorImplementation, template<class> class ... AdaptorImplementationPack, typename Arg, typename ...Args >
  decltype(auto) make_structure_manager(Arg arg, Args ...args) {
    auto manager_base = std::make_shared<ManagerImplementation>(ManagerImplementation{});

    auto test = AdaptorFactory<AdaptorImplementation, AdaptorImplementationPack... >(manager_base, args...);

    auto manager = test.get_manager();

    std::cout  <<std::endl;
    std::cout << manager->get_name() <<std::endl;

    manager->update(arg);

    return manager;
  }


}  // namespace rascal

#endif /* MAKE_STRUCTURE_MANAGER_H */