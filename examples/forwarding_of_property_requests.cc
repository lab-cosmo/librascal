/**
 * @file  forwarding_of_property_requests.cc
 *
 * @author Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   Oct 17 2019
 *
 * @brief example file to introduce default implementations dependent on traits.
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
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

#include <iostream>
#include <memory>

/**
 * This class simulates the framework of the structure manager to explain the
 * functionality of default implementations transparantly without the noise of
 * all the other functionalities of the structure manager classes.
 *
 * Framework translation to structure manager classes
 * Interface = StructureManager
 * RootImplementation = A root structure manager as StructureManagerCenters
 * AdaptorImplementation = An adaptor as AdaptorNeighbourList
 */

template <class Impl>
struct traits;

template <class Impl>
class Interface {
 public:
  typedef typename traits<Impl>::underlying_manager_t underlying_manager_t;
  Interface() = default;

  int get_id() { return this->implementation().get_id(); }

  void get_has_property() {
    if (!this->has_property) {
      std::cout << "ID " << this->get_id() << " has no property" << std::endl;
      if (typeid(underlying_manager_t) != typeid(Impl)) {
        std::cout << "ID " << this->get_id()
                  << " looks in underlying level for property" << std::endl;
        this->get_underlying_manager()->get_has_property();
      }
      return;
    }
    std::cout << "ID " << this->get_id() << " has property" << std::endl;
  }
  std::shared_ptr<underlying_manager_t> get_underlying_manager() {
    return this->implementation().get_underlying_manager_impl();
  }
  void set_has_property(bool has_property) {
    this->has_property = has_property;
  }

 protected:
  Impl & implementation() { return static_cast<Impl &>(*this); }

  int id{0};
  bool has_property{false};
};

class RootImplementation;

template <>
struct traits<RootImplementation> {
  typedef RootImplementation underlying_manager_t;
};

// using Interface's id with constructor
class RootImplementation
    : public Interface<RootImplementation>,
      public std::enable_shared_from_this<RootImplementation> {
 public:
  using this_traits = traits<RootImplementation>;
  typedef typename this_traits::underlying_manager_t underlying_manager_t;
  RootImplementation(int id, bool has_property) {
    this->id = id;
    this->has_property = has_property;
  }

  int get_id() {
    return this->id;  // this gets Interface's id
  }
  std::shared_ptr<underlying_manager_t> get_underlying_manager_impl() {
    // this line wrongly first created a new instance and then wrapped it in a
    // shared ptr instead of only wrapping this in a shared ptr return
    // std::make_shared<underlying_manager_t>(*this);
    return shared_from_this();
  }
};

template <class Impl>
class AdaptorImplementation;

template <class Impl>
struct traits<AdaptorImplementation<Impl>> {
  typedef Impl underlying_manager_t;
};

template <class Impl>
class AdaptorImplementation : public Interface<AdaptorImplementation<Impl>> {
 public:
  using this_traits = traits<AdaptorImplementation>;
  typedef typename this_traits::underlying_manager_t underlying_manager_t;
  AdaptorImplementation(int id, std::shared_ptr<Impl> underlying_manager,
                        bool has_property) {
    this->id = id;
    this->has_property = has_property;
    this->underlying_manager = underlying_manager;
  }

  int get_id() {
    return this->id;  // this gets Interface's id
  }
  std::shared_ptr<underlying_manager_t> get_underlying_manager_impl() {
    return this->underlying_manager;
  }
  std::shared_ptr<underlying_manager_t> underlying_manager{};
};

int main() {
  std::cout << "=====================" << std::endl;
  std::cout << "CASE: No one has property" << std::endl;
  std::cout << "=====================" << std::endl;
  auto impl{std::make_shared<RootImplementation>(0, false)};
  auto adaptor{std::make_shared<AdaptorImplementation<RootImplementation>>(
      1, impl, false)};
  auto adaptor2{std::make_shared<
      AdaptorImplementation<AdaptorImplementation<RootImplementation>>>(
      2, adaptor, false)};
  adaptor2->get_has_property();

  std::cout << "=====================" << std::endl;
  std::cout << "CASE: Lowest has property" << std::endl;
  std::cout << "=====================" << std::endl;
  impl->set_has_property(true);
  adaptor2->get_has_property();

  std::cout << "=====================" << std::endl;
  std::cout << "CASE: First adaptor has property" << std::endl;
  std::cout << "=====================" << std::endl;
  adaptor->set_has_property(true);
  adaptor2->get_has_property();

  std::cout << "=====================" << std::endl;
  std::cout << "CASE: Second adaptor has property" << std::endl;
  std::cout << "=====================" << std::endl;
  adaptor2->set_has_property(true);
  adaptor2->get_has_property();
  return 0;
}
