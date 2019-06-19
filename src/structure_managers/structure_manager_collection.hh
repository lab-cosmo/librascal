/**
 * file   structure_manager_collection.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *
 * @date   13 Jun 2019
 *
 * @brief Implementation of a container for structure managers
 *
 * Copyright 2019 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_COLLECTION_HH_
#define SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_COLLECTION_HH_

#include "structure_managers/structure_manager.hh"
#include "structure_managers/make_structure_manager.hh"
#include "structure_managers/property.hh"
#include "structure_managers/updateable_base.hh"
#include "math/math_utils.hh"
#include "rascal_utility.hh"
#include "json_io.hh"
#include "atomic_structure.hh"


namespace rascal {


  template<typename Manager,
            template <class> class... AdaptorImplementationPack>
  class ManagerCollection {
   public:
    using Self_t = ManagerCollection<Manager, AdaptorImplementationPack...>;
    using TypeHolder_t = StructureManagerTypeHolder<Manager, AdaptorImplementationPack...>;
    using Manager_t = typename TypeHolder_t::type;
    using ManagerPtr_t = std::shared_ptr<Manager_t>;
    using ManagerList_t = typename TypeHolder_t::type_list;
    using Hypers_t = typename Manager_t::Hypers_t;
    using traits = typename Manager_t::traits;
    using Data_t = std::vector<ManagerPtr_t>;
    using value_type = typename Data_t::value_type;

   protected:
    Data_t managers{};
    Hypers_t adaptor_inputs{};

   public:
    ManagerCollection() = default;

    explicit ManagerCollection(const Hypers_t& adaptor_inputs) {
      this->adaptor_inputs = adaptor_inputs;
    };


    //! Copy constructor
    ManagerCollection(const ManagerCollection & other) = delete;

    //! Move constructor
    ManagerCollection(ManagerCollection && other) = default;

    //! Destructor
    ~ManagerCollection() = default;

    //! Copy assignment operator
    ManagerCollection & operator=(const ManagerCollection & other) = delete;

    //! Move assignment operator
    ManagerCollection & operator=(ManagerCollection && other) = default;


    /**
     * Give the ManagerCollection the iterator functionality using Data_t
     * functionality
     */
    using iterator = typename Data_t::iterator;
    using const_iterator = typename Data_t::const_iterator;

    inline iterator begin() noexcept {
      return this->managers.begin();
    }
    inline const_iterator begin() const noexcept {
      return this->managers.begin();
    }

    inline iterator end() noexcept {
      return this->managers.end();
    }
    inline const_iterator end() const noexcept {
      return this->managers.end();
    }

    //! set the global inputs for the adaptors
    inline void set_adaptor_inputs(const Hypers_t& adaptor_inputs) {
      this->adaptor_inputs = adaptor_inputs;
    }

    inline const Hypers_t& get_adaptors_parameters() const {
      return this->adaptor_inputs;
    }

    /**
     * functions to add a(several) structures to the collection
     */
    inline void add_structure(const Hypers_t& structure, const Hypers_t& adaptor_inputs) {
      auto manager = make_structure_manager_stack<Manager, AdaptorImplementationPack...>(structure, adaptor_inputs);
      this->add_structure(manager);
    }

    inline void add_structure(const Hypers_t& structure) {
      auto manager = make_structure_manager_stack<Manager, AdaptorImplementationPack...>(structure, this->adaptor_inputs);
      this->add_structure(manager);
    }

    inline void add_structure(std::shared_ptr<Manager_t>& manager) {
      this->managers.emplace_back(manager);
    }

    /**
     * Function used from python. add empty structures to build the objects and
     * update them afterwards (small workaround because AtomicStructure<3>
     * can't be put in a json object).
     */
    void add_structures(const std::vector<AtomicStructure<3>>& atomic_structures) {
      Hypers_t structure = Hypers_t::object();
      for (const auto& atomic_structure : atomic_structures) {
        this->add_structure(structure);
        this->managers.back()->update(atomic_structure);
      }
    }

    void add_structures(const Hypers_t& structures, const Hypers_t& adaptors_inputs) {
      if (not structures.is_array()) {
        throw std::runtime_error(R"(Provide the structures as an array
        (or list) of json dictionary defining the structure)");
      }
      if (structures.size() != adaptors_inputs.size()) {
        throw std::runtime_error(R"(There should be as many structures as
        adaptors_inputs)");
      }

      for (int i_structure{0}; i_structure < structures.size(); ++i_structure) {
        this->add_structure(structures[i_structure], adaptors_inputs[i_structure]);
      }
    }

    void add_structures(const Hypers_t& structures) {
      if (not structures.is_array()) {
        throw std::runtime_error(R"(Provide the structures as an array
        (or list) of json dictionary defining the structure)");
      }

      for (auto& structure : structures) {
        this->add_structure(structure);
      }
    }

    /**
     * load structures from a json file
     *
     * @param filename path to the file containing the structures in
     * ase json format
     * @param start index of the first structure to include
     * @param length number of structure to include, -1 correspondons to all
     * after start
     *
     * Compatible file format are text and ubjson (binary).
     *
     * Note that start refers to 0 based indexing so 0 corresponds to the
     * first and 3 would corresponds to the 4th structure irrespective of the
     * actual indices in the file.
     */
    void add_structures(const std::string& filename, const int& start = 0,
                        int length = -1) {
      // important not to do brace initialization because it adds an extra
      // nesting layer
      json structures = json_io::load(filename);

      if (not structures.is_object()) {
        throw std::runtime_error(R"(The ase format's first level is a dictionary with indicies as keys to the structures)");
      }

      if (structures.count("ids") == 1) {
        // structures is in the ase format
        auto ids{structures["ids"].get<std::vector<int>>()};
        std::sort(ids.begin(), ids.end());
        ids.erase(ids.begin(), ids.begin() + start);
        if (length == -1) {
          length = ids.size();
        }
        ids.erase(ids.begin() + length, ids.end());

        for (auto& idx : ids) {
          this->add_structure(structures[std::to_string(idx)].get<Hypers_t>());
        }
      } else {
        throw std::runtime_error("The json structure format is not recognized");
      }
    }

    //! number of structure manager in the collection
    inline size_t size() const {
      return this->managers.size();
    }

  };

}

#endif  // SRC_STRUCTURE_MANAGERS_STRUCTURE_MANAGER_COLLECTION_HH_
