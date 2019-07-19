/**
 * file   adaptor_center_pairs.hh
 *
 * @author Felix Musil <felix.musil@epfl.ch>
 *         Markus Stricker <markus.stricker@epfl.ch>
 * @date   12 July 2019
 *
 * @brief implements an adaptor for structure_managers, adding the central atom
 * to the iteration over the neighbours
 *
 * Copyright  2019 Felix Musil COSMO (EPFL), LAMMM (EPFL)
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

#ifndef SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_
#define SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_

#include "structure_managers/adaptor_filter.hh"
#include "structure_managers/cluster_ref_key.hh"

namespace rascal {
  /**
   * Forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorCenterPairs;

  /**
   * Specialisation of traits for reset pair <code>MaxOrder</code> for adding
   * center-center pairs
   */
  template <class ManagerImplementation>
  struct StructureManager_traits<AdaptorCenterPairs<ManagerImplementation>> {
    constexpr static AdaptorTraits::Strict Strict{AdaptorTraits::Strict::no};
    constexpr static bool HasDistances{false};
    constexpr static bool HasDirectionVectors{
        ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    // New MaxOrder upon construction stays the same
    // todo(markus) reset PairOrder to zero since new pairs are added
    constexpr static size_t MaxOrder{ManagerImplementation::traits::MaxOrder};

    using LayerByOrder = typename LayerExtender<
        MaxOrder, typename ManagerImplementation::traits::LayerByOrder>::type;
  };

  /**
   * Adaptor that adds pairs of centers, i.e. `ii` pairs to the list of
   * neighbours.
   */
  template <class ManagerImplementation>
  class AdaptorCenterPairs
      : public AdaptorFilter<ManagerImplementation,
                             ManagerImplementation::traits::MaxOrder> {
   public:
    using Manager_t = AdaptorCenterPairs<ManagerImplementation>;
    // using Parent = StructureManager<Manager_t>;
    using Parent = AdaptorFilter<ManagerImplementation,
                                 ManagerImplementation::traits::MaxOrder>;
    using traits = StructureManager_traits<AdaptorCenterPairs>;
    using ImplementationPtr_t = std::shared_ptr<ManagerImplementation>;
    using Hypers_t = typename Parent::Hypers_t;

    //! Default constructor
    AdaptorCenterPairs() = delete;

    AdaptorCenterPairs(ImplementationPtr_t manager) : Parent{manager} {}

    AdaptorCenterPairs(ImplementationPtr_t manager, std::tuple<> /*val*/)
        : AdaptorCenterPairs{manager} {}

    AdaptorCenterPairs(ImplementationPtr_t manager,
                       const Hypers_t & /*adaptor_hypers*/)
        : AdaptorCenterPairs{manager} {}

    //! Copy constructor
    AdaptorCenterPairs(const AdaptorCenterPairs & other) = delete;

    //! Move constructor
    AdaptorCenterPairs(AdaptorCenterPairs && other) = delete;

    //! Destructor
    virtual ~AdaptorCenterPairs() = default;

    //! Copy assignment operator
    AdaptorCenterPairs & operator=(const AdaptorCenterPairs & other) = delete;

    //! Move assignment operator
    AdaptorCenterPairs & operator=(AdaptorCenterPairs && other) = delete;

    //! This is where the magic happens and the center-pairs are added along
    //! with all previous pairs
    void perform_filtering() final {
      for (auto && atom : this->manager) {
        // construct and add ii-pair

        std::array<int, 2> atom_tag_list{atom.get_atom_tag(),
                                         atom.get_atom_tag()};
        Eigen::Map<const Eigen::Matrix<size_t, 1, 1>> index_array{0};

        ClusterRefKey<2, 0> ii_cluster(atom_tag_list, index_array);

        // this->add_cluster(ii_cluster);

        auto && ii_pair{*(atom.begin())};
        ii_pair.set_atom_tag(1, atom.get_atom_tag());
        this->add_cluster(ii_pair);

        // auto new_size{atom.size() + 1};
        // std::cout << "atom.size(), new_size " << atom.size() << ", " <<
        // new_size
        //           << std::endl;
        // std::vector<int> atom_tag_list_new;
        // atom_tag_list_new.resize(new_size);
        // atom_tag_list_new[0] = atom.get_atom_tag();
        // auto i{0};
        // for (auto && pair : atom) {
        //   atom_tag_list_new[++i] = pair.get_atom_tag();
        // }

        // auto atom_tag_list{atom.get_atom_tag_list()};
        // std::cout << "atom_tag_list.size() " << atom_tag_list.size() <<
        // std::endl;
        // //std::array<int,

        // ii_pair.back() = ii_pair.front();
        // ii_pair.get_atom_type() = atom.get_atom_type();
        // ii_pair.get_position() = atom.get_position();

        // this->add_cluster(ii_pair);

        // add all other pairs in neighbour list
        for (auto && pair : atom) {
          this->add_cluster(pair);
        }
      }
    }
  };
}  // namespace rascal

#endif  // SRC_STRUCTURE_MANAGERS_ADAPTOR_CENTER_PAIRS_HH_
