#include "neighbourhood_managers/neighbourhood_manager_base.hh"
#include "neighbourhood_managers/property.hh"


#ifndef ADAPTOR_STRICT_H
#define ADAPTOR_STRICT_H

namespace rascal {
  /**
   * forward declaration for traits
   */
  template <class ManagerImplementation>
  class AdaptorStrict;

  /**
   * specialisation of traits for strict adaptor
   */
  template <class ManagerImplementation>
  struct NeighbourhoodManager_traits<AdaptorStrict<ManagerImplementation>> {

    constexpr static Strict Strict{AdaptorTraits::Strict::yes};
    constexpr static bool HasDistances{true};
    constexpr static bool HasDirectionVectors{ManagerImplementation::traits::HasDirectionVectors};
    constexpr static int Dim{ManagerImplementation::traits::Dim};
    constexpr static int MaxLevel{ManagerImplementation::traits::MaxLevel};
    using UpdateStruct = typename ManagerImplementation::traits::UpdateStruct;
  };

  /**
   * Adaptor that guarantees that only neighbours within the cutoff
   * are present. This interface should be implemented by all managers
   * with the trait AdaptorTraits::Strict::yes
   */
  template <class ManagerImplementation>
  class AdaptorStrict: public
  NeighbourhoodManagerBase<AdaptorStrict<ManagerImplementation>>
  {
  public:
    using Parent = ManagerImplementation;
    using traits = NeighbourhoodManager_traits<AdaptorStrict>;
    using AtomRef_t = typename Parent::AtomRef_t;
    template <int Level, int MaxLevel>
    using ClusterRef_t = typename Parent::template <Level, MaxLevel> ClusterRef_t;
    using PairRef_t = ClusterRef_t<2, traits::MaxLevel>;

    static_assert(traits::MaxLevel > 1>,
                  "ManagerImlementation needs to handle pairs");

    //! Default constructor
    AdaptorStrict() = delete;

    /**
     * construct a strict neighbourhood list from a given manager and cut-off radius
     */
    AdaptorStrict(ManagerImplementation& manager, double cut_off);

    //! Copy constructor
    AdaptorStrict(const AdaptorStrict &other) = delete;

    //! Move constructor
    AdaptorStrict(AdaptorStrict &&other) = default;

    //! Destructor
    virtual ~AdaptorStrict() = default;

    //! Copy assignment operator
    AdaptorStrict& operator=(const AdaptorStrict &other) = delete;

    //! Move assignment operator
    AdaptorStrict& operator=(AdaptorStrict &&other) = default;

    template<class ... Args>
    void update(Args&&... arguments) {
      this->manager.update(std::forward<Args>(arguments)...);
      do my stuff;
    }

    inline double get_cutoff() const;

    inline double get_distance(const PairRef_t & pair) const;

    inline size_t get_nb_clusters(int cluster_size) const;

    inline size_t get_size() const;

    inline Vector_ref get_position(const AtomRef_t & atom);

    // return the global id of an atom
    inline size_t get_atom_id(const Parent& , int i_atom_id) const;

    // return the global id of an atom
    template<int Level, int MaxLevel>
    inline size_t get_atom_id(const ClusterRef_t<Level, MaxLevel>& cluster,
                              int j_atom_id) const;

      inline Vector_ref get_neighbour_position(const AtomRef_t& atom, const AtomRef_t& center_atom,const int& j_linear_id);

    //! return atom type
    inline int get_atom_type(const AtomRef_t& atom);

    /**
     * return the linear index of cluster (i.e., the count at which
     * this cluster appears in an iteration
     */
    template<int Level, int MaxLevel>
    inline int get_offset_impl(const ClusterRef_t<Level, MaxLevel>& cluster) const;

  protected:
    ManagerImplementation & manager;
    Property<AdaptorStrict, double, 2> distance;
    double cut_off;
  private:
  };




  namespace internal {

    template<class 
    struct 

  }  // internal
  //----------------------------------------------------------------------------//
  template < class ManagerImplementation>
  AdaptorStrict<ManagerImplementation>::
  AdaptorStrict(ManagerImplementation & manager, double cut_off):
    manager{manager},
    distance{*this},
    cut_off{cut_off}
  {
    if (not internal::check_cut_off(manager, cut_off)) {
      throw std::runtime_error("underlying manager already has a smaller cut_off");
    }
  }
}  // rascal

#endif /* ADAPTOR_STRICT_H */
