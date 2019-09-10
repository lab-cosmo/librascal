/**
 * file   behler_feature.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   15 Jan 2019
 *
 * @brief utilities for the definition of input nodes (combinations of symmetry
 * functions with cutoff functions
 *
 * Copyright © 2019 Till Junge, COSMO (EPFL), LAMMM (EPFL)
 *
 * librascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * librascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with librascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef SRC_REPRESENTATIONS_BEHLER_FEATURE_HH_
#define SRC_REPRESENTATIONS_BEHLER_FEATURE_HH_

#include "representations/cutoff_functions.hh"
#include "representations/symmetry_functions.hh"

namespace rascal {

  template <class StructureManager>
  class BehlerFeatureBase {
   public:
    template <size_t Order>
    using ClusterRef_t = typename StructureManager::template ClusterRef<Order>;

    using StdSpecies =
        TupleStandardisation<int, StructureManager::traits::MaxOrder>;

    //! Default constructor
    BehlerFeatureBase() = delete;

    //! Constructor with symmetry function type
    BehlerFeatureBase(const SymmetryFunType sym_fun_type, const size_t order)
        : sym_fun_type{sym_fun_type}, order{order} {}

    //! Copy constructor
    BehlerFeatureBase(const BehlerFeatureBase & other);

    //! Move constructor
    BehlerFeatureBase(BehlerFeatureBase && other) noexcept;

    //! Destructor
    virtual ~BehlerFeatureBase() noexcept;

    //! Copy assignment operator
    BehlerFeatureBase & operator=(const BehlerFeatureBase & other);

    //! Move assignment operator
    BehlerFeatureBase & operator=(BehlerFeatureBase && other) noexcept;

    //! needs to be called after reading the input file and prior to the first
    //! evaluation. Attaches all necessary properties for precalculated values
    //! to the manager
    virtual void init(const UnitStyle & units) = 0;

    //! needs to be called in the beginning of every evaluation step, refreshes
    //! the precalculated properties if necessary
    virtual void prepare(StructureManager & manager) = 0;

    //! Main worker (raison d'être)
    virtual void apply(StructureManager & manager) const = 0;

    //! insert a parameter (sub-)json
    void add_params(const json & params) {
      const auto & type{params.at("type").get<std::string>()};
      if (type != get_name(this->sym_fun_type)) {
        std::stringstream error{};
        error << "Parameter set for function type '" << type
              << "' assigned to function of type '"
              << get_name(this->sym_fun_type) << "'.";
        throw std::runtime_error(error.str());
      }
      this->raw_params.push_back(params);
    }

   protected:
    std::vector<json> raw_params{};
    const SymmetryFunType sym_fun_type;
    const size_t order;
  };

  /* ---------------------------------------------------------------------- */
  template <SymmetryFunType SymFunType, internal::CutoffFunctionType CutFunType,
            class StructureManager>
  class BehlerFeature final : public BehlerFeatureBase<StructureManager> {
   public:
    using Parent = BehlerFeatureBase<StructureManager>;
    using SymmetryFunction = SymmetryFun<SymFunType>;
    using CutoffFunction = internal::CutoffFunction<CutFunType>;
    using StdSpecies = typename Parent::StdSpecies;
    template <size_t Order>
    using ClusterRef_t = typename Parent::ClusterRef_t;

    // stores parameter packs ordered by cutoff radius
    using ParamStorage =
        std::map<double,
                 Eigen::Matrix<double, SymmetryFun<SymFunType>::NbParams,
                               Eigen::Dynamic>>;

    constexpr static size_t MaxOrder{Parent::traits::MaxOrder};

    //! Default constructor
    BehlerFeature() : Parent(SymFunType) {}

    //! Copy constructor
    BehlerFeature(const BehlerFeature & other) = delete;

    //! Move constructor
    BehlerFeature(BehlerFeature && other) = default;

    //! Destructor
    ~BehlerFeature() = default;

    //! Copy assignment operator
    BehlerFeature & operator=(const BehlerFeature & other) = delete;

    //! Move assignment operator
    BehlerFeature & operator=(BehlerFeature && other) = default;

    void init(const UnitStyle & units) final;
    void prepare(StructureManager & manager) final;

    void apply(StructureManager & manager) const;

    StdSpecies get_species_combo() const;  // to implement
    size_t get_index() const;              // to implement

   protected:
    template <size_t Order>
    inline void eval_cluster(StructureManager & manager,
                             const ClusterRef_t<Order> & cluster);
    static constexpr size_t AtomOrder{1};
    static constexpr size_t PairOrder{2};
    static constexpr size_t AtomLayer{
        std::get<0>(StructureManager::traits::LayerByOrder::type)};
    static constexpr size_t PairLayer{
        std::get<1>(StructureManager::traits::LayerByOrder::type)};

    ParamStorage params{};
  };

}  // namespace rascal
#include "behler_feature_impl.hh"

#endif  // SRC_REPRESENTATIONS_BEHLER_FEATURE_HH_
