/**
 * @file   behler_fixtures.hh
 *
 * @author Till Junge <till.junge@epfl.ch>
 *
 * @date   18 Dec 2019
 *
 * @brief  common fixtures for testing Behler-Parinello descriptors
 *
 * Copyright © 2019 Till Junge
 *
 * rascal is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3, or (at
 * your option) any later version.
 *
 * rascal is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rascal; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

#ifndef TESTS_BEHLER_FIXTURES_HH_
#define TESTS_BEHLER_FIXTURES_HH_

#include "rascal/representations/cutoff_functions_inlineable.hh"
#include "rascal/representations/symmetry_functions.hh"

namespace rascal {

  template <InlCutoffFunctionType FunType>
  struct InlCutoffFunFixture {};

  template <>
  struct InlCutoffFunFixture<InlCutoffFunctionType::Cosine> {
    using CutoffFun = CutoffFunction<InlCutoffFunctionType::Cosine>;
    InlCutoffFunFixture() : cut_fun{unit_style, correct_input} {}
    const units::UnitStyle unit_style{units::metal};

    const double r_cut{1.1};
    CutoffFunctionBase::Hypers_t correct_input{
        {"params", {}}, {"r_cut", {{"value", r_cut}, {"unit", "Å"}}}};

    CutoffFunctionBase::Hypers_t incorrect_put{
        {"params", {}}, {"r_cut", {{"value", r_cut}, {"unit", "J"}}}};
    CutoffFun cut_fun;
  };

  template <SymmetryFunctionType FunType>
  struct SymmetryFunFixture {};

  template <>
  struct SymmetryFunFixture<SymmetryFunctionType::Gaussian> {
    using SymFun = SymmetryFunction<SymmetryFunctionType::Gaussian>;
    SymmetryFunFixture() : sym_fun{unit_style, correct_input} {}
    const units::UnitStyle unit_style{units::metal};

    const double r_ij{1.1};
    json correct_input{{"eta", {{"value", .1}, {"unit", "(Å)^(-1)"}}},
                       {"r_s", {{"value", 5.6}, {"unit", "Å"}}}};

    json incorrect_put{{"eta", {{"value", .1}, {"unit", "(Å)^(-2)"}}},
                       {"r_s", {{"value", 5.6}, {"unit", "Å"}}}};
    SymFun sym_fun;
  };

}  // namespace rascal

#endif  // TESTS_BEHLER_FIXTURES_HH_
