/**
 * file  benchmarks.hh
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief contains structures and functions to allow usage of a simplified
 *        usage of flexible benchmarks
 *
 * Copyright  2019 Alexander Goscinski, COSMO (EPFL), LAMMM (EPFL)
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
 * along with this software; see the file LICENSE. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <benchmark/benchmark.h>
#include <functional>
#include <map>
#include <iostream>

#include <sys/stat.h>
#include <unistd.h>
#include <fstream>

#include "json_io.hh"

namespace rascal {
  /* Helper function to create all tuple combinations between the vectors of
   * integers in arglists. 
   * ranges. The beginning of the range is always 0, therefore only the ending
   * of the range is given as argument.
   *
   * @param arglists -> vector of arguments for different parameters 
   * args_of_parameter1 x args_of_parameter2 x .... x args_of_parametern
   * 
   * @return -> given [[0,1],[2],[5,6]] output [[0,2,5],[0,2,6],[1,2,5],[1,2,6]]
   *
   * This is an adapted version of the `Ranges` function of
   * Gbenchmark.
   *
   * Reference:
   * https://github.com/google/benchmark/blob/master/src/benchmark_register.cc#L300-L332
   */
  std::vector<std::vector<int64_t>> combinations(std::vector<std::vector<int64_t>> arglists) {

    std::size_t total = 1;
    for (std::size_t i = 0; i < arglists.size(); i++) {
      total *= arglists[i].size();
    }
    std::vector<std::vector<int64_t>> args_;

    std::vector<std::size_t> ctr(arglists.size(), 0);

    for (std::size_t i = 0; i < total; i++) {
      std::vector<int64_t> tmp;
      tmp.reserve(arglists.size());

      for (std::size_t j = 0; j < arglists.size(); j++) {
        tmp.push_back(arglists[j].at(ctr[j]));
      }

      args_.push_back(std::move(tmp));

      for (std::size_t j = 0; j < arglists.size(); j++) {
        if (ctr[j] + 1 < arglists[j].size()) {
          ++ctr[j];
          break;
        }
        ctr[j] = 0;
      }
    }
    return args_;
  }

  /**
   * Takes the arguments for each parameter given by the template parameter Dataset and sets the benchmark arguments with an iterator of tuples of indices of each combination. These tuples of indices are used later to access the data in the static `Dataset` template parameter class.
   *
   * Note:
   * GBenchmark offers a functionality doing this with benchmarks using `Ranges`, but the arguments can only be given as a range with exponential step function. We would need a combinatiorial `DenseRange`, which has a linear step function, but is only offered for argument `Range`
   */
  template <class Dataset>
  void AllCombinationsArguments(benchmark::internal::Benchmark * b) {
    json data = Dataset::data();
    std::vector<std::vector<int64_t>> nbs_args;
    nbs_args.resize(data.size());
    // We have to cound backwards to keep the order in the benchmark as in the
    // json string, because the function creating the combinations, creates
    // them backwards
    size_t i{0};
    for (json::iterator it = data.begin(); it != data.end(); ++it) {
      std::vector<int64_t> args(it.value().size());
      std::iota(std::begin(args), std::end(args),
                0); 
      nbs_args.at(i) = args;
      i++;
    }
    auto args_comb = combinations(nbs_args);
    for (auto it = args_comb.begin(); it != args_comb.end(); ++it) {
      std::vector<int64_t> ve = *it;
      for (auto it2 = ve.begin(); it2 != ve.end(); ++it2) {
      }
      b->Args(*it);
    }
  }

  /**
   * Gives utilities to access the json data in the `Dataset` structure with indices given by the benchmark state. Since a json is an undordered collection of data, when it is iterated through the json to produce the indices for the state object, the indices do not necessary correspond to the order in the json string e.g. for a data structure of the form
   *
      static json data = {
        {"x", {1, 2, 3}},
        {"name", {"crystal.json", "cube.json"}}};
   * 
   * There is no guarantee what at position "x" is at the first position. Therefore, when generating the indices for the state, the resulting indices could be
   *
   * state.range = {{0,1,2},{0,1}} or {{0,1},{0,1,2}}
   *
   * This base fixture handles this by creating a map from json string to index.
   *
   * For {{0,1,2},{0,1}} state_range_index_to_key = {("x",0),("name",1)}
   * For {{0,1},{0,1,2}} state_range_index_to_key = {("x",1),("name",0)}
   */
  template <class Dataset>
  class BaseBFixture {
   protected:
    BaseBFixture() {
      json data = Dataset::data();
      this->build_state_range_index_to_key(data);
    }
    // We have do this with an iterator of the json file because it does not
    // perserve the insert order.
    void build_state_range_index_to_key(json & data) {
      int64_t i{0};
      for (json::iterator it = data.begin(); it != data.end(); ++it) {
        this->state_range_index_to_key.insert(
            std::pair<std::string, int64_t>(it.key(), i));
        i++;
      }
    }

    /**
     * @param state contaning a combination
     */
    template <typename T>
    T lookup(const json & data, std::string name,
             const ::benchmark::State & state) const {
      return this->template lookup<T>(
          data, name, state.range(this->state_range_index_to_key.at(name)));
    }

    template <typename T>
    T lookup(const json & data, std::string name, const int64_t index) const {
      return data[name].at(index).get<T>();
    }

    // json string -> position in json string
    std::map<std::string, int64_t> state_range_index_to_key;
  };

  inline bool file_exists(const char * name) {
    struct stat buffer;
    return (stat(name, &buffer) == 0);
  }
}  // namespace rascal
