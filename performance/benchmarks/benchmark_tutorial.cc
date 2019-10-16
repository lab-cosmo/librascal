/**
 * file  benchmark_tutorial.cc
 *
 * @author  Alexander Goscinski <alexander.goscinski@epfl.ch>
 *
 * @date   22 August 2019
 *
 * @brief tutorial how to write a benchmark file
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

#include "benchmarks.hh"

namespace rascal {
  /**
   * Structure holding the data, the data is created inside a static function,
   * because non primitives cannot be static as member variable.
   */
  struct SampleData {
    static const json data() {
      static json data = {{"x", {1, 2, 3}},
                          {"name", {"crystal.json", "cube.json"}}};
      return data;
    }
  };

  /**
   * We have the naming convention to call fixtures in the benchmarks BFixture
   * to potential prevent collision with the tests. All new BFixtures should
   * inherit from `BaseBFixture` to inherit a lookup table for the json
   * keys in the data json in the Dataset template parameter e.g. "x" and "name"
   * in the `SampleData`. According to the json standard a json does have no
   * deterministic order, therefore a map has to be built, mapping the json
   * keys to the correct position when accesing the json with an index.
   */
  template <typename Dataset>
  class MyBFixture : public BaseBFixture<Dataset> {
   public:
    using Parent = BaseBFixture<Dataset>;
    MyBFixture<Dataset>() : Parent() {}

    /**
     * By convention we use a setup function which handles the initialization of
     * new parameters, it should be invoked at the start of the benchmark.
     */
    void setup(benchmark::State & state) {
      std::string old_name = this->name;

      // initialize
      json data = Dataset::data();
      // lookup takes the indices inside the state and gets parameter x from
      // dataset at these indices
      this->x = this->template lookup<int>(data, "x", state);
      this->name = this->template lookup<std::string>(data, "name", state);

      // this check tests if the data has actually changed, if it has changed
      // then we do the costly initialization.
      if (this->name != old_name) {
        // Some very costly initialization
        std::cout << "Costly initialization." << std::endl;
      }
    }
    int x{};
    std::string name{};
  };

  template <class BFixture>
  void bm_example(benchmark::State & state, BFixture & fix) {
    fix.setup(state);
    std::cout << "bm_example invoked " << fix.x << ", " << fix.name
              << std::endl;
    //
    for (auto _ : state) {
      // do some tests for speed
    }
  }

  /**
   * The `all_combinations_of_arguments` function is used to produce all
   * combinations from the parameters included in a `Dataset` structure. The
   * parameters are accessible in the benchmark function with the benchmark
   * state. However, the google benchmark library only allows to give indices as
   * parameters. Therefore only combinations of indices are accessible from the
   * benchmark state. Furthermore, since a json does not have a deterministic
   * order by standard, we cannot directly use the indices to access the data in
   * the json. A map has to be built from the json key to the correct position
   * when accessing the json with an index. This is done in the `BaseBFixture`.
   * Inside the `BaseBFixture class the `lookup` function can be used with the
   * data json, the state and the targeted property to access it at the correct
   * position.
   */
  auto myfix{MyBFixture<SampleData>()};
  BENCHMARK_CAPTURE(bm_example, /* some_name_if_wanted */, myfix)
      ->Apply(all_combinations_of_arguments<SampleData>);

  /**
   * This benchmark should make clear how the different depths of iterations
   * work in google benchmark there is the `benchmark_repetitions` flag which
   * can be given when executing the benchmark. This parameter controls how
   * often each benchmark function is invoked. There are two more non-accessible
   * parameters controlling the number of invokes of a benchmark function. One
   * controls how often a benchmark is executed per repetition. Simply put, this
   * parameter is an additional repetiton of the benchmark function. The other
   * parameter controls the number of iterations inside the for-loop on the
   * state. These two hidden parameters can be manipulated implicitly with the
   * `benchmark_min_time` flag, but google benchmark itself handles how the time
   * changes the two parameters. To understand this better, please execute this
   * benchmark multiple times with different values.
   */
  void BM_StringCompare(benchmark::State & state) {
    std::cout << "Preiteration " << state.range(0) << std::endl;
    std::string s1(state.range(0), '-');
    std::string s2(state.range(0), '-');
    /* We do not print inside the iterations because this value can be seen when
     * the benchmarks are run.
     */
    // BEGIN ITERATION
    for (auto _ : state) {
      benchmark::DoNotOptimize(s1.compare(s2));
    }  // END ITERATION
    state.SetComplexityN(state.range(0));
  }

  /**
   * `RangeMultiplier` is the exponential base for all ranges given later by
   * `Range` or `Ranges`. The first parameter 1 << 10 results in parameters
   * (1,2,4,8), the second (1,2,4,8,16)
   */
  BENCHMARK(BM_StringCompare)->RangeMultiplier(2)->Range(1 << 10, 1 << 18);

  /**
   * Repetition repeats everything in this case 3 times and calculates
   * averages
   */
  BENCHMARK(BM_StringCompare)
      ->RangeMultiplier(2)
      ->Range(1 << 10, 1 << 18)
      ->Repetitions(3);

  /**
   * Google benchmark can estimate the complexity. To estimate the complexity
   * inside the benchmark the parameter determing the size is required.
   */
  BENCHMARK(BM_StringCompare)
      ->RangeMultiplier(2)
      ->Range(1 << 10, 1 << 18)
      ->Complexity();
  // alternatively an own complexity function can be given
  //->Range(1<<10, 1<<18)->Complexity([](auto n)->double{return n;});

}  // namespace rascal
BENCHMARK_MAIN();
