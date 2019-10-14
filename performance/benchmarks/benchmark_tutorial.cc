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
   * We have the name convention to call fixtures in the benchmarks BFixture to
   * potential prevent collision with the tests. All new BFixture should inherit
   * from BaseBFixture to eget a lookup table for the data json string in the
   * Dataset template parameter. Since a json does have no deterministic order,
   * a map has to build from json index
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
    int x;
    std::string name;
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
   * The `all_combinations_of_arguments` function is used to produce all combinations
   * from the parameters included in a `Dataset` structure. The parameters are
   * accessible in the benchmark function with the benchmark state. However, the
   * google benchmark library only allows to give indices as parameters.
   * Therefore only combinations of indices are accessible from the benchmark
   * state. Furthermore, since a json does have no deterministic order, we
   * cannot directly use the indices to access the data in the json string. A
   * map has to build from json string, to the index inside the json object.
   * This is done in the BaseBFixture. Inside a Fixture class the lookup
   * function can be used with the json data, the state and the targeted
   * property to access it at the right index.
   */
  auto myfix{MyBFixture<SampleData>()};
  BENCHMARK_CAPTURE(bm_example, , myfix)
      ->Apply(all_combinations_of_arguments<SampleData>);
  // we skip naming the benchmark here, to give the benchmark some_name use
  // BENCHMARK_CAPTURE(bm_example, some_name, myfix)

  /**
   * This benchmark should make clear how the different depth of iterations work
   * in google benchmark there are benchmark_repetitions which can be given when
   * executing the executable. This parameter controls how often the benchmark
   * is executed. In addition there are two additional parameters controlling
   * the number of execution. There is a hidden parameter which controls how
   * often a benchmark is executed for per repitition and there is the number of
   * iterations controlling how often the for loop on the state is executed.
   * These two hidden parameters can be manipulated implicitly with the
   * benchmark_min_time flag, but google benchmark itself handles how the time
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
   * RangeMultiplier is the exponential base for all Ranges given later. The
   * first parameter 1 << 10 results in parameters (1,2,4,8), the second
   * (1,2,4,8,16)
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
