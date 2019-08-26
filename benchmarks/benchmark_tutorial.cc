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
  namespace internal {
    class SampleData {
     public:
      static const json data() {
        return {
          {"x", {1,2,3}},
          {"name", {"crystal.json","cube.json"}}
        };
      }
    };

    template <typename Dataset>
    class MyFix : public BaseFixture<Dataset> {
     public:
      using Parent = BaseFixture<Dataset>;
      MyFix<Dataset>() : Parent() {}

      void SetUp(benchmark::State& state) {
        std::string old_name = this->name;

        // initialize
        json data = Dataset::data();
        this->x = this->template lookup<int>(data, "x", state);
        this->name = this->template lookup<std::string>(data, "name", state);

        if (this->name != old_name) {
          // Some very costly initialization
          std::cout << "Initialize file." << std::endl;
        }
      }
      int x;
      std::string name;
    };

    template <class Fix>
    void BM_Example(benchmark::State& state, Fix & fix) {
      fix.SetUp(state);
      std::cout << "BM_Example invoked " << fix.x << ", " << fix.name << std::endl;
      for (auto _ : state) {
        // do some tests for speed 
      }
    }

    auto myfix{MyFix<SampleData>()};
    BENCHMARK_CAPTURE(BM_Example, , myfix)->Apply(AllCombinationsArguments<SampleData>);
    // we skip name here BENCHMARK_CAPTURE(BM_Example, some_name, myfix)

    void BM_StringCompare(benchmark::State& state) {
      /*  
       / The benchmark function is executed around 10 times depending one the min time given.
       / There is a tradeoff between the number of preiterations and the number
       / of iterations. The total per test can be controlled with the 
       / benchmark_min_time flag.
      */
      std::cout << "Preiteration " << state.range(0) << std::endl;
      std::string s1(state.range(0), '-');
      std::string s2(state.range(0), '-');
      // BEGIN ITERATION
      for (auto _ : state) {
        benchmark::DoNotOptimize(s1.compare(s2));
      } // END ITERATION
      state.SetComplexityN(state.range(0));
    }

    BENCHMARK(BM_StringCompare)->RangeMultiplier(2)
        ->Range(1<<10, 1<<18);


    // Repetition repeats everything in this case 5 times and calculates
    // averages
    //BENCHMARK_REGISTER_F(MyFixture, BarTest)->Repetitions(5);

    // Complexity
    //BENCHMARK(BM_StringCompare)
    //    ->RangeMultiplier(2)->Range(1<<10, 1<<18)->Complexity();
    //     ->Range(1<<10, 1<<18)->Complexity([](auto n)->double{return n;});


    //TODO(alex)
    //- The benchmark does have not the option to initialize ressources only once for a benchmark. There is an issue for this which has not moved since end of 2018 https://github.com/google/benchmark/issues/743
  }
}
BENCHMARK_MAIN();
