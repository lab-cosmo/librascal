#include "benchmark_interpolator.hh"  

//class D {static const int x{2};}
//
//
//template <class Data>
//class MyFixture : public benchmark::Fixture {
//public:
//  void SetUp(const ::benchmark::State& ) {
//    std::cout << "SetUp invoked." << std::endl;
//    if (not(this->init)) {
//      this->x=Data::x;
//      this->init=true;
//    } else {
//      this->x++;
//    }
//  }
//
//  void TearDown(const ::benchmark::State&) {
//  }
//  int x;
//  bool init;
//};

//BENCHMARK_F(MyFixture, FooTest)(benchmark::State& state) {
//  std::cout << "FooTest invoked." << std::endl;
//  for (auto _ : state) {
//  }
//}

// this sadly does not work because GBench uses macros to create the name for the benchmark, and the template brackets disturb
// I have not found a method to be able to have Fixtures with input parameters and automatic SetUp und TearDown methods
//template <class Data>
//BENCHMARK_DEFINE_F(MyFixture<Data>, BarTest)(benchmark::State& state) {
//  std::cout << "BarTest invoked " << x << std::endl;
//  for (auto _ : state) {
//  }
//}
//BENCHMARK_DEFINE_F(MyFixture, BarTest2)(benchmark::State& state) {
//  std::cout << "BarTest2 invoked " << x << std::endl;
//  for (auto _ : state) {
//  }
//}

class MyClass {
public:
  void SetUp() {
  std::cout << "SetUp invoked." << std::endl;
  if (not(this->init)) {
    this->x=0;
    this->init=true;
  } else {
    this->x++;
  }
  }
  int x;
  bool init;
};

MyClass myclass;
void BM_Some(benchmark::State& state) {
  myclass.SetUp();
  std::cout << "BM_Some invoked " << myclass.x << std::endl;
  for (auto _ : state) {
  }
}


void BM_StringCompare(benchmark::State& state) {
  std::string s1(state.range(0), '-');
  std::string s2(state.range(0), '-');
  for (auto _ : state) {
    benchmark::DoNotOptimize(s1.compare(s2));
  }
  state.SetComplexityN(state.range(0));
}


//BENCHMARK_REGISTER_F(MyFixture, BarTest)->Repetitions(5);

//BENCHMARK(BM_Some);
//BENCHMARK(BM_Some);

// Complexity
BENCHMARK(BM_StringCompare)->RangeMultiplier(2)
    ->Range(1<<10, 1<<18)->Complexity([](auto n)->double{return n;});
//BENCHMARK(BM_StringCompare)
//    ->RangeMultiplier(2)->Range(1<<10, 1<<18)->Complexity(benchmark::oN);


BENCHMARK_MAIN();
//TODO(alex)
//- The benchmark does have not the option to initialize ressources only once for a benchmark. There is an issue for this which has not moved since end of 2018 https://github.com/google/benchmark/issues/743

