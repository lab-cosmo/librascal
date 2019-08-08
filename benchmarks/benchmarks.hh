#include <benchmark/benchmark.h>
#include <functional>
#include <map>
#include <iostream>

#include "json_io.hh"

namespace rascal {
  namespace internal {

    // Helper function to create all combinations from two arguments, because google only offers combinations of multiple Ranges with exponential incrementation. It is a modified copy of the Ranges function in google benchmark https://github.com/google/benchmark/blob/master/src/benchmark_register.cc#L300-L332
    std::vector<std::vector<int64_t>> combinations(std::vector<std::vector<int64_t>> arglists) {
      //arglists = {{0,1},{0,1,2}};//(ranges.size());
     
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


    // Because Ranges jumps in exponential
    template <class Dataset>
    void AllCombinationsArguments(benchmark::internal::Benchmark* b) {
      json data = Dataset::data();
      std::vector<std::vector<int64_t>> nbs_args;
      nbs_args.resize(data.size());
      // We have to cound backwards to keep the order in the benchmark as in the json string, because the function creating the combinations, creates them backwards
      size_t i{0};
      for (json::iterator it = data.begin(); it != data.end(); ++it) {
        //std::cout << it.key() << "=" << it.value().size() << std::endl;
        std::vector<int64_t> args(it.value().size()) ; // vector with 100 ints.
        std::iota (std::begin(args), std::end(args), 0); // Fill with 0, 1, ..., 99.
        nbs_args.at(i) = args;
        i++;
      }
      auto args_comb = combinations(nbs_args);
      for (auto it = args_comb.begin(); it != args_comb.end(); ++it) {
        std::vector<int64_t> ve = *it;
        for (auto it2 = ve.begin(); it2 != ve.end(); ++it2) {
          //std::cout << *it2 << " ";
        }
        //std::cout << std::endl;
        b->Args(*it);
      }
      //std::cout << std::endl;
    }


    template<class Dataset>
    class BaseFixture {
     protected:
      BaseFixture() {
        json data = Dataset::data();
        this->build_state_range_index_to_key(data);
      }
      // We have do this with an iterator of the json file because it does not
      // perserve the insert order.
      void build_state_range_index_to_key(json & data) {
        //state_range_index_to_key.resize(data.size()); // TODO(alex) change to capacity
        int64_t i{0};
        for (json::iterator it = data.begin(); it != data.end(); ++it) {
          this->state_range_index_to_key.insert(std::pair<std::string, int64_t>(it.key(),i));
          i++;
        }
      }
      
      template <typename T>
      T lookup(const json & data, std::string name, const ::benchmark::State& state) const {
        return this->template lookup<T>(data, name, state.range(this->state_range_index_to_key.at(name))); 
      }

      template <typename T>
      T lookup(const json & data, std::string name, const int64_t index) const {
        return data[name].at(index).get<T>();
      }

      std::map<std::string,int64_t> state_range_index_to_key;
    };
  } // internal
} // rascal
