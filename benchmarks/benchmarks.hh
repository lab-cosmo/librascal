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


    // TODO(alex) rename BaseBFixture to avoid conflict with tests
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
    
    inline bool file_exists(const char* name) {
      struct stat buffer;   
      return (stat (name, &buffer) == 0); 
    }

    template<class Matrix>
    void write_binary(const char* filename, const Matrix& matrix){
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
        out.write((char*) (&rows), sizeof(typename Matrix::Index));
        out.write((char*) (&cols), sizeof(typename Matrix::Index));
        out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
        out.close();
    }
    template<class Matrix>
    void read_binary(const char* filename, Matrix& matrix){
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        typename Matrix::Index rows=0, cols=0;
        in.read((char*) (&rows),sizeof(typename Matrix::Index));
        in.read((char*) (&cols),sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
        in.close();
    }

  } // internal
} // rascal

