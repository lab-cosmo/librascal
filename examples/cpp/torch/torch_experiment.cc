#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/structure_managers/make_structure_manager.hh"
#include "rascal/structure_managers/structure_manager_centers.hh"
#include "rascal/utils/json_io.hh"
#include "rascal/utils/basic_types.hh"
#include "rascal/math/utils.hh"

#include <torch/torch.h>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <chrono>

#include <torch/script.h> // One-stop header.
#include <memory>
#include <filesystem>

using namespace rascal;  // NOLINT

using Matrix_t = math::Matrix_t;

using Representation_t = CalculatorSphericalExpansion;
using Manager_t = AdaptorStrict<
    AdaptorCenterContribution<AdaptorNeighbourList<StructureManagerCenters>>>;
using Prop_t = typename CalculatorSphericalExpansion::Property_t<Manager_t>;
using PropGrad_t =
    typename CalculatorSphericalExpansion::PropertyGradient_t<Manager_t>;

// adapted from https://discuss.pytorch.org/t/data-transfer-between-libtorch-c-and-eigen/54156/6
torch::Tensor eigen_matrix_to_torch_tensor(Matrix_t eigen_mat) {
    // one can play more around with the options, not sure at the moment
    auto options =
      torch::TensorOptions()
        .dtype(torch::kDouble)
        .layout(torch::kStrided)
        .device(torch::kCPU) // kCUDA, 1
        .requires_grad(true);
     
    auto tensor = torch::empty({eigen_mat.rows(), eigen_mat.cols()}, options);
    double* data = tensor.data_ptr<double>();

    // TODO check from_blop 
    Eigen::Map<Matrix_t> eigen_mat_cast(data, tensor.size(0), tensor.size(1));
    eigen_mat_cast = eigen_mat.cast<double>();
    return tensor;
}

// copied code from sergeys script
std::ifstream get_file_input(std::string file_path) {
    std::ifstream ifs (file_path, std::ifstream::in);
    //std::ios_base::sync_with_stdio(false);
    ifs.tie(nullptr);
    return ifs;
}

// copied code from sergeys script
auto read_tensor(std::string file_path) {  
    std::ifstream input = get_file_input(file_path); 
    int dim;
    input >> dim;

    std::vector<long int> shape;
    for (int i = 0; i < dim; ++i) {
        long int now;
        input >> now;
        shape.push_back(now);
    }

    long int total = 1;
    for (int i = 0; i < dim; ++i) {
        total = total * shape[i];
    }

    std::vector<float> data;
    for (int i = 0; i < total; ++i) {
        float now;
        input >> now;
        data.push_back(now);
    }
    auto tensor = torch::from_blob(data.data(), total).clone();

    auto shape_proper_type = at::IntArrayRef(shape.data(), dim);
    return tensor.view(shape_proper_type);

}

int main(int argc, char * argv[]) {
  if (argc < 2) {
    std::cerr << "Must provide {experiment json filename}";
    std::cerr << std::endl;
    return -1;
  }

  std::string experiment_filename{argv[1]};
  std::string output_filename{argv[1]};

  // removes .json from filename
  const std::string ext(".json");
  if ( output_filename != ext &&
       output_filename.size() > ext.size() &&
       output_filename.substr(output_filename.size() - ext.size()) == ".json" )
  {
     // if so then strip them off
     output_filename = output_filename.substr(0, output_filename.size() - ext.size());
  }
  output_filename.append("_results.json");

  json experiment_hypers = json_io::load(experiment_filename);

  json sph_hypers;
  if (experiment_hypers.count("spherical_expansion_hypers")) {
    sph_hypers = experiment_hypers.at("spherical_expansion_hypers").get<json>();
  } else {
    throw std::runtime_error(std::string("Could not find spherical expansion hypers."));
    return 1;
  }

  json structures_absolute_path;
  if (experiment_hypers.count("structures_absolute_path")) {
    structures_absolute_path = experiment_hypers.at("structures_absolute_path").get<json>();
  } else {
    throw std::runtime_error(std::string("Could not find structures absolute path."));
    return 1;
  }

  std::string torch_model_absolute_path;
  if (experiment_hypers.count("torch_model_absolute_path")) {
    torch_model_absolute_path = experiment_hypers.at("torch_model_absolute_path").get<std::string>();
  } else {
    throw std::runtime_error(std::string("Could not find torch model absolute path."));
    return 1;
  }

  // the neighbourlist needs to also know the cutoff so we have to extract it here
  double cutoff;
  if (sph_hypers.count("cutoff_function")) {
    if (sph_hypers.at("cutoff_function").count("cutoff")) {
      if (sph_hypers.at("cutoff_function").at("cutoff").count("value")) {
        cutoff = sph_hypers.at("cutoff_function").at("cutoff").at("value").get<double>();
      } else {
        throw std::runtime_error(std::string("Could not find value in value."));
        return 1;
      }
    } else {
      throw std::runtime_error(std::string("Could not find cutoff in cutoff function hypers."));
      return 1;
    }
  } else {
    throw std::runtime_error(std::string("Could not find cutoff function in spherical expansion hypers."));
    return 1;
  }

  json adaptors;
  json ad1{{"name", "AdaptorNeighbourList"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  json ad1b{{"name", "AdaptorCenterContribution"},
            {"initialization_arguments", {}}};
  json ad2{{"name", "AdaptorStrict"},
           {"initialization_arguments", {{"cutoff", cutoff}}}};
  adaptors.emplace_back(ad1);
  adaptors.emplace_back(ad1b);
  adaptors.emplace_back(ad2);

  std::cout << "Loading structures...";
  json structures{{"filename", structures_absolute_path}};
  auto manager =
      make_structure_manager_stack<StructureManagerCenters,
                                   AdaptorNeighbourList,
                                   AdaptorCenterContribution, AdaptorStrict>(
          structures, adaptors);
  std::cout << " finished." << std::endl;


  // load torch model 
  std::cout << "Loading model...";
  torch::jit::script::Module module;
  module = torch::jit::load(torch_model_absolute_path);
  std::cout << " finished." << std::endl;

  //// DUMMY DATA
  // TODO(sergey) code is just to have some functional input for testing but in the end
  // the remapped representation should be used
  std::cout << "Loading dummy...";
  std::vector<std::string> names;

  for (const auto & entry : std::filesystem::directory_iterator("/ssd/local/code/pytorch_prototype/meeting/inputs/")) {
      std::string path_now = entry.path();
      char first_symbol = 0;
      for (int i = path_now.size() - 1; i >= 0; --i) {
          if (path_now[i] == '/') {
              first_symbol = path_now[i + 1];
              break;
          }
      }
      if (first_symbol != '.') { //filter out hidden files
          names.push_back(path_now);
      }
  }
  std::vector<torch::jit::IValue> remapped_sph_exp;

  std::sort(names.begin(), names.end());
  for (unsigned int i = 0 ; i < names.size(); ++i) {
      remapped_sph_exp.push_back(read_tensor(names[i]));
      //std::cout << names[i] << '\n';
  }
  std::cout << " finished." << std::endl;
  //// DUMMY DATA

  auto start = std::chrono::high_resolution_clock::now();

  //std::cout << "Loading comp rep...";
  Representation_t representation{sph_hypers};
  representation.compute(manager);
  //std::cout << " finished." << std::endl;

  auto && expansions_coefficients{
      *manager->template get_property<Prop_t>(representation.get_name())};

  auto finish_sph_exp_computation = std::chrono::high_resolution_clock::now();
  // this should do a memcopy, transforms sparse object to non sparse
  Matrix_t feature_matrix = expansions_coefficients.get_features();
  // at this point the expansion_coefficents could be freed, not sure if needed for benchmarks

  // this should *not* do a memcopy, TODO(alex) double check this 
  torch::Tensor feature_matrix_torch = eigen_matrix_to_torch_tensor(feature_matrix);
  auto finish_eigen_to_torch_format = std::chrono::high_resolution_clock::now();
  
  // TODO(sergey)
  // remapped_sph_exp = ...
  // reshape it the way you want

  auto finish_sph_exp_remapping = std::chrono::high_resolution_clock::now();

  //std::cout << "Computing forward...";
  at::Tensor output = module.forward(remapped_sph_exp).toTensor();
  //std::cout << " finished " << std::endl;

  auto finish = std::chrono::high_resolution_clock::now();
  // Some kind of output is needed, otherwise things could be optimized out 
  std::cout << "energies: " << output << '\n';

  std::chrono::duration<double> elapsed_sph_exp_computation{finish_sph_exp_computation - start};
  std::chrono::duration<double> elapsed_eigen_to_torch_format{finish_eigen_to_torch_format - finish_sph_exp_computation};
  std::chrono::duration<double> elapsed_sph_exp_remapping{finish_sph_exp_remapping - finish_sph_exp_computation};
  std::chrono::duration<double> elapsed_torch_prediction{finish - finish_sph_exp_remapping};
  std::chrono::duration<double> elapsed_overall{finish - start};

  std::ofstream output_file (output_filename);
  
  json output_json = {
      {"sph_exp_time", elapsed_sph_exp_computation.count()},
      {"eigen_to_torch_format_time", elapsed_eigen_to_torch_format.count()},
      {"sph_exp_remapping_time", elapsed_sph_exp_remapping.count()},
      {"torch_prediction_time", elapsed_torch_prediction.count()},
      {"overall_time", elapsed_overall.count()},
    };
  if (output_file.is_open()) {
    output_file << output_json << std::endl;
    output_file.close();
  } else {
    std::cout << "Unable to open file";
  }
}
