#include "rascal/representations/calculator_spherical_expansion.hh"
#include "rascal/structure_managers/adaptor_center_contribution.hh"
#include "rascal/structure_managers/adaptor_neighbour_list.hh"
#include "rascal/structure_managers/adaptor_strict.hh"
#include "rascal/structure_managers/atomic_structure.hh"
#include "rascal/structure_managers/structure_manager_collection.hh"
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
#include <math.h>
#include <Eigen/StdVector>

#include <torch/script.h> // One-stop header.
#include <memory>
#include <filesystem>

using namespace rascal;  // NOLINT

using Matrix_t = math::Matrix_t;
using Matrix_Ref = math::Matrix_Ref;

using Representation_t = CalculatorSphericalExpansion;
using Manager_t = AdaptorStrict<
    AdaptorCenterContribution<AdaptorNeighbourList<StructureManagerCenters>>>;
using ManagerCollection_t = AdaptorStrict<
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
  std::string timings_filename{argv[1]};
  std::string energies_filename{argv[1]};

  // removes .json from filename
  const std::string ext(".json");
  if ( timings_filename != ext &&
       timings_filename.size() > ext.size() &&
       timings_filename.substr(timings_filename.size() - ext.size()) == ".json" )
  {
     // if so then strip them off
     timings_filename = timings_filename.substr(0, timings_filename.size() - ext.size());
  }
  timings_filename.append("_timings.json");

  if ( energies_filename != ext &&
       energies_filename.size() > ext.size() &&
       energies_filename.substr(energies_filename.size() - ext.size()) == ".json" )
  {
     // if so then strip them off
     energies_filename = energies_filename.substr(0, energies_filename.size() - ext.size());
  }
  energies_filename.append("_energies.json");

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

  int max_radial;
  if (sph_hypers.count("max_radial")) {
    max_radial = sph_hypers.at("max_radial").get<int>();
  } else {
    throw std::runtime_error(std::string("Could not find max radial in spherical expansion hypers."));
    return 1;
  }

  int max_angular;
  if (sph_hypers.count("max_angular")) {
    max_angular = sph_hypers.at("max_angular").get<int>();
  } else {
    throw std::runtime_error(std::string("Could not find max angular in spherical expansion hypers."));
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
  auto managers{ManagerCollection<StructureManagerCenters, AdaptorNeighbourList, AdaptorCenterContribution, AdaptorStrict>(adaptors)};
  managers.add_structures(structures_absolute_path, 0, -1); // extendable with arguments: int start = 0, int length = -1)

  auto manager =
      make_structure_manager_stack<StructureManagerCenters, AdaptorNeighbourList, AdaptorCenterContribution, AdaptorStrict>(
          structures, adaptors);
  std::cout << " finished. Number of structures loaded: " << managers.size() << std::endl;

  // load torch model 
  std::cout << "Loading model...";
  torch::jit::script::Module module;
  module = torch::jit::load(torch_model_absolute_path);
  std::cout << " finished." << std::endl;

  std::vector<torch::jit::IValue> torch_model_inputs;
  // angular channels, species per center, structuer idx per center
  torch_model_inputs.reserve(max_angular+1 + 2);

  //// DUMMY DATA not used anymore
  //std::vector<torch::jit::IValue> remapped_sph_exp;
  //// code is just to have some functional input for testing but in the end
  //// the remapped representation should be used now
  //std::cout << "Loading dummy...";
  //std::vector<std::string> names;

  //for (const auto & entry : std::filesystem::directory_iterator("/home/alexgo//lib/pytorch_prototype/meeting/inputs/")) {
  //  std::string path_now = entry.path();
  //  char first_symbol = 0;
  //  for (int i = path_now.size() - 1; i >= 0; --i) {
  //      if (path_now[i] == '/') {
  //          first_symbol = path_now[i + 1];
  //          break;
  //      }
  //  }
  //  if (first_symbol != '.') { //filter out hidden files
  //      names.push_back(path_now);
  //  }
  //}

  //std::sort(names.begin(), names.end());
  //for (unsigned int i = 0 ; i < names.size(); ++i) {
  //  torch::Tensor t{read_tensor(names[i])};
  //  std::cout << "i: "<< i << " size: (";
  //  for (int64_t d = 0; d < t.dim(); ++d) {
  //    std::cout << t.size(d) << " ";
  //  }
  //  std::cout << ")" << std::endl;
  //  remapped_sph_exp.push_back(t);
  //}
  //std::cout << " finished." << std::endl;
  //// DUMMY DATA
  Representation_t representation{sph_hypers};

  auto start = std::chrono::high_resolution_clock::now();

  //std::cout << "Computing comp rep...";
  representation.compute(managers);
  auto finish_sph_exp_computation = std::chrono::high_resolution_clock::now();
  //std::cout << " finished." << std::endl;

  //// get_features for all structures and merge together like in pybind 
  auto property_name{managers.get_calculator_name(representation, false)};

  const auto & property_ =
      *managers[0]->template get_property<Prop_t>(property_name);
  // assume inner_size is consistent for all managers
  int inner_size{property_.get_nb_comp()};
  //std::cout << "rascal inner size: " << inner_size << std::endl;
  //std::cout << "should be equal to : " << max_radial * static_cast<int>(pow(max_angular + 1, 2)) << std::endl;

  auto all_keys{property_.get_keys()};

  // get features from rascal 
  auto n_rows{
      managers.template get_number_of_elements<Prop_t>(property_name)};
  size_t n_cols{all_keys.size() * inner_size};
  Matrix_t feature_matrix{Matrix_t(n_rows, n_cols)};
  int i_row{0};
  for (auto & manager : managers) {
    const auto & property =
        *manager->template get_property<Prop_t>(property_name);
    auto n_rows_manager = property.size();
    property.fill_dense_feature_matrix(
        feature_matrix.block(i_row, 0, n_rows_manager, n_cols), all_keys);
    i_row += n_rows_manager;
  }
  auto finish_densify = std::chrono::high_resolution_clock::now();
  //std::cout << "Densify finished" << std::endl;

  // remapping
  // only required for fixed-size vectorizable Eigen types
  // https://eigen.tuxfamily.org/dox/group__TopicFixedSizeVectorizable.html
  //std::vector<Matrix_t, Eigen::aligned_allocator<Matrix_t> > feature_matrices{};
  std::vector<Matrix_t> feature_matrices{};
  feature_matrices.reserve(max_angular + 1);
  int nb_species{static_cast<int>(all_keys.size())};
  int max_angular_size{static_cast<int>(pow(max_angular + 1, 2))};
  for (int angular_l = 0;  angular_l < max_angular + 1; angular_l++) {
    //std::cout << "size: " << feature_matrices.size() << std::endl;
    //std::cout << "capacity: " << feature_matrices.capacity() << std::endl;
    int n_cols{nb_species * max_radial * (2 * angular_l + 1)};
    int angular_size{2 * angular_l + 1};
    //std::cout << "n_cols " << n_cols << std::endl;
    feature_matrices.emplace_back(Matrix_t(n_rows, n_cols));
    //std::cout << "angular_l " << angular_l << std::endl;
    //std::cout << "feature_matrices shape " << feature_matrices.back().rows() << ", " << feature_matrices.back().cols() << std::endl;
    //std::cout << "feature_matrix shape " << feature_matrix.rows() << ", " << feature_matrix.cols() << std::endl;
    for (int i_species = 0; i_species <  static_cast<int>(all_keys.size()); i_species++) {
      for (int radial_n = 0; radial_n < max_radial; radial_n++) {
        int col_offset{i_species * max_radial * max_angular_size + radial_n * max_angular_size};
        //std::cout << "i_species * max_radial * max_angular_size + radial_n * max_angular_size " << i_species << " * " << max_radial << " * " << max_angular_size << " + " << radial_n << " * " << max_angular_size << std::endl;
        int matrices_col_offset{i_species * max_radial + radial_n};
        //std::cout << "angular_l, i_species, radial_n " << angular_l << ", " << i_species << ", " << radial_n << std::endl;
        //std::cout << "matrices_col_offset: " << matrices_col_offset << std::endl;
        //std::cout << "col_offset: " << col_offset << std::endl;
        //std::cout << "block feature_matrices " <<
        //feature_matrices.back().block(0, matrices_col_offset, feature_matrix.rows(), angular_size).rows()
        //<< " " <<
        //feature_matrices.back().block(0, matrices_col_offset, feature_matrix.rows(), angular_size).cols()
        //<< std::endl;

        //std::cout << "block " <<
        //feature_matrix.block(0, col_offset, feature_matrix.rows(), angular_size).rows()
        //<< " " <<
        //feature_matrix.block(0, col_offset, feature_matrix.rows(), angular_size).cols()
        //<< std::endl;
        feature_matrices.back().block(0, matrices_col_offset, feature_matrix.rows(), angular_size) = feature_matrix.block(0, col_offset, feature_matrix.rows(), angular_size);
      }
    }
    //std::cout << "shape: (" << feature_matrices.back().rows() << ", " << feature_matrices.back().cols() << ")" << std::endl;
    //std::cout << std::endl;
  }
  auto finish_remapping_angular = std::chrono::high_resolution_clock::now();
  //std::cout << "Remapping finished" << std::endl;
  ////


  auto options_features =
    torch::TensorOptions()
      .dtype(torch::kDouble)
      .layout(torch::kStrided)
      .device(torch::kCPU) // kCUDA, 1
      .requires_grad(true);

  for (int angular_l = 0;  angular_l < max_angular + 1; angular_l++) {
    at::Tensor tensor = torch::from_blob(feature_matrices.at(angular_l).data(), {feature_matrices.at(angular_l).rows(), feature_matrices.at(angular_l).cols()}, options_features).reshape({static_cast<long int>(n_rows), static_cast<long int>(all_keys.size() * max_radial), 2 * angular_l + 1});
    //std::cout << "tensor.dim()" << tensor.dim() << std::endl;
    torch_model_inputs.push_back(tensor);
  }

  auto finish_eigen_to_torch_format = std::chrono::high_resolution_clock::now();

  auto options =
    torch::TensorOptions()
      .dtype(torch::kInt32)
      .layout(torch::kStrided)
      .device(torch::kCPU) // kCUDA, 1
      .requires_grad(false); // ? don't know, I would think no
  torch::Tensor species_vec{torch::empty({static_cast<long int>(n_rows)}, options)};
  torch::Tensor structure_idx_vec{torch::empty({static_cast<long int>(n_rows)}, options)};
  int i_manager{0};
  i_row = 0;
  for (auto & manager : managers) {
    for (auto center : manager) {
      species_vec[i_row] = center.get_atom_type();
      structure_idx_vec[i_row] = i_manager;
      i_row++;
    }
    i_manager++;
  }
  torch_model_inputs.push_back(species_vec);
  torch_model_inputs.push_back(structure_idx_vec);

  //for (unsigned int i = 0 ; i < torch_model_inputs.size(); ++i) {
  //  std::cout << "i: "<< i << " size: (";
  //  for (int64_t d = 0; d < torch_model_inputs.at(i).toTensor().dim(); ++d) {

  //    std::cout << torch_model_inputs.at(i).toTensor().size(d) << " ";
  //  }
  //  std::cout << ")" << std::endl;
  //}

  auto finish_sph_exp_remapping = std::chrono::high_resolution_clock::now();

  //std::cout << "Computing forward...";
  at::Tensor energies = module.forward(torch_model_inputs).toTensor();
  //std::cout << " finished " << std::endl;

  auto finish = std::chrono::high_resolution_clock::now();
  // Output should be printed or stored, otherwise some steps could be optimized out 
  std::ofstream energies_file (energies_filename);
  if (energies_file.is_open()) {
    energies_file << energies << std::endl;
    energies_file.close();
  } else {
    std::cout << "Unable to open file";
  }

  std::chrono::duration<double> elapsed_sph_exp_computation{finish_sph_exp_computation - start};
  std::chrono::duration<double> elapsed_densify{finish_densify - finish_sph_exp_computation};
  std::chrono::duration<double> elapsed_remapping_angular{finish_remapping_angular - finish_densify};
  std::chrono::duration<double> elapsed_eigen_to_torch_format{finish_eigen_to_torch_format - finish_remapping_angular};
  std::chrono::duration<double> elapsed_sph_exp_remapping{finish_sph_exp_remapping - finish_sph_exp_computation};
  std::chrono::duration<double> elapsed_torch_prediction{finish - finish_sph_exp_remapping};
  std::chrono::duration<double> elapsed_overall{finish - start};

  
  json timings_json = {
      {"sph_exp_time", elapsed_sph_exp_computation.count()},
      {"eigen_to_torch_format_time", elapsed_eigen_to_torch_format.count()},
      {"sph_exp_remapping_time", elapsed_sph_exp_remapping.count()},
      {"torch_prediction_time", elapsed_torch_prediction.count()},
      {"overall_time", elapsed_overall.count()},
    };

  std::cout << timings_json << std::endl;
  std::cout << "Write output file...";
  std::ofstream timings_file (timings_filename);
  if (timings_file.is_open()) {
    timings_file << timings_json << std::endl;
    timings_file.close();
  } else {
    std::cout << "Unable to open file";
  }
  std::cout << " finished." << std::endl;
  return 0;
}
