#include <torch/script.h> // One-stop header.

#include <iostream>
#include <memory>
#include <filesystem>

std::ifstream get_file_input(std::string file_path) {
    std::ifstream ifs (file_path, std::ifstream::in);
    //std::ios_base::sync_with_stdio(false);
    ifs.tie(nullptr);
    return ifs;
}

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
int main(int argc, const char* argv[]) {
    if (argc != 3) {
        std::cerr << "usage: run <path-to-exported-script-module> <path-to-inputs-dir>\n";
        return -1;
    }
    auto inputs_path = std::string(argv[2]);

    std::vector<std::string> names;

    for (const auto & entry : std::filesystem::directory_iterator(inputs_path)) {
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
    std::vector<torch::jit::IValue> inputs;

    std::sort(names.begin(), names.end());
    for (int i = 0 ; i < names.size(); ++i) {
        inputs.push_back(read_tensor(names[i]));
        //std::cout << names[i] << '\n';
    }


    torch::jit::script::Module module;
    module = torch::jit::load(argv[1]);
    at::Tensor output = module.forward(inputs).toTensor();
    std::cout << "energies: " << output << '\n';

    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module = torch::jit::load(argv[1]);
        //at::Tensor output = module.forward(inputs).toTensor();
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model\n";
        return -1;
    }

    std::cout << "ok\n";
}
