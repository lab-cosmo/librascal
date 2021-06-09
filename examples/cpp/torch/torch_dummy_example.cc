#include "rascal/representations/calculator_spherical_expansion.hh"


#include <torch/torch.h>
#include <torch/script.h> // One-stop header.

#include <iostream>
#include <memory>
#include <filesystem>

int main() {
  torch::Tensor tensor = torch::rand({2, 3});
  std::cout << tensor << std::endl;
}
