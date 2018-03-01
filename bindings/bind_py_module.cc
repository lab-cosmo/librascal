#include <pybind11/pybind11.h>

using namespace pybind11::literals;
namespace py=pybind11;

PYBIND11_MODULE(_proteus, mod) {
  mod.doc() = "Hello, World!";
}
