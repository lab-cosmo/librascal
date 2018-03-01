#include <pybind11/pybind11.h>

using namespace pybind11::literals;
namespace py=pybind11;

extern void prot_dist_mat(py::module&);

PYBIND11_MODULE(_proteus, mod) {
  mod.doc() = "Hello, World!";
  // py::add_ostream_redirect(m, "ostream_redirect");
  prot_dist_mat(mod);
}
