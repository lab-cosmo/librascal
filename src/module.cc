#include <pybind11/iostream.h>
#include <cdist.cc>

namespace py = pybind11;

PYBIND11_MODULE(proteus, m){
	py::add_ostream_redirect(m, "ostream_redirect");
	prot_dist_mat(m);
}
