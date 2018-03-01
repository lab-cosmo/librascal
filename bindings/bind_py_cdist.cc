#include <pybind11/pybind11.h>
#include "cdist.h"

namespace py=pybind11;

void prot_dist_mat(py::module& m)
{
	m.def("cdist",&cdist);
}

