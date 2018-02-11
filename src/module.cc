#include "proteusMATH.h"
#include "basic_types.h"
#include <pybind11/iostream.h>

namespace py = pybind11;

PYBIND11_MODULE(proteus, m){

    py::add_ostream_redirect(m, "ostream_redirect");

    export_proteus_math_function(m);
    export_proteus_basic_types(m);

    }
