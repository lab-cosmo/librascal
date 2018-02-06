#include <pybind11/pybind11.h>


int add(int i, int j) {
    return i + j;
};

int subtract(int i, int j) {
    return i - j;
}


namespace py = pybind11;

void export_proteus_math_function(py::module& m){

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", &substract, R"pbdoc(
        subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    }



