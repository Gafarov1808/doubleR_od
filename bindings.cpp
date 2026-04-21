#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "doubleR.h"

namespace py = pybind11;

PYBIND11_MODULE(doubleR_module, m){
    m.doc() = "Double R Iteration Orbit Determination Module";

    py::class_<DoubleRIteration>(m, "DoubleRIteration")
        .def(py::init<
        std::array<double, 3>,
        std::array<double, 3>,
        std::array<double, 3>,
        std::array<double, 3>,
        std::array<double, 3>,
        std::array<double, 3>,
        double,
        double
    >())
    .def("calc_func", &DoubleRIteration::calc_func)
    .def("calc_der", &DoubleRIteration::calc_der)
    .def("solver", &DoubleRIteration::solver)
    .def("get_state", &DoubleRIteration::get_state)
    .def("set_state", [](DoubleRIteration& self, std::array<double, 6> state){
        self = {state[0], state[1], state[2], state[3], state[4], state[5]};
    })
    .def("print_state", &DoubleRIteration::print_state)
    .def("print_elements", &DoubleRIteration::print_elements)
    ;
}