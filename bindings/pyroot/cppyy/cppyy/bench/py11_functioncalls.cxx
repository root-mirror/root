#include <pybind11/pybind11.h>
#include "functioncalls.h"


namespace py = pybind11;

PYBIND11_MODULE(py11_functioncalls, m) {
//- group: empty-free --------------------------------------------------------
    m.def("empty_call", &empty_call);

//- group: empty-inst --------------------------------------------------------
    py::class_<EmptyCall>(m, "EmptyCall")
        .def(py::init<>())
        .def("empty_call", &EmptyCall::empty_call);


//- group: builtin-args-free -------------------------------------------------
    py::class_<Value>(m, "Value")
        .def(py::init<>());

    m.def("take_an_int",   &take_an_int);
    m.def("take_a_double", &take_a_double);
    m.def("take_a_struct", &take_a_struct);

//- group: builtin-args-free -------------------------------------------------
    py::class_<TakeAValue>(m, "TakeAValue")
        .def(py::init<>())
        .def("take_an_int",   &TakeAValue::take_an_int)
        .def("take_a_double", &TakeAValue::take_a_double)
        .def("take_a_struct", &TakeAValue::take_a_struct)

//- group: builtin-args-pass -------------------------------------------------
        .def("pass_int",      &TakeAValue::pass_int);


//- group: do-work-free ------------------------------------------------------
    m.def("do_work", &do_work);

//- group: do-work-inst ------------------------------------------------------
    py::class_<DoWork>(m, "DoWork")
        .def(py::init<>())
        .def("do_work", &DoWork::do_work);


//- group: overload-inst -----------------------------------------------------
    py::class_<OverloadedCall>(m, "OverloadedCall")
        .def(py::init<>())
        .def("add_it", (double (OverloadedCall::*)(int, int))     &OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(short))        &OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(long))         &OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(int, int, int))&OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(double))       &OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(float))        &OverloadedCall::add_it)
        .def("add_it", (double (OverloadedCall::*)(int))          &OverloadedCall::add_it);
}
