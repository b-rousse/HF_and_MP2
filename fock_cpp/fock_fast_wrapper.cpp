#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "fock_fast.hpp"

PYBIND11_MODULE(fock_fast, m)
{
    m.doc() = "Module for C++ implementations of the QM Project for the MolSSI SSS 2019 (c)";
    m.def("calculate_fock_matrix_fast", calculate_fock_matrix_fast, "Faster generation of the Fock matrix without einsum");
    m.def("test_fock_matrix_fast", test_fock_matrix_fast, "Same loops but no calculations as the Fock matrix generation for the purpose of timing.");
}
