#pragma once
#include <Eigen/Dense>
#include <string>
#include <vector>

int atom(int ao_index, int orbitals_per_atom);
std::string orb(int ao_index, std::vector<std::string> orbital_types, int orbitals_per_atom);
int ao_index(int atom_p, std::string orb_p, std::vector<std::string> orbital_types, int orbitals_per_atom);
double chi_on_atom(std::string orb1, std::string orb2, std::string orb3, std::map<std::string, double> model_parameters);
Eigen::MatrixXd calculate_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix, std::map<std::string, double> model_parameters, int orbitals_per_atom);
Eigen::MatrixXd test_fock_matrix_fast(Eigen::MatrixXd hamiltonian_matrix, Eigen::MatrixXd interaction_matrix, Eigen::MatrixXd density_matrix, std::map<std::string, double> model_parameters, int orbitals_per_atom);
