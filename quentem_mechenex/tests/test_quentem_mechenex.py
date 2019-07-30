"""
Unit and regression test for the quentem_mechenex package.
"""

# Import package, test suite, and other packages as needed
import quentem_mechenex as qm
import pytest
import sys
import numpy as np

def test_calculate_energy_ion():
    noble_gas_instance = qm.NobleGasModel()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    hartree_fock_instance = qm.HartreeFock(noble_gas_instance, atomic_coordinates)
    #hartree_fock_instance.density_matrix = hartree_fock_instance.calculate_atomic_density_matrix(noble_gas_instance)
    #hartree_fock_instance.density_matrix, hartree_fock_instance.fock_matrix = hartree_fock_instance.scf_cycle(noble_gas_instance)
    #hartree_fock_instance.energy_scf = hartree_fock_instance.calculate_energy_scf()
    #hartree_fock_instance.energy_ion = hartree_fock_instance.calculate_energy_ion(noble_gas_instance)
    #print(F'The SCF energy is  {hartree_fock_instance.energy_scf} and the ion energy is {hartree_fock_instance.energy_ion} ')
    #mp2_instance = qm.MP2(NobleGasModel,atomic_coordinates)
    #mp2_instance.mp2_energy = mp2_instance.calculate_energy_mp2()
    #print(F'The MP2 energy is {mp2_instance.mp2_energy}')
    #mp2_standalone_instance = qm.MP2NoHF(hartree_fock_instance, NobleGasModel)
    #mp2_iso_energy = mp2_standalone_instance.calculate_energy_mp2()
    expected_energy_ion = 5.091168824543142
    calculated_energy_ion = hartree_fock_instance.calculate_energy_ion (noble_gas_instance)
    assert np.isclose (expected_energy_ion, calculated_energy_ion)

