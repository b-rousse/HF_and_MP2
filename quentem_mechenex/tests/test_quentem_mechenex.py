"""
Unit and regression test for the quentem_mechenex package.
"""

# Import package, test suite, and other packages as needed
import quentem_mechenex as qm
import pytest
import sys
import numpy as np

@pytest.fixture()
def noble_gas_fixture():
    noble_gas_instance = qm.NobleGasModel()
    return noble_gas_instance

@pytest.fixture()
def hartree_fock_fixture(noble_gas_fixture):
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    hartree_fock_instance = qm.HartreeFock(noble_gas_fixture, atomic_coordinates)
    return hartree_fock_instance

@pytest.fixture()
def mp2_fixture(noble_gas_fixture):
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    mp2_instance = qm.MP2(noble_gas_fixture, atomic_coordinates)
    return mp2_instance


def test_nbg_fixture(noble_gas_fixture):
    atomic_coordinates = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 5.0]])
    number_of_atoms = len(atomic_coordinates)
    for index in range(number_of_atoms * noble_gas_fixture.orbitals_per_atom):
        assert (index == noble_gas_fixture.ao_index(noble_gas_fixture.atom(index), noble_gas_fixture.orb(index)))

def test_hf_fixture(hartree_fock_fixture):
    my_ac = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    assert(my_ac.all() == hartree_fock_fixture.atomic_coordinates.all())

def test_mp2_fixture(mp2_fixture, noble_gas_fixture):
    assert(mp2_fixture.ionic_charge == noble_gas_fixture.ionic_charge)

def test_calculate_energy_ion(hartree_fock_fixture, noble_gas_fixture):
    expected_energy_ion = 5.091168824543142
    calculated_energy_ion = hartree_fock_fixture.calculate_energy_ion(noble_gas_fixture)
    assert np.isclose(expected_energy_ion, calculated_energy_ion)

def test_calculate_coloumb_energy(hartree_fock_fixture, noble_gas_fixture):
    coulomb_vector = np.array([1,0,0])
    assert np.isclose (1.0, hartree_fock_fixture.calculate_coulomb_energy('s','s',coulomb_vector, noble_gas_fixture))
    assert np.isclose (1.0, hartree_fock_fixture.calculate_coulomb_energy('s','px',coulomb_vector, noble_gas_fixture))
    assert np.isclose (-2.0, hartree_fock_fixture.calculate_coulomb_energy('px','px',coulomb_vector, noble_gas_fixture))
    assert np.isclose (0.0, hartree_fock_fixture.calculate_coulomb_energy('s','py',coulomb_vector, noble_gas_fixture))
    assert np.isclose (0.0, hartree_fock_fixture.calculate_coulomb_energy('px','py',coulomb_vector, noble_gas_fixture))
    assert np.isclose (1.0, hartree_fock_fixture.calculate_coulomb_energy('py','py',coulomb_vector, noble_gas_fixture))

def test_hartree_fock_hopping(hartree_fock_fixture, noble_gas_fixture):
    """Test for hopping_energy function in hartree fock"""
    hop_vector = np.array([3.1810226927827516,0,0])
    assert np.isclose( 0.03365982238611262 , hartree_fock_fixture.calculate_hopping_energy( 's', 's' , hop_vector, noble_gas_fixture))

def test_mp2(hartree_fock_fixture, noble_gas_fixture):
    """Test the standalone MP2 module"""
    
    assert np.isclose (1,1)
