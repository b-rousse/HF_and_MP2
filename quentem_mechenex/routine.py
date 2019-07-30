"""
routine.py
HF and MP2 in python, with optimized Fock build in C++ wrapped with PYBIND11

Handles the primary functions
"""
import numpy as np
import quentem_mechenex.mp2 as mp2
import quentem_mechenex.mp2_no_hf as mp2_no_hf
import quentem_mechenex.hartree_fock as hf
import quentem_mechenex.noble_gas_model as noble_gas_model

MP2inheritsHF = True


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    NobleGasModel = noble_gas_model.NobleGasModel()
    atomic_coordinates = np.array([[0.0,0.0,0.0], [3.0,4.0,5.0]])
    if MP2inheritsHF:
        mp2_instance = mp2.MP2(NobleGasModel,atomic_coordinates)
        mp2_instance.density_matrix = mp2_instance.calculate_atomic_density_matrix(NobleGasModel)
        mp2_instance.density_matrix, mp2_instance.fock_matrix = mp2_instance.scf_cycle(NobleGasModel)
        mp2_instance.energy_scf = mp2_instance.calculate_energy_scf()
        mp2_instance.energy_ion = mp2_instance.calculate_energy_ion(NobleGasModel)
        print(F'The SCF energy is  {mp2_instance.energy_scf} and the ion energy is {mp2_instance.energy_ion} ')
        mp2_instance.mp2_energy = mp2_instance.calculate_energy_mp2()
        print(F'The MP2 energy is {mp2_instance.mp2_energy}')

    else:
        hartree_fock_instance = hf.HartreeFock(NobleGasModel, atomic_coordinates)
        hartree_fock_instance.density_matrix = hartree_fock_instance.calculate_atomic_density_matrix(NobleGasModel)
        hartree_fock_instance.density_matrix, hartree_fock_instance.fock_matrix = hartree_fock_instance.scf_cycle(NobleGasModel)
        hartree_fock_instance.energy_scf = hartree_fock_instance.calculate_energy_scf()
        hartree_fock_instance.energy_ion = hartree_fock_instance.calculate_energy_ion(NobleGasModel)
        print(F'The SCF energy is  {hartree_fock_instance.energy_scf} and the ion energy is {hartree_fock_instance.energy_ion} ')
        mp2_instance = mp2_no_hf.MP2NoHF(hartree_fock_instance, NobleGasModel)
        mp2_instance.mp2_energy = mp2_instance.calculate_energy_mp2()
        print(F'The MP2 energy is {mp2_instance.mp2_energy}')

