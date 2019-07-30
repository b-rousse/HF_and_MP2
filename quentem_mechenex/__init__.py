"""
quentem_mechenex
HF and MP2 in python, with optimized Fock build in C++ wrapped with PYBIND11
"""

# Add imports here
from .routine import *
from .mp2 import *
from .mp2_no_hf import *
from .noble_gas_model import * 
from .hartree_fock import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
