"""
quentem_mechenex
HF and MP2 in python, with optimized Fock build in C++ wrapped with PYBIND11
"""

# Add imports here
from .routine import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
