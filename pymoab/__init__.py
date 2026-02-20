""" 
Name

PyMOAB: A Python interface to Argonne National Lab's Mesh Oriented dAtaBase (MOAB)
====

Description
-----------

PyMOAB provides a means of interactively interrogating, modifying, and generating 
MOAB mesh files.

Much of the core functionality of MOAB has been implemented in the core.Core
module. Interaction with this is intended to be largely analogous to interaction
with MOAB via its native C++ API though some modifications have been made to
allow for interaction with the interface using the various native Python
constructs such as lists, tuples, etc.

MPI Support
-----------
If PyMOAB was built with MPI support (MOAB compiled with -DENABLE_MPI=ON and 
mpi4py available), the parallelcomm module provides parallel communication
functionality. Check pymoab.config.is_mpi_enabled() to determine if MPI
support is available.
"""

from importlib.metadata import version, PackageNotFoundError
from .paths import *

try:
    __version__ = version('MOAB')
except PackageNotFoundError:
    __version__ = "unknown"

# Import config for MPI availability check
try:
    from . import config
    def is_mpi_enabled():
        return config.is_mpi_enabled()
except ImportError:
    def is_mpi_enabled():
        return False

