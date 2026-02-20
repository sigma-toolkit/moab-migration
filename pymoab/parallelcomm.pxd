#cython: language_level=3

"""Implements parallel communication functionality for MOAB."""

from pymoab cimport moab
from mpi4py.MPI cimport Comm as MPIComm
from .rng cimport Range
from .core cimport Core

# Python class declaration
cdef class ParallelComm:
    """Parallel communication class for MOAB."""
    cdef moab.ParallelComm* inst
    cdef moab.Core* core
    cdef MPIComm comm
