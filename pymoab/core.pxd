#cython: language_level=3
""" Core Cython Header """

from pymoab cimport moab

cdef class Core:
    cdef moab.Core* inst
