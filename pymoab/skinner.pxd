#cython: language_level=3
"""Implements the skinner functionality to find the geometric skin entities."""

from pymoab cimport moab

cdef class Skinner:

    cdef moab.Skinner * inst
    cdef moab.Interface* interface
