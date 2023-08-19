#cython: language_level=3
"""Implements range functionality."""

from pymoab cimport moab

cdef class Range:

    cdef moab.Range * inst
