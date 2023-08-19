#cython: language_level=3
"""MOAB Tag Class"""
from pymoab cimport moab


cdef class Tag:
    cdef moab.TagInfo * inst
    cdef moab.Tag * ptr


cdef class _TagArray:
    cdef moab.TagInfo ** inst
    cdef moab.Tag * ptr
