#cython: language_level=3
"""MOAB Structured Mesh Interface"""

from pymoab cimport moab

cdef class ScdParData:
    cdef moab.ScdParData *inst

cdef class ScdInterface:

    cdef moab.ScdInterface * inst
    cdef moab.Interface * interface
    cdef object _core_ref  # prevent Core GC before ScdInterface

cdef class ScdBox:
    cdef moab.ScdBox* inst
    cdef object _scd_ref   # prevent ScdInterface GC before ScdBox
