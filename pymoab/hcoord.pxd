#cython: language_level=3
"""MOAB Homogenous Coordinate class"""

from pymoab cimport moab

cdef class HomCoord:
    cdef moab.HomCoord *inst
