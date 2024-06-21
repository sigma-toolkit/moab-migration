#cython: language_level=3
""" MOAB MeshTopoUtil """

from pymoab cimport moab

cdef class MeshTopoUtil:

    cdef moab.MeshTopoUtil* inst
    cdef moab.Interface* interface
