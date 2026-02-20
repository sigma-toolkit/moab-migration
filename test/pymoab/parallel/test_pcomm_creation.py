#!/usr/bin/env python3
"""
Test ParallelComm creation and basic properties.
"""

import sys
import os
import numpy as np

# Check if MPI is available first
try:
    from mpi4py import MPI
except ImportError:
    print("MPI not available - skipping parallel tests")
    print("To run parallel tests, ensure:")
    print("  1. MOAB was built with -DENABLE_MPI=ON")
    print("  2. mpi4py is installed")
    sys.exit(0)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pymoab'))

from pymoab import core, parallelcomm, config
from parallel_driver import run_parallel_tests, CHECK_EQ, CHECK, CHECK_PARALLEL, CHECK_PARALLEL_EQ

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Check capabilities
def test_capabilities():
    """Test MPI capability detection."""
    # Test basic MPI capability detection
    if config.MOAB_MPI_ENABLED:
        CHECK(config.has_basic_mpi())
        
        # Test MPI I/O capability if available
        if config.has_mpi_io():
            print(f"Process {rank}: MPI I/O support detected")
        else:
            print(f"Process {rank}: MPI I/O support not detected")

    # Test runtime capability detection
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    # Test capability properties
    if config.MOAB_MPI_ENABLED:
        CHECK(pcomm.has_basic_mpi)
        if config.has_mpi_io():
            CHECK(pcomm.has_mpi_io)
        CHECK(pcomm.has_full_mpi == (pcomm.has_basic_mpi and pcomm.has_mpi_io))
    else:
        CHECK(not pcomm.has_basic_mpi)
        CHECK(not pcomm.has_mpi_io)
        CHECK(not pcomm.has_full_mpi)

    print(f"Process {rank}: Capability detection successful")

def test_pcomm_creation():
    """Test ParallelComm creation."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    CHECK(pcomm is not None)
    
    # Test capability properties after creation
    if config.MOAB_MPI_ENABLED:
        CHECK(pcomm.has_basic_mpi)
        if config.has_mpi_io():
            CHECK(pcomm.has_mpi_io)
        CHECK(pcomm.has_full_mpi == (pcomm.has_basic_mpi and pcomm.has_mpi_io))

    print(f"Process {rank}: ParallelComm creation successful")


def test_pcomm_rank_size():
    """Test rank and size properties."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    CHECK_EQ(pcomm.rank, rank)
    CHECK_EQ(pcomm.size, size)


def test_pcomm_comm():
    """Test that communicator is stored correctly."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    CHECK(pcomm.comm is not None)
    CHECK_EQ(pcomm.comm.Get_rank(), rank)
    CHECK_EQ(pcomm.comm.Get_size(), size)


def test_get_comm_procs():
    """Test getting all processor ranks."""
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    
    procs = pcomm.get_comm_procs()
    CHECK_EQ(len(procs), size)
    CHECK_PARALLEL(len(procs) == size, "get_comm_procs size mismatch")


def test_mpi_comm_passing():
    """Test that MPI communicator works."""
    py_rank, py_size = parallelcomm.test_mpi_comm_passing(comm)
    CHECK_EQ(py_rank, rank)
    CHECK_EQ(py_size, size)


def main():
    tests = [
        test_capabilities,
        test_pcomm_creation,
        test_pcomm_rank_size,
        test_pcomm_comm,
        test_get_comm_procs,
        test_mpi_comm_passing,
    ]
    return run_parallel_tests(tests)

if __name__ == "__main__":
    sys.exit(main())
