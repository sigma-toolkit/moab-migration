#!/usr/bin/env python3
"""
Test script for PyMOAB parallel interface.
This tests the fixed parallel interface that uses the correct MPI pattern.
"""

import os
import sys
import unittest

try:
    from mpi4py import MPI
    from pymoab import config, core, parallelcomm
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure PyMOAB is properly installed and compiled.")
    sys.exit(1)

def test_capability_detection():
    """Test capability detection."""
    try:
        # Test Python-level capability detection
        from pymoab import config
        if config.MOAB_MPI_ENABLED:
            print(f"Python: Basic MPI support detected")
            if config.has_mpi_io():
                print(f"Python: MPI I/O support detected")
            else:
                print(f"Python: MPI I/O support not detected")
        else:
            print(f"Python: MPI support not detected")
        
        # Test runtime capability detection
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, MPI.COMM_WORLD)
        
        if config.MOAB_MPI_ENABLED:
            print(f"Python: Runtime basic MPI support: {pcomm.has_basic_mpi}")
            if config.has_mpi_io():
                print(f"Python: Runtime MPI I/O support: {pcomm.has_mpi_io}")
            print(f"Python: Runtime full MPI support: {pcomm.has_full_mpi}")
        
        return True
    except Exception as e:
        print(f"✗ Capability detection failed: {e}")
        return False

def test_parallel_comm_creation():
    """Test creating a ParallelComm object."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"Process {rank}/{size}: Testing ParallelComm creation")

    try:
        # Create MOAB Core
        mb = core.Core()

        # Create ParallelComm - this should work without MPI initialization errors
        pcomm = parallelcomm.ParallelComm(mb, comm)

        print(f"Process {rank}: ✓ ParallelComm created successfully")

        # Test basic properties
        pcomm_rank = pcomm.rank
        pcomm_size = pcomm.size

        print(f"Process {rank}: PyMOAB rank {pcomm_rank}, size {pcomm_size}")
        print(f"Process {rank}: Python rank {rank}, size {size}")

        # Verify rank and size match
        if pcomm_rank == rank and pcomm_size == size:
            print(f"Process {rank}: ✓ Rank and size match correctly")
            return True
        else:
            print(f"Process {rank}: ✗ Rank and size mismatch")
            return False

    except Exception as e:
        print(f"Process {rank}: ✗ Error creating ParallelComm: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_mpi_comm_passing():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"Process {rank}: Testing MPI communicator passing")

    try:
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(mb, comm)
        wrapped_comm = pcomm.comm

        print(f"Process {rank}: ✓ MPI communicator stored on ParallelComm")
        print(f"Process {rank}: Wrapped rank {wrapped_comm.Get_rank()}, size {wrapped_comm.Get_size()}")

        if wrapped_comm.Get_rank() == rank and wrapped_comm.Get_size() == size:
            print(f"Process {rank}: ✓ Results match correctly")
            return True

        print(f"Process {rank}: ✗ Results mismatch")
        return False

    except Exception as e:
        print(f"Process {rank}: ✗ Error in MPI communicator passing test: {e}")
        import traceback
        traceback.print_exc()
        return False

class ParallelIOTestCase(unittest.TestCase):
    def test_capability_detection_unittest(self):
        self.assertTrue(test_capability_detection())

    def test_mpi_comm_passing_unittest(self):
        self.assertTrue(test_mpi_comm_passing())

    def test_parallel_comm_creation_unittest(self):
        self.assertTrue(test_parallel_comm_creation())


def main():
    """Main test function."""
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ParallelIOTestCase)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
