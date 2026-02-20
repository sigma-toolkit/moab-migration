#!/usr/bin/env python3
"""
Test script for PyMOAB parallel interface.
This tests the fixed parallel interface that uses the correct MPI pattern.
"""

import sys
import os

# Add the pymoab directory to Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pymoab'))

try:
    from mpi4py import MPI
    from pymoab import core, parallelcomm
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
        from pymoab import core, parallelcomm
        mb = core.Core()
        pcomm = parallelcomm.ParallelComm(MPI.COMM_WORLD)
        
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
    """Test the MPI communicator passing functionality."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    print(f"Process {rank}: Testing MPI communicator passing")

    try:
        # Test the MPI Comm passing function
        py_rank, py_size = parallelcomm.test_mpi_comm_passing(comm)

        print(f"Process {rank}: ✓ MPI communicator passing test successful")
        print(f"Process {rank}: Python rank {py_rank}, size {py_size}")

        # Verify results
        if py_rank == rank and py_size == comm.Get_size():
            print(f"Process {rank}: ✓ Results match correctly")
            return True
        else:
            print(f"Process {rank}: ✗ Results mismatch")
            return False

    except Exception as e:
        print(f"Process {rank}: ✗ Error in MPI communicator passing test: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main test function."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"=== PyMOAB Parallel Interface Test ===")
    print(f"Process {rank}/{size} starting tests")

    # Test 1: Capability detection
    print(f"\n--- Test 1: Capability Detection ---")
    test1_success = test_capability_detection()

    # Test 2: MPI communicator passing
    print(f"\n--- Test 2: MPI Comm Passing ---")
    test2_success = test_mpi_comm_passing()

    # Test 3: ParallelComm creation
    print(f"\n--- Test 3: ParallelComm Creation ---")
    test3_success = test_parallel_comm_creation()

    # Final verification
    if test1_success and test2_success and test3_success:
        print(f"Process {rank}: ✓ ALL TESTS PASSED")
        print(f"Process {rank}: ✓ PyMOAB parallel interface is working correctly!")
        return 0
    else:
        print(f"Process {rank}: ✗ SOME TESTS FAILED")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)