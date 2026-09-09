#!/usr/bin/env python3

import sys

try:
    from mpi4py import MPI
except ImportError:
    print("MPI not available - skipping parallel tests")
    sys.exit(0)

from pymoab import core, parallelcomm
from parallel_driver import run_parallel_tests, CHECK_EQ

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def test_parallelcomm_lifecycle():
    mb = core.Core()
    pcomm = parallelcomm.ParallelComm(mb, comm)
    CHECK_EQ(pcomm.rank, rank)
    CHECK_EQ(pcomm.size, size)
    wrapped_comm = pcomm.comm
    CHECK_EQ(wrapped_comm.Get_rank(), rank)
    CHECK_EQ(wrapped_comm.Get_size(), size)



def main():
    return run_parallel_tests([
        test_parallelcomm_lifecycle,
    ])


if __name__ == "__main__":
    sys.exit(main())
