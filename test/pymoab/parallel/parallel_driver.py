import sys
import traceback
import numpy as np

# Check if MPI is available
try:
    from mpi4py import MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False

if MPI_AVAILABLE:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
else:
    rank = 0
    size = 1
    comm = None

class colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def check_mpi_enabled():
    """Check if MPI is enabled, exit gracefully if not."""
    if not MPI_AVAILABLE:
        print("MPI not available - skipping parallel tests")
        print("To run parallel tests, ensure:")
        print("  1. MOAB was built with -DENABLE_MPI=ON")
        print("  2. mpi4py is installed")
        sys.exit(0)

def run_parallel_tests(test_list):
    """Run tests with MPI wrapper - results only from rank 0."""
    check_mpi_enabled()
    
    ret_val = 0
    results = {}
    
    for test in test_list:
        local_pass = False
        local_error = None
        try:
            test()
            local_pass = True
        except:
            local_error = traceback.format_exc()
        
        results[test.__name__] = (local_pass, local_error)
    
    all_results = comm.gather(results, root=0)
    
    if rank == 0:
        print("\n" + "="*60)
        print("PyMOAB Parallel Test Results")
        print("="*60)
        print(f"Running on {size} processes\n")
        
        for proc_results in all_results:
            for test_name, (passed, error) in proc_results.items():
                if passed:
                    print(f"{colors.OKGREEN}PASS{colors.ENDC}: {test_name}")
                else:
                    print(f"{colors.FAIL}FAIL{colors.ENDC}: {test_name}")
                    if error:
                        print(f"  Error: {error.splitlines()[-1]}")
                    ret_val += 1
        
        print("="*60)
        if ret_val == 0:
            print(f"{colors.OKGREEN}ALL TESTS PASSED{colors.ENDC}")
        else:
            print(f"{colors.FAIL}{ret_val} TEST(S) FAILED{colors.ENDC}")
        print("="*60 + "\n")
    
    comm.Barrier()
    return ret_val

def CHECK(actual_value):
    CHECK_EQ(actual_value, True)

def CHECK_EQ(actual_value, expected_value):
    err_msg = "Expected value: {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    if isinstance(actual_value, np.ndarray) and isinstance(expected_value, np.ndarray):
        result = np.array_equal(actual_value, expected_value)
    else:
        result = actual_value == expected_value
    assert result, err_msg

def CHECK_NOT_EQ(actual_value, expected_value):
    err_msg = "Expected value: not {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    result = expected_value != actual_value
    assert result, err_msg

def CHECK_ITER_EQ(actual_value, expected_value):
    CHECK_EQ(len(actual_value), len(expected_value))
    for a,e in zip(actual_value, expected_value):
        if isinstance(a, str) and (e, str):
            CHECK_EQ(a,e)
            continue
        if hasattr(a, '__iter__') and hasattr(e, '__iter__'):
            CHECK_ITER_EQ(a,e)
        else:
            CHECK_EQ(a,e)

def CHECK_PARALLEL(condition, msg=""):
    """Check that passes on all ranks."""
    local_pass = condition
    all_pass = comm.allreduce(local_pass, op=MPI.LAND)
    if not all_pass:
        raise AssertionError(f"Parallel check failed: {msg}")

def CHECK_PARALLEL_EQ(actual_value, expected_value):
    """Check that value is equal across all ranks."""
    local_pass = (actual_value == expected_value)
    all_pass = comm.allreduce(local_pass, op=MPI.LAND)
    if not all_pass:
        raise AssertionError(f"Value {actual_value} differs across ranks")

def get_all_values(value):
    """Gather a value from all processes to rank 0."""
    return comm.allgather(value)

def get_root_value(value, root=0):
    """Gather value to root process."""
    return comm.bcast(value, root=root)
