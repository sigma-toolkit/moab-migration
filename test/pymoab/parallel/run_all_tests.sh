#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

PYTHON_BIN=${PYTHON_BIN:-/usr/bin/python3}
MPIEXEC_BIN=${MPIEXEC_BIN:-/opt/mpich/bin/mpiexec}

echo "Running PyMOAB Parallel Tests"
echo "=============================="

echo ""
echo "--- Test 1: ParallelComm Creation ---"
"$MPIEXEC_BIN" -n 2 "$PYTHON_BIN" test_pcomm_creation.py

echo ""
echo "--- Test 2: Ghost Exchange ---"
"$MPIEXEC_BIN" -n 2 "$PYTHON_BIN" test_ghost_exchange.py

echo ""
echo "--- Test 3: Shared Entities ---"
"$MPIEXEC_BIN" -n 2 "$PYTHON_BIN" test_shared_entities.py

echo ""
echo "--- Test 4: ParallelComm Shutdown ---"
"$MPIEXEC_BIN" -n 2 "$PYTHON_BIN" test_parallel_shutdown.py

echo ""
echo "=============================="
echo "All parallel tests completed!"
