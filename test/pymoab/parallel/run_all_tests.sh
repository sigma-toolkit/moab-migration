#!/bin/bash
# Run all parallel tests

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Running PyMOAB Parallel Tests"
echo "=============================="

echo ""
echo "--- Test 1: ParallelComm Creation ---"
mpiexec -n 2 python test_pcomm_creation.py

echo ""
echo "--- Test 2: Ghost Exchange ---"
mpiexec -n 2 python test_ghost_exchange.py

echo ""
echo "--- Test 3: Shared Entities ---"
mpiexec -n 2 python test_shared_entities.py

echo ""
echo "=============================="
echo "All parallel tests completed!"
