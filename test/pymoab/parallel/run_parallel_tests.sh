#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

PYTHON_BIN=${PYTHON_BIN:-python3}
MPIEXEC_BIN=${MPIEXEC_BIN:-mpiexec}

"$MPIEXEC_BIN" -n 2 "$PYTHON_BIN" test_pcomm_creation.py

if [ $? -eq 0 ]; then
    echo "Parallel tests passed!"
else
    echo "Parallel tests failed!"
    exit 1
fi
