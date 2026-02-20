#!/bin/bash

# Get the directory of this script
#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run the parallel tests with 2 processes
#cd "$SCRIPT_DIR"
mpiexec -n 2 python -m unittest test_parallel_io.py

# Check the exit status
if [ $? -eq 0 ]; then
    echo "Parallel tests passed!"
else
    echo "Parallel tests failed!"
    exit 1
fi
