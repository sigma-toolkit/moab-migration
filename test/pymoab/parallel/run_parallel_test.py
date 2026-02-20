#!/usr/bin/env python3
"""
Simple test runner for PyMOAB parallel interface.
"""

import sys
import os

# Add the test directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import test_parallel_io
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

if __name__ == "__main__":
    exit_code = test_parallel_io.main()
    sys.exit(exit_code)