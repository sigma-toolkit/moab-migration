#!/usr/bin/env python3
"""
Simple test runner for PyMOAB parallel interface.
Runs the pcomm creation tests as a quick smoke test.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import test_pcomm_creation
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

if __name__ == "__main__":
    sys.exit(test_pcomm_creation.main())
