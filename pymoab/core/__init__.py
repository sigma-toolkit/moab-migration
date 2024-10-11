"""
This module provides paths needed for the core MOAB interface.
It includes:
- Include directory path for MOAB headers
- Library directory paths for MOAB libraries
- List of all the shared libraries of MOAB
"""

import os
import glob
import sysconfig
from pymoab import __path__ as pymoab_path

def get_include_path():
    """Return the include directory path for MOAB headers."""
    return [os.path.join(pymoab_path[0], "core", "include")]

def get_library_path():
    """Return the library directory path, considering both 'lib' and 'lib64'."""
    lib_dir_candidates = ["lib", "lib64"]
    lib_paths = []
    
    for candidate in lib_dir_candidates:
        lib_path = os.path.join(pymoab_path[0], "core", candidate)
        if os.path.exists(lib_path):
            lib_paths.append(lib_path)

    return lib_paths

def get_shared_libraries():
    """List all the shared libraries of MOAB."""
    libs_path = os.path.join(sysconfig.get_paths()['platlib'], "moab.libs")
    if os.path.exists(libs_path):
        shared_libs = glob.glob(os.path.join(libs_path, "*"))
        return shared_libs
    return []

# Export variables for easy access
include = get_include_path()
lib = get_library_path()
shared_libs = get_shared_libraries()

# Export variables for easy access
__all__ = ["include", "lib", "shared_libs"]
