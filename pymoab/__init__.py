""" 
Name

PyMOAB: A Python interface to Argonne National Lab's Mesh Oriented dAtaBase (MOAB)
====

Description
-----------

PyMOAB provides a means of interactively interrogating, modifying, and generating 
MOAB mesh files.

Much of the core functionality of MOAB has been implemented in the core.Core
module. Interaction with this is intended to be largely analogous to interaction
with MOAB via its native C++ API though some modifications have been made to
allow for interaction with the interface using the various native Python
constructs such as lists, tuples, etc.
"""

import os
import sys
import glob
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version('MOAB')
except PackageNotFoundError:
    __version__ = "unknown"

def get_path(subdir, pattern="*"):
    """Helper function to return paths that match a given pattern within a subdirectory."""
    path = os.path.join(__path__[0], "core", subdir)
    return glob.glob(os.path.join(path, pattern)) if os.path.exists(path) else []

def get_include_path():
    """Return the include directory path for MOAB headers."""
    return os.path.join(__path__[0], "core", "include")

def get_core_libraries():
    """Return library paths and library directory paths."""
    lib_paths = [os.path.join(__path__[0], "core", subdir) for subdir in ["lib", "lib64"]]
    libs = [lib for subdir in lib_paths for lib in get_path(subdir)]
    return libs, [path for path in lib_paths if os.path.exists(path)][0]

def get_extra_libraries():
    """List all the extra libraries of MOAB."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.join(__path__[0], "..", "moab.libs")
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]
