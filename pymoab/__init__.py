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

def get_core_path(subdir, pattern="*"):
    """Helper function to return paths that match a given pattern within a subdirectory."""
    path = os.path.join(__path__[0], "core", subdir)
    return glob.glob(os.path.join(path, pattern)) if os.path.exists(path) else []

def get_include_path():
    """Return includes and include path for MOAB headers."""
    include = get_core_path("include")
    include_path = get_core_path("include", "")
    return include, include_path

def get_core_libraries():
    """Return libraries and library paths for MOAB."""
    lib = [libs for lib in ["lib", "lib64"] for libs in get_core_path(lib)]
    lib_paths = [libs for lib in ["lib", "lib64"] for libs in get_core_path(lib, "")]
    return lib, lib_paths

def get_extra_libraries():
    """Return the extra libraries installed by auditwheel or delocate."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.normpath(os.path.join(__path__[0], "..", "moab.libs"))
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include, include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include", "include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]
