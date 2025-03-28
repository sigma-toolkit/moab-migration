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

try:
    MOAB_CORE_BASE_PATH = os.path.join(__path__[0], "core")
except NameError:
    MOAB_CORE_BASE_PATH = None

if not MOAB_CORE_BASE_PATH or not os.path.exists(MOAB_CORE_BASE_PATH):
    import sysconfig
    MOAB_CORE_BASE_PATH = os.path.join(sysconfig.get_path("platlib"), "pymoab", "core")
    if not os.path.exists(MOAB_CORE_BASE_PATH):
        raise ImportError("MOAB is not installed. Please run 'pip install MOAB'.")
    warnings.warn(
        "It seems that PyMOAB is being run from its source directory. "
        "This setup is not recommended as it may lead to unexpected behavior, "
        "such as conflicts between source and installed versions. "
        "Please run your script from outside the MOAB source tree.",
        RuntimeWarning
    )

def get_core_path(subdir, pattern="*", recursive=False):
    """
    Helper function to return paths that match a given pattern within a subdirectory.

    Args:
        subdir (str): The subdirectory within the 'core' directory.
        pattern (str): The pattern to match files or directories.
        recursive (bool): Whether to search recursively in subdirectories.

    Returns:
        list: A list of matched paths.
    """
    search_pattern = os.path.join(MOAB_CORE_BASE_PATH, "**", pattern) if recursive else os.path.join(path, pattern)
    return glob.glob(search_pattern, recursive=recursive) if os.path.exists(path) else []

def get_include_path():
    """Return includes and include path for MOAB headers."""
    include = get_core_path("include", "*", recursive=True)
    include_path = get_core_path("include", "", recursive=False)
    return include, include_path

def get_core_libraries():
    """Return libraries and library paths for MOAB."""
    lib = [lib_file for lib in ["lib", "lib64"] for lib_file in get_core_path(lib, "libMOAB*", recursive=True)]
    lib_path = [lib_file for lib in ["lib", "lib64"] for lib_file in get_core_path(lib, "", recursive=False)]
    return lib, lib_path

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
