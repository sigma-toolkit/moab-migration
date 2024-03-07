#!/usr/bin/env python
import os
import sys
from setuptools import find_packages
from skbuild import setup

# Check for Windows
IS_NT = os.name == "nt"

# Check for Mac
IS_MAC = sys.platform == "darwin"

# Cmake args
cmake_args = [
    "-GNinja",
    "-DENABLE_PYMOAB:BOOL=ON",
    "-DPYTHON_EXECUTABLE:FILEPATH=" + sys.executable,
    "-DCMAKE_BUILD_TYPE:STRING=Release",
]

# Specify GCC as the compiler for Windows
if IS_NT:
    cmake_args.append("-GMinGW Makefiles")
    cmake_args.append("-DCMAKE_C_COMPILER:FILEPATH=gcc")
    cmake_args.append("-DCMAKE_CXX_COMPILER:FILEPATH=g++")


# Collect extension
extension = ["*.dll", "*.so", "*.dylib", "*.pyd", "*.pyo"]

# Setup configuration
setup(
    packages=find_packages(),
    package_data={
        "lib": extension,
        "pymoab": [
            "*.pxd",
        ] + extension,
    },
    cmake_args=cmake_args,
    cmake_install_dir=".",
)
