#!/usr/bin/env python
import sys
from setuptools import find_packages
from skbuild import setup

# Cmake args
cmake_args = [
    "-DPYTHON_EXECUTABLE:FILEPATH=" + sys.executable,
    "-DCMAKE_BUILD_TYPE:STRING=Release",
]

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
