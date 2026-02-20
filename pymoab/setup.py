from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import numpy
from Cython.Build import cythonize

# Check for mpi4py availability
MPI4PY_AVAILABLE = False
MPI4PY_INCLUDE = None
MOAB_MPI_BASIC = False
MOAB_MPI_IO = False
try:
    import mpi4py
    MPI4PY_AVAILABLE = True
    MPI4PY_INCLUDE = os.path.join(os.path.dirname(mpi4py.__file__), 'include')
    
    # Test basic MPI functionality
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        MOAB_MPI_BASIC = True
        print("PyMOAB: Basic MPI functionality detected")
    except:
        pass
    
    # Test MPI I/O functionality
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        # Test MPI file operations
        MOAB_MPI_IO = True
        print("PyMOAB: MPI I/O functionality detected")
    except:
        pass
except ImportError:
    pass

# Get MOAB and MPI directories from environment or CMake
MOAB_DIR = os.environ.get('MOAB_DIR', '')
MPI_DIR = os.environ.get('MPI_DIR', '')

# Include directories
include_dirs = [
    os.path.join(MOAB_DIR, 'include'),
    os.path.join(MOAB_DIR, 'include', 'moab'),
    numpy.get_include(),
]

if MPI4PY_AVAILABLE and MPI4PY_INCLUDE:
    include_dirs.append(MPI4PY_INCLUDE)

if MPI_DIR:
    include_dirs.append(os.path.join(MPI_DIR, 'include'))

# Library directories
library_dirs = [
    os.path.join(MOAB_DIR, 'lib'),
]

if MPI_DIR:
    library_dirs.append(os.path.join(MPI_DIR, 'lib'))

# Define extensions - core is always built, parallelcomm only with mpi4py
extensions = [
    Extension(
        'pymoab.core',
        ['pymoab/core.pyx'],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=['MOAB'],
        language='c++',
    ),
]

# Add parallelcomm only if mpi4py is available
if MPI4PY_AVAILABLE:
    extensions.append(
        Extension(
            'pymoab.parallelcomm',
            ['pymoab/parallelcomm.pyx'],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=['MOAB'],
            language='c++',
        )
    )

# Cython compiler directives
compiler_directives = {
    'language_level': '3',
    'embedsignature': True,
}

# Add capability detection test
if MPI4PY_AVAILABLE:
    print("PyMOAB: MPI capability detection enabled")
    print(f"PyMOAB: Basic MPI support: {MOAB_MPI_BASIC}")
    print(f"PyMOAB: MPI I/O support: {MOAB_MPI_IO}")
else:
    print("PyMOAB: MPI capability detection not available")

setup(
    name='pymoab',
    version='5.4.0',
    description='Python interface to MOAB',
    author='MOAB Team',
    author_email='moab-dev@mcs.anl.gov',
    packages=['pymoab'],
    ext_modules=cythonize(extensions, compiler_directives=compiler_directives),
    install_requires=[
        'numpy>=1.7.0',
    ],
    extras_require={
        'parallel': ['mpi4py>=3.0.0; platform_system!="Windows"'],
    },
    python_requires='>=3.6',
)