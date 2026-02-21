from setuptools import setup, Extension
import sys
import os
import numpy
from Cython.Build import cythonize

# Check for mpi4py availability
MPI4PY_AVAILABLE = False
MPI4PY_INCLUDE = None
try:
    import mpi4py
    MPI4PY_AVAILABLE = True
    MPI4PY_INCLUDE = os.path.join(os.path.dirname(mpi4py.__file__), 'include')
except ImportError:
    pass

# Get MOAB and MPI directories from environment
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

# Common extension kwargs
common_ext_kwargs = dict(
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    libraries=['MOAB'],
    language='c++',
)

# All serial modules
serial_modules = ['core', 'rng', 'scd', 'tag', 'hcoord', 'skinner', 'topo_util', 'types']

extensions = [
    Extension(f'pymoab.{mod}', [f'pymoab/{mod}.pyx'], **common_ext_kwargs)
    for mod in serial_modules
]

# Add parallelcomm only if mpi4py is available
if MPI4PY_AVAILABLE:
    extensions.append(
        Extension(
            'pymoab.parallelcomm',
            ['pymoab/parallelcomm.pyx'],
            **common_ext_kwargs,
        )
    )
    print(f"PyMOAB setup.py: mpi4py found, parallelcomm module will be built")
else:
    print(f"PyMOAB setup.py: mpi4py not found, parallelcomm module will NOT be built")

# Cython compiler directives
compiler_directives = {
    'language_level': '3',
    'embedsignature': True,
}

setup(
    name='pymoab',
    version='5.6.0',
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
    python_requires='>=3.8',
)
