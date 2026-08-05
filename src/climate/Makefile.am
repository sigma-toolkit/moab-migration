# Climate/E3SM-specific MOAB sources.
#
# This subtree isolates functionality that exists to support climate
# modeling workflows (E3SM component-coupler mesh migration, conservative
# remapping on the sphere, and the climate NetCDF grid readers/writers) from
# the general-purpose MOAB core and its domain-agnostic features.  Keeping it
# separate makes the climate stack easier to refactor and maintain without
# perturbing the core mesh database.
#
# Each component builds its own convenience library, which src/Makefile.am
# links directly into libMOAB.

AUTOMAKE_OPTIONS = foreign subdir-objects

SUBDIRS = io

EXTRA_DIST = CMakeLists.txt
