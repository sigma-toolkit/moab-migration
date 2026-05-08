#!/bin/bash
# DEPRECATED: this script previously printed a suggested MOAB ./configure
# command based on a hardcoded hostname-keyed database. The database had
# decayed (vesta, mira, blogin, theta, cori, edison are all retired), and
# the duplicated machine knowledge between this script and the install
# scripts in install/ was a maintenance hazard.
#
# Replacement:
#
#   install/install-moab.sh --list-machines        # show registered machines
#   install/install-moab.sh --machine=NAME         # build for that machine
#
# Push 2 will add a --print-only mode that emits just the resolved configure
# command (for users who want to inspect or edit before running). When that
# lands, this shim will be updated to delegate via:
#
#   exec install/install-moab.sh --print --profile=standalone "$@"
#
# Until then, this file remains as a placeholder so that EXTRA_DIST in
# Makefile.am keeps working.
echo "suggest_configuration.sh is deprecated." >&2
echo "Use: install/install-moab.sh --list-machines  (or --machine=NAME, see install/INSTALL-MOAB.md)" >&2
exit 1
