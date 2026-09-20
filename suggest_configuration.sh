#!/usr/bin/env bash
# DEPRECATED: machine knowledge has moved to install/install-moab.sh.
#
# This file used to print a hardcoded ./configure command keyed by hostname,
# with a database that decayed (vesta, mira, blogin, theta, cori, edison were
# all retired but still listed). The replacement integrates with the unified
# install-moab.sh's machine database and standalone profile:
#
#   install/install-moab.sh --print --profile=standalone [other flags]
#
# This shim forwards everything you pass to it through that pipeline so
# legacy invocations continue to print a suggested configure command,
# but now sourced from a maintained registry instead of stale paths.
#
# Examples:
#   ./suggest_configuration.sh                            # uses auto-detected machine
#   ./suggest_configuration.sh --machine=bebop            # for a specific machine
#   ./suggest_configuration.sh --build-system=cmake       # cmake variant
#
# Note: --profile=standalone deliberately does NOT touch your environment;
# it trusts whatever modules / env vars you have loaded. For --profile=e3sm
# (which auto-loads E3SM's modules + env), use install-moab.sh directly.

exec "$(dirname "${BASH_SOURCE[0]}")/install/install-moab.sh" \
    --print --profile=standalone "$@"
