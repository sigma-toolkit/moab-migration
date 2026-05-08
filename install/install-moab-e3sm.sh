#!/usr/bin/env bash
# Backwards-compat shim. The unified entry point is install-moab.sh.
# This shim exists so existing E3SM workflows / docs that reference
# install-moab-e3sm.sh by name continue to work. New users should call
# install-moab.sh directly.
exec "$(dirname "${BASH_SOURCE[0]}")/install-moab.sh" "$@"
