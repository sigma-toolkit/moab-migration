#!/usr/bin/env bash
# DEPRECATED: machine-specific invocations now live in the machine database
# inside install-moab.sh. Use one of:
#
#   ./install-moab.sh --machine=auto         # auto-detect via hostname / NERSC_HOST
#   ./install-moab.sh --machine=bebop        # force a specific entry
#   ./install-moab.sh --list-machines        # show all registered entries
#
# Push 2 will additionally support --profile=e3sm to source modules + env
# vars from $E3SM_ROOT/cime_config/machines/config_machines.xml so that
# machine-specific TPL paths come from E3SM's source of truth rather than
# being hardcoded here.
#
# This shim forwards any extra args to install-moab.sh with --machine=auto
# so that legacy callers keep working.
exec "$(dirname "${BASH_SOURCE[0]}")/install-moab.sh" --machine=auto "$@"
