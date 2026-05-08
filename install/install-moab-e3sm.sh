#!/usr/bin/env bash
# Backwards-compat shim for the original install-moab-e3sm.sh name. The
# unified entry point is install-moab.sh; this shim preserves the old
# contract by selecting --profile=e3sm (the unified script's default
# profile is now 'standalone' to suit MOAB downstream users; only callers
# of the legacy *-e3sm.sh name get the auto-loaded E3SM env behavior).
#
# --profile=e3sm is prepended so any user-supplied --profile= override
# (which appears later in $@) wins via last-flag-wins parsing.
exec "$(dirname "${BASH_SOURCE[0]}")/install-moab.sh" --profile=e3sm "$@"
