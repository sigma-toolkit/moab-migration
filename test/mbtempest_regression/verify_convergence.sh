#!/bin/bash
# ==============================================================================
# verify_convergence.sh
# Verify that remapping error decreases between two mbtempest runs.
#
# Usage:
#   verify_convergence.sh <map_coarse.nc> <map_fine.nc> <workdir>
#
# The script reads the L_inf error from stdout captured during the mbtempest
# runs (stored in .err files alongside the map files) and verifies:
#   error_fine < error_coarse
#
# Alternatively, if two error values are passed directly:
#   verify_convergence.sh --values <error_coarse> <error_fine> [<label>]
# ==============================================================================

set -euo pipefail

if [ "${1:-}" = "--values" ]; then
    # Direct value comparison mode
    err_coarse="$2"
    err_fine="$3"
    label="${4:-convergence}"

    # Use awk for floating-point comparison
    result=$(awk "BEGIN { print ($err_fine < $err_coarse) ? \"PASS\" : \"FAIL\" }")
    ratio=$(awk "BEGIN { printf \"%.4f\", $err_fine / $err_coarse }")

    if [ "$result" = "PASS" ]; then
        echo "PASS: ${label}: error decreased from ${err_coarse} to ${err_fine} (ratio=${ratio})"
        exit 0
    else
        echo "FAIL: ${label}: error did NOT decrease: coarse=${err_coarse}, fine=${err_fine} (ratio=${ratio})"
        exit 1
    fi
fi

echo "Usage: $0 --values <error_coarse> <error_fine> [label]"
exit 1
