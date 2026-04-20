#!/bin/bash
# run_convergence_test.sh
# Run two mbtempest commands and verify that L_inf error decreases.
#
# Usage:
#   run_convergence_test.sh <test_name> <mbtempest> <mpi_prefix> -- <coarse_args...> -- <fine_args...>
#
# The script runs mbtempest twice, extracts L_inf error from each,
# and exits 0 if fine error < coarse error, 1 otherwise.

set -uo pipefail

TEST_NAME="$1"; shift
MBTEMPEST="$1"; shift
MPI_PREFIX="$1"; shift

# Parse coarse and fine args separated by --
COARSE_ARGS=()
FINE_ARGS=()
phase="skip_first"
for arg in "$@"; do
    if [ "$arg" = "--" ]; then
        if [ "$phase" = "skip_first" ]; then
            phase="coarse"
        elif [ "$phase" = "coarse" ]; then
            phase="fine"
        fi
    else
        if [ "$phase" = "coarse" ]; then
            COARSE_ARGS+=("$arg")
        elif [ "$phase" = "fine" ]; then
            FINE_ARGS+=("$arg")
        fi
    fi
done

extract_linf() {
    echo "$1" | grep "L_inf error" | tail -1 | sed 's/.*= *//'
}

# Run coarse case
echo "${TEST_NAME}: Running coarse case..."
if [ -n "$MPI_PREFIX" ]; then
    COARSE_OUT=$(${MPI_PREFIX} "${MBTEMPEST}" "${COARSE_ARGS[@]}" 2>&1)
else
    COARSE_OUT=$("${MBTEMPEST}" "${COARSE_ARGS[@]}" 2>&1)
fi
COARSE_RC=$?
if [ $COARSE_RC -ne 0 ]; then
    echo "FAIL: coarse run exited with code $COARSE_RC"
    echo "$COARSE_OUT" | tail -5
    exit 1
fi

# Run fine case
echo "${TEST_NAME}: Running fine case..."
if [ -n "$MPI_PREFIX" ]; then
    FINE_OUT=$(${MPI_PREFIX} "${MBTEMPEST}" "${FINE_ARGS[@]}" 2>&1)
else
    FINE_OUT=$("${MBTEMPEST}" "${FINE_ARGS[@]}" 2>&1)
fi
FINE_RC=$?
if [ $FINE_RC -ne 0 ]; then
    echo "FAIL: fine run exited with code $FINE_RC"
    echo "$FINE_OUT" | tail -5
    exit 1
fi

# Extract L_inf errors
ERR_COARSE=$(extract_linf "$COARSE_OUT")
ERR_FINE=$(extract_linf "$FINE_OUT")

if [ -z "$ERR_COARSE" ] || [ -z "$ERR_FINE" ]; then
    echo "FAIL: could not extract L_inf errors"
    echo "  Coarse output: $(echo "$COARSE_OUT" | grep L_inf)"
    echo "  Fine output: $(echo "$FINE_OUT" | grep L_inf)"
    exit 1
fi

# Compare
RESULT=$(awk "BEGIN { print ($ERR_FINE < $ERR_COARSE) ? \"PASS\" : \"FAIL\" }")
RATIO=$(awk "BEGIN { printf \"%.4f\", $ERR_FINE / $ERR_COARSE }")

echo "${TEST_NAME}: L_inf coarse=${ERR_COARSE}, fine=${ERR_FINE}, ratio=${RATIO}"

if [ "$RESULT" = "PASS" ]; then
    echo "PASS: error decreased (ratio=${RATIO})"
    exit 0
else
    echo "FAIL: error did NOT decrease"
    exit 1
fi
