#!/bin/bash
# ==============================================================================
# generate_baselines.sh
# Generate golden reference files for mbtempest regression tests
#
# Usage:
#   cd <build_dir>
#   cmake -DMOAB_BUILD_MBTEMPEST_TESTS=ON -DMOAB_MBTEMPEST_GENERATE_BASELINES=ON ..
#   make -j$(nproc)
#   ctest -L mbtempest  # or run this script directly
#
# This script is the manual alternative to the CMake baseline generation mode.
# It runs mbtempest with all test configurations and copies the map outputs
# to the baselines directory.
#
# After generation, upload baselines to the external server and set
# MBTEMPEST_BASELINE_URL in CMake.
# ==============================================================================

set -euo pipefail

# Paths — adjust these for your build
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${1:-$(pwd)}"
MBTEMPEST="${BUILD_DIR}/bin/mbtempest"
WORK_DIR="${BUILD_DIR}/test/mbtempest_regression/work"
BASELINE_DIR="${BUILD_DIR}/test/mbtempest_regression/baselines"

if [ ! -x "${MBTEMPEST}" ]; then
  echo "ERROR: mbtempest not found at ${MBTEMPEST}"
  echo "Usage: $0 [build_dir]"
  exit 1
fi

mkdir -p "${WORK_DIR}" "${BASELINE_DIR}"
cd "${WORK_DIR}"

echo "============================================"
echo "Generating mbtempest regression baselines"
echo "Working directory: ${WORK_DIR}"
echo "Baseline directory: ${BASELINE_DIR}"
echo "============================================"

# Helper: run mbtempest and copy output map as baseline
run_and_save() {
  local name="$1"
  local output="$2"
  shift 2
  echo ""
  echo "--- ${name} ---"
  echo "Command: ${MBTEMPEST} $@"
  if "${MBTEMPEST}" "$@"; then
    if [ -f "${output}" ]; then
      cp "${output}" "${BASELINE_DIR}/gold_${output}"
      echo "  -> Saved baseline: gold_${output}"
    else
      echo "  WARNING: Expected output ${output} not found"
    fi
  else
    echo "  FAILED: mbtempest returned non-zero exit code"
    return 1
  fi
}

# ==============================================================================
# Step 1: Generate meshes
# ==============================================================================

echo ""
echo "=== Generating meshes ==="

"${MBTEMPEST}" -t 0 -r 5  -f cs5.h5m
"${MBTEMPEST}" -t 0 -r 10 -f cs10.h5m
"${MBTEMPEST}" -t 0 -r 25 -f cs25.h5m
"${MBTEMPEST}" -t 1 -r 5  -f rll5.h5m
"${MBTEMPEST}" -t 1 -r 10 -f rll10.h5m
"${MBTEMPEST}" -t 1 -r 25 -f rll25.h5m
"${MBTEMPEST}" -t 2 -r 5  -f ico5.h5m
"${MBTEMPEST}" -t 2 -r 5  -d -f icod5.h5m
"${MBTEMPEST}" -t 2 -r 10 -f ico10.h5m
"${MBTEMPEST}" -t 2 -r 10 -d -f icod10.h5m

echo "Meshes generated."

# ==============================================================================
# Step 2: FV-FV remapping baselines
# ==============================================================================

echo ""
echo "=== FV-FV remapping baselines ==="

run_and_save "CS5→RLL5 FV-FV o1" map_cs5_rll5_fv_fv_o1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_cs5_rll5_fv_fv_o1.nc

run_and_save "CS5→ICO5 FV-FV o1" map_cs5_ico5_fv_fv_o1.nc \
  -t 5 -l cs5.h5m -l ico5.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_cs5_ico5_fv_fv_o1.nc

run_and_save "CS5→ICOD5 FV-FV o1" map_cs5_icod5_fv_fv_o1.nc \
  -t 5 -l cs5.h5m -l icod5.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_cs5_icod5_fv_fv_o1.nc

run_and_save "RLL5→CS5 FV-FV o1" map_rll5_cs5_fv_fv_o1.nc \
  -t 5 -l rll5.h5m -l cs5.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_rll5_cs5_fv_fv_o1.nc

run_and_save "ICO5→CS5 FV-FV o1" map_ico5_cs5_fv_fv_o1.nc \
  -t 5 -l ico5.h5m -l cs5.h5m -w -m fv -m fv -o 1 -o 1 --verify --var SH \
  -f map_ico5_cs5_fv_fv_o1.nc

run_and_save "CS10→RLL10 FV-FV o2" map_cs10_rll10_fv_fv_o2.nc \
  -t 5 -l cs10.h5m -l rll10.h5m -w -m fv -m fv -o 2 -o 2 --verify --var SH \
  -f map_cs10_rll10_fv_fv_o2.nc

run_and_save "CS10→RLL10 FV-FV o3" map_cs10_rll10_fv_fv_o3.nc \
  -t 5 -l cs10.h5m -l rll10.h5m -w -m fv -m fv -o 3 -o 3 --verify --var SH \
  -f map_cs10_rll10_fv_fv_o3.nc

# ==============================================================================
# Step 3: SE-FV baselines
# ==============================================================================

echo ""
echo "=== SE-FV remapping baselines ==="

run_and_save "CS5→RLL5 CGLL-FV o1-o1" map_cs5_rll5_cgll_fv_o1_o1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m fv -o 1 -o 1 --verify --var SH \
  -f map_cs5_rll5_cgll_fv_o1_o1.nc

run_and_save "CS5→RLL5 CGLL-FV o4-o1" map_cs5_rll5_cgll_fv_o4_o1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m fv -o 4 -o 1 --verify --var SV \
  -f map_cs5_rll5_cgll_fv_o4_o1.nc

run_and_save "CS5→RLL5 DGLL-FV o1-o1" map_cs5_rll5_dgll_fv_o1_o1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m dgll -m fv -o 1 -o 1 --verify --var SH \
  -f map_cs5_rll5_dgll_fv_o1_o1.nc

run_and_save "CS5→RLL5 DGLL-FV o4-o1" map_cs5_rll5_dgll_fv_o4_o1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m dgll -m fv -o 4 -o 1 --verify --var SV \
  -f map_cs5_rll5_dgll_fv_o4_o1.nc

# ==============================================================================
# Step 4: FV-SE baselines (serial only)
# ==============================================================================

echo ""
echo "=== FV-SE remapping baselines ==="

run_and_save "CS5→RLL5 FV-CGLL o1-o4" map_cs5_rll5_fv_cgll_o1_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m cgll -o 1 -o 4 --verify --var SV \
  -f map_cs5_rll5_fv_cgll_o1_o4.nc

run_and_save "CS5→RLL5 FV-DGLL o1-o4" map_cs5_rll5_fv_dgll_o1_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m dgll -o 1 -o 4 --verify --var SV \
  -f map_cs5_rll5_fv_dgll_o1_o4.nc

# ==============================================================================
# Step 5: SE-SE baselines (serial only)
# ==============================================================================

echo ""
echo "=== SE-SE remapping baselines ==="

run_and_save "CS5→RLL5 CGLL-CGLL o4-o4" map_cs5_rll5_cgll_cgll_o4_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m cgll -o 4 -o 4 --verify --var SV \
  -f map_cs5_rll5_cgll_cgll_o4_o4.nc

run_and_save "CS5→RLL5 DGLL-DGLL o4-o4" map_cs5_rll5_dgll_dgll_o4_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m dgll -m dgll -o 4 -o 4 --verify --var SV \
  -f map_cs5_rll5_dgll_dgll_o4_o4.nc

run_and_save "CS5→RLL5 CGLL-DGLL o4-o4" map_cs5_rll5_cgll_dgll_o4_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m dgll -o 4 -o 4 --verify --var SV \
  -f map_cs5_rll5_cgll_dgll_o4_o4.nc

run_and_save "CS5→RLL5 DGLL-CGLL o4-o4" map_cs5_rll5_dgll_cgll_o4_o4.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m dgll -m cgll -o 4 -o 4 --verify --var SV \
  -f map_cs5_rll5_dgll_cgll_o4_o4.nc

# ==============================================================================
# Step 6: FV sub-method baselines
# ==============================================================================

echo ""
echo "=== FV sub-method baselines ==="

run_and_save "CS5→RLL5 bilinear" map_cs5_rll5_fv_fv_bilin.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod bilin \
  --verify --var SH -f map_cs5_rll5_fv_fv_bilin.nc

run_and_save "CS5→RLL5 delaunay" map_cs5_rll5_fv_fv_delaunay.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod delaunay \
  --verify --var SH -f map_cs5_rll5_fv_fv_delaunay.nc

run_and_save "CS5→RLL5 intbilin" map_cs5_rll5_fv_fv_intbilin.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod intbilin \
  --verify --var SH -f map_cs5_rll5_fv_fv_intbilin.nc

run_and_save "CS5→RLL5 invdist" map_cs5_rll5_fv_fv_invdist.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --fvmethod invdist \
  --verify --var SH -f map_cs5_rll5_fv_fv_invdist.nc

# ==============================================================================
# Step 7: Monotonicity baselines
# ==============================================================================

echo ""
echo "=== Monotonicity baselines ==="

for mono in 1 2 3; do
  run_and_save "CS5→RLL5 FV-FV mono${mono}" "map_cs5_rll5_fv_fv_mono${mono}.nc" \
    -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --monotonicity ${mono} \
    --verify --var SH -f "map_cs5_rll5_fv_fv_mono${mono}.nc"
done

run_and_save "CS5→RLL5 CGLL-FV mono1" map_cs5_rll5_cgll_fv_mono1.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m fv -o 4 -o 1 --monotonicity 1 \
  --verify --var SV -f map_cs5_rll5_cgll_fv_mono1.nc

# ==============================================================================
# Step 8: CAAS limiter baselines
# ==============================================================================

echo ""
echo "=== CAAS limiter baselines ==="

for caas in 1 2 3 4; do
  run_and_save "CS5→RLL5 FV-FV CAAS${caas}" "map_cs5_rll5_fv_fv_caas${caas}.nc" \
    -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --limiter ${caas} \
    --verify --var SH -f "map_cs5_rll5_fv_fv_caas${caas}.nc"
done

# ==============================================================================
# Step 9: Advanced option baselines
# ==============================================================================

echo ""
echo "=== Advanced option baselines ==="

run_and_save "CS5→RLL5 gnomonic" map_cs5_rll5_gnomonic.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --gnomonic \
  --verify --var SH -f map_cs5_rll5_gnomonic.nc

run_and_save "CS5→RLL5 advfront" map_cs5_rll5_advfront.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --advfront \
  --verify --var SH -f map_cs5_rll5_advfront.nc

run_and_save "CS5→RLL5 noconserve" map_cs5_rll5_noconserve.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --noconserve \
  --verify --var SH -f map_cs5_rll5_noconserve.nc

run_and_save "CS5→RLL5 convexity" map_cs5_rll5_convexity.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --enforce_convexity \
  --verify --var SH -f map_cs5_rll5_convexity.nc

run_and_save "CS5→RLL5 nobubble" map_cs5_rll5_nobubble.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m cgll -m fv -o 4 -o 1 --nobubble \
  --verify --var SV -f map_cs5_rll5_nobubble.nc

run_and_save "CS5→RLL5 sparse" map_cs5_rll5_sparse.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --sparseconstraints \
  --verify --var SH -f map_cs5_rll5_sparse.nc

run_and_save "CS5→RLL5 volumetric" map_cs5_rll5_volumetric.nc \
  -t 5 -l cs5.h5m -l rll5.h5m -w -m fv -m fv -o 1 -o 1 --volumetric \
  --verify --var SH -f map_cs5_rll5_volumetric.nc

# ==============================================================================
# Step 10: Convergence baselines
# ==============================================================================

echo ""
echo "=== Convergence baselines ==="

run_and_save "CS10→RLL10 CGLL-FV o4-o1" map_cs10_rll10_cgll_fv_o4_o1.nc \
  -t 5 -l cs10.h5m -l rll10.h5m -w -m cgll -m fv -o 4 -o 1 --verify --var SV \
  -f map_cs10_rll10_cgll_fv_o4_o1.nc

# ==============================================================================
# Done
# ==============================================================================

echo ""
echo "============================================"
echo "Baseline generation complete."
echo "Files saved to: ${BASELINE_DIR}"
echo ""
echo "Baseline files generated:"
ls -la "${BASELINE_DIR}"/gold_*.nc 2>/dev/null | wc -l
echo " files total"
echo ""
echo "Next steps:"
echo "  1. Upload baselines to external server"
echo "  2. Set MBTEMPEST_BASELINE_URL in CMake"
echo "  3. Rebuild with MOAB_MBTEMPEST_GENERATE_BASELINES=OFF"
echo "============================================"
