# CDash submission target for the MOAB project.
#
# Read by `include(CTest)` (which is invoked from the top-level CMakeLists.txt)
# and by the `ctest` command-line driver when submitting dashboards.
#
# Submission flow:
#   ctest -D Experimental         # one-shot configure/build/test/submit
#   ctest -D Nightly              # for scheduled runs (uses NIGHTLY_START_TIME)
#   ctest -D Continuous           # for CI loops
#
# Per-machine identity (CTEST_SITE / CTEST_BUILD_NAME) is auto-derived in
# config/CTestBuildName.cmake from the resolved feature flags.
set(CTEST_PROJECT_NAME       "MOAB")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD        "https")
set(CTEST_DROP_SITE          "my.cdash.org")
set(CTEST_DROP_LOCATION      "/submit.php?project=MOAB")
set(CTEST_DROP_SITE_CDASH    TRUE)
