# run_convergence_test.cmake
# CMake script to run two mbtempest invocations and verify convergence.
#
# Inputs (set via -D on the cmake command line):
#   MBTEMPEST       - path to mbtempest binary
#   COARSE_ARGS     - semicolon-separated args for the coarse/low-order run
#   FINE_ARGS       - semicolon-separated args for the fine/high-order run
#   VERIFY_SCRIPT   - path to verify_convergence.sh
#   TEST_NAME       - test name for reporting
#   MPI_PREFIX      - mpi launch prefix (may be empty)
#   MPI_SUFFIX      - mpi launch suffix (may be empty)

# Build command lists — args arrive as semicolon-separated CMake lists
set(COARSE_ARGS_LIST ${COARSE_ARGS})
set(FINE_ARGS_LIST ${FINE_ARGS})
separate_arguments(MPI_PREFIX_LIST UNIX_COMMAND "${MPI_PREFIX}")
separate_arguments(MPI_SUFFIX_LIST UNIX_COMMAND "${MPI_SUFFIX}")

# Run coarse case
message(STATUS "${TEST_NAME}: Running coarse case...")
execute_process(
  COMMAND ${MPI_PREFIX_LIST} ${MBTEMPEST} ${COARSE_ARGS_LIST} ${MPI_SUFFIX_LIST}
  OUTPUT_VARIABLE COARSE_OUTPUT
  ERROR_VARIABLE COARSE_ERROR
  RESULT_VARIABLE COARSE_RESULT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT COARSE_RESULT EQUAL 0)
  message(FATAL_ERROR "${TEST_NAME}: Coarse run failed with exit code ${COARSE_RESULT}\n${COARSE_ERROR}")
endif()

# Run fine case
message(STATUS "${TEST_NAME}: Running fine case...")
execute_process(
  COMMAND ${MPI_PREFIX_LIST} ${MBTEMPEST} ${FINE_ARGS_LIST} ${MPI_SUFFIX_LIST}
  OUTPUT_VARIABLE FINE_OUTPUT
  ERROR_VARIABLE FINE_ERROR
  RESULT_VARIABLE FINE_RESULT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT FINE_RESULT EQUAL 0)
  message(FATAL_ERROR "${TEST_NAME}: Fine run failed with exit code ${FINE_RESULT}\n${FINE_ERROR}")
endif()

# Extract L_inf errors from output
# Combined stdout+stderr since mbtempest prints to both
set(ALL_COARSE "${COARSE_OUTPUT}\n${COARSE_ERROR}")
set(ALL_FINE "${FINE_OUTPUT}\n${FINE_ERROR}")

string(REGEX MATCH "L_inf error *= *([0-9.e+\\-]+)" _match_coarse "${ALL_COARSE}")
set(ERR_COARSE "${CMAKE_MATCH_1}")
string(REGEX MATCH "L_inf error *= *([0-9.e+\\-]+)" _match_fine "${ALL_FINE}")
set(ERR_FINE "${CMAKE_MATCH_1}")

if(NOT ERR_COARSE OR NOT ERR_FINE)
  message(FATAL_ERROR "${TEST_NAME}: Could not extract L_inf errors.\n  Coarse output: ${ALL_COARSE}\n  Fine output: ${ALL_FINE}")
endif()

message(STATUS "${TEST_NAME}: L_inf coarse=${ERR_COARSE}, fine=${ERR_FINE}")

# Compare using the verify script
execute_process(
  COMMAND bash ${VERIFY_SCRIPT} --values ${ERR_COARSE} ${ERR_FINE} ${TEST_NAME}
  RESULT_VARIABLE VERIFY_RESULT
  OUTPUT_VARIABLE VERIFY_OUTPUT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "${VERIFY_OUTPUT}")
if(NOT VERIFY_RESULT EQUAL 0)
  message(FATAL_ERROR "${TEST_NAME}: Convergence check FAILED")
endif()
