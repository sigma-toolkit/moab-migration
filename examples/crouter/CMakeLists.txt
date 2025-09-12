project(test_distributed_reduction)

find_package(MOAB REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

include(${MOAB_DIR}/MOABConfig.cmake)

# Add test program
add_executable(test_distributed_reduction
    test_distributed_reduction.cpp
    determine_communication_pattern.cpp
)

target_link_libraries(test_distributed_reduction
    PRIVATE
    ${MOAB_LIBRARIES}
    Eigen3::Eigen
    ${MPI_CXX_LIBRARIES}
)

add_executable(test_dof_communication
    test_dof_communication.cpp
    determine_communication_pattern.cpp
)

target_link_libraries(test_dof_communication
    PRIVATE
    ${MOAB_LIBRARIES}
    Eigen3::Eigen
    ${MPI_CXX_LIBRARIES}
)


# Add test
add_test(NAME test_distributed_reduction
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
    $<TARGET_FILE:test_distributed_reduction>
)
add_test(NAME test_dof_communication
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
    $<TARGET_FILE:test_dof_communication>
)
