/**
 * @file StructuredMeshSimple.cpp
 * @brief Example demonstrating creation and query of structured meshes in MOAB
 *
 * This example shows how to:
 * - Create structured meshes in 1D, 2D, or 3D
 * - Work with structured mesh interfaces (ScdInterface)
 * - Handle both serial and parallel structured mesh creation
 * - Query mesh entities and their connectivity
 * - Access element coordinates and connectivity
 * - Use parametric indexing for structured meshes
 *
 * In serial mode, a single N*N*N block of elements is created.
 * In parallel mode, each processor gets an N*N*N block arranged
 * in a 1D column, sharing vertices and faces at interfaces.
 *
 * @author MOAB Development Team
 * @date 2024
 *

 * \brief Show creation and query of structured mesh, serial or parallel, through MOAB's structured
 * mesh interface. This is an example showing creation and query of a 3D structured mesh.  In
 * serial, a single N*N*N block of elements is created; in parallel, each proc gets an N*N*N block,
 * with blocks arranged in a 1d column, sharing vertices and faces at their interfaces (proc 0 has
 * no left neighbor and proc P-1 no right neighbor). Each square block of hex elements is then
 * referenced by its ijk parameterization. 1D and 2D examples could be made simply by changing the
 * dimension parameter passed into the MOAB functions. \n
 *
 * <b>This example </b>:
 *    -# Instantiate MOAB and get the structured mesh interface
 *    -# Decide what the local parameters of the mesh will be, based on parallel/serial and rank.
 *    -# Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
 *    -# Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d,
 * resp.
 *    -# Loop over elements in d nested loops over i, j, k; for each (i,j,k):
 *      -# Get the element corresponding to (i,j,k)
 *      -# Get the connectivity of the element
 *      -# Get the coordinates of the vertices comprising that element
 *    -# Release the structured mesh interface and destroy the MOAB instance
 *
 * <b> To run: </b> ./StructuredMeshSimple [d [N] ] \n
 * (default values so can run w/ no user interaction)
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */

#include "moab/Core.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CN.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include <mpi.h>
#endif

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

// Using declarations for cleaner code
using moab::Core;
using moab::EntityHandle;
using moab::ErrorCode;
using moab::Range;
using moab::ScdBox;
using moab::ScdInterface;

// Constants
namespace {
    constexpr int DEFAULT_DIMENSION = 3;
    constexpr int DEFAULT_ELEMENTS_PER_SIDE = 10;
} // namespace

int main(int argc, char** argv) {
    // Initialize MPI if available
#ifdef MOAB_HAVE_MPI
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed" << std::endl;
        return 1;
    }

    // Ensure MPI is properly finalized when we exit
    struct MPIFinalizer {
        void operator()(int*) const { MPI_Finalize(); }
    };
    std::unique_ptr<int, MPIFinalizer> mpi_guard(nullptr);
#endif

    // Parse command line options
    int dimension = DEFAULT_DIMENSION;
    int elements_per_side = DEFAULT_ELEMENTS_PER_SIDE;

    try {
        ProgOptions opts;
        opts.addOpt<int>("dim,d", "Dimension of mesh (default=3)", &dimension);
        opts.addOpt<int>(",n", "Number of elements on a side (default=10)", &elements_per_side);
        opts.parseCommandLine(argc, argv);

        // Validate input
        if (dimension < 1 || dimension > 3) {
            throw std::invalid_argument("Dimension must be 1, 2, or 3");
        }
        if (elements_per_side < 1) {
            throw std::invalid_argument("Number of elements must be positive");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing command line: " << e.what() << std::endl;
        return 1;
    }

    // Initialize MOAB
    std::unique_ptr<Core> moab_instance = std::make_unique<Core>();
    if (!moab_instance) {
        std::cerr << "Failed to create MOAB instance" << std::endl;
        return 1;
    }

    // Get the structured mesh interface
    ScdInterface* scd_interface = nullptr;
    MB_CHK_SET_ERR(
        moab_instance->query_interface(scd_interface),
        "Failed to get ScdInterface"
    );
    if (!scd_interface) {
        std::cerr << "Failed to get ScdInterface: null pointer returned" << std::endl;
        return 1;
    }

    // 1. Determine local parameters based on parallel/serial and rank
    int rank = 0;
    int ilow = 0;
    int ihigh = elements_per_side;

#ifdef MOAB_HAVE_MPI
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Calculate local range for this process
    const int elements_per_proc = elements_per_side / nprocs;
    const int remainder = elements_per_side % nprocs;

    // Distribute elements among processes
    ilow = rank * elements_per_proc + std::min(rank, remainder);
    ihigh = ilow + elements_per_proc + (rank < remainder ? 1 : 0);
#endif

    std::string write_options = "";
    if (nprocs > 1) write_options = "PARALLEL=WRITE_PART";

    // 2. Create a N^d structured mesh, which includes (N+1)^d vertices and N^d elements.
    ScdBox* box = nullptr;

    // Calculate box bounds based on dimension
    const moab::HomCoord low(
        ilow,
        (dimension > 1) ? 0 : -1,
        (dimension > 2) ? 0 : -1
    );

    const moab::HomCoord high(
        ihigh,
        (dimension > 1) ? elements_per_side : -1,
        (dimension > 2) ? elements_per_side : -1
    );

    // Create the structured mesh box
    MB_CHK_SET_ERR(
        scd_interface->construct_box(
            low, high,
            nullptr,    // No coordinates array
            0,          // No coordinates
            box,        // Output parameter
            nullptr,    // Periodicity
            nullptr,    // Parallel data
            true,       // No coordinates
            0 // Resolve dimensionality
        ),
        "Failed to construct structured box"
    );

    if (!box) {
        std::cerr << "Failed to construct structured box: null box returned" << std::endl;
        return 1;
    }

    // 3. Get the vertices and elements from moab and check their numbers against (N+1)^d and N^d
    Range vertices, elements;

    // Get all vertices (dimension = 0) and elements (dimension = mesh dimension)
    MB_CHK_SET_ERR(moab_instance->get_entities_by_dimension(0, 0, vertices), "Failed to get vertices");
    MB_CHK_SET_ERR(moab_instance->get_entities_by_dimension(0, dimension, elements), "Failed to get elements");

#ifdef MOAB_HAVE_MPI
    std::size_t local_entities_size[2] = {vertices.size(), elements.size()}, global_entities_size[2] = {0, 0};
    MPI_Allreduce(local_entities_size, global_entities_size, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

    // Calculate expected number of vertices and elements
    const std::size_t expected_vertices = static_cast<std::size_t>(std::pow(elements_per_side + 1, dimension));
    const std::size_t expected_elements = static_cast<std::size_t>(std::pow(elements_per_side, dimension));

    // Only print from rank 0 in parallel
    auto print_message = [rank](const auto& message) {
        if (rank == 0) {
            std::cout << message << std::endl;
        }
    };

    // Verify the number of created elements and vertices
    if (expected_elements == global_entities_size[1]) {
        const auto element_type = moab_instance->type_from_handle(*elements.begin());
        std::ostringstream msg;
        msg << "Created " << global_entities_size[1] << " "
            << moab::CN::EntityTypeName(element_type) << " elements and "
            << global_entities_size[0] << " vertices.";
        print_message(msg.str());
    } else {
        std::ostringstream err_msg;
        err_msg << "Error: Expected " << expected_elements << " elements and "
                << expected_vertices << " vertices, but got "
                << global_entities_size[1] << " elements and "
                << global_entities_size[0] << " vertices.";
        print_message(err_msg.str());
        return 1;
    }

    // 4. Loop over elements in nested loops over i, j, k based on dimension
    // const int i_max = elements_per_side - 1;  // 0-based indexing
    const int j_max = (dimension > 1) ? elements_per_side - 1 : 0;
    const int k_max = (dimension > 2) ? elements_per_side - 1 : 0;

    // Pre-allocate vectors to avoid reallocation

    std::vector<double> coordinates;
    for (int k = 0; k <= k_max; ++k) {
        for (int j = 0; j <= j_max; ++j) {
            for (int i = ilow; i < ihigh; ++i) {
                // 4a. Get the element corresponding to (i,j,k)
                const EntityHandle element = box->get_element(i, j, k);
                if (0 == element) {
                    std::cerr << "Failed to get element at (" << i << ", " << j << ", " << k << ")" << std::endl;
                    // return 1;
                    continue;
                }

                // 4b. Get the connectivity of the element
                std::vector< EntityHandle > connectivity;
                MB_CHK_SET_ERR(
                    moab_instance->get_connectivity(&element, 1, connectivity),
                    "Failed to get connectivity for element at (" << i << ", " << j << ", " << k << ") with handle " << element
                );

                // 4c. Get the coordinates of the vertices comprising that element
                coordinates.resize(3 * connectivity.size());
                MB_CHK_SET_ERR(
                    moab_instance->get_coords(connectivity.data(), connectivity.size(), coordinates.data()),
                    "Failed to get coordinates for element at (" << i << ", " << j << ", " << k << ") with handle " << element
                );
            }
        }
    }

    // 5. Write the mesh file to disk
    MB_CHK_SET_ERR(moab_instance->write_file("structured_mesh.h5m", "h5m", write_options.c_str() ), "Failed to write mesh file");

    // 6. Clean up
    MB_CHK_SET_ERR(moab_instance->release_interface(scd_interface), "Failed to release interface");

    // MPI_Finalize is handled by the mpi_guard destructor if MPI is enabled
    return 0;
}
