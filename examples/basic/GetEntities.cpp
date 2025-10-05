/** @example GetEntities.cpp
 * @brief Example demonstrating entity querying and connectivity access.
 *
 * This example demonstrates how to:
 * - Get all entities in the mesh database
 * - Access entity connectivity (vertex connectivity for elements)
 * - Query vertex adjacencies (elements connected to vertices)
 * - Iterate through entities and examine their properties
 *
 * The program reads a mesh file and reports connectivity information
 * for all non-vertex entities and adjacency information for vertices.
 *
 * To run: ./GetEntities [meshfile]
 * (default: uses hex01.vtk in the current directory)
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

// Using declarations for cleaner code
using moab::Core;
using moab::EntityHandle;
using moab::ErrorCode;
using moab::Range;

// Constants
namespace {
    const char* const DEFAULT_MESH_FILE = "hex01.vtk";
    const int DIMENSION = 3;  // For adjacency queries
} // namespace

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [meshfile]\n"
              << "  meshfile - Path to the mesh file (default: " << MESH_DIR << "/" << DEFAULT_MESH_FILE << ")\n";
}

int main(int argc, char** argv) {
    try {
        // Parse command line arguments
        if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
            print_usage(argv[0]);
            return 0;
        }

        // Get mesh file path
        const std::string mesh_file = (argc > 1)
            ? argv[1]
            : std::string(MESH_DIR) + "/" + std::string(DEFAULT_MESH_FILE);

        // Initialize MOAB
        auto moab = std::make_unique<Core>();
        if (!moab) {
            std::cerr << "Error: Failed to create MOAB instance\n";
            return 1;
        }

        // Load the mesh
        MB_CHK_SET_ERR(
            moab->load_mesh(mesh_file.c_str()),
            "Failed to load mesh file: " << mesh_file
        );

        // Get all entities in the mesh
        Range entities;
        MB_CHK_SET_ERR(
            moab->get_entities_by_handle(0, entities),
            "Failed to retrieve entities from mesh"
        );

        // Process each entity
        for (const auto& entity : entities) {
            const auto entity_type = moab->type_from_handle(entity);
            const auto entity_id = moab->id_from_handle(entity);

            if (entity_type == moab::MBVERTEX) {
                // For vertices, get adjacent elements
                Range adjacent_elements;
                MB_CHK_SET_ERR(
                    moab->get_adjacencies(&entity, 1, DIMENSION, false, adjacent_elements),
                    "Failed to get adjacencies for vertex " << entity_id
                );

                std::cout << "Vertex " << entity_id << " is connected to "
                          << adjacent_elements.size() << " elements\n";
                if (!adjacent_elements.empty()) {
                    std::cout << "  Adjacent elements: ";
                    for (const auto& elem : adjacent_elements) {
                        std::cout << moab->id_from_handle(elem) << " ";
                    }
                    std::cout << "\n";
                }
            }
            else if (entity_type < moab::MBENTITYSET) {
                // For elements, get connectivity
                const EntityHandle* connectivity = nullptr;
                int num_vertices = 0;

                MB_CHK_SET_ERR(
                    moab->get_connectivity(entity, connectivity, num_vertices),
                    "Failed to get connectivity for " << moab::CN::EntityTypeName(entity_type)
                    << " " << entity_id
                );

                std::cout << moab::CN::EntityTypeName(entity_type) << " " << entity_id
                          << " has " << num_vertices << " vertices: ";

                for (int i = 0; i < num_vertices; ++i) {
                    std::cout << moab->id_from_handle(connectivity[i]);
                    if (i < num_vertices - 1) std::cout << ", ";
                }
                std::cout << "\n";
            }
        }

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "Error: Unknown exception occurred\n";
        return 1;
    }
}
