/** @example DeleteEdges.cpp
 * This example demonstrates edge deletion from a mesh.
 * It shows how to load a mesh from a file,
 * retrieve all edges (1D entities) from the mesh,
 * delete all edges from the mesh,
 * and write the modified mesh to a new file.
 *
 * The resulting mesh will have vertices and higher-dimensional entities
 * (faces, volumes) but no edges, which can be useful for certain
 * mesh processing workflows.
 *
 * To run: ./DeleteEdges [meshfile] [outfile]
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"
#include <iostream>
#include <memory>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string( MESH_DIR ) + string( "/hex01.vtk" );
string out_file       = string( "outFile.h5m" );

int main( int argc, char** argv )
{
    if( argc == 1 )
    {
        std::cout << "Usage: " << argv[0] << " [meshfile] [outfile]\n";
        return 0;
    }
    std::string mesh_file = ( argc > 1 ) ? argv[1] : std::string( MESH_DIR ) + "/hex01.vtk";
    std::string out_file  = ( argc > 2 ) ? argv[2] : "outFile.h5m";

    auto mb = std::make_unique< Core >();
    if( !mb )
    {
        std::cerr << "Error: Could not allocate MOAB Core instance.\n";
        return 1;
    }
    if( MB_SUCCESS != mb->load_mesh( mesh_file.c_str() ) )
    {
        std::cerr << "Error: Could not load mesh file '" << mesh_file << "'\n";
        return 1;
    }

    Range edges;
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 1, edges ) );
    MB_CHK_ERR( mb->delete_entities( edges ) );

    MB_CHK_SET_ERR( mb->write_file( out_file.c_str() ), "Error: Could not write output mesh file '" + out_file + "'" );
    std::cout << "Deleted " << edges.size() << " edges. Output written to '" << out_file << "'.\n";
    return 0;
}
