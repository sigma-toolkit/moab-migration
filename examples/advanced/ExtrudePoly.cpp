/**
 * @example ExtrudePoly.cpp
 * Example demonstrating extrusion of 2D polygons to create 3D polyhedra
 *
 * This example shows how to:
 * - Load a 2D mesh from a file
 * - Extrude 2D polygons to create 3D polyhedra (prisms)
 * - Create multiple layers with specified thickness
 * - Generate lateral faces connecting the layers
 * - Handle different polygon types (triangles, quads, etc.)
 * - Write the extruded mesh to a new file
 *
 * The extrusion process creates layers of vertices above the base mesh,
 * then connects them to form polyhedra with top, bottom, and lateral faces.
 *
 * \param[in] [meshfile] The 2D mesh file to load (default: io/poly8-10.vtk)
 * \param[in] [outfile]  The output, extruded 3D mesh file name (default: polyhedra.vtk)
 * \param[in] [nlayers]  The number of axial layers (default: 1)
 * \param[in] [thickness_layer]  The thickness of axial layers (default: 1.0)
 *
 * \return 0 on success, 1 on failure
 *
 * \par Usage:
 * \code
 * ./ExtrudePoly [meshfile] [outfile] [nlayers] [thickness_layer]
 * $> ./ExtrudePoly poly2d.h5m poly3d.h5m 5 1.0
 * \endcode
 *
 * @author MOAB Development Team
 * @date Last updated: 2025
 */

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// at most 20 edges per polygon
#define MAXEDGES 20
// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string( MESH_DIR ) + string( "/io/poly8-10.vtk" );
string output         = string( "polyhedra.vtk" );
int layers            = 1;
double layer_thick    = 1.0;

int main( int argc, char** argv )
{
    // Get MOAB instance
    Interface* mb = new( std::nothrow ) Core;
    if( NULL == mb ) return 1;

    // Need option handling here for input filename
    if( argc > 1 )
    {
        // User has input a mesh file
        test_file_name = argv[1];
    }
    if( argc > 2 ) output = argv[2];

    if( argc > 3 ) layers = atoi( argv[3] );

    if( argc > 4 ) layer_thick = atof( argv[4] );

    std::cout << "Run: " << argv[0] << " " << test_file_name << " " << output << " " << layers << " " << layer_thick
              << "\n";
    // Load the mesh from vtk file
    MB_CHK_ERR( mb->load_mesh( test_file_name.c_str() ) );

    // Get verts entities, by type
    Range verts;
    MB_CHK_ERR( mb->get_entities_by_type( 0, MBVERTEX, verts ) );

    // Get faces, by dimension, so we stay generic to entity type
    Range faces;
    MB_CHK_ERR( mb->get_entities_by_dimension( 0, 2, faces ) );
    cout << "Number of vertices is " << verts.size() << endl;
    cout << "Number of faces is " << faces.size() << endl;

    Range edges;
    // create all edges
    MB_CHK_ERR( mb->get_adjacencies( faces, 1, true, edges, Interface::UNION ) );

    cout << "Number of edges is " << edges.size() << endl;

    std::vector< double > coords;
    int nvPerLayer = (int)verts.size();
    coords.resize( 3 * nvPerLayer );

    MB_CHK_ERR( mb->get_coords( verts, &coords[0] ) );
    // create first vertices
    Range* newVerts = new Range[layers + 1];
    newVerts[0]     = verts;  // just for convenience
    for( int ii = 0; ii < layers; ii++ )
    {
        for( int i = 0; i < nvPerLayer; i++ )
            coords[3 * i + 2] += layer_thick;

        MB_CHK_ERR( mb->create_vertices( &coords[0], nvPerLayer, newVerts[ii + 1] ) );
    }
    // for each edge, we will create layers quads
    int nquads = edges.size() * layers;
    ReadUtilIface* read_iface;
    MB_CHK_SET_ERR( mb->query_interface( read_iface ), "Error in query_interface" );

    EntityHandle start_elem, *connect;
    // Create quads
    MB_CHK_SET_ERR( read_iface->get_element_connect( nquads, 4, MBQUAD, 0, start_elem, connect ),
                    "Error in get_element_connect" );
    int nedges = (int)edges.size();

    int indexConn = 0;
    for( int j = 0; j < nedges; j++ )
    {
        EntityHandle edge = edges[j];

        const EntityHandle* conn2 = NULL;
        int num_nodes;
        MB_CHK_ERR( mb->get_connectivity( edge, conn2, num_nodes ) );
        if( 2 != num_nodes ) MB_CHK_ERR( MB_FAILURE );

        int i0 = verts.index( conn2[0] );
        int i1 = verts.index( conn2[1] );
        for( int ii = 0; ii < layers; ii++ )
        {
            connect[indexConn++] = newVerts[ii][i0];
            connect[indexConn++] = newVerts[ii][i1];
            connect[indexConn++] = newVerts[ii + 1][i1];
            connect[indexConn++] = newVerts[ii + 1][i0];
        }
    }

    int nfaces                = (int)faces.size();
    EntityHandle* allPolygons = new EntityHandle[nfaces * ( layers + 1 )];
    for( int i = 0; i < nfaces; i++ )
        allPolygons[i] = faces[i];

    // vertices are parallel to the base vertices
    int indexVerts[MAXEDGES]                  = { 0 };  // polygons with at most MAXEDGES edges
    EntityHandle newConn[MAXEDGES]            = { 0 };
    EntityHandle polyhedronConn[MAXEDGES + 2] = { 0 };

    // edges will be used to determine the lateral faces of polyhedra (prisms)
    int indexEdges[MAXEDGES] = { 0 };  // index of edges in base polygon
    for( int j = 0; j < nfaces; j++ )
    {
        EntityHandle polyg = faces[j];

        const EntityHandle* connp = NULL;
        int num_nodes;
        MB_CHK_ERR( mb->get_connectivity( polyg, connp, num_nodes ) );

        for( int i = 0; i < num_nodes; i++ )
        {
            indexVerts[i]             = verts.index( connp[i] );
            int i1                    = ( i + 1 ) % num_nodes;
            EntityHandle edgeVerts[2] = { connp[i], connp[i1] };
            // get edge adjacent to these vertices
            Range adjEdges;
            MB_CHK_ERR( mb->get_adjacencies( edgeVerts, 2, 1, false, adjEdges ) );
            if( adjEdges.size() < 1 ) MB_CHK_SET_ERR( MB_FAILURE, " did not find edge " );
            indexEdges[i] = edges.index( adjEdges[0] );
            if( indexEdges[i] < 0 ) MB_CHK_SET_ERR( MB_FAILURE, "did not find edge in range" );
        }

        for( int ii = 0; ii < layers; ii++ )
        {
            // create a polygon on each layer
            for( int i = 0; i < num_nodes; i++ )
                newConn[i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1

            MB_CHK_ERR( mb->create_element( MBPOLYGON, newConn, num_nodes, allPolygons[nfaces * ( ii + 1 ) + j] ) );

            // now create a polyhedra with top, bottom and lateral swept faces
            // first face is the bottom
            polyhedronConn[0] = allPolygons[nfaces * ii + j];
            // second face is the top
            polyhedronConn[1] = allPolygons[nfaces * ( ii + 1 ) + j];
            // next add lateral quads, in order of edges, using the start_elem
            // first layer of quads has EntityHandle from start_elem to start_elem+nedges-1
            // second layer of quads has EntityHandle from start_elem + nedges to
            // start_elem+2*nedges-1 ,etc
            for( int i = 0; i < num_nodes; i++ )
            {
                polyhedronConn[2 + i] = start_elem + ii + layers * indexEdges[i];
            }
            // Create polyhedron
            EntityHandle polyhedron;
            MB_CHK_ERR( mb->create_element( MBPOLYHEDRON, polyhedronConn, 2 + num_nodes, polyhedron ) );
        }
    }
    MB_CHK_ERR( mb->write_file( output.c_str() ) );

    delete mb;

    return 0;
}
