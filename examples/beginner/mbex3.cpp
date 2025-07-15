/** @example mbex3.cpp
 * \brief beginner tutorial, example 3: Demonstrates creating sets and tags
 *
 * In this example, we read in the VTK file (mbex2.vtk) generated
 * in example 2, and demonstrate how to create sets and tags.
 * We will create a set for each hexahedron, and add a tag to
 * each set with some data.
 *
 * \author Milad Fatenejad
 */

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

// The moab/ScdInterface.hpp contains code which defines the moab
// structured mesh interface
#include "moab/ScdInterface.hpp"

#include <iostream>

// ****************
// *              *
// *     main     *
// *              *
// ****************
int main()
{
    moab::ErrorCode rval;
    moab::Core mbint;

    // **********************************
    // *   Create the Structured Mesh   *
    // **********************************

    // As before, we have to create an array defining the coordinate of
    // each vertex:
    const unsigned NUMVTX                  = 27;
    const double vertex_coords[3 * NUMVTX] = {
        0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 0, 0, 2, 0, 1, 2, 0, 2, 2, 0,

        0, 0, 1, 1, 0, 1, 2, 0, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 0, 2, 1, 1, 2, 1, 2, 2, 1,

        0, 0, 2, 1, 0, 2, 2, 0, 2, 0, 1, 2, 1, 1, 2, 2, 1, 2, 0, 2, 2, 1, 2, 2, 2, 2, 2 };

    // moab::ScdInterface is the structured mesh interface class for
    // MOAB.
    moab::ScdInterface* scdint;

    // The query_interface method will create a structured mesh instance
    // for mbcore and will point scdint to it. This is how you tell moab
    // that our moab::Core instance is going to represent a structured
    // mesh.
    MB_CHK_SET_ERR( mbint.query_interface( scdint ), "mbint.query_interface failed" );

    // Structured meshes a divided into "boxes". Each box represents a
    // little structured mesh. A single mesh, for example a
    // block-structured mesh, can contain many individual boxes. In this
    // example, we want to create a single, 2x2x2 box. The construct_box
    // method will do this for us.
    moab::ScdBox* scdbox = NULL;
    MB_CHK_SET_ERR( scdint->construct_box( moab::HomCoord( 0, 0, 0 ), moab::HomCoord( 2, 2, 2 ), vertex_coords, NUMVTX,
                                           scdbox ),
                    "scdint->construct_box failed" );

    // moab::HomCoord is a little class that is used to represent a
    // coordinate in logical space. Above, we told MOAB that we want
    // indexes which extend from 0,0,0 to 2,2,2 for vertexes. Since the
    // i,j,k start/end indexes are all different, MOAB knows that our
    // mesh consists of hexes. If we had gone from 0,0,0 to 1,1,0 then
    // MOAB would construct a 2x2x1 mesh of quadrilaterals.

    // ***************
    // *   Inspect   *
    // ***************

    // Our mesh now exists and contains a single box! Lets loop over
    // that box and manually print out the handle of each vertex/hex:
    unsigned i, j, k;

    // Print out the entity handle associated with each hex:
    for( i = 0; i < 2; ++i )
        for( j = 0; j < 2; ++j )
            for( k = 0; k < 2; ++k )
            {
                std::cout << "Hex (" << i << "," << j << "," << k << ") "
                          << "has handle: " << scdbox->get_element( i, j, k ) << std::endl;
            }

    // Print out the entity handle associated with each vertex:
    for( i = 0; i < 3; ++i )
        for( j = 0; j < 3; ++j )
            for( k = 0; k < 3; ++k )
            {
                std::cout << "Vertex (" << i << "," << j << "," << k << ") "
                          << "has handle: " << scdbox->get_vertex( i, j, k ) << std::endl;
            }

    // Of course, you can still use all of the functionality defined in
    // mbint to get coordinates, connectivity, etc...

    // ***************************
    // *   Write Mesh to Files   *
    // ***************************
    MB_CHK_SET_ERR( mbint.write_file( "mbex3.vtk" ), "write_file(mbex3.vtk) failed" );

    return 0;
}
