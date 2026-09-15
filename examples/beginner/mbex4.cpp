/** @example mbex4.cpp
 * \brief beginner tutorial, example 4: Demonstrates creating a structured mesh
 *
 * In this example, we create a 2x2x2 mesh that is identical to the
 * previous example. However, in this case we will use the structured
 * mesh interface since the mesh we created is logically
 * structured. There are many advantages to using the structured mesh
 * interface...such as memory savings, speed, ease-of-use...
 *
 * \author Milad Fatenejad
 */

// The moab/Core.hpp header file is needed for all MOAB work...
#include "moab/Core.hpp"

// The moab/ScdInterface.hpp contains code which defines the moab
// structured mesh interface
#include "moab/ScdInterface.hpp"

#include <iostream>
#include <cmath>

// ****************
// *              *
// *     main     *
// *              *
// ****************
int main()
{
    moab::ErrorCode rval;
    moab::Core mbint;

    // ***********************
    // *   Create the Mesh   *
    // ***********************

    // First, lets make the mesh. It will be a 100 by 100 uniform grid
    // (there will be 100x100 quads, 101x101 vertexes) with dx = dy =
    // 0.1. Unlike the previous example, we will first make the mesh,
    // then set the coordinates one at a time.

    const unsigned NI = 100;
    const unsigned NJ = 100;

    // moab::ScdInterface is the structured mesh interface class for
    // MOAB.
    moab::ScdInterface* scdint;

    // Tell MOAB that our mesh is structured:
    MB_CHK_SET_ERR( mbint.query_interface( scdint ), "mbint.query_interface failed" );

    // Create the mesh:
    moab::ScdBox* scdbox = NULL;
    MB_CHK_SET_ERR( scdint->construct_box( moab::HomCoord( 0, 0, 0 ), moab::HomCoord( NI, NJ, 0 ), NULL, 0, scdbox ),
                    "scdint->construct_box failed" );

    // MOAB knows to make quads instead of hexes because the last start
    // and end indexes are the same (0). Note that it is still a "3D"
    // mesh because each vertex coordinate is still defined using three
    // numbers although every element in the mesh is a quadrilateral.

    // ******************************
    // *   Set Vertex Coordinates   *
    // ******************************

    // The "NULL" and "0" arguments in the call to construct_box are
    // where we could specify the vertex coordinates. Since we didn't
    // give any coordinates, every vertex is given a position of 0,0,0
    // by default. Now we will set the vertex coordinates...

    const double DX = 0.1;
    const double DY = 0.1;

    for( unsigned i = 0; i < NI + 1; i++ )
        for( unsigned j = 0; j < NJ + 1; j++ )
        {
            // First, get the entity handle:
            moab::EntityHandle handle = scdbox->get_vertex( i, j );

            // Compute the coordinate:
            double coord[3] = { DX * i, DY * j, 0.0 };

            // Change the coordinate of the vertex:
            mbint.set_coords( &handle, 1, coord );
        }

    // *******************
    // *   Attach Tags   *
    // *******************

    // The vertex coordinates have been defined, now let's attach some
    // data to the mesh. In MOAB this is done using "tags". Tags are
    // little bits of information that can be attached to any mesh
    // entity. In our example, we want to create two tags. The
    // "temperature" tag will be attached to each quad and will be 1
    // double. The "velocity" tag will be attached to each vertex and
    // will be an array of two doubles.

    // zero and twozeros represent the initial tag

    // Create the tags:
    moab::Tag temp_tag;
    double temp_default_value = 0.0;
    rval                      = mbint.tag_get_handle( "temperature", 1, moab::MB_TYPE_DOUBLE, temp_tag,
                                                      moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &temp_default_value );MB_CHK_SET_ERR( rval, "mbint.tag_get_handle(temperature) failed" );

    moab::Tag vel_tag;
    double vel_default_value[2] = { 0.0, 0.0 };
    rval = mbint.tag_get_handle( "velocity", 2, moab::MB_TYPE_DOUBLE, vel_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT,
                                 vel_default_value );MB_CHK_SET_ERR( rval, "mbint.tag_get_handle(velocity) failed" );

    // Note that when we created each tag, we specified two flags:
    //
    // The moab::MB_TAG_DENSE flag tells MOAB that this is a dense
    // tag. Dense tags will get automatically assigned to entities which
    // have continuous handles. For this example, this means that once
    // we set a tag on one vertex, memory will be allocated for
    // assigning the tag to all vertexes. The same is true of
    // quads. Dense tags are much more efficient when assigning tags to
    // lots of entities. If you only want to assign a tag to a few
    // entities, it is more efficient to use sparse tags
    // (moab::MB_TAG_SPARSE).
    //
    // The moab::MB_TAG_CREAT flag tells MOAB to create the tag if it
    // doesn't already exist.

    // The tags have now been created, now we have to attach them to
    // entities and set their values. NOTE: I am going to do this in a
    // manner which emphasizes clarity - this is not the most efficient
    // approach - that will be saved for later tutorial examples.

    // Loop through each quad and set the temperature:
    for( unsigned i = 0; i < NI; i++ )
        for( unsigned j = 0; j < NJ; j++ )
        {
            // Get the handle for this quad:
            moab::EntityHandle handle = scdbox->get_element( i, j );

            // Compute the temperature...
            double xc          = DX * ( i + 0.5 );
            double yc          = DY * ( j + 0.5 );
            double r           = std::sqrt( xc * xc + yc * yc );
            double temperature = std::exp( -0.5 * r );

            // Set the temperature on a single quad:
            MB_CHK_SET_ERR( mbint.tag_set_data( temp_tag, &handle, 1, &temperature ),
                            "mbint.tag_set_data(temp_tag) failed" );
        }

    // Loop through each vertex and set the velocity:
    for( unsigned i = 0; i < NI + 1; i++ )
        for( unsigned j = 0; j < NJ + 1; j++ )
        {
            // Get the handle for this vertex:
            moab::EntityHandle handle = scdbox->get_vertex( i, j );
            double velocity[2]        = { static_cast< double >( i ), static_cast< double >( j ) };

            // Set the velocity on a vertex:
            MB_CHK_SET_ERR( mbint.tag_set_data( vel_tag, &handle, 1, velocity ), "mbint.tag_set_data(vel_tag) failed" );
        }

    // ***************************
    // *   Write Mesh to Files   *
    // ***************************

    // NOTE: Some visualization software (such as VisIt) may not
    // interpret the velocity tag as a vector and you may not be able to
    // plot it. But you should be able to plot the temperature on top of
    // the mesh.

    MB_CHK_SET_ERR( mbint.write_file( "mbex4.vtk" ), "write_file(mbex4.vtk) failed" );

    return 0;
}
