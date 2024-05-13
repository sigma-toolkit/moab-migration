#ifndef __mesh_utilities_hpp__
#define __mesh_utilities_hpp__

#include "RemapMPASROMS.hpp"

// Remapping related includes
#include "FiniteVolumeTools.h"

/**
 * @brief Scales the coordinates of a set of nodes.
 *
 * This function scales the coordinates of a set of nodes by a given factor.
 *
 * @param mb The MOAB interface.
 * @param nodes The vector of node handles.
 * @param R The scaling factor.
 * @param is_cartesian Flag indicating whether the coordinates are in Cartesian or spherical system.
 * @return The MOAB error code.
 */
moab::ErrorCode ScaleCoords( moab::Interface* mb,
                             std::vector< moab::EntityHandle >& nodes,
                             double R,
                             bool is_cartesian )
{
    moab::ErrorCode rval;
    double posi[3], posf[3], len = 0;

    // one by one, get the node and project it on the sphere, with a radius given
    // the center of the sphere is at 0,0,0
    for( auto nit = nodes.begin(); nit != nodes.end(); ++nit )
    {
        moab::EntityHandle nd = *nit;

        if( !is_cartesian )
        {
            rval = mb->get_coords( &nd, 1, posi );MB_CHK_ERR( rval );
            const double lat = posi[1] * 3.14159265358979323846 / 180;
            const double lon = posi[0] * 3.14159265358979323846 / 180;
            posf[0]          = cos( lat ) * cos( lon );  // x coordinate
            posf[1]          = cos( lat ) * sin( lon );  // y
            posf[2]          = sin( lat );               // z
            // spherical_to_cart( posi[1], posi[0], R, posf );
            // dbgprint( nd << " lat=" << posi[1] << ", lon=" << posi[0] << "; Cartesian = [" << posf[0] << ", " << posf[1] << ", " << posf[2] << "]" );
            // std::cout << "ERROR: We do not know how to scale z-coordinate if not cartesian\n";
            // exit(1);
        }
        else
        {
            rval = mb->get_coords( &nd, 1, posf );MB_CHK_ERR( rval );
        }

        len = std::sqrt( posf[0] * posf[0] + posf[1] * posf[1] + posf[2] * posf[2] );
        // len = std::sqrt( posf[0] * posf[0] + posf[1] * posf[1] );
        if( len < 1e-12 )
        {
            std::cout << nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2]
                      << ": Failed with length == 0." << std::endl;
            return moab::MB_FAILURE;
        }

        if( nit == nodes.begin() ) std::cout << "Length of the MPAS node magnitude: " << len << std::endl;

        // rescale to radius
        posf[0] *= R / len;
        posf[1] *= R / len;
        posf[2] *= R / len;
        // posf[2] = 0.0;
        // if( is_cartesian ) posf[2] *= R / len;
        // else
        //     posf[2] /= len;

        // if (is_threed)
        //     dbgprint( nd << " X=" << posf[0] << ", Y=" << posf[1] << ", Z = " << posf[2]  );
        rval = mb->set_coords( &nd, 1, posf );MB_CHK_ERR( rval );
    }
    return moab::MB_SUCCESS;
}

moab::ErrorCode ExtrudePolygonsToPolyhedra( RuntimeContext& context,
                                            std::vector< double >& layer_thickness,
                                            moab::EntityHandle& poly2dset,
                                            moab::EntityHandle& outputset,
                                            const bool is_mpas,
                                            int zlayers )
{
    using namespace moab;
    using namespace std;
    ErrorCode rval;
    Interface* mb = context.moab_interface;
    if( outputset == 0 )
    {
        rval = mb->create_meshset( moab::MESHSET_SET, outputset );MB_CHK_SET_ERR( rval, "Can't create new set" );
    }

    // Get verts entities, by type
    Range verts, edges, faces;
    // rval = mb->get_entities_by_type( poly2dset, MBVERTEX, verts );MB_CHK_ERR( rval );
    rval = mb->get_entities_by_dimension( poly2dset, 0, verts );MB_CHK_ERR( rval );

    // Get faces, by dimension, so we stay generic to entity type
    rval = mb->get_entities_by_dimension( poly2dset, 2, faces );MB_CHK_ERR( rval );
    // std::cout << "Number of 2D entities: vertices = " << verts.size() << ", faces = " << faces.size() << endl;

    // add the initial faces to the first set
    if( is_mpas )
    {
        rval = mb->add_entities( outputset, verts );MB_CHK_ERR( rval );
        rval = mb->add_entities( outputset, faces );MB_CHK_ERR( rval );
    }

    // Create all edges
    rval = mb->get_adjacencies( faces, 1, true, edges, Interface::UNION );MB_CHK_ERR( rval );

    const size_t nverts = verts.size();
    const size_t nedges = edges.size();
    const size_t nfaces = faces.size();
    const int nlayers   = static_cast< int >( layer_thickness.size() / nfaces );
    const size_t nquads = nedges * nlayers;
    std::vector< double > coords( 3 * nverts );

    // output some information
    if( context.proc_id == 0 )
    {
        std::cout << " Input 2D " << ( is_mpas ? "MPAS" : "ROMS" ) << " Mesh details ::" << std::endl;
        std::cout << "\tNumber of Vertices = " << nverts << std::endl;
        std::cout << "\t          Edges    = " << edges.size() << std::endl;
        std::cout << "\t          Faces    = " << nfaces << std::endl;
    }

    // get the vertex coordinates for the polygonal mesh
    rval = mb->get_coords( verts, &coords[0] );MB_CHK_ERR( rval );

    Tag gidTag = mb->globalId_tag();

    Tag parentTag;
    rval = mb->tag_get_handle( "ColumnParent", 1, moab::MB_TYPE_INTEGER, parentTag,
                               moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( rval );

    std::vector< int > gidData( nverts ), gidParentVertexData( nverts ), gidParentFaceData( nfaces );
    rval = mb->tag_get_data( gidTag, verts, gidParentVertexData.data() );MB_CHK_ERR( rval );
    rval = mb->tag_get_data( gidTag, faces, gidParentFaceData.data() );MB_CHK_ERR( rval );

    std::vector< int > layerchildren;
    // create first vertices
    Range* newVerts = new Range[nlayers + 1];
    newVerts[0]     = verts;  // just for convenience
    for( int ii = 0; ii < nlayers; ii++ )
    {
        for( size_t i = 0; i < nverts; i++ )
        {
            // coords[3 * i + 2] -= 0.5;
            // if( false )
            {
                Range eladjs;
                const EntityHandle vtx = verts[i];
                rval                   = mb->get_adjacencies( &vtx, 1, 2, false, eladjs, Interface::UNION );MB_CHK_ERR( rval );

                if( eladjs.size() )
                {
                    double thickness = 0.0;
                    double invweight = 0.0;
                    for( size_t k = 0; k < eladjs.size(); ++k )
                    {
                        int il = faces.index( eladjs[k] );
                        if( il < 0 ) continue;

                        if( layer_thickness[il * nlayers + ii] < 0 )
                        {
                            // printf( "Thickness value for layer %d: element %zu  = %f\n", ii, k,
                            //        layer_thickness[il * nlayers + ii] );
                            // exit( 1 );
                            continue;
                        }
                        invweight += 1.0;
                        thickness += layer_thickness[il * nlayers + ii];
                    }
                    thickness /= invweight;

                    if( !is_mpas && thickness < 1e-10 )
                    {
                        printf( "Thickness value for layer %d: vertex %zu, adj = %zu = %f\n", ii, i, eladjs.size(),
                                thickness );
                        // exit( 1 );
                        thickness = 0.0;
                    }

                    // Subtract or Add depending on the z-Direction to extrude
                    coords[3 * i + 2] -= thickness;
                }
                else
                    printf( "Vertex %zu has no adjacencies, coord = %f\n", i, coords[3 * i + 2] );
            }
        }

        rval = mb->create_vertices( &coords[0], nverts, newVerts[ii + 1] );MB_CHK_ERR( rval );
        std::iota( gidData.begin(), gidData.end(), nverts * ( 1 + ii ) );
        rval = mb->tag_set_data( gidTag, newVerts[ii + 1], gidData.data() );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( parentTag, newVerts[ii + 1], gidParentVertexData.data() );MB_CHK_ERR( rval );

        rval = mb->add_entities( outputset, newVerts[ii + 1] );MB_CHK_ERR( rval );
    }

    EntityHandle start_elem;
    std::vector< EntityHandle > allPolygons;
    if( is_mpas )
    {
        // for each edge, we will create nlayers quads
        ReadUtilIface* read_iface;
        rval = mb->query_interface( read_iface );MB_CHK_SET_ERR( rval, "Error in query_interface" );

        // Create quads
        EntityHandle* connect;
        rval = read_iface->get_element_connect( nquads, 4, MBQUAD, 0, start_elem, connect );MB_CHK_SET_ERR( rval, "Error in get_element_connect" );
        Range quads( start_elem, start_elem + nquads );

        // ---------------------------------------------------------------------------
        int indexConn = 0;
        for( size_t j = 0; j < nedges; j++ )
        {
            EntityHandle edge = edges[j];

            const EntityHandle* conn2 = NULL;
            int nnodes;
            rval = mb->get_connectivity( edge, conn2, nnodes );MB_CHK_ERR( rval );
            if( 2 != nnodes ) MB_CHK_ERR( MB_FAILURE );

            int i0 = verts.index( conn2[0] );
            int i1 = verts.index( conn2[1] );
            for( int ii = 0; ii < nlayers; ii++ )
            {
                connect[indexConn++] = newVerts[ii][i0];
                connect[indexConn++] = newVerts[ii][i1];
                connect[indexConn++] = newVerts[ii + 1][i1];
                connect[indexConn++] = newVerts[ii + 1][i0];
            }
        }

        // rval = mb->add_entities( outputset, quads );MB_CHK_ERR( rval );

        std::vector< int > allPolygonsGID( nquads );
        // allPolygonsGID.resize( nfaces * ( nlayers + 1 ) + nquads );

        // std::vector< int > allPolygonsGID;
        // GIDS for lateral quads will be at the end of the list
        // std::iota( allPolygonsGID.begin(), allPolygonsGID.end(), static_cast< int >( nfaces * ( nlayers + 1 ) ) + 1 );
        // TODO: this fails. Need to fix
        // rval = mb->tag_set_data( gidTag, quads, allPolygonsGID.data() );MB_CHK_ERR( rval );

        // next allocate for the x-y extruded faces
        allPolygons.resize( nfaces * ( nlayers + 1 ) );
        // rval = mb->tag_get_data( gidTag, faces, allPolygonsGID.data() );MB_CHK_ERR( rval );
        for( size_t i = 0; i < nfaces; i++ )
        {
            allPolygons[i] = faces[i];
        }
    }

    // vertices are parallel to the base vertices
    int gidElem                               = 1;
    int ipolygon                              = nfaces;
    int indexVerts[MAXEDGES]                  = { 0 };  // polygons with at most MAXEDGES edges
    EntityHandle polyhedronConn[MAXEDGES + 2] = { 0 };
    EntityHandle vertexConn[MAXEDGES * 2]     = { 0 };
    // edges will be used to determine the lateral faces of polyhedra (prisms)
    int indexEdges[MAXEDGES] = { 0 };  // index of edges in base polygon
    std::vector< int > minlevelFace, maxlevelFace;
    moab::Tag minlvlTag, maxlvlTag;
    moab::Tag mpas_soltags[nvars], mpas_soltags_new[nvars];
    std::vector< double > src_data( mpas_zreflevels * nvars );
    if( is_mpas )
    {
        rval = mb->tag_get_handle( "minLevelCell", 1, moab::MB_TYPE_INTEGER, minlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( rval );
        minlevelFace.resize( nfaces );

        rval = mb->tag_get_data( minlvlTag, faces, minlevelFace.data() );MB_CHK_ERR( rval );

        rval = mb->tag_get_handle( "maxLevelCell", 1, moab::MB_TYPE_INTEGER, maxlvlTag, moab::MB_TAG_DENSE );MB_CHK_ERR( rval );
        maxlevelFace.resize( nfaces );
        rval = mb->tag_get_data( maxlvlTag, faces, maxlevelFace.data() );MB_CHK_ERR( rval );

        for( size_t it = 0; it < nvars; ++it )
        {
            rval = mb->tag_get_handle( mpas_tagnames[it], mpas_zreflevels, moab::MB_TYPE_DOUBLE, mpas_soltags[it],
                                       moab::MB_TAG_DENSE );MB_CHK_ERR( rval );

            rval = mb->tag_get_handle( mpas_ele_tagnames[it], 1, moab::MB_TYPE_DOUBLE, mpas_soltags_new[it],
                                       moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( rval );
        }
    }

    for( int ii = 0; ii < nlayers; ii++ )
    {
        for( size_t j = 0; j < nfaces; j++ )
        {
            // only add this extruded MPAS element if it is within the accepted layer mask
            if( is_mpas && ( ii + 1 < minlevelFace[j] || ii + 1 > maxlevelFace[j] ) ) continue;

            const EntityHandle polyg = faces[j];
            const EntityType etype   = mb->type_from_handle( polyg );
            const int polyGID        = gidParentFaceData[j];

            // printf( "Polygon %d has %d nodes\n", j, nnodes );

            const EntityHandle* connp = nullptr;
            int nnodes;
            rval = mb->get_connectivity( polyg, connp, nnodes );MB_CHK_ERR( rval );

            std::vector< int > vecents( nnodes + 2, polyGID );
            if( is_mpas )
            {
                const int orig_nodes = nnodes;

                // account for padded polygons
                while( connp[nnodes - 2] == connp[nnodes - 1] && nnodes > 3 )
                    nnodes--;

                std::copy( connp, connp + nnodes, vertexConn );
                // we had padded entities
                if( orig_nodes != nnodes )
                {
                    rval = mb->set_connectivity( polyg, vertexConn, nnodes );MB_CHK_ERR( rval );
                }

                for( int i = 0; i < nnodes; i++ )
                {
                    indexVerts[i]             = verts.index( connp[i] );
                    int i1                    = ( i + 1 ) % nnodes;
                    EntityHandle edgeVerts[2] = { connp[i], connp[i1] };
                    // get edge adjacent to these vertices
                    Range adjEdges;
                    rval = mb->get_adjacencies( edgeVerts, 2, 1, false, adjEdges );MB_CHK_ERR( rval );
                    if( adjEdges.size() < 1 ) MB_CHK_SET_ERR( MB_FAILURE, " did not find edge " );
                    indexEdges[i] = edges.index( adjEdges[0] );
                    if( indexEdges[i] < 0 ) MB_CHK_SET_ERR( MB_FAILURE, "did not find edge in range" );
                }

                for( size_t it = 0; it < nvars; ++it )
                {
                    // get the source data from tag
                    rval = mb->tag_get_data( mpas_soltags[it], &polyg, 1, src_data.data() + it * zlayers );MB_CHK_ERR( rval );
                }
            }
            else
            {
                for( int i = 0; i < nnodes; i++ )
                {
                    // vertexConn[nnodes + i] = connp[i];
                    indexVerts[i] = verts.index( connp[i] );
                }
            }

            {
                // only add this extruded MPAS element if it is within the accepted layer mask
                // if( is_mpas && ( ii + 1 < minlevelFace[j] || ii + 1 > maxlevelFace[j] ) ) continue;

                // create a polygon on each layer
                if( is_mpas )
                {
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[nnodes + i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1
                }
                else
                {
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[i] = newVerts[ii + 1][indexVerts[i]];  // vertices in layer ii+1
                    for( int i = 0; i < nnodes; i++ )
                        vertexConn[nnodes + i] = newVerts[ii][indexVerts[i]];  // vertices in layer ii+1
                }

                EntityHandle polyhedron;
                if( is_mpas )
                {
                    rval =
                        mb->create_element( etype, &vertexConn[nnodes], nnodes, allPolygons[nfaces * ( ii + 1 ) + j] );MB_CHK_ERR( rval );
                    // allPolygonsGID[nfaces * ( ii + 1 ) + j] = ipolygon++;

                    // now create a polyhedra with top, bottom and lateral swept faces
                    // first face is the bottom
                    polyhedronConn[0] = allPolygons[nfaces * ii + j];
                    // next add lateral quads, in order of edges, using the start_elem
                    // first layer of quads has EntityHandle from start_elem to start_elem+nedges-1
                    // second layer of quads has EntityHandle from start_elem + nedges to
                    // start_elem+2*nedges-1 ,etc
                    for( int i = 0; i < nnodes; i++ )
                    {
                        polyhedronConn[1 + i] = start_elem + ii + nlayers * indexEdges[i];
                    }
                    // second face is the top
                    polyhedronConn[1 + nnodes] = allPolygons[nfaces * ( ii + 1 ) + j];

                    // Create polyhedron
                    rval = mb->create_element( MBPOLYHEDRON, polyhedronConn, 2 + nnodes, polyhedron );MB_CHK_ERR( rval );

                    rval = mb->add_entities( outputset, polyhedronConn, nnodes + 2 );MB_CHK_ERR( rval );

                    rval = mb->tag_set_data( parentTag, polyhedronConn, nnodes + 2, vecents.data() );MB_CHK_ERR( rval );

                    for( size_t it = 0; it < nvars; ++it )
                    {
                        // get the source data from tag
                        rval = mb->tag_set_data( mpas_soltags_new[it], &polyhedron, 1,
                                                 src_data.data() + it * zlayers + ii );MB_CHK_ERR( rval );
                    }

                    rval = mb->tag_set_data( gidTag, &allPolygons[nfaces * ( ii + 1 ) + j], 1, &ipolygon );MB_CHK_ERR( rval );
                    ipolygon++;
                }
                else
                {
                    EntityType extrudedType;
                    switch( etype )
                    {
                        case MBTRI:
                            extrudedType = MBPRISM;
                            break;
                        case MBQUAD:
                            extrudedType = MBHEX;
                            break;
                        default:
                            MB_CHK_SET_ERR( MB_FAILURE, "Unsupported standard 2D element type found for extrusion." );
                    }
                    rval = mb->create_element( extrudedType, vertexConn, 2 * nnodes, polyhedron );MB_CHK_ERR( rval );
                }

                rval = mb->add_entities( outputset, &polyhedron, 1 );MB_CHK_ERR( rval );
                rval = mb->tag_set_data( gidTag, &polyhedron, 1, &gidElem );MB_CHK_ERR( rval );
                gidElem++;

                rval = mb->tag_set_data( parentTag, &polyhedron, 1, &polyGID );MB_CHK_ERR( rval );
            }
        }
    }

    // output some information
    if( context.proc_id == 0 )
    {
        std::cout << " Output 3D " << ( is_mpas ? "MPAS" : "ROMS" ) << " Mesh details ::" << std::endl;
        std::cout << "\tNumber of Vertices = " << nverts * ( nlayers + 1 ) << std::endl;
        std::cout << "\t          Elements = " << gidElem - 1 << std::endl;
    }

    return moab::MB_SUCCESS;
}

#endif  // __mesh_utilities_hpp__
