/**
 * @example iMOABApplyMap.cpp
 * @brief iMOAB example: define source/target/intersection apps, load meshes and a map file,
 *        associate map to vertex/edge/face, set analytical data, project to target, write for visualization.
 *
 * Call sequence demonstrated:
 * - iMOAB_Initialize / RegisterApplication (source, target, intersection)
 * - iMOAB_LoadMesh (source and target)
 * - Either: iMOAB_LoadMapFile (map from disk) with source/target entity types (vertex/edge/face)
 *   Or: iMOAB_ComputeMeshIntersectionOnSphere + iMOAB_ComputeScalarProjectionWeights
 * - iMOAB_DefineTagStorage (source and target, entity type: vertex, edge, or face)
 * - iMOAB_SetDoubleTagStorage (analytical field on source)
 * - iMOAB_ApplyScalarProjectionWeights
 * - iMOAB_WriteMesh (target with projected data)
 *
 * @note Requires MOAB built with TempestRemap and NetCDF (for map file I/O).
 *
 * @par Usage:
 * @code
 * # With precomputed map file:
 * iMOABApplyMap -s source.h5m -t target.h5m -m map.nc -o target_with_field.h5m
 *
 * # Without map file (compute intersection and weights on the fly):
 * iMOABApplyMap -s source.h5m -t target.h5m -o target_with_field.h5m
 *
 * # Entity types: --source_entity face --target_entity face (default). Options: vertex, edge, face.
 * @endcode
 */

#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP

#include "moab/iMOAB.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CHECKIERR( ierr, message ) \
    do {                           \
        if( 0 != ( ierr ) )         \
        {                           \
            std::cerr << ( message ) << std::endl; \
            return 1;               \
        }                           \
    } while( 0 )

static double analytical_field_vertex( double x, double y, double z )
{
    double r = std::sqrt( x * x + y * y + z * z );
    if( r < 1.e-12 ) return 0.0;
    double lat = std::asin( z / r );
    double lon = std::atan2( y, x );
    return std::sin( 2 * lat ) * std::cos( lon );
}

static double analytical_field_element( int index, int total )
{
    (void)total;
    return 0.5 + 0.5 * std::sin( index * 0.1 );
}

int main( int argc, char* argv[] )
{
    std::string sourceMeshFile, targetMeshFile, mapFile, outputFile;
    std::string sourceEntityStr = "face", targetEntityStr = "face";
    bool haveMapFile = false;

    for( int i = 1; i < argc; ++i )
    {
        std::string arg = argv[i];
        if( arg == "-s" && i + 1 < argc )
            sourceMeshFile = argv[++i];
        else if( arg == "-t" && i + 1 < argc )
            targetMeshFile = argv[++i];
        else if( arg == "-m" && i + 1 < argc )
        {
            mapFile = argv[++i];
            haveMapFile = true;
        }
        else if( arg == "-o" && i + 1 < argc )
            outputFile = argv[++i];
        else if( arg == "--source_entity" && i + 1 < argc )
            sourceEntityStr = argv[++i];
        else if( arg == "--target_entity" && i + 1 < argc )
            targetEntityStr = argv[++i];
        else if( arg == "-h" || arg == "--help" )
        {
            std::cout << "Usage: " << argv[0]
                      << " -s <source.h5m> -t <target.h5m> [-m <map.nc>] -o <output.h5m>\n"
                      << "      [--source_entity vertex|edge|face] [--target_entity vertex|edge|face]\n"
                      << "Without -m: compute intersection and weights (sphere) then apply.\n";
            return 0;
        }
    }

    // Apply built-in defaults so the binary can run as a ctest without arguments.
    if( sourceMeshFile.empty() ) sourceMeshFile = std::string( MOAB_MESH_DIR ) + "unittest/srcWithSolnTag.h5m";
    if( targetMeshFile.empty() ) targetMeshFile = std::string( MOAB_MESH_DIR ) + "unittest/outTri15_8.h5m";
    if( outputFile.empty() ) outputFile = "/tmp/imoab_applymap_test_output.h5m";
    if( !haveMapFile )
    {
        mapFile    = std::string( MOAB_MESH_DIR ) + "unittest/mapNE20_FV15.nc";
        haveMapFile = true;
    }

    if( sourceMeshFile.empty() || targetMeshFile.empty() || outputFile.empty() )
    {
        std::cerr << "Required: -s source_mesh -t target_mesh -o output_mesh. Optional: -m map_file\n";
        return 1;
    }

    int srcEntityType = IMOAB_FACE_ENTITY;
    int tgtEntityType = IMOAB_FACE_ENTITY;
    if( sourceEntityStr == "vertex" )
        srcEntityType = IMOAB_VERTEX_ENTITY;
    else if( sourceEntityStr == "edge" )
        srcEntityType = IMOAB_EDGE_ENTITY;
    else if( sourceEntityStr == "face" )
        srcEntityType = IMOAB_FACE_ENTITY;
    else
        std::cerr << "Unknown --source_entity " << sourceEntityStr << ", using face\n";

    if( targetEntityStr == "vertex" )
        tgtEntityType = IMOAB_VERTEX_ENTITY;
    else if( targetEntityStr == "edge" )
        tgtEntityType = IMOAB_EDGE_ENTITY;
    else if( targetEntityStr == "face" )
        tgtEntityType = IMOAB_FACE_ENTITY;
    else
        std::cerr << "Unknown --target_entity " << targetEntityStr << ", using face\n";

#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
    MPI_Comm comm = MPI_COMM_WORLD;
#endif

    ErrCode ierr = iMOAB_Initialize( argc, argv );
    CHECKIERR( ierr, "iMOAB_Initialize failed" );

    int srcAppID, tgtAppID, intxAppID;
    iMOAB_AppID srcPID  = &srcAppID;
    iMOAB_AppID tgtPID  = &tgtAppID;
    iMOAB_AppID intxPID = &intxAppID;

    int srcCompID = 1, tgtCompID = 2, intxCompID = 10;
    ierr = iMOAB_RegisterApplication( "SOURCE",
#ifdef MOAB_HAVE_MPI
                                      &comm,
#else
                                      nullptr,
#endif
                                      &srcCompID, srcPID );
    CHECKIERR( ierr, "Register source app failed" );

    ierr = iMOAB_RegisterApplication( "TARGET",
#ifdef MOAB_HAVE_MPI
                                       &comm,
#else
                                       nullptr,
#endif
                                       &tgtCompID, tgtPID );
    CHECKIERR( ierr, "Register target app failed" );

    ierr = iMOAB_RegisterApplication( "INTERSECTION",
#ifdef MOAB_HAVE_MPI
                                      &comm,
#else
                                      nullptr,
#endif
                                      &intxCompID, intxPID );
    CHECKIERR( ierr, "Register intersection app failed" );

    const char* read_opts = "";
    int num_ghost_layers  = 0;
    CHECKIERR( iMOAB_LoadMesh( srcPID, sourceMeshFile.c_str(), read_opts, &num_ghost_layers ), "Load source mesh failed" );
    CHECKIERR( iMOAB_LoadMesh( tgtPID, targetMeshFile.c_str(), read_opts, &num_ghost_layers ), "Load target mesh failed" );
    CHECKIERR( ierr, "Load target mesh failed" );

    int nverts[3], nelem[3], nedges[3], nfaces[3];
    CHECKIERR( iMOAB_GetMeshInfo( srcPID, nverts, nelem, nullptr, nullptr, nullptr, nedges, nfaces ), "Get source mesh info failed" );
    CHECKIERR( ierr, "Get source mesh info failed" );
    std::cout << "Source: " << nverts[0] << " vertices, " << nelem[0] << " elements, " << nedges[0] << " edges, " << nfaces[0] << " faces (local).\n";

    CHECKIERR( iMOAB_GetMeshInfo( tgtPID, nverts, nelem, nullptr, nullptr, nullptr, nedges, nfaces ), "Get target mesh info failed" );
    std::cout << "Target: " << nverts[0] << " vertices, " << nelem[0] << " elements, " << nedges[0] << " edges, " << nfaces[0] << " faces (local).\n";

    const std::string weights_id = "scalar";
    const std::string source_tag  = "analytical_field";
    const std::string target_tag = "projected_field";

    if( haveMapFile )
    {
#ifdef MOAB_HAVE_NETCDF
        int src_disc = 3, tgt_disc = 3;
        int arearead = 0;
        int src_ent = srcEntityType, tgt_ent = tgtEntityType;
        ierr = iMOAB_LoadMapFile( srcPID, tgtPID, intxPID, &src_disc, &tgt_disc, &arearead,
                                 weights_id.c_str(), mapFile.c_str(), &src_ent, &tgt_ent );
        CHECKIERR( ierr, "Load map file failed" );
        std::cout << "Loaded map from " << mapFile << " (source_entity=" << sourceEntityStr
                  << ", target_entity=" << targetEntityStr << ").\n";
#else
        std::cerr << "Map file requires NetCDF; rebuild MOAB with NetCDF.\n";
        return 1;
#endif
    }
    else
    {
        ierr = iMOAB_ComputeMeshIntersectionOnSphere( srcPID, tgtPID, intxPID );
        CHECKIERR( ierr, "Compute mesh intersection failed" );

        const char* src_disc_method = "fv", * tgt_disc_method = "fv";
        const char* src_dof_tag = "GLOBAL_ID", * tgt_dof_tag = "GLOBAL_ID";
        int src_order = 1, tgt_order = 1;
        int fNoBubble = 1, fMonotone = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fInverseDist = 0;
        ierr = iMOAB_ComputeScalarProjectionWeights(
            intxPID, weights_id.c_str(), src_disc_method, &src_order, tgt_disc_method, &tgt_order, nullptr,
            &fNoBubble, &fMonotone, &fVolumetric, &fInverseDist, &fNoConserve, &fValidate, src_dof_tag, tgt_dof_tag );
        CHECKIERR( ierr, "Compute scalar projection weights failed" );
        std::cout << "Computed intersection and weights on sphere.\n";
    }

    int tag_type      = IMOAB_DENSE_DOUBLE_TAG;
    int ncomp         = 1;
    int tag_index_src = 0, tag_index_tgt = 0;
    ierr = iMOAB_DefineTagStorage( srcPID, source_tag.c_str(), &tag_type, &ncomp, &tag_index_src );
    CHECKIERR( ierr, "Define source tag failed" );
    ierr = iMOAB_DefineTagStorage( tgtPID, target_tag.c_str(), &tag_type, &ncomp, &tag_index_tgt );
    CHECKIERR( ierr, "Define target tag failed" );

    ierr = iMOAB_GetMeshInfo( srcPID, nverts, nelem, nullptr, nullptr, nullptr, nedges, nfaces );
    CHECKIERR( ierr, "Get source mesh info (2) failed" );

    int nent = 0;
    if( srcEntityType == IMOAB_VERTEX_ENTITY )
        nent = nverts[2];
    else if( srcEntityType == IMOAB_EDGE_ENTITY )
        nent = nedges[2];
    else if( srcEntityType == IMOAB_FACE_ENTITY )
    {
        nent = nfaces[2];
        if( nent == 0 ) nent = nelem[2];
    }
    else
        nent = nelem[2];

    if( srcEntityType == IMOAB_VERTEX_ENTITY && nent > 0 )
    {
        int coords_len = 3 * nent;
        std::vector< double > coords( coords_len );
        ierr = iMOAB_GetVisibleVerticesCoordinates( srcPID, &coords_len, coords.data() );
        CHECKIERR( ierr, "Get source vertex coordinates failed" );
        std::vector< double > vals( nent );
        for( int i = 0; i < nent; ++i )
            vals[i] = analytical_field_vertex( coords[3 * i], coords[3 * i + 1], coords[3 * i + 2] );
        int ent_type = IMOAB_VERTEX_ENTITY;
        ierr = iMOAB_SetDoubleTagStorage( srcPID, source_tag.c_str(), &nent, &ent_type, vals.data() );
        CHECKIERR( ierr, "Set source vertex tag failed" );
    }
    else if( nent > 0 )
    {
        std::vector< double > vals( nent );
        for( int i = 0; i < nent; ++i )
            vals[i] = analytical_field_element( i, nent );
        int ent_type = srcEntityType;
        ierr = iMOAB_SetDoubleTagStorage( srcPID, source_tag.c_str(), &nent, &ent_type, vals.data() );
        CHECKIERR( ierr, "Set source entity tag failed" );
    }

    int filter_weights = 0;
    ierr = iMOAB_ApplyScalarProjectionWeights( intxPID, &filter_weights, weights_id.c_str(),
                                              source_tag.c_str(), target_tag.c_str() );
    CHECKIERR( ierr, "Apply scalar projection weights failed" );

    const char* write_opts = "";
    ierr = iMOAB_WriteMesh( tgtPID, outputFile.c_str(), write_opts );
    CHECKIERR( ierr, "Write target mesh failed" );
    std::cout << "Wrote " << outputFile << " with projected field '" << target_tag << "' (visualize with Visit/ParaView).\n";

    ierr = iMOAB_DeregisterApplication( intxPID );
    CHECKIERR( ierr, "Deregister intersection failed" );
    ierr = iMOAB_DeregisterApplication( tgtPID );
    CHECKIERR( ierr, "Deregister target failed" );
    ierr = iMOAB_DeregisterApplication( srcPID );
    CHECKIERR( ierr, "Deregister source failed" );

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "iMOAB_Finalize failed" );

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

#else

#include <iostream>
int main( int, char** )
{
    std::cerr << "This example requires MOAB built with TempestRemap (MOAB_HAVE_TEMPESTREMAP).\n";
    return 1;
}

#endif
