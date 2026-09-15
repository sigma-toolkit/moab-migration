/*
 * edge_maps_test.cpp
 * will take source and target spherical meshes, and create edge maps using intersection file
 *  between them
 *  First, every edge in source and target will be decomposed in edges that are adjacent to
 *   the intersection polygons
 *
 */
#include <iostream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "moab/Core.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif
#include "moab/earthsystem/intx_mesh/Intx2MeshEdges.hpp"
#include "moab/earthsystem/intx_mesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"
#include <cmath>

using namespace moab;

int main( int argc, char* argv[] )
{

    std::string sourceFile, targetFile, intersectionFile, mapEdgeTargetFile;
    //sourceFile =
    //        "../sandbox/MeshFiles/e3sm/edge_maps/source_1.h5m";  // it also has data associated to edges
    //targetFile =
    //        "../sandbox/MeshFiles/e3sm/edge_maps/target_1.h5m";  //
    intersectionFile  = "intx_edges.h5m";
    mapEdgeTargetFile = "target_edge.nc";

    ProgOptions opts;
    opts.addOpt< std::string >( "source,s", "first mesh filename (source)", &sourceFile );
    opts.addOpt< std::string >( "target,t", "second mesh filename (target)", &targetFile );
    opts.addOpt< std::string >( "intersectionFile,i", "output intersection file", &intersectionFile );
    opts.addOpt< std::string >( "edgeTarget,p", "output map edge target file", &mapEdgeTargetFile );

    double R             = 1.;  // input
    double epsrel        = 1.e-12;
    double boxeps        = 1.e-4;
    double areaTolerance = 5.e-12;
    intersectionFile     = "intx.h5m";
    opts.addOpt< double >( "radius,R", "radius for model intx", &R );
    opts.addOpt< double >( "epsilon,e", "relative error in intx", &epsrel );
    opts.addOpt< double >( "boxerror,b", "relative error for box boundaries", &boxeps );
    opts.addOpt< double >( "areaTol,a", "area recovery tolerance", &areaTolerance );

    opts.addOpt< void >( "writeFiles,w", "write files of interest" );
    opts.addOpt< void >( "kdtreeOption,k", "use kd tree for intersection" );

    opts.parseCommandLine( argc, argv );

    bool write_files_rank = opts.numOptSet( "writeFiles" ) > 0;
    bool brute_force      = opts.numOptSet( "kdtreeOption" ) > 0;

    int rank = 0, size = 1;
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    std::string optsRead = ( size == 1 ? ""
                                       : std::string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" ) +
                                             std::string( ";PARALLEL_RESOLVE_SHARED_ENTS" ) );
#else
    std::string optsRead;
#endif
    // read meshes in 2 file sets
    ErrorCode rval;
    Core* mb = new Core();
    EntityHandle sf1, sf2, outputSet;

    // create meshsets and load files

    rval = mb->create_meshset( MESHSET_SET, sf1 );MB_CHK_ERR( rval );
    rval = mb->create_meshset( MESHSET_SET, sf2 );MB_CHK_ERR( rval );
    if( 0 == rank ) std::cout << "Loading mesh file " << sourceFile << "\n";
    rval = mb->load_file( sourceFile.c_str(), &sf1, optsRead.c_str() );MB_CHK_ERR( rval );
    if( 0 == rank ) std::cout << "Loading mesh file " << targetFile << "\n";
    rval = mb->load_file( targetFile.c_str(), &sf2, optsRead.c_str() );MB_CHK_ERR( rval );

    if( 0 == rank )
    {
        std::cout << "Radius:  " << R << "\n";
        std::cout << "relative eps:  " << epsrel << "\n";
        std::cout << "box eps:  " << boxeps << "\n";
        if( brute_force )
            std::cout << " use kd tree for intersection \n";
        else
            std::cout << " use advancing front for intersection \n";

        std::cout << " area tolerance:" << areaTolerance << "\n";
        std::cout << " target edge file: " << mapEdgeTargetFile << "\n";
    }
    rval = mb->create_meshset( MESHSET_SET, outputSet );MB_CHK_ERR( rval );

    // fix radius of both meshes, to be consistent with input R
    rval = moab::IntxUtils::ScaleToRadius( mb, sf1, R );MB_CHK_ERR( rval );
    rval = moab::IntxUtils::ScaleToRadius( mb, sf2, R );MB_CHK_ERR( rval );

#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = ParallelComm::get_pcomm( mb, 0 );
#endif
    Intx2MeshEdges worker( mb );
    IntxAreaUtils areaAdaptor;

    worker.set_error_tolerance( R * epsrel );
    worker.set_box_error( boxeps );
#ifdef MOAB_HAVE_MPI
    worker.set_parallel_comm( pcomm );
#endif
    // worker.SetEntityType(moab::MBQUAD);
    worker.set_radius_source_mesh( R );
    worker.set_radius_destination_mesh( R );
    // worker.enable_debug();

    rval = worker.FindMaxEdges( sf1, sf2 );MB_CHK_ERR( rval );

    EntityHandle covering_set;
#ifdef MOAB_HAVE_MPI
    if( size > 1 )
    {
        Range local_verts;
        rval = worker.build_processor_euler_boxes( sf2, local_verts );MB_CHK_ERR( rval );  // output also the local_verts
        if( write_files_rank )
        {
            std::stringstream outf;
            outf << "second_mesh" << rank << ".h5m";
            rval = mb->write_file( outf.str().c_str(), 0, 0, &sf2, 1 );MB_CHK_ERR( rval );
        }
    }
    if( size > 1 )
    {
        double elapsed = MPI_Wtime();
        rval           = mb->create_meshset( moab::MESHSET_SET, covering_set );MB_CHK_SET_ERR( rval, "Can't create new set" );
        bool gnomonic      = true;
        int order          = 0;     // we should not need ghost layers here
        bool include_edges = true;  // this is by default false; make it true for this case, for edge maps computation
        rval               = worker.construct_covering_set( sf1, covering_set, gnomonic, order, include_edges );MB_CHK_ERR( rval );  // lots of communication if mesh is distributed very differently
        elapsed = MPI_Wtime() - elapsed;
        if( 0 == rank ) std::cout << "\nTime to communicate the mesh = " << elapsed << std::endl;
        if( write_files_rank )
        {
            std::stringstream cof;
            cof << "covering_mesh" << rank << ".h5m";
            rval = mb->write_file( cof.str().c_str(), 0, 0, &covering_set, 1 );MB_CHK_ERR( rval );
        }
    }
    else
#endif
        covering_set = sf1;

    if( 0 == rank ) std::cout << "Computing intersections ..\n";
#ifdef MOAB_HAVE_MPI
    double elapsed = MPI_Wtime();
#endif
    if( brute_force )
    {
        rval = worker.intersect_meshes_kdtree( covering_set, sf2, outputSet );MB_CHK_SET_ERR( rval, "failed to intersect meshes with slow method" );
    }
    else
    {
        rval = worker.intersect_meshes( covering_set, sf2, outputSet );MB_CHK_SET_ERR( rval, "failed to intersect meshes" );
    }
#ifdef MOAB_HAVE_MPI
    elapsed = MPI_Wtime() - elapsed;
    if( 0 == rank ) std::cout << "\nTime to compute the intersection between meshes = " << elapsed << std::endl;
#endif
    if( write_files_rank )
    {
        std::stringstream outf;
        outf << "intersect" << rank << ".h5m";
        rval = mb->write_file( outf.str().c_str(), 0, 0, &outputSet, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
    }
    double intx_area    = areaAdaptor.area_on_sphere( mb, outputSet, R );
    double arrival_area = areaAdaptor.area_on_sphere( mb, sf2, R );
    std::cout << "On rank : " << rank << " arrival area: " << arrival_area << "  intersection area:" << intx_area
              << " rel error: " << fabs( ( intx_area - arrival_area ) / arrival_area ) << "\n";

#ifdef MOAB_HAVE_MPI
#ifdef MOAB_HAVE_HDF5_PARALLEL
    std::ostringstream intx_str;
    intx_str << "p" << pcomm->size() << "_" << intersectionFile;
    rval = mb->write_file( intx_str.str().c_str(), 0, "PARALLEL=WRITE_PART", &outputSet, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
    if( 0 == rank ) std::cout << " Wrote intx file: " << intx_str.str() << "\n";
#else
    // write intx set on rank 0, in serial; we cannot write in parallel
    if( 0 == rank )
    {
        rval = mb->write_file( intersectionFile.c_str(), 0, 0, &outputSet, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
    }
#endif

#else
    rval = mb->write_file( intersectionFile.c_str(), 0, 0, &outputSet, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
#endif
    moab::Tag fractionTag;
    moab::Tag numSubTag;
    moab::Tag areaDiffTag;
    moab::Tag areaTag;
    double defVal = 0.;
    rval = mb->tag_get_handle( "EdgeRecoveryFraction", 1, MB_TYPE_DOUBLE, fractionTag, MB_TAG_DENSE | MB_TAG_CREAT,
                               &defVal );MB_CHK_SET_ERR( rval, "can't create edge recovery fraction tag" );
    rval = mb->tag_get_handle( "NumSubEnts", 1, MB_TYPE_DOUBLE, numSubTag, MB_TAG_DENSE | MB_TAG_CREAT, &defVal );MB_CHK_SET_ERR( rval, "can't create NumSubEnts tag" );
    rval = mb->tag_get_handle( "AreaDiff", 1, MB_TYPE_DOUBLE, areaDiffTag, MB_TAG_DENSE | MB_TAG_CREAT, &defVal );MB_CHK_SET_ERR( rval, "can't create AreaDiff tag" );
    rval = mb->tag_get_handle( "Area", 1, MB_TYPE_DOUBLE, areaTag, MB_TAG_DENSE | MB_TAG_CREAT, &defVal );MB_CHK_SET_ERR( rval, "can't create Area tag" );
    /*    std::map<EntityHandle, std::vector<EntityHandle>> edgeVertices; // for each recovered edge, the chain of vertices that form subedges
    std::map<EntityHandle, std::vector<int>> edgePolygons; // for each recovered edge, the list of intersected polygons;
    moab::Range recoveredPolys;*/

    //bool sourceEdgeMap = false;
    rval = worker.EdgeSplits( areaTolerance );MB_CHK_SET_ERR( rval, "failed to compute edge splits for target" );
#ifdef MOAB_HAVE_MPI
#ifdef MOAB_HAVE_HDF5_PARALLEL
    std::ostringstream h5mFile;
    h5mFile << "p" << pcomm->size() << "_targetWithEdges.h5m";
    rval = mb->write_file( h5mFile.str().c_str(), 0, "PARALLEL=WRITE_PART", &sf2, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
    if( 0 == rank ) std::cout << " Wrote file with edge mapping info: " << h5mFile.str() << "\n";
#endif
#endif

#ifdef MOAB_HAVE_PNETCDF
#ifdef MOAB_HAVE_MPI
    std::ostringstream file_str;
    file_str << "p" << pcomm->size() << "_" << mapEdgeTargetFile;
    rval = worker.write_edge_map_parallel( file_str.str().c_str() );MB_CHK_SET_ERR( rval, "failed to write edge map for target" );
    if( 0 == rank ) std::cout << " Wrote netcdf file with edge mapping info: " << file_str.str() << "\n";
#else
#ifdef MOAB_HAVE_NETCDF
    rval = worker.write_edge_map( mapEdgeTargetFile.c_str() );MB_CHK_SET_ERR( rval, "failed to write edge map for target" );
#endif
#endif
#endif

    delete mb;
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
