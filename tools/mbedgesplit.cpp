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
#include "moab/IntxMesh/Intx2MeshOnSphere.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"
#include <cmath>

using namespace moab;

int main( int argc, char* argv[] )
{

    std::string sourceFile, targetFile, intersectionFile, edgeFile;
    //sourceFile =
    //        "../sandbox/MeshFiles/e3sm/edge_maps/source_1.h5m";  // it also has data associated to edges
    //targetFile =
    //        "../sandbox/MeshFiles/e3sm/edge_maps/target_1.h5m";  //
    intersectionFile = "intx_edges.h5m";

    ProgOptions opts;
    opts.addOpt< std::string >( "source,s", "first mesh filename (source)", &sourceFile );
    opts.addOpt< std::string >( "target,t", "second mesh filename (target)", &targetFile );
    opts.addOpt< std::string >( "intersectionFile,i", "output intersection file", &intersectionFile );

    double R      = 1.;  // input
    double epsrel = 1.e-12;
    double boxeps = 1.e-4;
    double areaTolerance = 5.e-12;
    intersectionFile    = "intx.h5m";
    opts.addOpt< double >( "radius,R", "radius for model intx", &R );
    opts.addOpt< double >( "epsilon,e", "relative error in intx", &epsrel );
    opts.addOpt< double >( "boxerror,b", "relative error for box boundaries", &boxeps );
    opts.addOpt< double >( "areaTol,a", "area recovery tolerance", &areaTolerance);

    opts.addOpt<void>( "outputFraction,f", "output fraction of areas" );
    opts.addOpt<void>( "writeFiles,w", "write files of interest" );
    opts.addOpt<void>( "kdtreeOption,k", "use kd tree for intersection" );

    opts.parseCommandLine( argc, argv );

    bool output_fraction = opts.numOptSet("outputFraction") > 0;
    bool write_files_rank = opts.numOptSet("writeFiles") > 0;
    bool brute_force = opts.numOptSet("kdtreeOption") > 0;

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
    Core moab;
    Interface* mb = &moab;  // global
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
        std::cout << " use kd tree for intersection: " << brute_force << "\n";
        std::cout << " area tolerance:" << areaTolerance <<"\n";
    }
    rval = mb->create_meshset( MESHSET_SET, outputSet );MB_CHK_ERR( rval );

    // fix radius of both meshes, to be consistent with input R
    rval = moab::IntxUtils::ScaleToRadius( mb, sf1, R );MB_CHK_ERR( rval );
    rval = moab::IntxUtils::ScaleToRadius( mb, sf2, R );MB_CHK_ERR( rval );


#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = ParallelComm::get_pcomm( mb, 0 );
#endif
    Intx2MeshOnSphere worker( mb );
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
        bool gnomonic = true;
        int order = 0; // we should not need ghost layers here
        bool include_edges = true; // this is by default false; make it true for this case, for edge maps computation
        rval          = worker.construct_covering_set( sf1, covering_set, gnomonic, order, include_edges );MB_CHK_ERR( rval );  // lots of communication if mesh is distributed very differently
        elapsed = MPI_Wtime() - elapsed;
        if( 0 == rank ) std::cout << "\nTime to communicate the mesh = " << elapsed << std::endl;
        // area fraction of the covering set that needed to be communicated from other processors
        // number of elements in the covering set communicated, compared to total number of elements
        // in the covering set
        if( output_fraction )
        {
            EntityHandle comm_set;  // set with elements communicated from other tasks
            rval = mb->create_meshset( MESHSET_SET, comm_set );MB_CHK_ERR( rval );
            // see how much more different is compared to sf1
            rval = mb->unite_meshset( comm_set, covering_set );MB_CHK_ERR( rval );  // will have to subtract from covering set, initial set
            // subtract
            rval = mb->subtract_meshset( comm_set, sf1 );MB_CHK_ERR( rval );
            // compute fractions
            double area_cov_set = areaAdaptor.area_on_sphere( mb, covering_set, R );
            assert( area_cov_set > 0 );
            double comm_area = areaAdaptor.area_on_sphere( mb, comm_set, R );
            // more important is actually the number of elements communicated
            int num_cov_cells, num_comm_cells;
            rval = mb->get_number_entities_by_dimension( covering_set, 2, num_cov_cells );MB_CHK_ERR( rval );
            rval = mb->get_number_entities_by_dimension( comm_set, 2, num_comm_cells );MB_CHK_ERR( rval );
            double fraction_area = comm_area / area_cov_set;
            double fraction_num_cells =
                (double)num_comm_cells / num_cov_cells;  // determine min, max, average of these fractions

            double max_fraction_area, max_fraction_num_cells, min_fraction_area, min_fraction_num_cells;
            double average_fraction_area, average_fraction_num_cells;
            MPI_Reduce( &fraction_area, &max_fraction_area, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
            MPI_Reduce( &fraction_num_cells, &max_fraction_num_cells, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
            MPI_Reduce( &fraction_area, &min_fraction_area, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
            MPI_Reduce( &fraction_num_cells, &min_fraction_num_cells, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
            MPI_Reduce( &fraction_area, &average_fraction_area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            MPI_Reduce( &fraction_num_cells, &average_fraction_num_cells, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
            average_fraction_area /= size;
            average_fraction_num_cells /= size;
            if( rank == 0 )
            {
                std::cout << " fraction area:      min: " << min_fraction_area << " max: " << max_fraction_area
                          << " average :" << average_fraction_area << " \n";
                std::cout << " fraction num cells: min: " << min_fraction_num_cells
                          << " max: " << max_fraction_num_cells << " average: " << average_fraction_num_cells << " \n";
            }
        }
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
    // the output set does not have the intx vertices on the boundary shared, so they will be


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
    rval = mb->write_file( intersectionFile.c_str(), 0, "PARALLEL=WRITE_PART", &outputSet, 1 );MB_CHK_SET_ERR( rval, "failed to write intx file" );
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
    bool sourceEdgeMap = true;
    moab::Tag fractionTag;
    moab::Tag numSubTag;
    moab::Tag areaDiffTag;
    moab::Tag areaTag;
    double defVal = 0.;
    rval = mb->tag_get_handle( "EdgeRecoveryFraction", 1, MB_TYPE_DOUBLE, fractionTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                   &defVal );MB_CHK_SET_ERR( rval, "can't create edge recovery fraction tag" );
    rval = mb->tag_get_handle( "NumSubEnts", 1, MB_TYPE_DOUBLE, numSubTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                       &defVal );MB_CHK_SET_ERR( rval, "can't create NumSubEnts tag" );
    rval = mb->tag_get_handle( "AreaDiff", 1, MB_TYPE_DOUBLE, areaDiffTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                           &defVal );MB_CHK_SET_ERR( rval, "can't create AreaDiff tag" );
    rval = mb->tag_get_handle( "Area", 1, MB_TYPE_DOUBLE, areaTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                               &defVal );MB_CHK_SET_ERR( rval, "can't create Area tag" );
    std::map<EntityHandle, std::vector<EntityHandle>> edgeVertices; // for each recovered edge, the chain of vertices that form subedges
    std::map<EntityHandle, std::vector<int>> edgePolygons; // for each recovered edge, the list of intersected polygons;
    moab::Range recoveredPolys;
    rval = moab::IntxUtils::EdgeMap(mb, sf1, outputSet, sourceEdgeMap,
        edgeVertices, edgePolygons, recoveredPolys, areaTolerance );MB_CHK_SET_ERR( rval, "failed to compute edge map for source" );
#ifdef MOAB_HAVE_PNETCDF
#ifdef MOAB_HAVE_MPI
    if (size > 1)
    {
       std::ostringstream file_str;
       file_str << "source_edge_p" << pcomm->size() << ".nc";
       rval = moab::IntxUtils::write_edge_map_parallel(file_str.str().c_str(), pcomm, mb, sf1, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for source file" );
    }
    else
    {
#ifdef MOAB_HAVE_NETCDF
       rval = moab::IntxUtils::write_edge_map("source_edge.nc", mb, sf1, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for source file" );
#endif
    }
#endif
#else
#ifdef MOAB_HAVE_NETCDF
    rval = moab::IntxUtils::write_edge_map("source_edge.nc", mb, sf1, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for source file" );
#endif
#endif
    if (1==size)
    {
        rval = mb->write_file("source_withEdges.h5m", 0, 0, &sf1, 1);MB_CHK_SET_ERR( rval, "failed rewrite initial source" );
    }

    recoveredPolys.clear();
    edgeVertices.clear();
    edgePolygons.clear();
    sourceEdgeMap = false;
    rval = moab::IntxUtils::EdgeMap(mb, sf2, outputSet, sourceEdgeMap,
        edgeVertices, edgePolygons, recoveredPolys, areaTolerance );MB_CHK_SET_ERR( rval, "failed to compute edge map for target" );
#ifdef MOAB_HAVE_PNETCDF
#ifdef MOAB_HAVE_MPI
    if (size > 1)
    {
        std::ostringstream file_str;
        file_str << "target_edge_p" << pcomm->size() << ".nc";
        rval = moab::IntxUtils::write_edge_map_parallel(file_str.str().c_str(), pcomm, mb, sf2, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for target" );
    }
    else
    {
#ifdef MOAB_HAVE_NETCDF
        rval = moab::IntxUtils::write_edge_map("target_edge.nc", mb, sf2, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for target" );
#endif
    }
#endif
#else
#ifdef MOAB_HAVE_NETCDF
    rval = moab::IntxUtils::write_edge_map("target_edge.nc", mb, sf2, edgeVertices, edgePolygons, recoveredPolys);MB_CHK_SET_ERR( rval, "failed to write edge map for target" );
#endif
#endif
    if (1==size)
    {
        rval = mb->write_file("target_withEdges.h5m", 0, 0, &sf2, 1);MB_CHK_SET_ERR( rval, "failed rewrite initial target" );
    }
    return 0;
}
