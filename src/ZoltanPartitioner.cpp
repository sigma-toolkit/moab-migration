/* $Id$
 *
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

/** Get a mesh from MOAB and write a Zoltan partition set back into MOAB and to
 *  a file
 *
 *
 */

#include <iostream>
#include <cassert>
#include <sstream>
#include <algorithm>

#include "moab/ZoltanPartitioner.hpp"
#include "moab/Interface.hpp"
#include "Internals.hpp"
#include "moab/Range.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
// used for gnomonic projection
#include "moab/IntxMesh/IntxUtils.hpp"

using namespace moab;

#define RR \
    if( MB_SUCCESS != result ) return result

static double* Points      = NULL;
static int* GlobalIds      = NULL;
static int NumPoints       = 0;
static int* NumEdges       = NULL;
static int* NborGlobalId   = NULL;
static int* NborProcs      = NULL;
static double* ObjWeights  = NULL;
static double* EdgeWeights = NULL;
static int* Parts          = NULL;

const bool debug = false;

ZoltanPartitioner::ZoltanPartitioner( Interface* impl,
#ifdef MOAB_HAVE_MPI
                                      ParallelComm* parcomm,
#endif
                                      const bool use_coords,
                                      int argc,
                                      char** argv
                                      )
    : PartitionerBase< int >( impl,
                              use_coords
#ifdef MOAB_HAVE_MPI
                              ,
                              parcomm
#endif
                              ),
      myZZ( NULL ), myNumPts( 0 ), argcArg( argc ), argvArg( argv )
{
}

ZoltanPartitioner::~ZoltanPartitioner()
{
    if( NULL != myZZ ) delete myZZ;
}

ErrorCode ZoltanPartitioner::balance_mesh( const char* zmethod,
                                           const char* other_method,
                                           const bool write_as_sets,
                                           const bool write_as_tags )
{
    if( !strcmp( zmethod, "RR" ) && !strcmp( zmethod, "RCB" ) && !strcmp( zmethod, "RIB" ) &&
        !strcmp( zmethod, "HSFC" ) && !strcmp( zmethod, "Hypergraph" ) && !strcmp( zmethod, "PHG" ) &&
        !strcmp( zmethod, "PARMETIS" ) && !strcmp( zmethod, "OCTPART" ) )
    {
        std::cout << "ERROR node " << mbpc->proc_config().proc_rank() << ": Method must be "
                  << "RR, RCB, RIB, HSFC, Hypergraph (PHG), PARMETIS, or OCTPART" << std::endl;
        return MB_FAILURE;
    }

    std::vector< double > pts;  // x[0], y[0], z[0], ... from MOAB
    std::vector< int > ids;     // point ids from MOAB
    std::vector< int > adjs, length;
    Range elems;

    // Get a mesh from MOAB and divide it across processors.

    ErrorCode result;

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        result = assemble_graph( 3, pts, ids, adjs, length, elems );RR;
    }

    myNumPts = mbInitializePoints( (int)ids.size(), &pts[0], &ids[0], &adjs[0], &length[0] );

    // Initialize Zoltan.  This is a C call.  The simple C++ code
    // that creates Zoltan objects does not keep track of whether
    // Zoltan_Initialize has been called.

    float version;

    Zoltan_Initialize( argcArg, argvArg, &version );

    // Create Zoltan object.  This calls Zoltan_Create.
    if( NULL == myZZ ) myZZ = new Zoltan( mbpc->comm() );

    if( NULL == zmethod || !strcmp( zmethod, "RCB" ) )
        SetRCB_Parameters();
    else if( !strcmp( zmethod, "RIB" ) )
        SetRIB_Parameters();
    else if( !strcmp( zmethod, "HSFC" ) )
        SetHSFC_Parameters();
    else if( !strcmp( zmethod, "Hypergraph" ) || !strcmp( zmethod, "PHG" ) )
        if( NULL == other_method )
            SetHypergraph_Parameters( "auto" );
        else
            SetHypergraph_Parameters( other_method );
    else if( !strcmp( zmethod, "PARMETIS" ) )
    {
        if( NULL == other_method )
            SetPARMETIS_Parameters( "RepartGDiffusion" );
        else
            SetPARMETIS_Parameters( other_method );
    }
    else if( !strcmp( zmethod, "OCTPART" ) )
    {
        if( NULL == other_method )
            SetOCTPART_Parameters( "2" );
        else
            SetOCTPART_Parameters( other_method );
    }

    // Call backs:

    myZZ->Set_Num_Obj_Fn( mbGetNumberOfAssignedObjects, NULL );
    myZZ->Set_Obj_List_Fn( mbGetObjectList, NULL );
    myZZ->Set_Num_Geom_Fn( mbGetObjectSize, NULL );
    myZZ->Set_Geom_Multi_Fn( mbGetObject, NULL );
    myZZ->Set_Num_Edges_Multi_Fn( mbGetNumberOfEdges, NULL );
    myZZ->Set_Edge_List_Multi_Fn( mbGetEdgeList, NULL );

    // Perform the load balancing partitioning

    int changes;
    int numGidEntries;
    int numLidEntries;
    int numImport;
    ZOLTAN_ID_PTR importGlobalIds;
    ZOLTAN_ID_PTR importLocalIds;
    int* importProcs;
    int* importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalIds;
    ZOLTAN_ID_PTR exportLocalIds;
    int* exportProcs;
    int* exportToPart;

    int rc = myZZ->LB_Partition( changes, numGidEntries, numLidEntries, numImport, importGlobalIds, importLocalIds,
                                 importProcs, importToPart, numExport, exportGlobalIds, exportLocalIds, exportProcs,
                                 exportToPart );

    rc = mbGlobalSuccess( rc );

    if( !rc )
    {
        mbPrintGlobalResult( "==============Result==============", myNumPts, numImport, numExport, changes );
    }
    else
    {
        return MB_FAILURE;
    }

    // take results & write onto MOAB partition sets

    int* assignment;

    mbFinalizePoints( (int)ids.size(), numExport, exportLocalIds, exportProcs, &assignment );

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        result = write_partition( mbpc->proc_config().proc_size(), elems, assignment, write_as_sets, write_as_tags );

        if( MB_SUCCESS != result ) return result;

        free( (int*)assignment );
    }

    // Free the memory allocated for lists returned by LB_Parition()

    myZZ->LB_Free_Part( &importGlobalIds, &importLocalIds, &importProcs, &importToPart );
    myZZ->LB_Free_Part( &exportGlobalIds, &exportLocalIds, &exportProcs, &exportToPart );

    // Implementation note:  A Zoltan object contains an MPI communicator.
    //   When the Zoltan object is destroyed, it uses it's MPI communicator.
    //   So it is important that the Zoltan object is destroyed before
    //   the MPI communicator is destroyed.  To ensure this, dynamically
    //   allocate the Zoltan object, so you can explicitly destroy it.
    //   If you create a Zoltan object on the stack, it's destructor will
    //   be invoked atexit, possibly after the communicator's
    //   destructor.

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::repartition( std::vector< double >& x,
                                          std::vector< double >& y,
                                          std::vector< double >& z,
                                          int StartID,
                                          const char* zmethod,
                                          Range& localGIDs )
{
    //
    int nprocs = mbpc->proc_config().proc_size();
    int rank   = mbpc->proc_config().proc_rank();
    clock_t t  = clock();
    // form pts and ids, as in assemble_graph
    std::vector< double > pts;  // x[0], y[0], z[0], ... from MOAB
    pts.resize( x.size() * 3 );
    std::vector< int > ids;  // point ids from MOAB
    ids.resize( x.size() );
    for( size_t i = 0; i < x.size(); i++ )
    {
        pts[3 * i]     = x[i];
        pts[3 * i + 1] = y[i];
        pts[3 * i + 2] = z[i];
        ids[i]         = StartID + (int)i;
    }
    // do not get out until done!

    Points       = &pts[0];
    GlobalIds    = &ids[0];
    NumPoints    = (int)x.size();
    NumEdges     = NULL;
    NborGlobalId = NULL;
    NborProcs    = NULL;
    ObjWeights   = NULL;
    EdgeWeights  = NULL;
    Parts        = NULL;

    float version;
    if( rank == 0 ) std::cout << "Initializing zoltan..." << std::endl;

    Zoltan_Initialize( argcArg, argvArg, &version );

    // Create Zoltan object.  This calls Zoltan_Create.
    if( NULL == myZZ ) myZZ = new Zoltan( mbpc->comm() );

    if( NULL == zmethod || !strcmp( zmethod, "RCB" ) )
        SetRCB_Parameters();
    else if( !strcmp( zmethod, "RIB" ) )
        SetRIB_Parameters();
    else if( !strcmp( zmethod, "HSFC" ) )
        SetHSFC_Parameters();

    // set # requested partitions
    char buff[10];
    snprintf( buff, 10, "%d", nprocs );
    int retval = myZZ->Set_Param( "NUM_GLOBAL_PARTITIONS", buff );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    // request all, import and export
    retval = myZZ->Set_Param( "RETURN_LISTS", "ALL" );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    myZZ->Set_Num_Obj_Fn( mbGetNumberOfAssignedObjects, NULL );
    myZZ->Set_Obj_List_Fn( mbGetObjectList, NULL );
    myZZ->Set_Num_Geom_Fn( mbGetObjectSize, NULL );
    myZZ->Set_Geom_Multi_Fn( mbGetObject, NULL );
    myZZ->Set_Num_Edges_Multi_Fn( mbGetNumberOfEdges, NULL );
    myZZ->Set_Edge_List_Multi_Fn( mbGetEdgeList, NULL );

    // Perform the load balancing partitioning

    int changes;
    int numGidEntries;
    int numLidEntries;
    int num_import;
    ZOLTAN_ID_PTR import_global_ids, import_local_ids;
    int* import_procs;
    int* import_to_part;
    int num_export;
    ZOLTAN_ID_PTR export_global_ids, export_local_ids;
    int *assign_procs, *assign_parts;

    if( rank == 0 )
        std::cout << "Computing partition using " << ( zmethod ? zmethod : "RCB" ) << " method for " << nprocs
                  << " processors..." << std::endl;

    retval = myZZ->LB_Partition( changes, numGidEntries, numLidEntries, num_import, import_global_ids, import_local_ids,
                                 import_procs, import_to_part, num_export, export_global_ids, export_local_ids,
                                 assign_procs, assign_parts );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    if( rank == 0 )
    {
        std::cout << " time to LB_partition " << ( clock() - t ) / (double)CLOCKS_PER_SEC << "s. \n";
        t = clock();
    }

    std::sort( import_global_ids, import_global_ids + num_import, std::greater< int >() );
    std::sort( export_global_ids, export_global_ids + num_export, std::greater< int >() );

    Range iniGids( (EntityHandle)StartID, (EntityHandle)StartID + x.size() - 1 );
    Range imported, exported;
    std::copy( import_global_ids, import_global_ids + num_import, range_inserter( imported ) );
    std::copy( export_global_ids, export_global_ids + num_export, range_inserter( exported ) );
    localGIDs = subtract( iniGids, exported );
    localGIDs = unite( localGIDs, imported );
    // Free data structures allocated by Zoltan::LB_Partition
    retval = myZZ->LB_Free_Part( &import_global_ids, &import_local_ids, &import_procs, &import_to_part );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;
    retval = myZZ->LB_Free_Part( &export_global_ids, &export_local_ids, &assign_procs, &assign_parts );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::partition_inferred_mesh( EntityHandle sfileset,
                                                      size_t num_parts,
                                                      int part_dim,
                                                      const bool write_as_sets,
                                                      int projection_type )
{
    ErrorCode result;

    moab::Range elverts;
    result = mbImpl->get_entities_by_dimension( sfileset, part_dim, elverts );RR;

    std::vector< double > elcoords( elverts.size() * 3 );
    result = mbImpl->get_coords( elverts, &elcoords[0] );RR;

    std::vector< std::vector< EntityHandle > > part_assignments( num_parts );
    int part, proc;

    // Loop over coordinates
    for( size_t iel = 0; iel < elverts.size(); iel++ )
    {
        // Gather coordinates into temporary array
        double* ecoords = &elcoords[iel * 3];

        // do a projection if needed
        if( projection_type > 0 ) IntxUtils::transform_coordinates( ecoords, projection_type );

        // Compute the coordinate's part assignment
        myZZ->LB_Point_PP_Assign( ecoords, proc, part );

        // Store the part assignment in the return array
        part_assignments[part].push_back( elverts[iel] );
    }

    // get the partition set tag
    Tag part_set_tag = mbpc->partition_tag();

    // If some entity sets exist, let us delete those first.
    if( write_as_sets )
    {
        moab::Range oldpsets;
        result = mbImpl->get_entities_by_type_and_tag( sfileset, MBENTITYSET, &part_set_tag, NULL, 1, oldpsets,
                                                       Interface::UNION );RR;
        result = mbImpl->remove_entities( sfileset, oldpsets );RR;
    }

    size_t nparts_assigned = 0;
    for( size_t partnum = 0; partnum < num_parts; ++partnum )
    {
        std::vector< EntityHandle >& partvec = part_assignments[partnum];

        nparts_assigned += ( partvec.size() ? 1 : 0 );

        if( write_as_sets )
        {
            EntityHandle partNset;
            result = mbImpl->create_meshset( moab::MESHSET_SET, partNset );RR;
            if( partvec.size() )
            {
                result = mbImpl->add_entities( partNset, &partvec[0], partvec.size() );RR;
            }
            result = mbImpl->add_parent_child( sfileset, partNset );RR;

            // part numbering should start from 0
            int ipartnum = (int)partnum;

            // assign partitions as a sparse tag by grouping elements under sets
            result = mbImpl->tag_set_data( part_set_tag, &partNset, 1, &ipartnum );RR;
        }
        else
        {
            /* assign as a dense vector to all elements
               allocate integer-size partitions */
            std::vector< int > assignment( partvec.size(),
                                           partnum );  // initialize all values to partnum
            result = mbImpl->tag_set_data( part_set_tag, &partvec[0], partvec.size(), &assignment[0] );RR;
        }
    }

    if( nparts_assigned != num_parts )
    {
        std::cout << "WARNING: The inference yielded lesser number of parts (" << nparts_assigned
                  << ") than requested by user (" << num_parts << ").\n";
    }

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::partition_mesh_and_geometry( const double part_geom_mesh_size,
                                                          const int nparts,
                                                          const char* zmethod,
                                                          const char* other_method,
                                                          double imbal_tol,
                                                          const int part_dim,
                                                          const bool write_as_sets,
                                                          const bool write_as_tags,
                                                          const int obj_weight,
                                                          const int edge_weight,
                                                          const int projection_type,
                                                          const bool recompute_rcb_box,
                                                          const bool print_time )
{
    // should only be called in serial
    if( mbpc->proc_config().proc_size() != 1 )
    {
        std::cout << "ZoltanPartitioner::partition_mesh_and_geometry must be called in serial." << std::endl;
        return MB_FAILURE;
    }
    clock_t t = clock();
    if( NULL != zmethod && strcmp( zmethod, "RR" ) && strcmp( zmethod, "RCB" ) && strcmp( zmethod, "RIB" ) &&
        strcmp( zmethod, "HSFC" ) && strcmp( zmethod, "Hypergraph" ) && strcmp( zmethod, "PHG" ) &&
        strcmp( zmethod, "PARMETIS" ) && strcmp( zmethod, "OCTPART" ) )
    {
        std::cout << "ERROR node " << mbpc->proc_config().proc_rank() << ": Method must be "
                  << "RCB, RIB, HSFC, Hypergraph (PHG), PARMETIS, or OCTPART" << std::endl;
        return MB_FAILURE;
    }

    bool part_geom = false;
    if( 0 == strcmp( zmethod, "RR" ) || 0 == strcmp( zmethod, "RCB" ) || 0 == strcmp( zmethod, "RIB" ) ||
        0 == strcmp( zmethod, "HSFC" ) )
        part_geom = true;       // so no adjacency / edges needed
    std::vector< double > pts;  // x[0], y[0], z[0], ... from MOAB
    std::vector< int > ids;     // point ids from MOAB
    std::vector< int > adjs, length, parts;
    std::vector< double > obj_weights, edge_weights;
    Range elems;

    // Get a mesh from MOAB and diide it across processors.

    ErrorCode result = MB_SUCCESS;

    // short-circuit everything if RR partition is requested
    if( !strcmp( zmethod, "RR" ) )
    {
        if( part_geom_mesh_size < 0. )
        {
            // get all elements
            result = mbImpl->get_entities_by_dimension( 0, part_dim, elems );RR;
            if( elems.empty() ) return MB_FAILURE;
            // make a trivial assignment vector
            std::vector< int > assign_vec( elems.size() );
            int num_per = elems.size() / nparts;
            int extra   = elems.size() % nparts;
            if( extra ) num_per++;
            int nstart = 0;
            for( int i = 0; i < nparts; i++ )
            {
                if( i == extra ) num_per--;
                std::fill( &assign_vec[nstart], &assign_vec[nstart + num_per], i );
                nstart += num_per;
            }

            result = write_partition( nparts, elems, &assign_vec[0], write_as_sets, write_as_tags );RR;
        }
        else
        {
        }

        return result;
    }

    std::cout << "Assembling graph..." << std::endl;

    if( part_geom_mesh_size < 0. )
    {
        // if (!part_geom) {
        result = assemble_graph( part_dim, pts, ids, adjs, length, elems, part_geom, projection_type );RR;
    }
    else
    {
      MB_CHK_SET_ERR( MB_FAILURE, "Geometry partitions not supported.\n" );
    }
    if( print_time )
    {
        std::cout << " time to assemble graph: " << ( clock() - t ) / (double)CLOCKS_PER_SEC << "s. \n";
        t = clock();
    }
    double* o_wgt = NULL;
    double* e_wgt = NULL;
    if( obj_weights.size() > 0 ) o_wgt = &obj_weights[0];
    if( edge_weights.size() > 0 ) e_wgt = &edge_weights[0];

    myNumPts = mbInitializePoints( (int)ids.size(), &pts[0], &ids[0], &adjs[0], &length[0], o_wgt, e_wgt, &parts[0],
                                   part_geom );

    // Initialize Zoltan.  This is a C call.  The simple C++ code
    // that creates Zoltan objects does not keep track of whether
    // Zoltan_Initialize has been called.

    if( print_time )
    {
        std::cout << " time to initialize points: " << ( clock() - t ) / (double)CLOCKS_PER_SEC << "s. \n";
        t = clock();
    }
    float version;

    std::cout << "Initializing zoltan..." << std::endl;

    Zoltan_Initialize( argcArg, argvArg, &version );

    // Create Zoltan object.  This calls Zoltan_Create.
    if( NULL == myZZ ) myZZ = new Zoltan( mbpc->comm() );

    if( NULL == zmethod || !strcmp( zmethod, "RCB" ) )
        SetRCB_Parameters( recompute_rcb_box );
    else if( !strcmp( zmethod, "RIB" ) )
        SetRIB_Parameters();
    else if( !strcmp( zmethod, "HSFC" ) )
        SetHSFC_Parameters();
    else if( !strcmp( zmethod, "Hypergraph" ) || !strcmp( zmethod, "PHG" ) )
    {
        if( NULL == other_method || ( other_method[0] == '\0' ) )
            SetHypergraph_Parameters( "auto" );
        else
            SetHypergraph_Parameters( other_method );

        if( imbal_tol )
        {
            std::ostringstream str;
            str << imbal_tol;
            myZZ->Set_Param( "IMBALANCE_TOL", str.str().c_str() );  // no debug messages
        }
    }
    else if( !strcmp( zmethod, "PARMETIS" ) )
    {
        if( NULL == other_method )
            SetPARMETIS_Parameters( "RepartGDiffusion" );
        else
            SetPARMETIS_Parameters( other_method );
    }
    else if( !strcmp( zmethod, "OCTPART" ) )
    {
        if( NULL == other_method )
            SetOCTPART_Parameters( "2" );
        else
            SetOCTPART_Parameters( other_method );
    }

    // set # requested partitions
    char buff[10];
    snprintf( buff, 10, "%d", nparts );
    int retval = myZZ->Set_Param( "NUM_GLOBAL_PARTITIONS", buff );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    // request only partition assignments
    retval = myZZ->Set_Param( "RETURN_LISTS", "PARTITION ASSIGNMENTS" );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    if( obj_weight > 0 )
    {
        std::ostringstream str;
        str << obj_weight;
        retval = myZZ->Set_Param( "OBJ_WEIGHT_DIM", str.str().c_str() );
        if( ZOLTAN_OK != retval ) return MB_FAILURE;
    }

    if( edge_weight > 0 )
    {
        std::ostringstream str;
        str << edge_weight;
        retval = myZZ->Set_Param( "EDGE_WEIGHT_DIM", str.str().c_str() );
        if( ZOLTAN_OK != retval ) return MB_FAILURE;
    }

    // Call backs:

    myZZ->Set_Num_Obj_Fn( mbGetNumberOfAssignedObjects, NULL );
    myZZ->Set_Obj_List_Fn( mbGetObjectList, NULL );
    myZZ->Set_Num_Geom_Fn( mbGetObjectSize, NULL );
    myZZ->Set_Geom_Multi_Fn( mbGetObject, NULL );
    myZZ->Set_Num_Edges_Multi_Fn( mbGetNumberOfEdges, NULL );
    myZZ->Set_Edge_List_Multi_Fn( mbGetEdgeList, NULL );
    if( part_geom_mesh_size > 0. )
    {
        myZZ->Set_Part_Multi_Fn( mbGetPart, NULL );
    }

    // Perform the load balancing partitioning

    int changes;
    int numGidEntries;
    int numLidEntries;
    int dumnum1;
    ZOLTAN_ID_PTR dum_local, dum_global;
    int *dum1, *dum2;
    int num_assign;
    ZOLTAN_ID_PTR assign_gid, assign_lid;
    int *assign_procs, *assign_parts;

    std::cout << "Computing partition using " << ( zmethod ? zmethod : "RCB" ) << " method for " << nparts
              << " processors..." << std::endl;
#ifndef NDEBUG
#if 0
  if (NULL == zmethod || !strcmp(zmethod, "RCB"))
    Zoltan_Generate_Files(myZZ->Get_C_Handle(), (char*)zmethod, 1, 1, 0, 0);

  if ( !strcmp(zmethod, "PHG"))
      Zoltan_Generate_Files(myZZ->Get_C_Handle(), (char*)zmethod, 1, 0, 1, 1);
#endif
#endif
    retval = myZZ->LB_Partition( changes, numGidEntries, numLidEntries, dumnum1, dum_global, dum_local, dum1, dum2,
                                 num_assign, assign_gid, assign_lid, assign_procs, assign_parts );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    if( print_time )
    {
        std::cout << " time to LB_partition " << ( clock() - t ) / (double)CLOCKS_PER_SEC << "s. \n";
        t = clock();
    }
    // take results & write onto MOAB partition sets
    std::cout << "Saving partition information to MOAB..." << std::endl;

    if( part_geom_mesh_size < 0. )
    {
        // if (!part_geom) {
        result = write_partition( nparts, elems, assign_parts, write_as_sets, write_as_tags );
    }
    else
    {
      MB_CHK_SET_ERR( MB_FAILURE, "Geometry partitions not supported.\n" );
    }

    if( print_time )
    {
        std::cout << " time to write partition in memory " << ( clock() - t ) / (double)CLOCKS_PER_SEC << "s. \n";
        t = clock();
    }

    if( MB_SUCCESS != result ) return result;

    // Free the memory allocated for lists returned by LB_Parition()
    myZZ->LB_Free_Part( &assign_gid, &assign_lid, &assign_procs, &assign_parts );

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::include_closure()
{
    ErrorCode result;
    Range ents;
    Range adjs;
    std::cout << "Adding closure..." << std::endl;

    for( Range::iterator rit = partSets.begin(); rit != partSets.end(); ++rit )
    {

        // get the top-dimensional entities in the part
        result = mbImpl->get_entities_by_handle( *rit, ents, true );RR;

        if( ents.empty() ) continue;

        // get intermediate-dimensional adjs and add to set
        for( int d = mbImpl->dimension_from_handle( *ents.begin() ) - 1; d >= 1; d-- )
        {
            adjs.clear();
            result = mbImpl->get_adjacencies( ents, d, false, adjs, Interface::UNION );RR;
            result = mbImpl->add_entities( *rit, adjs );RR;
        }

        // now get vertices and add to set; only need to do for ents, not for adjs
        adjs.clear();
        result = mbImpl->get_adjacencies( ents, 0, false, adjs, Interface::UNION );RR;
        result = mbImpl->add_entities( *rit, adjs );RR;

        ents.clear();
    }

    // now go over non-part entity sets, looking for contained entities
    Range sets, part_ents;
    result = mbImpl->get_entities_by_type( 0, MBENTITYSET, sets );RR;
    for( Range::iterator rit = sets.begin(); rit != sets.end(); ++rit )
    {
        // skip parts
        if( partSets.find( *rit ) != partSets.end() ) continue;

        // get entities in this set, recursively
        ents.clear();
        result = mbImpl->get_entities_by_handle( *rit, ents, true );RR;

        // now check over all parts
        for( Range::iterator rit2 = partSets.begin(); rit2 != partSets.end(); ++rit2 )
        {
            part_ents.clear();
            result = mbImpl->get_entities_by_handle( *rit2, part_ents, false );RR;
            Range int_range = intersect( ents, part_ents );
            if( !int_range.empty() )
            {
                // non-empty intersection, add to part set
                result = mbImpl->add_entities( *rit2, &( *rit ), 1 );RR;
            }
        }
    }

    // finally, mark all the part sets as having the closure
    Tag closure_tag;
    result =
        mbImpl->tag_get_handle( "INCLUDES_CLOSURE", 1, MB_TYPE_INTEGER, closure_tag, MB_TAG_SPARSE | MB_TAG_CREAT );RR;

    std::vector< int > closure_vals( partSets.size(), 1 );
    result = mbImpl->tag_set_data( closure_tag, partSets, &closure_vals[0] );RR;

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::assemble_graph( const int dimension,
                                             std::vector< double >& coords,
                                             std::vector< int >& moab_ids,
                                             std::vector< int >& adjacencies,
                                             std::vector< int >& length,
                                             Range& elems,
                                             bool part_geom,
                                             int projection_type )
{
    // assemble a graph with vertices equal to elements of specified dimension, edges
    // signified by list of other elements to which an element is connected

    // get the elements of that dimension
    ErrorCode result = mbImpl->get_entities_by_dimension( 0, dimension, elems );
    if( MB_SUCCESS != result || elems.empty() ) return result;

    // assign global ids
    if( assign_global_ids )
    {
        EntityHandle rootset = 0;
        result               = mbpc->assign_global_ids( rootset, dimension, 1, true, true );RR;
    }

    // now assemble the graph, calling MeshTopoUtil to get bridge adjacencies through d-1
    // dimensional neighbors
    MeshTopoUtil mtu( mbImpl );
    Range adjs;
    // can use a fixed-size array 'cuz the number of lower-dimensional neighbors is limited
    // by MBCN
    int neighbors[5 * MAX_SUB_ENTITIES];
    double avg_position[3];
    int moab_id;

    // get the global id tag handle
    Tag gid = mbImpl->globalId_tag();

    for( Range::iterator rit = elems.begin(); rit != elems.end(); ++rit )
    {

        if( !part_geom )
        {
            // get bridge adjacencies
            adjs.clear();
            result = mtu.get_bridge_adjacencies( *rit, ( dimension > 0 ? dimension - 1 : 3 ), dimension, adjs );RR;

            // get the graph vertex ids of those
            if( !adjs.empty() )
            {
                assert( adjs.size() < 5 * MAX_SUB_ENTITIES );
                result = mbImpl->tag_get_data( gid, adjs, neighbors );RR;
            }

            // copy those into adjacencies vector
            length.push_back( (int)adjs.size() );
            std::copy( neighbors, neighbors + adjs.size(), std::back_inserter( adjacencies ) );
        }

        // get average position of vertices
        result = mtu.get_average_position( *rit, avg_position );RR;

        // get the graph vertex id for this element
        result = mbImpl->tag_get_data( gid, &( *rit ), 1, &moab_id );RR;

        // copy those into coords vector
        moab_ids.push_back( moab_id );
        // transform coordinates to spherical coordinates, if requested
        if( projection_type > 0 ) IntxUtils::transform_coordinates( avg_position, projection_type );

        std::copy( avg_position, avg_position + 3, std::back_inserter( coords ) );
    }

    if( debug )
    {
        std::cout << "Length vector: " << std::endl;
        std::copy( length.begin(), length.end(), std::ostream_iterator< int >( std::cout, ", " ) );
        std::cout << std::endl;
        std::cout << "Adjacencies vector: " << std::endl;
        std::copy( adjacencies.begin(), adjacencies.end(), std::ostream_iterator< int >( std::cout, ", " ) );
        std::cout << std::endl;
        std::cout << "Moab_ids vector: " << std::endl;
        std::copy( moab_ids.begin(), moab_ids.end(), std::ostream_iterator< int >( std::cout, ", " ) );
        std::cout << std::endl;
        std::cout << "Coords vector: " << std::endl;
        std::copy( coords.begin(), coords.end(), std::ostream_iterator< double >( std::cout, ", " ) );
        std::cout << std::endl;
    }

    return MB_SUCCESS;
}

ErrorCode ZoltanPartitioner::write_partition( const int nparts,
                                              Range& elems,
                                              const int* assignment,
                                              const bool write_as_sets,
                                              const bool write_as_tags )
{
    ErrorCode result;

    // get the partition set tag
    Tag part_set_tag;
    int dum_id = -1, i;
    result     = mbImpl->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_set_tag,
                                         MB_TAG_SPARSE | MB_TAG_CREAT, &dum_id );RR;

    // get any sets already with this tag, and clear them
    Range tagged_sets;
    result =
        mbImpl->get_entities_by_type_and_tag( 0, MBENTITYSET, &part_set_tag, NULL, 1, tagged_sets, Interface::UNION );RR;
    if( !tagged_sets.empty() )
    {
        result = mbImpl->clear_meshset( tagged_sets );RR;
        if( !write_as_sets )
        {
            result = mbImpl->tag_delete_data( part_set_tag, tagged_sets );RR;
        }
    }

    if( write_as_sets )
    {
        // first, create partition sets and store in vector
        partSets.clear();

        if( nparts > (int)tagged_sets.size() )
        {
            // too few partition sets - create missing ones
            int num_new = nparts - tagged_sets.size();
            for( i = 0; i < num_new; i++ )
            {
                EntityHandle new_set;
                result = mbImpl->create_meshset( MESHSET_SET, new_set );RR;
                tagged_sets.insert( new_set );
            }
        }
        else if( nparts < (int)tagged_sets.size() )
        {
            // too many partition sets - delete extras
            int num_del = tagged_sets.size() - nparts;
            for( i = 0; i < num_del; i++ )
            {
                EntityHandle old_set = tagged_sets.pop_back();
                result               = mbImpl->delete_entities( &old_set, 1 );RR;
            }
        }
        // assign partition sets to vector
        partSets.swap( tagged_sets );

        // write a tag to those sets denoting they're partition sets, with a value of the
        // proc number
        int* dum_ids = new int[nparts];
        for( i = 0; i < nparts; i++ )
            dum_ids[i] = i;

        result = mbImpl->tag_set_data( part_set_tag, partSets, dum_ids );RR;
        // found out by valgrind when we run mbpart
        delete[] dum_ids;
        dum_ids = NULL;

        // assign entities to the relevant sets
        std::vector< EntityHandle > tmp_part_sets;
        // int N = (int)elems.size();
        std::copy( partSets.begin(), partSets.end(), std::back_inserter( tmp_part_sets ) );
        /*Range::reverse_iterator riter;
        for (i = N-1, riter = elems.rbegin(); riter != elems.rend(); ++riter, i--) {
          int assigned_part = assignment[i];
          part_ranges[assigned_part].insert(*riter);
          //result = mbImpl->add_entities(tmp_part_sets[assignment[i]], &(*rit), 1); RR;
        }*/

        Range::iterator rit;
        for( i = 0, rit = elems.begin(); rit != elems.end(); ++rit, i++ )
        {
            result = mbImpl->add_entities( tmp_part_sets[assignment[i]], &( *rit ), 1 );RR;
        }
        /*for (i=0; i<nparts; i++)
        {
          result = mbImpl->add_entities(tmp_part_sets[i], part_ranges[i]); RR;
        }
        delete [] part_ranges;*/
        // check for empty sets, warn if there are any
        Range empty_sets;
        for( rit = partSets.begin(); rit != partSets.end(); ++rit )
        {
            int num_ents = 0;
            result       = mbImpl->get_number_entities_by_handle( *rit, num_ents );
            if( MB_SUCCESS != result || !num_ents ) empty_sets.insert( *rit );
        }
        if( !empty_sets.empty() )
        {
            std::cout << "WARNING: " << empty_sets.size() << " empty sets in partition: ";
            for( rit = empty_sets.begin(); rit != empty_sets.end(); ++rit )
                std::cout << *rit << " ";
            std::cout << std::endl;
        }
    }

    if( write_as_tags )
    {
        // allocate integer-size partitions
        result = mbImpl->tag_set_data( part_set_tag, elems, assignment );RR;
    }

    return MB_SUCCESS;
}

void ZoltanPartitioner::SetRCB_Parameters( const bool recompute_rcb_box )
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nRecursive Coordinate Bisection" << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );  // no debug messages
    myZZ->Set_Param( "LB_METHOD", "RCB" );  // recursive coordinate bisection
    // myZZ->Set_Param( "RCB_RECOMPUTE_BOX", "1" );  // recompute RCB box if needed ?

    // RCB parameters:

    myZZ->Set_Param( "RCB_OUTPUT_LEVEL", "1" );
    myZZ->Set_Param( "KEEP_CUTS", "1" );  // save decomposition so that we can infer partitions
    // myZZ->Set_Param("RCB_RECTILINEAR_BLOCKS", "1"); // don't split point on boundary
    if( recompute_rcb_box ) myZZ->Set_Param( "RCB_RECOMPUTE_BOX", "1" );
}

void ZoltanPartitioner::SetRIB_Parameters()
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nRecursive Inertial Bisection" << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );  // no debug messages
    myZZ->Set_Param( "LB_METHOD", "RIB" );  // Recursive Inertial Bisection

    // RIB parameters:

    myZZ->Set_Param( "KEEP_CUTS", "1" );  // save decomposition
    myZZ->Set_Param( "AVERAGE_CUTS", "1" );
}

void ZoltanPartitioner::SetHSFC_Parameters()
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nHilbert Space Filling Curve" << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );   // no debug messages
    myZZ->Set_Param( "LB_METHOD", "HSFC" );  // perform Hilbert space filling curve

    // HSFC parameters:

    myZZ->Set_Param( "KEEP_CUTS", "1" );  // save decomposition
}

void ZoltanPartitioner::SetHypergraph_Parameters( const char* phg_method )
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nHypergraph (or PHG): " << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );         // no debug messages
    myZZ->Set_Param( "LB_METHOD", "Hypergraph" );  // Hypergraph (or PHG)

    // Hypergraph (or PHG) parameters:
    myZZ->Set_Param( "PHG_COARSEPARTITION_METHOD", phg_method );  // CoarsePartitionMethod
}

void ZoltanPartitioner::SetPARMETIS_Parameters( const char* parmetis_method )
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nPARMETIS: " << parmetis_method << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );       // no debug messages
    myZZ->Set_Param( "LB_METHOD", "PARMETIS" );  // the ParMETIS library

    // PARMETIS parameters:

    myZZ->Set_Param( "PARMETIS_METHOD", parmetis_method );  // method in the library
}

void ZoltanPartitioner::SetOCTPART_Parameters( const char* oct_method )
{
    if( mbpc->proc_config().proc_rank() == 0 ) std::cout << "\nOctree Partitioning: " << oct_method << std::endl;
    // General parameters:

    myZZ->Set_Param( "DEBUG_LEVEL", "0" );      // no debug messages
    myZZ->Set_Param( "LB_METHOD", "OCTPART" );  // octree partitioning

    // OCTPART parameters:

    myZZ->Set_Param( "OCT_METHOD", oct_method );  // the SFC to be used
    myZZ->Set_Param( "OCT_OUTPUT_LEVEL", "3" );
}

int ZoltanPartitioner::mbInitializePoints( int npts,
                                           double* pts,
                                           int* ids,
                                           int* adjs,
                                           int* length,
                                           double* obj_weights,
                                           double* edge_weights,
                                           int* parts,
                                           bool part_geom )
{
    unsigned int i;
    int j;
    int *numPts, *nborProcs = NULL;
    int sum, ptsPerProc, ptsAssigned, mySize;
    MPI_Status stat;
    double* sendPts;
    int* sendIds;
    int* sendEdges  = NULL;
    int* sendNborId = NULL;
    int* sendProcs;

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        /* divide pts to start */

        numPts      = (int*)malloc( sizeof( int ) * mbpc->proc_config().proc_size() );
        ptsPerProc  = npts / mbpc->proc_config().proc_size();
        ptsAssigned = 0;

        for( i = 0; i < mbpc->proc_config().proc_size() - 1; i++ )
        {
            numPts[i] = ptsPerProc;
            ptsAssigned += ptsPerProc;
        }

        numPts[mbpc->proc_config().proc_size() - 1] = npts - ptsAssigned;

        mySize  = numPts[mbpc->proc_config().proc_rank()];
        sendPts = pts + ( 3 * numPts[0] );
        sendIds = ids + numPts[0];
        sum     = 0;  // possible no adjacency sent
        if( !part_geom )
        {
            sendEdges = length + numPts[0];

            for( j = 0; j < numPts[0]; j++ )
                sum += length[j];

            sendNborId = adjs + sum;

            for( j = numPts[0]; j < npts; j++ )
                sum += length[j];

            nborProcs = (int*)malloc( sizeof( int ) * sum );
        }
        for( j = 0; j < sum; j++ )
            if( ( i = adjs[j] / ptsPerProc ) < mbpc->proc_config().proc_size() )
                nborProcs[j] = i;
            else
                nborProcs[j] = mbpc->proc_config().proc_size() - 1;

        sendProcs = nborProcs + ( sendNborId - adjs );

        for( i = 1; i < mbpc->proc_config().proc_size(); i++ )
        {
            MPI_Send( &numPts[i], 1, MPI_INT, i, 0x00, mbpc->comm() );
            MPI_Send( sendPts, 3 * numPts[i], MPI_DOUBLE, i, 0x01, mbpc->comm() );
            MPI_Send( sendIds, numPts[i], MPI_INT, i, 0x03, mbpc->comm() );
            MPI_Send( sendEdges, numPts[i], MPI_INT, i, 0x06, mbpc->comm() );
            sum = 0;

            for( j = 0; j < numPts[i]; j++ )
                sum += sendEdges[j];

            MPI_Send( sendNborId, sum, MPI_INT, i, 0x07, mbpc->comm() );
            MPI_Send( sendProcs, sum, MPI_INT, i, 0x08, mbpc->comm() );
            sendPts += ( 3 * numPts[i] );
            sendIds += numPts[i];
            sendEdges += numPts[i];
            sendNborId += sum;
            sendProcs += sum;
        }

        free( numPts );
    }
    else
    {
        MPI_Recv( &mySize, 1, MPI_INT, 0, 0x00, mbpc->comm(), &stat );
        pts    = (double*)malloc( sizeof( double ) * 3 * mySize );
        ids    = (int*)malloc( sizeof( int ) * mySize );
        length = (int*)malloc( sizeof( int ) * mySize );
        if( obj_weights != NULL ) obj_weights = (double*)malloc( sizeof( double ) * mySize );
        MPI_Recv( pts, 3 * mySize, MPI_DOUBLE, 0, 0x01, mbpc->comm(), &stat );
        MPI_Recv( ids, mySize, MPI_INT, 0, 0x03, mbpc->comm(), &stat );
        MPI_Recv( length, mySize, MPI_INT, 0, 0x06, mbpc->comm(), &stat );
        sum = 0;

        for( j = 0; j < mySize; j++ )
            sum += length[j];

        adjs = (int*)malloc( sizeof( int ) * sum );
        if( edge_weights != NULL ) edge_weights = (double*)malloc( sizeof( double ) * sum );
        nborProcs = (int*)malloc( sizeof( int ) * sum );
        MPI_Recv( adjs, sum, MPI_INT, 0, 0x07, mbpc->comm(), &stat );
        MPI_Recv( nborProcs, sum, MPI_INT, 0, 0x08, mbpc->comm(), &stat );
    }

    Points       = pts;
    GlobalIds    = ids;
    NumPoints    = mySize;
    NumEdges     = length;
    NborGlobalId = adjs;
    NborProcs    = nborProcs;
    ObjWeights   = obj_weights;
    EdgeWeights  = edge_weights;
    Parts        = parts;

    return mySize;
}

void ZoltanPartitioner::mbFinalizePoints( int npts,
                                          int numExport,
                                          ZOLTAN_ID_PTR exportLocalIDs,
                                          int* exportProcs,
                                          int** assignment )
{
    int* MyAssignment;
    int i;
    int numPts;
    MPI_Status stat;
    int* recvA;

    /* assign pts to start */

    if( mbpc->proc_config().proc_rank() == 0 )
        MyAssignment = (int*)malloc( sizeof( int ) * npts );
    else
        MyAssignment = (int*)malloc( sizeof( int ) * NumPoints );

    for( i = 0; i < NumPoints; i++ )
        MyAssignment[i] = mbpc->proc_config().proc_rank();

    for( i = 0; i < numExport; i++ )
        MyAssignment[exportLocalIDs[i]] = exportProcs[i];

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        /* collect pts */
        recvA = MyAssignment + NumPoints;

        for( i = 1; i < (int)mbpc->proc_config().proc_size(); i++ )
        {
            MPI_Recv( &numPts, 1, MPI_INT, i, 0x04, mbpc->comm(), &stat );
            MPI_Recv( recvA, numPts, MPI_INT, i, 0x05, mbpc->comm(), &stat );
            recvA += numPts;
        }

        *assignment = MyAssignment;
    }
    else
    {
        MPI_Send( &NumPoints, 1, MPI_INT, 0, 0x04, mbpc->comm() );
        MPI_Send( MyAssignment, NumPoints, MPI_INT, 0, 0x05, mbpc->comm() );
        free( MyAssignment );
    }
}

int ZoltanPartitioner::mbGlobalSuccess( int rc )
{
    int fail = 0;
    unsigned int i;
    int* vals = (int*)malloc( mbpc->proc_config().proc_size() * sizeof( int ) );

    MPI_Allgather( &rc, 1, MPI_INT, vals, 1, MPI_INT, mbpc->comm() );

    for( i = 0; i < mbpc->proc_config().proc_size(); i++ )
    {
        if( vals[i] != ZOLTAN_OK )
        {
            if( 0 == mbpc->proc_config().proc_rank() )
            {
                mbShowError( vals[i], "Result on process " );
            }
            fail = 1;
        }
    }

    free( vals );
    return fail;
}

void ZoltanPartitioner::mbPrintGlobalResult( const char* s, int begin, int import, int exp, int change )
{
    unsigned int i;
    int* v1 = (int*)malloc( 4 * sizeof( int ) );
    int* v2 = NULL;
    int* v;

    v1[0] = begin;
    v1[1] = import;
    v1[2] = exp;
    v1[3] = change;

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        v2 = (int*)malloc( 4 * mbpc->proc_config().proc_size() * sizeof( int ) );
    }

    MPI_Gather( v1, 4, MPI_INT, v2, 4, MPI_INT, 0, mbpc->comm() );

    if( mbpc->proc_config().proc_rank() == 0 )
    {
        fprintf( stdout, "======%s======\n", s );
        for( i = 0, v = v2; i < mbpc->proc_config().proc_size(); i++, v += 4 )
        {
            fprintf( stdout, "%u: originally had %d, import %d, exp %d, %s\n", i, v[0], v[1], v[2],
                     v[3] ? "a change of partitioning" : "no change" );
        }
        fprintf( stdout, "==========================================\n" );
        fflush( stdout );

        free( v2 );
    }

    free( v1 );
}

void ZoltanPartitioner::mbShowError( int val, const char* s )
{
    if( s ) printf( "%s ", s );

    switch( val )
    {
        case ZOLTAN_OK:
            printf( "%d: SUCCESSFUL\n", mbpc->proc_config().proc_rank() );
            break;
        case ZOLTAN_WARN:
            printf( "%d: WARNING\n", mbpc->proc_config().proc_rank() );
            break;
        case ZOLTAN_FATAL:
            printf( "%d: FATAL ERROR\n", mbpc->proc_config().proc_rank() );
            break;
        case ZOLTAN_MEMERR:
            printf( "%d: MEMORY ALLOCATION ERROR\n", mbpc->proc_config().proc_rank() );
            break;
        default:
            printf( "%d: INVALID RETURN CODE\n", mbpc->proc_config().proc_rank() );
            break;
    }
}

/**********************
** call backs
**********************/

int mbGetNumberOfAssignedObjects( void* /* userDefinedData */, int* err )
{
    *err = 0;
    return NumPoints;
}

void mbGetObjectList( void* /* userDefinedData */,
                      int /* numGlobalIds */,
                      int /* numLids */,
                      ZOLTAN_ID_PTR gids,
                      ZOLTAN_ID_PTR lids,
                      int wgt_dim,
                      float* obj_wgts,
                      int* err )
{
    for( int i = 0; i < NumPoints; i++ )
    {
        gids[i] = GlobalIds[i];
        lids[i] = i;
        if( wgt_dim > 0 ) obj_wgts[i] = ObjWeights[i];
    }

    *err = 0;
}

int mbGetObjectSize( void* /* userDefinedData */, int* err )
{
    *err = 0;
    return 3;
}

void mbGetObject( void* /* userDefinedData */,
                  int /* numGlobalIds */,
                  int /* numLids */,
                  int numObjs,
                  ZOLTAN_ID_PTR /* gids */,
                  ZOLTAN_ID_PTR lids,
                  int numDim,
                  double* pts,
                  int* err )
{
    int i, id, id3;
    int next = 0;

    if( numDim != 3 )
    {
        *err = 1;
        return;
    }

    for( i = 0; i < numObjs; i++ )
    {
        id = lids[i];

        if( ( id < 0 ) || ( id >= NumPoints ) )
        {
            *err = 1;
            return;
        }

        id3 = lids[i] * 3;

        pts[next++] = Points[id3];
        pts[next++] = Points[id3 + 1];
        pts[next++] = Points[id3 + 2];
    }
}

void mbGetNumberOfEdges( void* /* userDefinedData */,
                         int /* numGlobalIds */,
                         int /* numLids */,
                         int numObjs,
                         ZOLTAN_ID_PTR /* gids */,
                         ZOLTAN_ID_PTR lids,
                         int* numEdges,
                         int* err )
{
    int i, id;
    int next = 0;

    for( i = 0; i < numObjs; i++ )
    {
        id = lids[i];

        if( ( id < 0 ) || ( id >= NumPoints ) )
        {
            *err = 1;
            return;
        }

        numEdges[next++] = NumEdges[id];
    }
}

void mbGetEdgeList( void* /* userDefinedData */,
                    int /* numGlobalIds */,
                    int /* numLids */,
                    int numObjs,
                    ZOLTAN_ID_PTR /* gids */,
                    ZOLTAN_ID_PTR lids,
                    int* /* numEdges */,
                    ZOLTAN_ID_PTR nborGlobalIds,
                    int* nborProcs,
                    int wgt_dim,
                    float* edge_wgts,
                    int* err )
{
    int i, id, idSum, j;
    int next = 0;

    for( i = 0; i < numObjs; i++ )
    {
        id = lids[i];

        if( ( id < 0 ) || ( id >= NumPoints ) )
        {
            *err = 1;
            return;
        }

        idSum = 0;

        for( j = 0; j < id; j++ )
            idSum += NumEdges[j];

        for( j = 0; j < NumEdges[id]; j++ )
        {
            nborGlobalIds[next] = NborGlobalId[idSum];
            nborProcs[next]     = NborProcs[idSum];
            if( wgt_dim > 0 ) edge_wgts[next] = EdgeWeights[idSum];
            next++;
            idSum++;
        }
    }
}

void mbGetPart( void* /* userDefinedData */,
                int /* numGlobalIds */,
                int /* numLids */,
                int numObjs,
                ZOLTAN_ID_PTR /* gids */,
                ZOLTAN_ID_PTR lids,
                int* part,
                int* err )
{
    int i, id;
    int next = 0;

    for( i = 0; i < numObjs; i++ )
    {
        id = lids[i];

        if( ( id < 0 ) || ( id >= NumPoints ) )
        {
            *err = 1;
            return;
        }

        part[next++] = Parts[id];
    }
}

// new methods for partition in parallel, used by migrate in iMOAB
ErrorCode ZoltanPartitioner::partition_owned_cells( Range& primary,
                                                    std::multimap< int, int >& extraGraphEdges,
                                                    std::map< int, int > procs,
                                                    int& numNewPartitions,
                                                    std::map< int, Range >& distribution,
                                                    int met )
{
    // start copy
    MeshTopoUtil mtu( mbImpl );
    Range adjs;

    std::vector< int > adjacencies;
    std::vector< int > ids;
    ids.resize( primary.size() );
    std::vector< int > length;
    std::vector< double > coords;
    std::vector< int > nbor_proc;
    // can use a fixed-size array 'cuz the number of lower-dimensional neighbors is limited
    // by MBCN
    int neighbors[5 * MAX_SUB_ENTITIES];
    int neib_proc[5 * MAX_SUB_ENTITIES];
    double avg_position[3];
    int moab_id;
    int primaryDim = mbImpl->dimension_from_handle( *primary.rbegin() );
    // get the global id tag handle
    Tag gid = mbImpl->globalId_tag();

    ErrorCode rval = mbImpl->tag_get_data( gid, primary, &ids[0] );MB_CHK_ERR( rval );

    // mbpc is member in base class, PartitionerBase
    int rank = mbpc->rank();  // current rank , will be put on regular neighbors
    int i    = 0;
    for( Range::iterator rit = primary.begin(); rit != primary.end(); ++rit, i++ )
    {
        EntityHandle cell = *rit;
        // get bridge adjacencies for each cell
        if( 1 == met )
        {
            adjs.clear();
            rval = mtu.get_bridge_adjacencies( cell, ( primaryDim > 0 ? primaryDim - 1 : 3 ), primaryDim, adjs );MB_CHK_ERR( rval );

            // get the graph vertex ids of those
            if( !adjs.empty() )
            {
                assert( adjs.size() < 5 * MAX_SUB_ENTITIES );
                rval = mbImpl->tag_get_data( gid, adjs, neighbors );MB_CHK_ERR( rval );
            }
            // if adjacent to neighbor partitions, add to the list
            int size_adjs = (int)adjs.size();
            moab_id       = ids[i];
            for( int k = 0; k < size_adjs; k++ )
                neib_proc[k] = rank;  // current rank
            if( extraGraphEdges.find( moab_id ) != extraGraphEdges.end() )
            {
                // it means that the current cell is adjacent to a cell in another partition ; maybe
                // a few
                std::pair< std::multimap< int, int >::iterator, std::multimap< int, int >::iterator > ret;
                ret = extraGraphEdges.equal_range( moab_id );
                for( std::multimap< int, int >::iterator it = ret.first; it != ret.second; ++it )
                {
                    int otherID          = it->second;
                    neighbors[size_adjs] = otherID;         // the id of the other cell, across partition
                    neib_proc[size_adjs] = procs[otherID];  // this is how we built this map, the cell id maps to
                                                            // what proc it came from
                    size_adjs++;
                }
            }
            // copy those into adjacencies vector
            length.push_back( size_adjs );
            std::copy( neighbors, neighbors + size_adjs, std::back_inserter( adjacencies ) );
            std::copy( neib_proc, neib_proc + size_adjs, std::back_inserter( nbor_proc ) );
        }
        else if( 2 == met )
        {
            if( TYPE_FROM_HANDLE( cell ) == MBVERTEX )
            {
                rval = mbImpl->get_coords( &cell, 1, avg_position );MB_CHK_ERR( rval );
            }
            else
            {
                rval = mtu.get_average_position( cell, avg_position );MB_CHK_ERR( rval );
            }
            std::copy( avg_position, avg_position + 3, std::back_inserter( coords ) );
        }
    }

#ifdef VERBOSE
    std::stringstream ff2;
    ff2 << "zoltanInput_" << mbpc->rank() << ".txt";
    std::ofstream ofs;
    ofs.open( ff2.str().c_str(), std::ofstream::out );
    ofs << "Length vector: " << std::endl;
    std::copy( length.begin(), length.end(), std::ostream_iterator< int >( ofs, ", " ) );
    ofs << std::endl;
    ofs << "Adjacencies vector: " << std::endl;
    std::copy( adjacencies.begin(), adjacencies.end(), std::ostream_iterator< int >( ofs, ", " ) );
    ofs << std::endl;
    ofs << "Moab_ids vector: " << std::endl;
    std::copy( ids.begin(), ids.end(), std::ostream_iterator< int >( ofs, ", " ) );
    ofs << std::endl;
    ofs << "Coords vector: " << std::endl;
    std::copy( coords.begin(), coords.end(), std::ostream_iterator< double >( ofs, ", " ) );
    ofs << std::endl;
    ofs.close();
#endif

    // these are static var in this file, and used in the callbacks
    Points = NULL;
    if( 1 != met ) Points = &coords[0];
    GlobalIds    = &ids[0];
    NumPoints    = (int)ids.size();
    NumEdges     = &length[0];
    NborGlobalId = &adjacencies[0];
    NborProcs    = &nbor_proc[0];
    ObjWeights   = NULL;
    EdgeWeights  = NULL;
    Parts        = NULL;

    float version;
    if( mbpc->rank() == 0 ) std::cout << "Initializing zoltan..." << std::endl;

    Zoltan_Initialize( argcArg, argvArg, &version );

    // Create Zoltan object.  This calls Zoltan_Create.
    if( NULL == myZZ ) myZZ = new Zoltan( mbpc->comm() );

    // set # requested partitions
    char buff[10];
    snprintf( buff, 10, "%d", numNewPartitions );
    int retval = myZZ->Set_Param( "NUM_GLOBAL_PARTITIONS", buff );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    // request parts assignment
    retval = myZZ->Set_Param( "RETURN_LISTS", "PARTS" );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

    myZZ->Set_Num_Obj_Fn( mbGetNumberOfAssignedObjects, NULL );
    myZZ->Set_Obj_List_Fn( mbGetObjectList, NULL );
    // due to a bug in zoltan, if method is graph partitioning, do not pass coordinates!!
    if( 2 == met )
    {
        myZZ->Set_Num_Geom_Fn( mbGetObjectSize, NULL );
        myZZ->Set_Geom_Multi_Fn( mbGetObject, NULL );
        SetRCB_Parameters();  // geometry
    }
    else if( 1 == met )
    {
        myZZ->Set_Num_Edges_Multi_Fn( mbGetNumberOfEdges, NULL );
        myZZ->Set_Edge_List_Multi_Fn( mbGetEdgeList, NULL );
        SetHypergraph_Parameters( "auto" );
    }

    // Perform the load balancing partitioning

    int changes;
    int numGidEntries;
    int numLidEntries;
    int num_import;
    ZOLTAN_ID_PTR import_global_ids, import_local_ids;
    int* import_procs;
    int* import_to_part;
    int num_export;
    ZOLTAN_ID_PTR export_global_ids, export_local_ids;
    int *assign_procs, *assign_parts;

    if( mbpc->rank() == 0 )
        std::cout << "Computing partition using method (1-graph, 2-geom):" << met << " for " << numNewPartitions
                  << " parts..." << std::endl;

#ifndef NDEBUG
#if 0
  static int counter=0; // it may be possible to call function multiple times in a simulation
  // give a way to not overwrite the files
  // it should work only with a modified version of Zoltan
  std::stringstream basename;
  if (1==met)
  {
    basename << "phg_" << counter++;
    Zoltan_Generate_Files(myZZ->Get_C_Handle(), (char*)(basename.str().c_str()), 1, 0, 1, 0);
  }
  else if (2==met)
  {
    basename << "rcb_" << counter++;
    Zoltan_Generate_Files(myZZ->Get_C_Handle(), (char*)(basename.str().c_str()), 1, 1, 0, 0);
  }
#endif
#endif
    retval = myZZ->LB_Partition( changes, numGidEntries, numLidEntries, num_import, import_global_ids, import_local_ids,
                                 import_procs, import_to_part, num_export, export_global_ids, export_local_ids,
                                 assign_procs, assign_parts );
    if( ZOLTAN_OK != retval ) return MB_FAILURE;

#ifdef VERBOSE
    std::stringstream ff3;
    ff3 << "zoltanOutput_" << mbpc->rank() << ".txt";
    std::ofstream ofs3;
    ofs3.open( ff3.str().c_str(), std::ofstream::out );
    ofs3 << " export elements on rank " << rank << " \n";
    ofs3 << "\t index \t gb_id \t local \t proc \t part \n";
    for( int k = 0; k < num_export; k++ )
    {
        ofs3 << "\t" << k << "\t" << export_global_ids[k] << "\t" << export_local_ids[k] << "\t" << assign_procs[k]
             << "\t" << assign_parts[k] << "\n";
    }
    ofs3.close();
#endif

    // basically each local cell is assigned to a part

    assert( num_export == (int)primary.size() );
    for( i = 0; i < num_export; i++ )
    {
        EntityHandle cell = primary[export_local_ids[i]];
        distribution[assign_parts[i]].insert( cell );
    }

    Zoltan::LB_Free_Part( &import_global_ids, &import_local_ids, &import_procs, &import_to_part );
    Zoltan::LB_Free_Part( &export_global_ids, &export_local_ids, &assign_procs, &assign_parts );

    delete myZZ;
    myZZ = NULL;

    // clear arrays that were resized locally, to free up local memory

    std::vector< int >().swap( adjacencies );
    std::vector< int >().swap( ids );
    std::vector< int >().swap( length );
    std::vector< int >().swap( nbor_proc );
    std::vector< double >().swap( coords );

    return MB_SUCCESS;
}
