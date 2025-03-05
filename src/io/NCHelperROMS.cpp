/*
 * NCHelperROMS.cpp
 *
 *  Created on: Mar. 4, 2025
 *      Author: iulian
 */

#include "NCHelperROMS.hpp"
#include "moab/FileOptions.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "AEntityFactory.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelMergeMesh.hpp"
#endif

#include <cmath>
#include <sstream>

namespace moab
{
bool NCHelperROMS::can_read_file( ReadNC* readNC, int /*fileId*/ )
{
    std::vector< std::string >& dimNames = readNC->dimNames;

    // If dimension names "n" AND "ni" AND "nj" AND "nv" exist then it should be the Domain grid
    // some files do not have n, which should be just ni*nj
    if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "xi_vert" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "eta_vert" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "xi_rho" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "eta_rho" ) ) != dimNames.end() ) )
    {
        // we do not test anymore if it is an actual CAM file
        return true;
    }

    return false;
}


ErrorCode NCHelperROMS::init_mesh_vals()
{
    Interface*& mbImpl                                = _readNC->mbImpl;
    std::vector< std::string >& dimNames              = _readNC->dimNames;
    std::vector< int >& dimLens                       = _readNC->dimLens;
    std::map< std::string, ReadNC::VarData >& varInfo = _readNC->varInfo;
    DebugOutput& dbgOut                               = _readNC->dbgOut;
#ifdef MOAB_HAVE_MPI
    bool& isParallel = _readNC->isParallel;
#endif
    int& partMethod     = _readNC->partMethod;
    ScdParData& parData = _readNC->parData;

    ErrorCode rval;

    // Look for names of i/j dimensions
    // First i
    std::vector< std::string >::iterator vit;
    unsigned int idx;
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "xi_vert" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        //MB_SET_ERR( MB_FAILURE, "Couldn't find 'xi_vert' variable" );
    }
    iDim     = idx;
    gDims[0] = 0;
    gDims[3] = dimLens[idx]-1;

    // Then j
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "eta_vert" ) ) != dimNames.end() )
    {
        idx = vit - dimNames.begin();
    }

    else
    {
        //MB_SET_ERR( MB_FAILURE, "Couldn't find 'eta_vert' variable" );
    }
    jDim     = idx;
    gDims[1] = 0;
    gDims[4] = dimLens[idx]-1;  // Add 2 for the pole points ? not needed

    if (iDim >= 0 and jDim >= 0)
        vertices_exist = true;
    // do not use gcdims ? or use only gcdims?

    // Try a truly 2D mesh
    gDims[2] = -1;
    gDims[5] = -1;

    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "xi_rho" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        //MB_SET_ERR( MB_FAILURE, "Couldn't find 'xi_vert' variable" );
    }
    iCDim     = idx;
    gCDims[0] = 0;
    gCDims[3] = dimLens[idx]-1;

    // Then j
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "eta_rho" ) ) != dimNames.end() )
    {
        idx = vit - dimNames.begin();
    }

    else
    {
        //MB_SET_ERR( MB_FAILURE, "Couldn't find 'eta_vert' variable" );
    }
    jCDim     = idx;
    gCDims[1] = 0;
    gCDims[4] = dimLens[idx]-1;  // Add 2 for the pole points ? not needed

    // Parse options to get subset
    int rank = 0, procs = 1;
#ifdef MOAB_HAVE_MPI
    if( isParallel )
    {
        ParallelComm*& myPcomm = _readNC->myPcomm;
        rank                   = myPcomm->proc_config().proc_rank();
        procs                  = myPcomm->proc_config().proc_size();
    }
#endif
    if( procs > 1 )
    {
        for( int i = 0; i < 6; i++ )
            parData.gDims[i] = gCDims[i];
        parData.partMethod = partMethod;
        int pdims[3];

        locallyPeriodic[0] = locallyPeriodic[1] = locallyPeriodic[2] = 0;
        rval = ScdInterface::compute_partition( procs, rank, parData, lDims, locallyPeriodic, pdims );MB_CHK_ERR( rval );
        for( int i = 0; i < 3; i++ )
            parData.pDims[i] = pdims[i];

        dbgOut.tprintf( 1, "Partition: %dx%d (out of %dx%d)\n", lDims[3] - lDims[0], lDims[4] - lDims[1],
                        gDims[3] - gDims[0], gDims[4] - gDims[1] );
        if( 0 == rank )
            dbgOut.tprintf( 1, "Contiguous chunks of size %d bytes.\n",
                            8 * ( lDims[3] - lDims[0] ) * ( lDims[4] - lDims[1] ) );
    }
    else
    {
        for( int i = 0; i < 6; i++ )
            lDims[i] = gCDims[i];
        locallyPeriodic[0] = globallyPeriodic[0];
    }

    // Now get actual coordinate values for vertices and cell centers
    lCDims[0] = lDims[0];

    lCDims[3] = lDims[3] ; //

    // will always be non-periodic
    lCDims[1] = lDims[1];
    lCDims[4] = lDims[4] ;


    dbgOut.tprintf( 1, "I=%d-%d, J=%d-%d\n", lDims[0], lDims[3], lDims[1], lDims[4] );
    dbgOut.tprintf( 1, "%d elements, %d vertices\n", ( lCDims[3] - lCDims[0] ) * ( lCDims[4] - lCDims[1] ),
                    ( lDims[3] - lDims[0] ) * ( lDims[4] - lDims[1] ) );

    // For each variable, determine the entity location type and number of levels
    std::map< std::string, ReadNC::VarData >::iterator mit;
    for( mit = varInfo.begin(); mit != varInfo.end(); ++mit )
    {
        ReadNC::VarData& vd = ( *mit ).second;

        // Default entLoc is ENTLOCSET
        if( std::find( vd.varDims.begin(), vd.varDims.end(), tDim ) != vd.varDims.end() )
        {
            if( ( std::find( vd.varDims.begin(), vd.varDims.end(), iCDim ) != vd.varDims.end() ) &&
                ( std::find( vd.varDims.begin(), vd.varDims.end(), jCDim ) != vd.varDims.end() ) )
                vd.entLoc = ReadNC::ENTLOCFACE;
            else if( ( std::find( vd.varDims.begin(), vd.varDims.end(), jDim ) != vd.varDims.end() ) &&
                     ( std::find( vd.varDims.begin(), vd.varDims.end(), iCDim ) != vd.varDims.end() ) )
                vd.entLoc = ReadNC::ENTLOCNSEDGE;
            else if( ( std::find( vd.varDims.begin(), vd.varDims.end(), jCDim ) != vd.varDims.end() ) &&
                     ( std::find( vd.varDims.begin(), vd.varDims.end(), iDim ) != vd.varDims.end() ) )
                vd.entLoc = ReadNC::ENTLOCEWEDGE;
        }

        // Default numLev is 0
        if( std::find( vd.varDims.begin(), vd.varDims.end(), levDim ) != vd.varDims.end() ) vd.numLev = nLevels;
    }

    std::vector< std::string > ijdimNames( 2 );

    ijdimNames[0] = "__xi_rho";
    ijdimNames[1] = "__eta_rho";
    std::string tag_name;
    Tag tagh;

    // __<dim_name>_LOC_MINMAX (for slon, slat, lon and lat)
    for( unsigned int i = 0; i != ijdimNames.size(); i++ )
    {
        std::vector< int > val( 2, 0 );
        if( ijdimNames[i] == "__xi_rho" )
        {
            val[0] = lDims[0];
            val[1] = lDims[3];
        }
        else if( ijdimNames[i] == "__eta_rho" )
        {
            val[0] = lDims[1];
            val[1] = lDims[4];
        }

        std::stringstream ss_tag_name;
        ss_tag_name << ijdimNames[i] << "_LOC_MINMAX";
        tag_name = ss_tag_name.str();
        rval     = mbImpl->tag_get_handle( tag_name.c_str(), 2, MB_TYPE_INTEGER, tagh, MB_TAG_SPARSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Trouble creating conventional tag " << tag_name );
        rval = mbImpl->tag_set_data( tagh, &_fileSet, 1, &val[0] );MB_CHK_SET_ERR( rval, "Trouble setting data to conventional tag " << tag_name );
        if( MB_SUCCESS == rval ) dbgOut.tprintf( 2, "Conventional tag %s is created.\n", tag_name.c_str() );
    }

    // __<dim_name>_LOC_VALS (for slon, slat, lon and lat)
    // Assume all have the same data type as lon (expected type is float or double)

    // do not need conventional tags
    Tag convTagsCreated = 0;
    int def_val         = 0;
    rval                = mbImpl->tag_get_handle( "__CONV_TAGS_CREATED", 1, MB_TYPE_INTEGER, convTagsCreated,
                                                  MB_TAG_SPARSE | MB_TAG_CREAT, &def_val );MB_CHK_SET_ERR( rval, "Trouble getting _CONV_TAGS_CREATED tag" );
    int create_conv_tags_flag = 1;
    rval                      = mbImpl->tag_set_data( convTagsCreated, &_fileSet, 1, &create_conv_tags_flag );MB_CHK_SET_ERR( rval, "Trouble setting _CONV_TAGS_CREATED tag" );

    return MB_SUCCESS;
}

ErrorCode NCHelperROMS::create_mesh( Range& faces )
{
    Interface*& mbImpl = _readNC->mbImpl;
    // std::string& fileName = _readNC->fileName;
    Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
    // const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
    DebugOutput& dbgOut = _readNC->dbgOut;


    ErrorCode rval;
    int success = 0;

    int local_elems = ( lCDims[4] - lCDims[1] ) * ( lCDims[3] - lCDims[0] );
    dbgOut.tprintf( 1, "local cells: %d \n", local_elems );
    int local_vertices = ( lDims[4] - lDims[1] +1) * ( lDims[3] - lDims[0] +1);
    dbgOut.tprintf( 1, "local vertices: %d \n", local_vertices );


    // count how many will be with mask 1 here
    // basically, read the mask variable on the local elements;
    std::string maskstr( "mask_rho" );
    ReadNC::VarData& vmask = _readNC->varInfo[maskstr];

    // mask is (nj, ni)
    vmask.readStarts.push_back( lCDims[1] );
    vmask.readStarts.push_back( lCDims[0] );
    vmask.readCounts.push_back( lCDims[4] - lCDims[1] );
    vmask.readCounts.push_back( lCDims[3] - lCDims[0] );
    std::vector< double > mask( local_elems );
    success = NCFUNCAG( _vara_double )( _fileId, vmask.varId, &vmask.readStarts[0], &vmask.readCounts[0], &mask[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read int data for mask_rho variable " );

    std::vector< int > gids( local_elems );
    int elem_index      = 0;
    int global_row_size = gCDims[3] - gCDims[0]+1;  // this is along first dimension in global decomposition
    // create global id array for cells, for all cells, including those with 0 mask; which will be not used eventually
    for( int j = lCDims[1]; j < lCDims[4]; j++ )
        for( int i = lCDims[0]; i < lCDims[3]; i++ )
        {
            gids[elem_index] = j * global_row_size + i + 1;
            elem_index++;
        }
    std::vector< int > gidv( local_vertices );
    int vertex_index      = 0;
    int global_row_size_vert = gCDims[3] - gCDims[0]+2;  // this is along first dimension in global decomposition
    // create global id array for vertices
    for( int j = lDims[1]; j < lDims[4]; j++ )
        for( int i = lDims[0]; i < lDims[3]; i++ )
        {
            gidv[vertex_index] = j * global_row_size_vert + i + 1;
            vertex_index++;
        }

    std::vector< NCDF_SIZE > startsv( 3 );
    startsv[0] = vmask.readStarts[0];
    startsv[1] = vmask.readStarts[1];
    startsv[2] = 0;
    std::vector< NCDF_SIZE > countsv( 3 );
    countsv[0] = vmask.readCounts[0]+1;
    countsv[1] = vmask.readCounts[1]+1; // more for vertices
    countsv[2] = 0;  // number of vertices per element

    // read xv and yv coords for vertices, and create elements;
    std::string xvstr( "lon_vert" );
    ReadNC::VarData& var_xv = _readNC->varInfo[xvstr];
    std::vector< double > lonv( local_vertices );

    std::string yvstr( "lat_vert" );
    ReadNC::VarData& var_yv = _readNC->varInfo[yvstr];
    std::vector< double > latv( local_vertices );

    success = NCFUNCAG( _vara_double )( _fileId, var_xv.varId, &startsv[0], &countsv[0], &lonv[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for lon_vert variable " );

    success = NCFUNCAG( _vara_double )( _fileId, var_yv.varId, &startsv[0], &countsv[0], &latv[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for lat_vert variable " );

    // create vertices

    std::vector< double* > arrays;
    EntityHandle start_vertex;
    rval = _readNC->readMeshIface->get_node_coords( 3, local_vertices, 0, start_vertex, arrays );MB_CHK_SET_ERR( rval, "Failed to create local vertices" );

    IntxUtils::SphereCoords sph;
    sph.R = 1;

    for( int i = 0; i < local_vertices ; ++i )
    {
        sph.lon = lonv[i]/180.*M_PI;
        sph.lat = latv[i]/180.*M_PI;
        CartVect pos = IntxUtils::spherical_to_cart(sph);
        arrays[0][i]=pos[0];
        arrays[1][i]=pos[1];
        arrays[2][i]=pos[2];
    }

    // Add local vertices to current file set
    Range local_verts_range( start_vertex, start_vertex + local_vertices - 1 );
    rval = _readNC->mbImpl->add_entities( _fileSet, local_verts_range );MB_CHK_SET_ERR( rval, "Failed to add local vertices to current file set" );

    rval = _readNC->mbImpl->tag_set_data(mGlobalIdTag, local_verts_range, &gidv[0]); MB_CHK_SET_ERR( rval, "Failed to set vertices global id" );

    int procs             = 1;
    int rank = 0;

#ifdef MOAB_HAVE_MPI

    bool& isParallel      = _readNC->isParallel;
    ParallelComm* myPcomm = NULL;
    if( isParallel )
    {
        myPcomm = _readNC->myPcomm;
        procs   = myPcomm->proc_config().proc_size();
        rank = myPcomm->proc_config().proc_rank();
    }

#endif

    int nb_with_mask1 = 0;
    for( int i = 0; i < local_elems; i++ )
        if( 1 == mask[i] ) nb_with_mask1++;
    dbgOut.tprintf( 1, "local cells with mask 1: %d \n", nb_with_mask1 );

    EntityHandle* conn_arr;
    EntityHandle vtx_handle;
    Range tmp_range;

    // set connectivity into that space

    EntityHandle start_cell;
    EntityType mdb_type = MBQUAD;

    int num_actual_cells = local_elems;

    if( num_actual_cells > 0 )
    {
        rval = _readNC->readMeshIface->get_element_connect( num_actual_cells, 4, mdb_type, 0, start_cell, conn_arr );MB_CHK_SET_ERR( rval, "Failed to create local cells" );
        tmp_range.insert( start_cell, start_cell + num_actual_cells - 1 );
    }


    // quad cells , explicit connectivity
    int row_cell = lCDims[3] - lCDims[0];
    int row_vert = row_cell + 1;
    for (int j=lCDims[1]; j< lCDims[4]; j++)
    {
        int j1 = j - lCDims[1];
        for (int i = lCDims[0]; i< lCDims[3]; i++)
        {
            int i1 = i - lCDims[0];
            int index_cells = j1* row_cell + i1;
            conn_arr[index_cells*4 ] = start_vertex + j1 * row_vert + i1;
            conn_arr[index_cells*4 + 1] = start_vertex + j1 * row_vert + i1 + 1;
            conn_arr[index_cells*4 + 2] = start_vertex + ( j1 +1 ) * row_vert + i1 + 1;
            conn_arr[index_cells*4 + 3] = start_vertex + ( j1 +1 ) * row_vert + i1 ;
        }
    }

    // Add local cells to current file set
    Range local_cells( start_cell, start_cell + local_elems - 1 );
    rval = _readNC->mbImpl->add_entities( _fileSet, local_cells );MB_CHK_SET_ERR( rval, "Failed to add local cells to current file set" );

    rval = _readNC->mbImpl->tag_set_data(mGlobalIdTag, local_cells, &gids[0]); MB_CHK_SET_ERR( rval, "Failed to set vertices global id" );

    Tag mask_tag;
    rval = _readNC->mbImpl->tag_get_handle("maskd", 1, MB_TYPE_DOUBLE, mask_tag,
            MB_TAG_DENSE | MB_TAG_EXCL | MB_TAG_CREAT ); MB_CHK_SET_ERR( rval, "Failed to create mask tag" );
    rval = _readNC->mbImpl->tag_set_data(mask_tag, local_cells, &mask[0]); MB_CHK_SET_ERR( rval, "Failed to set cells masks " );
    std::stringstream local_file_name;
    local_file_name << "roms_task." << procs <<"."<< rank <<".h5m";
    rval = _readNC->mbImpl->write_file(local_file_name.str().c_str(), NULL, "", &_fileSet, 1); MB_CHK_SET_ERR( rval, "Failed to write local file" );



#ifdef MOAB_HAVE_MPI
    if( isParallel )
    {
       rval   = _readNC->myPcomm->resolve_shared_ents( _fileSet, -1, -1, &mGlobalIdTag );MB_CHK_ERR( rval );
    }

#endif

    // set some global id on vertices, and also some cell global id
#if 0

        // create now vertex arrays, size vertex_map.size()
        for( elem_index = 0; elem_index < local_elems; elem_index++ )
        {
            if( culling && 0 == mask[elem_index] )
                continue;  // nothing to do, do not advance elem_index in actual moab arrays
            // set area and fraction on those elements too
            for( int k = 0; k < nv; k++ )
            {
                int index_v_arr = nv * elem_index + k;
                if( !nv_last ) index_v_arr = k * local_elems + elem_index;
                if( nv > 1 )
                {
                    double x      = xv[index_v_arr];
                    double y      = yv[index_v_arr];
                    double cosphi = cos( pideg * y );
                    double zmult  = sin( pideg * y );
                    double xmult  = cosphi * cos( x * pideg );
                    double ymult  = cosphi * sin( x * pideg );
                    Node3D pt( xmult, ymult, zmult );
                    conn_arr[index * nv + k] = vertex_map[pt];
                }
            }
            EntityHandle cell = start_vertex + index;
            if( nv > 1 ) cell = start_cell + index;
            // set other tags, like xc, yc, frac, area
            rval = mbImpl->tag_set_data( xcTag, &cell, 1, &xc[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set xc tag" );
            rval = mbImpl->tag_set_data( ycTag, &cell, 1, &yc[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set yc tag" );
            rval = mbImpl->tag_set_data( areaTag, &cell, 1, &area[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set area tag" );
            rval = mbImpl->tag_set_data( fracTag, &cell, 1, &frac[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set frac tag" );
            rval = mbImpl->tag_set_data( maskTag, &cell, 1, &mask[elem_index] );MB_CHK_SET_ERR( rval, "Failed to set mask tag" );

            // set the global id too:
            int globalId = gids[elem_index];
            rval         = mbImpl->tag_set_data( mGlobalIdTag, &cell, 1, &globalId );MB_CHK_SET_ERR( rval, "Failed to set global id tag" );
            index++;
        }

        rval = mbImpl->add_entities( _fileSet, tmp_range );MB_CHK_SET_ERR( rval, "Failed to add new cells to current file set" );

        // modify local file set, to merge coincident vertices, and to correct repeated vertices in elements
        std::vector< Tag > tagList;
        tagList.push_back( mGlobalIdTag );
        tagList.push_back( xcTag );
        tagList.push_back( ycTag );
        tagList.push_back( areaTag );
        tagList.push_back( fracTag );
        tagList.push_back( maskTag );  // not sure this is needed though ? on cells or on vertices?
        rval = IntxUtils::remove_padded_vertices( mbImpl, _fileSet, tagList );MB_CHK_SET_ERR( rval, "Failed to remove duplicate vertices" );

        rval = mbImpl->get_entities_by_dimension( _fileSet, 2, faces );MB_CHK_ERR( rval );
        Range all_verts;
        rval = mbImpl->get_connectivity( faces, all_verts );MB_CHK_ERR( rval );
        //printf(" range vert size :%ld \n", all_verts.size());
        rval = mbImpl->add_entities( _fileSet, all_verts );MB_CHK_ERR( rval );

        // need to add adjacencies; TODO: fix this for all nc readers
        // copy this logic from migrate mesh in par comm graph
        Core* mb                 = (Core*)mbImpl;
        AEntityFactory* adj_fact = mb->a_entity_factory();
        if( !adj_fact->vert_elem_adjacencies() )
            adj_fact->create_vert_elem_adjacencies();
        else
        {
            for( Range::iterator it = faces.begin(); it != faces.end(); ++it )
            {
                EntityHandle eh          = *it;
                const EntityHandle* conn = NULL;
                int num_nodes            = 0;
                rval                     = mb->get_connectivity( eh, conn, num_nodes );MB_CHK_ERR( rval );
                adj_fact->notify_create_entity( eh, conn, num_nodes );
            }
        }
    }

#ifdef MOAB_HAVE_MPI
    myPcomm = _readNC->myPcomm;  // we will have to set the global id on vertices in any case, even
    if( myPcomm && procs >= 2 )
    {
        double tol = 1.e-12;  // this is the same as static tolerance in NCHelper
        ParallelMergeMesh pmm( myPcomm, tol );
        rval = pmm.merge( _fileSet,
                          /* do not do local merge*/ false,
                          /*  2d cells*/ 2 );MB_CHK_SET_ERR( rval, "Failed to merge vertices in parallel" );
        // assign global ids only for vertices, cells have them fine
        rval = myPcomm->assign_global_ids( _fileSet, /*dim*/ 0 );MB_CHK_ERR( rval );
    }
#endif

#endif
    return MB_SUCCESS;
}
}
