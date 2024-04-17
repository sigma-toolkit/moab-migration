#include "NCHelperDomain.hpp"
#include "moab/FileOptions.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "AEntityFactory.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelMergeMesh.hpp"
#endif
#ifdef MOAB_HAVE_ZOLTAN
#include "moab/ZoltanPartitioner.hpp"
#include "moab/TupleList.hpp"
#endif

#include <cmath>
#include <sstream>

namespace moab
{

bool NCHelperDomain::can_read_file( ReadNC* readNC, int fileId )
{
    std::vector< std::string >& dimNames = readNC->dimNames;

    // If dimension names "n" AND "ni" AND "nj" AND "nv" exist then it should be the Domain grid
    if(  //(std::find( dimNames.begin(), dimNames.end(), std::string( "n" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "ni" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "nj" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "nv" ) ) != dimNames.end() ) )
    {
        // Make sure it is CAM grid
        /* std::map< std::string, ReadNC::AttData >::iterator attIt = readNC->globalAtts.find( "source" );
        if( attIt == readNC->globalAtts.end() ) return false;
        unsigned int sz = attIt->second.attLen;
        std::string att_data;
        att_data.resize( sz + 1 );
        att_data[sz] = '\000';
        int success =
            NCFUNC( get_att_text )( fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0] );
        if( success ) return false;*/
        /*if (att_data.find("CAM") == std::string::npos)
          return false;*/

        return true;
    }

    return false;
}

ErrorCode NCHelperDomain::init_mesh_vals()
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
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "ni" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'ni' variable" );
    }
    iDim     = idx;
    gDims[0] = 0;
    gDims[3] = dimLens[idx];

    // Then j
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "nj" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'nj' variable" );
    }
    jDim     = idx;
    gDims[1] = 0;
    gDims[4] = dimLens[idx];  // Add 2 for the pole points ? not needed

    // do not use gcdims ? or use only gcdims?

    // Try a truly 2D mesh
    gDims[2] = -1;
    gDims[5] = -1;

    // Get number of vertices per cell
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "nv" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'nv' dimension" );
    }
    nvDim = idx;
    nv    = dimLens[idx];

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
            parData.gDims[i] = gDims[i];
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
            lDims[i] = gDims[i];
        locallyPeriodic[0] = globallyPeriodic[0];
    }

    // Now get actual coordinate values for vertices and cell centers
    lCDims[0] = lDims[0];

    lCDims[3] = lDims[3];

    // For FV models, will always be non-periodic in j
    lCDims[1] = lDims[1];
    lCDims[4] = lDims[4];

#if 0
  // Resize vectors to store values later
  if (-1 != lDims[0])
    ilVals.resize(lDims[3] - lDims[0] + 1);
  if (-1 != lCDims[0])
    ilCVals.resize(lCDims[3] - lCDims[0] + 1);
  if (-1 != lDims[1])
    jlVals.resize(lDims[4] - lDims[1] + 1);
  if (-1 != lCDims[1])
    jlCVals.resize(lCDims[4] - lCDims[1] + 1);
  if (nTimeSteps > 0)
    tVals.resize(nTimeSteps);
#endif

    dbgOut.tprintf( 1, "I=%d-%d, J=%d-%d\n", lDims[0], lDims[3], lDims[1], lDims[4] );
    dbgOut.tprintf( 1, "%d elements, %d vertices\n", ( lDims[3] - lDims[0] ) * ( lDims[4] - lDims[1] ),
                    ( lDims[3] - lDims[0] ) * ( lDims[4] - lDims[1] ) * nv );

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
    ijdimNames[0] = "__ni";
    ijdimNames[1] = "__nj";
    std::string tag_name;
    Tag tagh;

    // __<dim_name>_LOC_MINMAX (for slon, slat, lon and lat)
    for( unsigned int i = 0; i != ijdimNames.size(); i++ )
    {
        std::vector< int > val( 2, 0 );
        if( ijdimNames[i] == "__ni" )
        {
            val[0] = lDims[0];
            val[1] = lDims[3];
        }
        else if( ijdimNames[i] == "__nj" )
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
    switch( varInfo["xc"].varDataType )
    {
        case NC_FLOAT:
        case NC_DOUBLE:
            break;
        default:
            MB_SET_ERR( MB_FAILURE, "Unexpected variable data type for 'lon'" );
    }

    // do not need conventional tags
    Tag convTagsCreated = 0;
    int def_val         = 0;
    rval                = mbImpl->tag_get_handle( "__CONV_TAGS_CREATED", 1, MB_TYPE_INTEGER, convTagsCreated,
                                                  MB_TAG_SPARSE | MB_TAG_CREAT, &def_val );MB_CHK_SET_ERR( rval, "Trouble getting _CONV_TAGS_CREATED tag" );
    int create_conv_tags_flag = 1;
    rval                      = mbImpl->tag_set_data( convTagsCreated, &_fileSet, 1, &create_conv_tags_flag );MB_CHK_SET_ERR( rval, "Trouble setting _CONV_TAGS_CREATED tag" );

    return MB_SUCCESS;
}

ErrorCode NCHelperDomain::create_mesh( Range& faces )
{
    Interface*& mbImpl = _readNC->mbImpl;
    // std::string& fileName = _readNC->fileName;
    Tag& mGlobalIdTag = _readNC->mGlobalIdTag;
    // const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
    DebugOutput& dbgOut = _readNC->dbgOut;

    bool& culling = _readNC->culling;
    bool& repartition = _readNC->repartition;
    /*int& gatherSetRank = _readNC->gatherSetRank;
    int& trivialPartitionShift = _readNC->trivialPartitionShift;*/
    /*

      int rank = 0;
      int procs = 1;
    #ifdef MOAB_HAVE_MPI
      bool& isParallel = _readNC->isParallel;
      if (isParallel) {
        ParallelComm*& myPcomm = _readNC->myPcomm;
        rank = myPcomm->proc_config().proc_rank();
        procs = myPcomm->proc_config().proc_size();
      }
    #endif
    */

    ErrorCode rval;
    int success = 0;

    /*
      bool create_gathers = false;
      if (rank == gatherSetRank)
        create_gathers = true;

      // Shift rank to obtain a rotated trivial partition
      int shifted_rank = rank;
      if (procs >= 2 && trivialPartitionShift > 0)
        shifted_rank = (rank + trivialPartitionShift) % procs;*/

    // how many will have mask 0 or 1
    // how many will have a fraction  ? we will not instantiate all elements; only those with mask 1
    // ? also, not all vertices, only those that belong to mask 1 elements ? we will not care about
    // duplicate vertices; maybe another time ? we will start reading masks, vertices
    int local_elems = ( lDims[4] - lDims[1] ) * ( lDims[3] - lDims[0] );
    dbgOut.tprintf( 1, "local cells: %d \n", local_elems );

    // count how many will be with mask 1 here
    // basically, read the mask variable on the local elements;
    std::string maskstr( "mask" );
    ReadNC::VarData& vmask = _readNC->varInfo[maskstr];

    // mask is (nj, ni)
    vmask.readStarts.push_back( lDims[1] );
    vmask.readStarts.push_back( lDims[0] );
    vmask.readCounts.push_back( lDims[4] - lDims[1] );
    vmask.readCounts.push_back( lDims[3] - lDims[0] );
    std::vector< int > mask( local_elems );
    success = NCFUNCAG( _vara_int )( _fileId, vmask.varId, &vmask.readStarts[0], &vmask.readCounts[0], &mask[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read int data for mask variable " );

    std::vector<int> gids(local_elems);
    int elem_index = 0;
    int global_row_size = gDims[3] - gDims[0];  // this is along first dimension in global decomposition
    // create global id array for cells, for all cells, including those with 0 mask; which will be not used eventually
    for( int j = lCDims[1]; j < lCDims[4]; j++ )
        for( int i = lCDims[0]; i < lCDims[3]; i++ )
        {
            gids[elem_index] = j * global_row_size + i + 1;
            elem_index++;
        }
    std::vector< NCDF_SIZE > startsv( 3 );
    startsv[0] = vmask.readStarts[0];
    startsv[1] = vmask.readStarts[1];
    startsv[2] = 0;
    std::vector< NCDF_SIZE > countsv( 3 );
    countsv[0] = vmask.readCounts[0];
    countsv[1] = vmask.readCounts[1];
    countsv[2] = nv;  // number of vertices per element

    // read xv and yv coords for vertices, and create elements;
    std::string xvstr( "xv" );
    ReadNC::VarData& var_xv = _readNC->varInfo[xvstr];
    std::vector< double > xv( local_elems * nv );

    std::vector< NCDF_SIZE > starts_vv, counts_vv;
    starts_vv.resize( 3 );
    counts_vv.resize( 3 );
    bool nv_last = true;
    if( _readNC->dimNames[var_xv.varDims[0]] == std::string( "nv" ) )  // it means that nv is the first dimension
    {
        starts_vv[0] = 0;
        starts_vv[1] = startsv[0];
        starts_vv[2] = startsv[1];
        counts_vv[0] = nv;
        counts_vv[1] = countsv[0];
        counts_vv[2] = countsv[1];
        nv_last      = false;
    }
    else
    {
        starts_vv = startsv;
        counts_vv = countsv;
    }
    dbgOut.tprintf( 1, " nv is the last dimension in xv array (0 or 1)  %d \n", (int)nv_last );
    success = NCFUNCAG( _vara_double )( _fileId, var_xv.varId, &starts_vv[0], &counts_vv[0], &xv[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for xv variable " );

    std::string yvstr( "yv" );
    ReadNC::VarData& var_yv = _readNC->varInfo[yvstr];
    std::vector< double > yv( local_elems * nv );
    success = NCFUNCAG( _vara_double )( _fileId, var_yv.varId, &starts_vv[0], &counts_vv[0], &yv[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for yv variable " );

    // read other variables, like xc, yc, frac, area
    std::string xcstr( "xc" );
    ReadNC::VarData& var_xc = _readNC->varInfo[xcstr];
    std::vector< double > xc( local_elems );
    success = NCFUNCAG( _vara_double )( _fileId, var_xc.varId, &vmask.readStarts[0], &vmask.readCounts[0], &xc[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for xc variable " );

    std::string ycstr( "yc" );
    ReadNC::VarData& var_yc = _readNC->varInfo[ycstr];
    std::vector< double > yc( local_elems );
    success = NCFUNCAG( _vara_double )( _fileId, var_yc.varId, &vmask.readStarts[0], &vmask.readCounts[0], &yc[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for yc variable " );

    std::string fracstr( "frac" );
    ReadNC::VarData& var_frac = _readNC->varInfo[fracstr];
    std::vector< double > frac( local_elems );
    if( var_frac.varId >= 0 )
    {
        success =
            NCFUNCAG( _vara_double )( _fileId, var_frac.varId, &vmask.readStarts[0], &vmask.readCounts[0], &frac[0] );
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for frac variable " );
    }
    std::string areastr( "area" );
    std::vector< double > area( local_elems );
    ReadNC::VarData& var_area = _readNC->varInfo[areastr];
    success = NCFUNCAG( _vara_double )( _fileId, var_area.varId, &vmask.readStarts[0], &vmask.readCounts[0], &area[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read double data for area variable " );
    // create tags for them
    Tag areaTag, fracTag, xcTag, ycTag;
    rval = mbImpl->tag_get_handle( "area", 1, MB_TYPE_DOUBLE, areaTag, MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Trouble creating area tag" );
    rval = mbImpl->tag_get_handle( "frac", 1, MB_TYPE_DOUBLE, fracTag, MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Trouble creating frac tag" );
    rval = mbImpl->tag_get_handle( "xc", 1, MB_TYPE_DOUBLE, xcTag, MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Trouble creating xc tag" );
    rval = mbImpl->tag_get_handle( "yc", 1, MB_TYPE_DOUBLE, ycTag, MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_SET_ERR( rval, "Trouble creating yc tag" );

    // create tags for GRID_IMASK, which will be the same name as the Scrip helper tag that holds the mask
    // create the maskTag GRID_IMASK, with default value of 1
    Tag maskTag;
    int def_val = 1;
    rval = mbImpl->tag_get_handle( "GRID_IMASK", 1, MB_TYPE_INTEGER, maskTag, MB_TAG_DENSE | MB_TAG_CREAT, &def_val );MB_CHK_SET_ERR( rval, "Trouble creating GRID_IMASK tag" );

    // will now look to repartition the cells, using zoltan, looking at the xc and yc coordinates in 2d, convert to 3d,
    // and decide based on those partitioning info to what task to send each cell, along with its vertices, and global id
#ifdef MOAB_HAVE_MPI
    int rank              = 0;
    int procs             = 1;
    bool& isParallel      = _readNC->isParallel;
    ParallelComm* myPcomm = NULL;
    if( isParallel )
    {
        myPcomm = _readNC->myPcomm;
        rank    = myPcomm->proc_config().proc_rank();
        procs   = myPcomm->proc_config().proc_size();
    }

    if( procs >= 2 && repartition )
    {
        // Redistribute local cells after trivial partition (e.g. apply Zoltan partition)
        rval = redistribute_cells( myPcomm, xc, yc, xv, yv, frac, mask, area, gids, nv, nv_last );MB_CHK_SET_ERR( rval, "Failed to redistribute local cells" );
        local_elems = (int) xc.size();
        dbgOut.tprintf( 1, "local cells after repartition: %d \n", local_elems );
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
    EntityType mdb_type = MBVERTEX;
    if( nv == 3 )
        mdb_type = MBTRI;
    else if( nv == 4 )
        mdb_type = MBQUAD;
    else if( nv > 4 )  // (nv > 4)
        mdb_type = MBPOLYGON;
    // for nv = 1 , type is vertex

    int num_actual_cells = local_elems;
    if( culling ) num_actual_cells = nb_with_mask1;

    if( nv > 1 && num_actual_cells > 0 )
    {
        rval = _readNC->readMeshIface->get_element_connect( num_actual_cells, nv, mdb_type, 0, start_cell, conn_arr );MB_CHK_SET_ERR( rval, "Failed to create local cells" );
        tmp_range.insert( start_cell, start_cell + num_actual_cells - 1 );
    }

    // Create vertices; first identify different ones, with a tolerance
    std::map< Node3D, EntityHandle > vertex_map;

    if( num_actual_cells > 0 )
    {
        // Set vertex coordinates
        // will read all xv, yv, but use only those with correct mask on

          // total index in netcdf arrays
        const double pideg = acos( -1.0 ) / 180.0;

        for(elem_index = 0; elem_index < local_elems; elem_index++ )
        {
            if( culling && 0 == mask[elem_index] )
                continue;  // nothing to do, do not advance elem_index in actual moab arrays
            // set area and fraction on those elements too
            dbgOut.tprintf( 3, "elem index  %d \n", elem_index );
            for( int k = 0; k < nv; k++ )
            {
                int index_v_arr = nv * elem_index + k;
                if( !nv_last ) index_v_arr = k * local_elems + elem_index;
                double x, y;
                if( nv > 1 )
                {
                    x             = xv[index_v_arr];
                    y             = yv[index_v_arr];
                    double cosphi = cos( pideg * y );
                    double zmult  = sin( pideg * y );
                    double xmult  = cosphi * cos( x * pideg );
                    double ymult  = cosphi * sin( x * pideg );
                    Node3D pt( xmult, ymult, zmult );
                    vertex_map[pt] = 0;
                    dbgOut.tprintf( 3, "   %d  x=%5.1f y=%5.1f  pt: %f %f %f \n", k, x, y, pt.coords[0], pt.coords[1],
                                    pt.coords[2] );
                }
                else
                {
                    x = xc[elem_index];
                    y = yc[elem_index];
                    Node3D pt( x, y, 0 );
                    vertex_map[pt] = 0;
                }
            }
        }
        int nLocalVertices = (int)vertex_map.size();
        std::vector< double* > arrays;
        EntityHandle start_vertex;
        rval = _readNC->readMeshIface->get_node_coords( 3, nLocalVertices, 0, start_vertex, arrays );MB_CHK_SET_ERR( rval, "Failed to create local vertices" );

        vtx_handle = start_vertex;
        // Copy vertex coordinates into entity sequence coordinate arrays
        // and copy handle into vertex_map.
        double *x = arrays[0], *y = arrays[1], *z = arrays[2];
        for( auto i = vertex_map.begin(); i != vertex_map.end(); ++i )
        {
            i->second = vtx_handle;
            ++vtx_handle;
            *x = i->first.coords[0];
            ++x;
            *y = i->first.coords[1];
            ++y;
            *z = i->first.coords[2];
            ++z;
        }

        // int nj = gDims[4]-gDims[1]; // is it about 1 in irregular cases

        // int local_row_size  = lCDims[3] - lCDims[0];
        int index           = 0;  // consider the mask for advancing in moab arrays;

        // create now vertex arrays, size vertex_map.size()
        for( elem_index = 0; elem_index < local_elems; elem_index++)
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
    myPcomm = _readNC->myPcomm; // we will have to set the global id on vertices in any case, even
    if( myPcomm && procs >= 2)
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

    return MB_SUCCESS;
}

#ifdef MOAB_HAVE_MPI
ErrorCode NCHelperDomain::redistribute_cells( ParallelComm * myPcomm,
                                  std::vector<double> & xc, // center x
                                  std::vector<double> & yc, // center y
                                  std::vector<double> & xv, // vertex coords
                                  std::vector<double> & yv,
                                  std::vector<double> & frac, // fractions
                                  std::vector<int> & masks,    // mask
                                  std::vector<double> & area, // area
                                  std::vector<int> & gids,    // global ids
                                  int nv,                     // number of vertices per cell
                                  bool nv_last)               // type of xv, yv, first or last
{

#ifdef MOAB_HAVE_ZOLTAN
    // use zoltan and
    int num_local_cells=(int)gids.size();
    std::vector<double> xi(num_local_cells), yi(num_local_cells), zi(num_local_cells);
    const double pideg = acos( -1.0 ) / 180.0;
    for (size_t i=0; i<xc.size(); i++)
    {
        double x = xc[i];
        double y = yc[i];
        double cosphi = cos( pideg * y );
        double zmult  = sin( pideg * y );
        double xmult  = cosphi * cos( x * pideg );
        double ymult  = cosphi * sin( x * pideg );
        xi[i] = xmult;
        yi[i] = ymult;
        zi[i] = zmult;
    }
    Interface*& mbImpl         = _readNC->mbImpl;
    ZoltanPartitioner* mbZTool = new ZoltanPartitioner( mbImpl, myPcomm, false, 0, NULL );
    int start_cell_idx = gids[0];
    std::vector<int> dest(num_local_cells);
    ErrorCode rval             = mbZTool->repartition_to_procs( xi, yi, zi, gids, "RCB", dest );MB_CHK_SET_ERR( rval, "Error in Zoltan partitioning" );
    delete mbZTool;
    // now use crystal router to send the arrays to the right places
    moab::TupleList tl;
    unsigned numr = 2 * nv + 4;         //  doubles: area, centerlon, centerlat, frac, xv, yv,
    tl.initialize( 3, 0, 0, numr, num_local_cells );  // to proc, dof, mask
    tl.enableWriteAccess();
    // populate
    for( unsigned i = 0; i < num_local_cells; i++ )
    {
        int gdof    = gids[i];
        int to_proc = dest[i];
        int mask    = masks[i];
        int n                  = tl.get_n();
        tl.vi_wr[3 * n]        = to_proc;
        tl.vi_wr[3 * n + 1]    = gdof;
        tl.vi_wr[3 * n + 2]    = mask;
        tl.vr_wr[n * numr]     = area[i];
        tl.vr_wr[n * numr + 1] = xc[i];
        tl.vr_wr[n * numr + 2] = yc[i];
        tl.vr_wr[n * numr + 3] = frac[i];
        for( int k = 0; k < nv; k++ )
        {
            int index_v_arr = nv * i + k;
            if( !nv_last ) index_v_arr = k * num_local_cells + n;
            tl.vr_wr[ n * numr + 4 + k ]       = xv[index_v_arr];
            tl.vr_wr[ n * numr + 4 + nv + k]   = yv[index_v_arr];
        }
        tl.inc_n();
    }

    // now do the heavy communication
    ( myPcomm->proc_config().crystal_router() )->gs_transfer( 1, tl, 0 );

    // after communication, on each processor we should have tuples coming in
    // rearrange the vectors by global id
    moab::TupleList::buffer sort_buffer;
    int N = tl.get_n();
    sort_buffer.buffer_init( N );
    tl.sort( 1, &sort_buffer ); // 1 is the index for global id
    xc.resize(N);
    yc.resize(N);
    xv.resize(N*nv);
    yv.resize(N*nv);
    frac.resize(N);
    masks.resize(N);
    area.resize(N);
    gids.resize(N);
    for (int n=0; n<N; n++)
    {

        gids[n]  = tl.vi_wr[3 * n + 1];
        masks[n] = tl.vi_wr[3 * n + 2];
        area[n]  = tl.vr_wr[n * numr];
        xc[n]    = tl.vr_wr[n * numr + 1];
        yc[n]    = tl.vr_wr[n * numr + 2];
        frac[n]  = tl.vr_wr[n * numr + 3] ;
        for( int k = 0; k < nv; k++ )
        {
            int index_v_arr = nv * n + k;
            if( !nv_last ) index_v_arr = k * N + n;
            xv[index_v_arr] = tl.vr_wr[ n * numr + 4 + k ];
            yv[index_v_arr] = tl.vr_wr[ n * numr + 4 + nv + k];
        }
    }

    return MB_SUCCESS;
#else
    return MB_FAILURE;
#endif

}

#endif


}  // namespace moab
