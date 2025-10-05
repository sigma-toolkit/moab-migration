/*
 * NCHelperESMF.cpp
 *
 *  Created on: Sep. 12, 2023
 */

#include "NCHelperESMF.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"

#ifdef MOAB_HAVE_ZOLTAN
#include "moab/ZoltanPartitioner.hpp"
#endif

#include <cmath>

namespace moab
{

const int DEFAULT_MAX_EDGES_PER_CELL = 10;
const double pideg                   = acos( -1.0 ) / 180.0;

NCHelperESMF::NCHelperESMF( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
    : UcdNCHelper( readNC, fileId, opts, fileSet ), maxEdgesPerCell( DEFAULT_MAX_EDGES_PER_CELL ), coordDim( 0 ),
      centerCoordsId( -1 ), degrees( true )
{
}

bool NCHelperESMF::can_read_file( ReadNC* readNC )
{
    std::vector< std::string >& dimNames = readNC->dimNames;
    if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "nodeCount" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "elementCount" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "maxNodePElement" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "coordDim" ) ) != dimNames.end() ) )
    {
        return true;
    }

    return false;
}

ErrorCode NCHelperESMF::init_mesh_vals()
{
    std::vector< std::string >& dimNames              = _readNC->dimNames;
    std::vector< int >& dimLens                       = _readNC->dimLens;
    std::map< std::string, ReadNC::VarData >& varInfo = _readNC->varInfo;

    ErrorCode rval;
    unsigned int idx;
    std::vector< std::string >::iterator vit;

    // Get max edges per cell reported in the ESMF file header
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "maxNodePElement" ) ) != dimNames.end() )
    {
        idx             = vit - dimNames.begin();
        maxEdgesPerCell = dimLens[idx];
        if( maxEdgesPerCell > DEFAULT_MAX_EDGES_PER_CELL )
        {
            MB_SET_ERR( MB_INVALID_SIZE,
                        "maxEdgesPerCell read from the ESMF file header has exceeded " << DEFAULT_MAX_EDGES_PER_CELL );
        }
    }

    // Get number of cells
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "elementCount" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'elementCount' dimension" );
    }
    cDim   = idx;
    nCells = dimLens[idx];

    // Get number of vertices per cell (max)
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "maxNodePElement" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'maxNodePElement' dimension" );
    }

    maxEdgesPerCell = dimLens[idx];

    // Get number of vertices
    vDim = -1;
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "nodeCount" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'nodeCount' dimension" );
    }
    vDim      = idx;
    nVertices = dimLens[idx];

    // Get number of coord dimensions
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "coordDim" ) ) != dimNames.end() )
        idx = vit - dimNames.begin();
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find 'coordDim' dimension" );
    }

    coordDim = dimLens[idx];

    int success = NCFUNC( inq_varid )( _fileId, "centerCoords", &centerCoordsId );
    if( success ) centerCoordsId = -1;  // no center coords variable

    // decide now the units, by looking at nodeCoords; they should always exist
    int nodeCoordsId;
    success = NCFUNC( inq_varid )( _fileId, "nodeCoords", &nodeCoordsId );
    if( success ) MB_CHK_SET_ERR( MB_FAILURE, "Trouble getting nodeCoords" );

    auto vmit = varInfo.find( "nodeCoords" );
    if( varInfo.end() == vmit )
        MB_SET_ERR( MB_FAILURE, "Couldn't find variable "
                                    << "nodeCoords" );
    ReadNC::VarData& glData = vmit->second;
    auto attIt              = glData.varAtts.find( "units" );
    if( attIt != glData.varAtts.end() )
    {
        unsigned int sz = attIt->second.attLen;
        std::string att_data;
        att_data.resize( sz + 1 );
        att_data[sz] = '\000';
        success =
            NCFUNC( get_att_text )( _fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0] );
        if( 0 == success && att_data.find( "radians" ) != std::string::npos ) degrees = false;
    }

    // Hack: create dummy variables for dimensions (like nCells) with no corresponding coordinate
    // variables
    rval = create_dummy_variables();MB_CHK_SET_ERR( rval, "Failed to create dummy variables" );

    return MB_SUCCESS;
}

ErrorCode NCHelperESMF::create_mesh( Range& faces )
{
    bool& noMixedElements = _readNC->noMixedElements;
    DebugOutput& dbgOut   = _readNC->dbgOut;

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

    if( procs >= 2 )
    {
        // Shift rank to obtain a rotated trivial partition
        int shifted_rank           = rank;
        int& trivialPartitionShift = _readNC->trivialPartitionShift;
        if( trivialPartitionShift > 0 ) shifted_rank = ( rank + trivialPartitionShift ) % procs;

        // Compute the number of local cells on this proc
        nLocalCells = int( std::floor( 1.0 * nCells / procs ) );

        // The starting global cell index in the ESMF file for this proc
        int start_cell_idx = shifted_rank * nLocalCells;

        // Number of extra cells after equal split over procs
        int iextra = nCells % procs;

        // Allocate extra cells over procs
        if( shifted_rank < iextra ) nLocalCells++;
        start_cell_idx += std::min( shifted_rank, iextra );

        start_cell_idx++;  // 0 based -> 1 based

        // Redistribute local cells after trivial partition (e.g. apply Zoltan partition)
        ErrorCode rval = redistribute_local_cells( start_cell_idx, myPcomm );MB_CHK_SET_ERR( rval, "Failed to redistribute local cells after trivial partition" );
    }
    else
    {
        nLocalCells = nCells;
        localGidCells.insert( 1, nLocalCells );
    }
#else
    nLocalCells = nCells;
    localGidCells.insert( 1, nLocalCells );
#endif
    dbgOut.tprintf( 1, " localGidCells.psize() = %d\n", (int)localGidCells.psize() );
    dbgOut.tprintf( 1, " localGidCells.size() = %d\n", (int)localGidCells.size() );

    // Read number of edges on each local cell, to calculate actual maxEdgesPerCell
    int nEdgesOnCellVarId;
    int success = NCFUNC( inq_varid )( _fileId, "numElementConn", &nEdgesOnCellVarId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of numElementConn" );
    std::vector< int > num_edges_on_local_cells( nLocalCells );

#ifdef MOAB_HAVE_PNETCDF
    size_t nb_reads = localGidCells.psize();
    std::vector< int > requests( nb_reads );
    std::vector< int > statuss( nb_reads );
    size_t idxReq = 0;
#endif
    size_t indexInArray = 0;
    for( Range::pair_iterator pair_iter = localGidCells.pair_begin(); pair_iter != localGidCells.pair_end();
         ++pair_iter )
    {
        EntityHandle starth  = pair_iter->first;
        EntityHandle endh    = pair_iter->second;
        NCDF_SIZE read_start = (NCDF_SIZE)( starth - 1 );
        NCDF_SIZE read_count = (NCDF_SIZE)( endh - starth + 1 );

        // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
        success = NCFUNCREQG( _vara_int )( _fileId, nEdgesOnCellVarId, &read_start, &read_count,
                                           &( num_edges_on_local_cells[indexInArray] ), &requests[idxReq++] );
#else
        success = NCFUNCAG( _vara_int )( _fileId, nEdgesOnCellVarId, &read_start, &read_count,
                                         &( num_edges_on_local_cells[indexInArray] ) );
#endif
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read nEdgesOnCell data in a loop" );

        // Increment the index for next subrange
        indexInArray += ( endh - starth + 1 );
    }

#ifdef MOAB_HAVE_PNETCDF
    // Wait outside the loop
    success = NCFUNC( wait_all )( _fileId, requests.size(), &requests[0], &statuss[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
#endif

    // Read vertices on each local cell, to get localGidVerts and cell connectivity later
    int verticesOnCellVarId;
    success = NCFUNC( inq_varid )( _fileId, "elementConn", &verticesOnCellVarId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of elementConn" );
    std::vector< int > vertices_on_local_cells( nLocalCells * maxEdgesPerCell );
#ifdef MOAB_HAVE_PNETCDF
    idxReq = 0;
#endif
    indexInArray = 0;
    for( Range::pair_iterator pair_iter = localGidCells.pair_begin(); pair_iter != localGidCells.pair_end();
         ++pair_iter )
    {
        EntityHandle starth      = pair_iter->first;
        EntityHandle endh        = pair_iter->second;
        NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( starth - 1 ), 0 };
        NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( endh - starth + 1 ),
                                     static_cast< NCDF_SIZE >( maxEdgesPerCell ) };

        // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
        success = NCFUNCREQG( _vara_int )( _fileId, verticesOnCellVarId, read_starts, read_counts,
                                           &( vertices_on_local_cells[indexInArray] ), &requests[idxReq++] );
#else
        success = NCFUNCAG( _vara_int )( _fileId, verticesOnCellVarId, read_starts, read_counts,
                                         &( vertices_on_local_cells[indexInArray] ) );
#endif
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read verticesOnCell data in a loop" );

        // Increment the index for next subrange
        indexInArray += ( endh - starth + 1 ) * maxEdgesPerCell;
    }

#ifdef MOAB_HAVE_PNETCDF
    // Wait outside the loop
    success = NCFUNC( wait_all )( _fileId, requests.size(), &requests[0], &statuss[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
#endif

    // Correct local cell vertices array, replace the padded vertices with the last vertices
    // in the corresponding cells; sometimes the padded vertices are 0, sometimes a large
    // vertex id. Make sure they are consistent to our padded option
    for( int local_cell_idx = 0; local_cell_idx < nLocalCells; local_cell_idx++ )
    {
        int num_edges             = num_edges_on_local_cells[local_cell_idx];
        int idx_in_local_vert_arr = local_cell_idx * maxEdgesPerCell;
        int last_vert_idx         = vertices_on_local_cells[idx_in_local_vert_arr + num_edges - 1];
        for( int i = num_edges; i < maxEdgesPerCell; i++ )
            vertices_on_local_cells[idx_in_local_vert_arr + i] = last_vert_idx;
    }

    // Create local vertices
    EntityHandle start_vertex;
    ErrorCode rval = create_local_vertices( vertices_on_local_cells, start_vertex );MB_CHK_SET_ERR( rval, "Failed to create local vertices for MPAS mesh" );

    // Create local cells, either unpadded or padded
    if( noMixedElements )
    {
        rval = create_padded_local_cells( vertices_on_local_cells, start_vertex, faces );MB_CHK_SET_ERR( rval, "Failed to create padded local cells for MPAS mesh" );
    }
    else
    {
        rval = create_local_cells( vertices_on_local_cells, num_edges_on_local_cells, start_vertex, faces );MB_CHK_SET_ERR( rval, "Failed to create local cells for MPAS mesh" );
    }

    return MB_SUCCESS;
}

#ifdef MOAB_HAVE_MPI
ErrorCode NCHelperESMF::redistribute_local_cells( int start_cell_idx, ParallelComm* pco )
{
    // If possible, apply Zoltan partition
    // will start from trivial partition in cell space
    // will read cell connectivities, coordinates of vertices in conn, and then compute centers of the cells
    //
#ifdef MOAB_HAVE_ZOLTAN
    if( ScdParData::RCBZOLTAN == _readNC->partMethod )
    {
        // if we have centerCoords, use them; if not, compute them for Zoltan to work
        std::vector< double > xverts, yverts, zverts;
        xverts.resize( nLocalCells );
        yverts.resize( nLocalCells );
        zverts.resize( nLocalCells );

        if( centerCoordsId >= 0 )
        {
            // read from file
            std::vector< double > coords( nLocalCells * coordDim );

            NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( start_cell_idx - 1 ), 0 };
            NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( nLocalCells ),
                                         static_cast< NCDF_SIZE >( coordDim ) };
            int success = NCFUNCAG( _vara_double )( _fileId, centerCoordsId, read_starts, read_counts, &( coords[0] ) );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed on reading center coordinates" );
            if( 2 == coordDim )
            {
                double factor = 1.;
                if( degrees ) factor = pideg;
                for( int i = 0; i < nLocalCells; i++ )
                {
                    double lon    = coords[i * 2] * factor;
                    double lat    = coords[i * 2 + 1] * factor;
                    double cosphi = cos( lat );
                    zverts[i]     = sin( lat );
                    xverts[i]     = cosphi * cos( lon );
                    yverts[i]     = cosphi * sin( lon );
                }
            }
        }
        else
        {
            // Read connectivities
            // Read vertices on each local cell, to get localGidVerts and cell connectivity later
            int verticesOnCellVarId;
            int success = NCFUNC( inq_varid )( _fileId, "elementConn", &verticesOnCellVarId );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of elementConn" );
            std::vector< int > vertices_on_local_cells( nLocalCells * maxEdgesPerCell );

            NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( start_cell_idx - 1 ), 0 };
            NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( nLocalCells ),
                                         static_cast< NCDF_SIZE >( maxEdgesPerCell ) };

            success = NCFUNCAG( _vara_int )( _fileId, verticesOnCellVarId, read_starts, read_counts,
                                             &( vertices_on_local_cells[0] ) );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed on reading elementConn variable" );
            std::vector< int > num_edges_on_local_cells( nLocalCells );

            NCDF_SIZE read_start = start_cell_idx - 1;
            NCDF_SIZE read_count = (NCDF_SIZE)( nLocalCells );

            int nEdgesOnCellVarId;
            success = NCFUNC( inq_varid )( _fileId, "numElementConn", &nEdgesOnCellVarId );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of numElementConn" );
            success = NCFUNCAG( _vara_int )( _fileId, nEdgesOnCellVarId, &read_start, &read_count,
                                             &( num_edges_on_local_cells[0] ) );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get numElementConn" );
            // Make a copy of vertices_on_local_cells for sorting (keep original one to set cell
            // connectivity later)

            // Correct local cell vertices array, replace the padded vertices with the last vertices
            // in the corresponding cells; sometimes the padded vertices are 0, sometimes a large
            // vertex id. Make sure they are consistent to our padded option
            for( int local_cell_idx = 0; local_cell_idx < nLocalCells; local_cell_idx++ )
            {
                int num_edges             = num_edges_on_local_cells[local_cell_idx];
                int idx_in_local_vert_arr = local_cell_idx * maxEdgesPerCell;
                int last_vert_idx         = vertices_on_local_cells[idx_in_local_vert_arr + num_edges - 1];
                for( int i = num_edges; i < maxEdgesPerCell; i++ )
                    vertices_on_local_cells[idx_in_local_vert_arr + i] = last_vert_idx;
            }

            std::vector< int > vertices_on_local_cells_sorted( vertices_on_local_cells );
            std::sort( vertices_on_local_cells_sorted.begin(), vertices_on_local_cells_sorted.end() );
            std::copy( vertices_on_local_cells_sorted.rbegin(), vertices_on_local_cells_sorted.rend(),
                       range_inserter( localGidVerts ) );
            nLocalVertices = localGidVerts.size();

            // Read  nodeCoords for local vertices
            double* coords = new double[localGidVerts.size() * coordDim];
            int nodeCoordVarId;
            success = NCFUNC( inq_varid )( _fileId, "nodeCoords", &nodeCoordVarId );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of nodeCoords" );
#ifdef MOAB_HAVE_PNETCDF
            size_t nb_reads = localGidVerts.psize();
            std::vector< int > requests( nb_reads );
            std::vector< int > statuss( nb_reads );
            size_t idxReq = 0;
#endif
            size_t indexInArray = 0;
            for( Range::pair_iterator pair_iter = localGidVerts.pair_begin(); pair_iter != localGidVerts.pair_end();
                 ++pair_iter )
            {
                EntityHandle starth      = pair_iter->first;
                EntityHandle endh        = pair_iter->second;
                NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( starth - 1 ), 0 };
                NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( endh - starth + 1 ),
                                             static_cast< NCDF_SIZE >( coordDim ) };

                // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
                success = NCFUNCREQG( _vara_double )( _fileId, nodeCoordVarId, read_starts, read_counts,
                                                      &coords[indexInArray], &requests[idxReq++] );
#else
                success = NCFUNCAG( _vara_double )( _fileId, nodeCoordVarId, read_starts, read_counts,
                                                    &coords[indexInArray] );
#endif
                if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read nodeCoords data in a loop" );

                // Increment the index for next subrange
                indexInArray += ( endh - starth + 1 ) * coordDim;
            }

#ifdef MOAB_HAVE_PNETCDF
            // Wait outside the loop
            success = NCFUNC( wait_all )( _fileId, requests.size(), &requests[0], &statuss[0] );
            if( success ) MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
#endif

            double* coords3d = new double[localGidVerts.size() * 3];
            // now convert from lat/lon to 3d
            if( 2 == coordDim )
            {
                // basically convert from degrees to xyz on a sphere
                double factor = 1.;
                if( degrees ) factor = pideg;
                for( int i = 0; i < (int)localGidVerts.size(); i++ )
                {
                    double lon          = coords[i * 2] * factor;
                    double lat          = coords[i * 2 + 1] * factor;
                    double cosphi       = cos( lat );
                    coords3d[3 * i + 2] = sin( lat );
                    coords3d[3 * i]     = cosphi * cos( lon );
                    coords3d[3 * i + 1] = cosphi * sin( lon );
                }
            }

            // find the center for each cell
            for( int i = 0; i < nLocalCells; i++ )
            {
                int nv   = num_edges_on_local_cells[i];
                double x = 0., y = 0., z = 0.;
                for( int j = 0; j < nv; j++ )
                {
                    int vertexId = vertices_on_local_cells[i * maxEdgesPerCell + j];
                    // now find what index is in gid
                    int index = localGidVerts.index( vertexId );  // it should be a binary search!
                    x += coords3d[3 * index];
                    y += coords3d[3 * index + 1];
                    z += coords3d[3 * index + 2];
                }
                if( nv != 0 )
                {
                    x /= nv;
                    y /= nv;
                    z /= nv;
                }
                xverts[i] = x;
                yverts[i] = y;
                zverts[i] = z;
            }
            delete[] coords3d;
        }
        // Zoltan partition using RCB; maybe more studies would be good, as to which partition
        // is better
        Interface*& mbImpl         = _readNC->mbImpl;
        DebugOutput& dbgOut        = _readNC->dbgOut;
        ZoltanPartitioner* mbZTool = new ZoltanPartitioner( mbImpl, pco, false, 0, NULL );
        ErrorCode rval = mbZTool->repartition( xverts, yverts, zverts, start_cell_idx, "RCB", localGidCells );MB_CHK_SET_ERR( rval, "Error in Zoltan partitioning" );
        delete mbZTool;

        dbgOut.tprintf( 1, "After Zoltan partitioning, localGidCells.psize() = %d\n", (int)localGidCells.psize() );
        dbgOut.tprintf( 1, "                           localGidCells.size() = %d\n", (int)localGidCells.size() );

        // This is important: local cells are now redistributed, so nLocalCells might be different!
        nLocalCells = localGidCells.size();
        localGidVerts.clear();
        return MB_SUCCESS;
    }
#else
    UNUSED(pco);
#endif

    // By default, apply trivial partition
    localGidCells.insert( start_cell_idx, start_cell_idx + nLocalCells - 1 );

    return MB_SUCCESS;
}
#endif

ErrorCode NCHelperESMF::create_local_vertices( const std::vector< int >& vertices_on_local_cells,
                                               EntityHandle& start_vertex )
{
    Interface*& mbImpl      = _readNC->mbImpl;
    Tag& mGlobalIdTag       = _readNC->mGlobalIdTag;
    const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
    DebugOutput& dbgOut     = _readNC->dbgOut;

    // Make a copy of vertices_on_local_cells for sorting (keep original one to set cell
    // connectivity later)
    std::vector< int > vertices_on_local_cells_sorted( vertices_on_local_cells );
    std::sort( vertices_on_local_cells_sorted.begin(), vertices_on_local_cells_sorted.end() );
    std::copy( vertices_on_local_cells_sorted.rbegin(), vertices_on_local_cells_sorted.rend(),
               range_inserter( localGidVerts ) );
    nLocalVertices = localGidVerts.size();

    dbgOut.tprintf( 1, "   localGidVerts.psize() = %d\n", (int)localGidVerts.psize() );
    dbgOut.tprintf( 1, "   localGidVerts.size() = %d\n", (int)localGidVerts.size() );

    // Create local vertices
    std::vector< double* > arrays;
    ErrorCode rval =
        _readNC->readMeshIface->get_node_coords( 3, nLocalVertices, 0, start_vertex, arrays, nLocalVertices );MB_CHK_SET_ERR( rval, "Failed to create local vertices" );

    // Add local vertices to current file set
    Range local_verts_range( start_vertex, start_vertex + nLocalVertices - 1 );
    rval = _readNC->mbImpl->add_entities( _fileSet, local_verts_range );MB_CHK_SET_ERR( rval, "Failed to add local vertices to current file set" );

    // Get ptr to GID memory for local vertices
    int count  = 0;
    void* data = NULL;
    rval       = mbImpl->tag_iterate( mGlobalIdTag, local_verts_range.begin(), local_verts_range.end(), count, data );MB_CHK_SET_ERR( rval, "Failed to iterate global id tag on local vertices" );
    assert( count == nLocalVertices );
    int* gid_data = (int*)data;
    std::copy( localGidVerts.begin(), localGidVerts.end(), gid_data );

    // Duplicate GID data, which will be used to resolve sharing
    if( mpFileIdTag )
    {
        rval = mbImpl->tag_iterate( *mpFileIdTag, local_verts_range.begin(), local_verts_range.end(), count, data );MB_CHK_SET_ERR( rval, "Failed to iterate file id tag on local vertices" );
        assert( count == nLocalVertices );
        int bytes_per_tag = 4;
        rval              = mbImpl->tag_get_bytes( *mpFileIdTag, bytes_per_tag );MB_CHK_SET_ERR( rval, "Can't get number of bytes for file id tag" );
        if( 4 == bytes_per_tag )
        {
            gid_data = (int*)data;
            std::copy( localGidVerts.begin(), localGidVerts.end(), gid_data );
        }
        else if( 8 == bytes_per_tag )
        {  // Should be a handle tag on 64 bit machine?
            long* handle_tag_data = (long*)data;
            std::copy( localGidVerts.begin(), localGidVerts.end(), handle_tag_data );
        }
    }

#ifdef MOAB_HAVE_PNETCDF
    size_t nb_reads = localGidVerts.psize();
    std::vector< int > requests( nb_reads );
    std::vector< int > statuss( nb_reads );
    size_t idxReq = 0;
#endif

    // Read  nodeCoords for local vertices
    double* coords = new double[localGidVerts.size() * coordDim];
    int nodeCoordVarId;
    int success = NCFUNC( inq_varid )( _fileId, "nodeCoords", &nodeCoordVarId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of nodeCoords" );
    size_t indexInArray = 0;
    for( Range::pair_iterator pair_iter = localGidVerts.pair_begin(); pair_iter != localGidVerts.pair_end();
         ++pair_iter )
    {
        EntityHandle starth      = pair_iter->first;
        EntityHandle endh        = pair_iter->second;
        NCDF_SIZE read_starts[2] = { static_cast< NCDF_SIZE >( starth - 1 ), 0 };
        NCDF_SIZE read_counts[2] = { static_cast< NCDF_SIZE >( endh - starth + 1 ),
                                     static_cast< NCDF_SIZE >( coordDim ) };

        // Do a partial read in each subrange
#ifdef MOAB_HAVE_PNETCDF
        success = NCFUNCREQG( _vara_double )( _fileId, nodeCoordVarId, read_starts, read_counts, &coords[indexInArray],
                                              &requests[idxReq++] );
#else
        success = NCFUNCAG( _vara_double )( _fileId, nodeCoordVarId, read_starts, read_counts, &coords[indexInArray] );
#endif
        if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read nodeCoords data in a loop" );

        // Increment the index for next subrange
        indexInArray += ( endh - starth + 1 ) * coordDim;
    }

#ifdef MOAB_HAVE_PNETCDF
    // Wait outside the loop
    success = NCFUNC( wait_all )( _fileId, requests.size(), &requests[0], &statuss[0] );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
#endif

    // now convert from lat/lon to 3d
    if( 2 == coordDim )
    {
        // basically convert from degrees to xyz on a sphere
        double factor = 1;
        if( degrees ) factor = pideg;
        for( int i = 0; i < (int)localGidVerts.size(); i++ )
        {
            double lon    = coords[i * 2] * factor;
            double lat    = coords[i * 2 + 1] * factor;
            double cosphi = cos( lat );
            arrays[2][i]  = sin( lat );
            arrays[0][i]  = cosphi * cos( lon );
            arrays[1][i]  = cosphi * sin( lon );
        }
    }
    delete[] coords;
    return MB_SUCCESS;
}

ErrorCode NCHelperESMF::create_local_cells( const std::vector< int >& vertices_on_local_cells,
                                            const std::vector< int >& num_edges_on_local_cells,
                                            EntityHandle start_vertex,
                                            Range& faces )
{
    Interface*& mbImpl = _readNC->mbImpl;
    Tag& mGlobalIdTag  = _readNC->mGlobalIdTag;

    // Divide local cells into groups based on the number of edges
    Range local_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
    // Insert larger values before smaller ones to increase efficiency
    for( int i = nLocalCells - 1; i >= 0; i-- )
    {
        int num_edges = num_edges_on_local_cells[i];
        local_cells_with_n_edges[num_edges].insert( localGidCells[i] );  // Global cell index
    }

    std::vector< int > num_edges_on_cell_groups;
    for( int i = 3; i <= maxEdgesPerCell; i++ )
    {
        if( local_cells_with_n_edges[i].size() > 0 ) num_edges_on_cell_groups.push_back( i );
    }
    int numCellGroups = (int)num_edges_on_cell_groups.size();
    EntityHandle* conn_arr_local_cells_with_n_edges[DEFAULT_MAX_EDGES_PER_CELL + 1];
    for( int i = 0; i < numCellGroups; i++ )
    {
        int num_edges_per_cell = num_edges_on_cell_groups[i];
        int num_group_cells    = (int)local_cells_with_n_edges[num_edges_per_cell].size();

        EntityType typeEl = MBTRI;
        if( num_edges_per_cell == 4 ) typeEl = MBQUAD;
        if( num_edges_per_cell > 4 ) typeEl = MBPOLYGON;
        // Create local cells for each non-empty cell group
        EntityHandle start_element;
        ErrorCode rval =
            _readNC->readMeshIface->get_element_connect( num_group_cells, num_edges_per_cell, typeEl, 0, start_element,
                                                         conn_arr_local_cells_with_n_edges[num_edges_per_cell],
                                                         num_group_cells );MB_CHK_SET_ERR( rval, "Failed to create local cells" );
        faces.insert( start_element, start_element + num_group_cells - 1 );

        // Add local cells to current file set
        Range local_cells_range( start_element, start_element + num_group_cells - 1 );
        rval = _readNC->mbImpl->add_entities( _fileSet, local_cells_range );MB_CHK_SET_ERR( rval, "Failed to add local cells to current file set" );

        // Get ptr to gid memory for local cells
        int count  = 0;
        void* data = NULL;
        rval = mbImpl->tag_iterate( mGlobalIdTag, local_cells_range.begin(), local_cells_range.end(), count, data );MB_CHK_SET_ERR( rval, "Failed to iterate global id tag on local cells" );
        assert( count == num_group_cells );
        int* gid_data = (int*)data;
        std::copy( local_cells_with_n_edges[num_edges_per_cell].begin(),
                   local_cells_with_n_edges[num_edges_per_cell].end(), gid_data );

        // Set connectivity array with proper local vertices handles
        for( int j = 0; j < num_group_cells; j++ )
        {
            EntityHandle global_cell_idx =
                local_cells_with_n_edges[num_edges_per_cell][j];          // Global cell index, 1 based
            int local_cell_idx = localGidCells.index( global_cell_idx );  // Local cell index, 0 based
            assert( local_cell_idx != -1 );

            for( int k = 0; k < num_edges_per_cell; k++ )
            {
                EntityHandle global_vert_idx =
                    vertices_on_local_cells[local_cell_idx * maxEdgesPerCell + k];  // Global vertex index, 1 based
                int local_vert_idx = localGidVerts.index( global_vert_idx );        // Local vertex index, 0 based
                assert( local_vert_idx != -1 );
                conn_arr_local_cells_with_n_edges[num_edges_per_cell][j * num_edges_per_cell + k] =
                    start_vertex + local_vert_idx;
            }
        }
    }

    return MB_SUCCESS;
}

ErrorCode NCHelperESMF::create_padded_local_cells( const std::vector< int >& vertices_on_local_cells,
                                                   EntityHandle start_vertex,
                                                   Range& faces )
{
    Interface*& mbImpl = _readNC->mbImpl;
    Tag& mGlobalIdTag  = _readNC->mGlobalIdTag;

    // Create cells for this cell group
    EntityHandle start_element;
    EntityHandle* conn_arr_local_cells = NULL;
    ErrorCode rval = _readNC->readMeshIface->get_element_connect( nLocalCells, maxEdgesPerCell, MBPOLYGON, 0,
                                                                  start_element, conn_arr_local_cells, nLocalCells );MB_CHK_SET_ERR( rval, "Failed to create local cells" );
    faces.insert( start_element, start_element + nLocalCells - 1 );

    // Add local cells to current file set
    Range local_cells_range( start_element, start_element + nLocalCells - 1 );
    rval = _readNC->mbImpl->add_entities( _fileSet, local_cells_range );MB_CHK_SET_ERR( rval, "Failed to add local cells to current file set" );

    // Get ptr to GID memory for local cells
    int count  = 0;
    void* data = NULL;
    rval       = mbImpl->tag_iterate( mGlobalIdTag, local_cells_range.begin(), local_cells_range.end(), count, data );MB_CHK_SET_ERR( rval, "Failed to iterate global id tag on local cells" );
    assert( count == nLocalCells );
    int* gid_data = (int*)data;
    std::copy( localGidCells.begin(), localGidCells.end(), gid_data );

    // Set connectivity array with proper local vertices handles
    // vertices_on_local_cells array was already corrected to have the last vertices padded
    // no need for extra checks considering
    for( int local_cell_idx = 0; local_cell_idx < nLocalCells; local_cell_idx++ )
    {
        for( int i = 0; i < maxEdgesPerCell; i++ )
        {
            EntityHandle global_vert_idx =
                vertices_on_local_cells[local_cell_idx * maxEdgesPerCell + i];  // Global vertex index, 1 based
            int local_vert_idx = localGidVerts.index( global_vert_idx );        // Local vertex index, 0 based
            assert( local_vert_idx != -1 );
            conn_arr_local_cells[local_cell_idx * maxEdgesPerCell + i] = start_vertex + local_vert_idx;
        }
    }

    return MB_SUCCESS;
}

}  // namespace moab
