/*
 * NCHelperTOPO.cpp
 *
 *  Created on: Sep. 12, 2023
 *
 *  Tested with: mpiexec -n 4 tools/mbconvert  ../../usgs-rawdata.nc usgs-rawdata.h5m -O "DEBUG_IO=3" -O "PARALLEL=READ_PART" -o "PARALLEL=WRITE_PART"
 */

#include "NCHelperTOPO.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"

#ifdef MOAB_HAVE_ZOLTAN
#include "moab/ZoltanPartitioner.hpp"
#endif

#include <cmath>

namespace moab
{

const double pideg                   = acos( -1.0 ) / 180.0;

NCHelperTOPO::NCHelperTOPO( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
    : UcdNCHelper( readNC, fileId, opts, fileSet ), coordDim( 0 ), degrees( true )
{
}

bool NCHelperTOPO::can_read_file( ReadNC* readNC )
{
    std::vector< std::string >& dimNames = readNC->dimNames;
    std::map< std::string, ReadNC::VarData >& varInfo = readNC->varInfo;

    // Check for USGS topo format: lat/lon dimensions and htopo variable
    if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "lon" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "lat" ) ) != dimNames.end() ) &&
        ( varInfo.find( "htopo" ) != varInfo.end() ) )
    {
        return true;
    }
    else if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "longitude" ) ) != dimNames.end() ) &&
        ( std::find( dimNames.begin(), dimNames.end(), std::string( "latitude" ) ) != dimNames.end() ) &&
        ( varInfo.find( "htopo" ) != varInfo.end() ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

ErrorCode NCHelperTOPO::init_mesh_vals()
{
    std::vector< std::string >& dimNames              = _readNC->dimNames;
    std::vector< int >& dimLens                       = _readNC->dimLens;
    std::map< std::string, ReadNC::VarData >& varInfo = _readNC->varInfo;

    ErrorCode rval;
    unsigned int idx;
    std::vector< std::string >::iterator vit;

    // Find lat/lon dimensions for USGS topo format
    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "lat" ) ) != dimNames.end() )
    {
        latDimName = "lat";
        idx = vit - dimNames.begin();
        nLatVals = dimLens[idx];
    }
    else if( ( vit = std::find( dimNames.begin(), dimNames.end(), "latitude" ) ) != dimNames.end() )
    {
        latDimName = "latitude";
        idx = vit - dimNames.begin();
        nLatVals = dimLens[idx];
    }
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find latitude dimension" );
    }

    if( ( vit = std::find( dimNames.begin(), dimNames.end(), "lon" ) ) != dimNames.end() )
    {
        lonDimName = "lon";
        idx = vit - dimNames.begin();
        nLonVals = dimLens[idx];
    }
    else if( ( vit = std::find( dimNames.begin(), dimNames.end(), "longitude" ) ) != dimNames.end() )
    {
        lonDimName = "longitude";
        idx = vit - dimNames.begin();
        nLonVals = dimLens[idx];
    }
    else
    {
        MB_SET_ERR( MB_FAILURE, "Couldn't find longitude dimension" );
    }

    // For USGS topo format, we create a point cloud (vertices only)
    nVertices = nLatVals * nLonVals;
    nCells = 0;  // No cells for point cloud
    coordDim = 2;  // lat/lon coordinates
    degrees = true;  // USGS data is typically in degrees

    // Check units of lat/lon variables to confirm degrees
    auto latVarIt = varInfo.find( latDimName );
    if( latVarIt != varInfo.end() )
    {
        auto attIt = latVarIt->second.varAtts.find( "units" );
        if( attIt != latVarIt->second.varAtts.end() )
        {
            unsigned int sz = attIt->second.attLen;
            std::string att_data;
            att_data.resize( sz + 1 );
            att_data[sz] = '\000';
            int success = NCFUNC( get_att_text )( _fileId, attIt->second.attVarId, attIt->second.attName.c_str(), &att_data[0] );
            if( 0 == success && att_data.find( "radians" ) != std::string::npos ) degrees = false;
        }
    }

    // Hack: create dummy variables for dimensions (like nCells) with no corresponding coordinate
    // variables
    rval = create_dummy_variables();MB_CHK_SET_ERR( rval, "Failed to create dummy variables" );

    return MB_SUCCESS;
}

ErrorCode NCHelperTOPO::create_mesh( Range& faces )
{
    Interface*& mbImpl      = _readNC->mbImpl;
    Tag& mGlobalIdTag       = _readNC->mGlobalIdTag;
    const Tag*& mpFileIdTag = _readNC->mpFileIdTag;
    DebugOutput& dbgOut     = _readNC->dbgOut;
    ErrorCode rval;

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

    dbgOut.tprintf( 2, "Rank %d: isParallel=%s, procs=%d\n", rank, isParallel ? "true" : "false", procs );

    // Determine local vertex range for parallel reading
    size_t nLocalVertices;
    size_t startLatIdx = 0, endLatIdx = nLatVals - 1;
    size_t startLonIdx = 0, endLonIdx = nLonVals - 1;

    if( isParallel && procs > 1 )
    {
        // Simple 1D decomposition along latitude dimension
        size_t shifted_rank           = rank;
        int& trivialPartitionShift = _readNC->trivialPartitionShift;
        if( trivialPartitionShift > 0 ) shifted_rank = ( rank + trivialPartitionShift ) % procs;

        dbgOut.tprintf( 2, "Rank %d: Computing parallel decomposition with %d processes\n", rank, procs );

        // Compute local latitude range
        size_t nLocalLats = size_t( std::floor( 1.0 * nLatVals / procs ) );
        startLatIdx = shifted_rank * nLocalLats;

        // Handle remainder
        size_t iextra = nLatVals % procs;
        if( shifted_rank < iextra ) nLocalLats++;
        startLatIdx += std::min( shifted_rank, iextra );

        endLatIdx = startLatIdx + nLocalLats - 1;
        nLocalVertices = nLocalLats * nLonVals;

        // dbgOut.tprintf( 2, "Rank %d: Before MPI_Scan: lat range [%zu,%zu] (%zu lats) x lon range [%zu,%zu] (%d lons)\n",
        //                rank, startLatIdx, endLatIdx, nLocalLats, startLonIdx, endLonIdx, nLonVals );
        // size_t offsetLatIdx = (rank == 0) ? 0 : startLatIdx;
        // MPI_Scan( &offsetLatIdx, &startLatIdx, 1, MPI_UNSIGNED_LONG, MPI_SUM, myPcomm->comm() );

        dbgOut.tprintf( 2, "Rank %d: After MPI_Scan: Assigned lat range [%zu,%zu] (%zu lats) x lon range [%zu,%zu] (%d lons)\n",
                       rank, startLatIdx, endLatIdx, nLocalLats, startLonIdx, endLonIdx, nLonVals );
    }
    else
    {
        // Even in serial mode, use decomposition for very large datasets to avoid PNetCDF limits
        if( nVertices > INT_MAX )
        {
            // Use a single "chunk" that fits within PNetCDF limits
            int maxChunkSize = INT_MAX / 2;  // Conservative limit
            if( nVertices > maxChunkSize )
            {
                // For very large datasets, just read a representative subset
                size_t nLocalLats = std::min( (size_t)1000, (size_t)nLatVals );  // Read up to 1000 latitude bands
                endLatIdx = startLatIdx + nLocalLats - 1;
                nLocalVertices = nLocalLats * nLonVals;

                dbgOut.tprintf( 1, "Warning: Dataset too large for single read. Reading subset: %zu vertices\n", nLocalVertices );
            }
            else
            {
                nLocalVertices = nVertices;
            }
        }
        else
        {
            nLocalVertices = nVertices;
        }
    }
#else
    // Same logic for non-MPI builds
    if( nVertices > INT_MAX )
    {
        size_t maxChunkSize = INT_MAX / 2;
        if( nVertices > maxChunkSize )
        {
            size_t nLocalLats = std::min( (size_t)1000, (size_t)nLatVals );
            endLatIdx = startLatIdx + nLocalLats - 1;
            nLocalVertices = nLocalLats * nLonVals;

            dbgOut.tprintf( 1, "Warning: Dataset too large for single read. Reading subset: %zu vertices\n", nLocalVertices );
        }
        else
        {
            nLocalVertices = nVertices;
        }
    }
    else
    {
        nLocalVertices = nVertices;
    }
#endif

    dbgOut.tprintf( 1, " Creating %zu local vertices from lat range [%zu,%zu] x lon range [%zu,%zu]\n",
                    nLocalVertices, startLatIdx, endLatIdx, startLonIdx, endLonIdx );

    // Create local vertices for the point cloud
    std::vector< double* > arrays;
    EntityHandle start_vertex;
    rval = _readNC->readMeshIface->get_node_coords( 3, nLocalVertices, 0, start_vertex, arrays, nLocalVertices );
    MB_CHK_SET_ERR( rval, "Failed to create local vertices" );

    // Add local vertices to current file set
    Range local_verts_range( start_vertex, start_vertex + nLocalVertices - 1 );
    rval = mbImpl->add_entities( _fileSet, local_verts_range );
    MB_CHK_SET_ERR( rval, "Failed to add local vertices to current file set" );

    // Set up global IDs for vertices - use long to handle >4 billion vertices
    int count = 0;
    void* data = NULL;
    rval = mbImpl->tag_iterate( mGlobalIdTag, local_verts_range.begin(), local_verts_range.end(), count, data );
    MB_CHK_SET_ERR( rval, "Failed to iterate global id tag on local vertices" );
    assert( count == (int)nLocalVertices );

    // Check the actual tag data type to handle properly
    int tag_size;
    rval = mbImpl->tag_get_bytes( mGlobalIdTag, tag_size );
    MB_CHK_SET_ERR( rval, "Failed to get global id tag size" );

    std::cout << "Global ID tag size: " << tag_size << std::endl;

    // Assign global IDs based on lat/lon indices
    size_t gid_idx = 0;
    for( size_t lat_idx = startLatIdx; lat_idx <= endLatIdx; lat_idx++ )
    {
        for( size_t lon_idx = startLonIdx; lon_idx <= endLonIdx; lon_idx++ )
        {
            size_t global_id = lat_idx * nLonVals + lon_idx + 1;  // 1-based global ID

            if( tag_size == 4 )
            {
                int* gid_data = reinterpret_cast<int*>(data);
                // Check for overflow
                if( global_id > INT_MAX )
                {
                    MB_SET_ERR( MB_FAILURE, "Global ID " << global_id << " exceeds 32-bit integer limit. Use 64-bit build or smaller mesh." );
                }
                gid_data[gid_idx] = static_cast<int>(global_id);
            }
            else if( tag_size == 8 )
            {
                size_t* gid_data = reinterpret_cast<size_t*>(data);
                gid_data[gid_idx] = global_id;
            }
            else
            {
                MB_SET_ERR( MB_FAILURE, "Unsupported global id tag size: " << tag_size );
            }
            gid_idx++;
        }
    }

    // Duplicate GID data for file ID tag if needed
    if( mpFileIdTag )
    {
        void* fid_data_ptr = nullptr;
        rval = mbImpl->tag_iterate( *mpFileIdTag, local_verts_range.begin(), local_verts_range.end(), count, fid_data_ptr );
        MB_CHK_SET_ERR( rval, "Failed to iterate file id tag on local vertices" );
        assert( count == (int)nLocalVertices );

        int fid_bytes_per_tag = 4;
        rval = mbImpl->tag_get_bytes( *mpFileIdTag, fid_bytes_per_tag );
        MB_CHK_SET_ERR( rval, "Can't get number of bytes for file id tag" );

        // Copy global IDs to file ID tag
        for( size_t i = 0; i < nLocalVertices; i++ )
        {
            size_t lat_idx = startLatIdx + i / (endLonIdx - startLonIdx + 1);
            size_t lon_idx = startLonIdx + i % (endLonIdx - startLonIdx + 1);
            mbGIDType global_id = static_cast<mbGIDType>(lat_idx * nLonVals + lon_idx + 1);

            if( fid_bytes_per_tag == 4 )
            {
                int* fid_data = reinterpret_cast<int*>(fid_data_ptr);
                if( global_id > INT_MAX )
                {
                    MB_SET_ERR( MB_FAILURE, "Global ID " << global_id << " exceeds 32-bit integer limit for file ID tag." );
                }
                fid_data[i] = static_cast<int>(global_id);
            }
            else if( fid_bytes_per_tag == 8 )
            {
                mbGIDType* fid_data = reinterpret_cast<mbGIDType*>(fid_data_ptr);
                fid_data[i] = global_id;
            }
        }
    }

    // Read lat/lon coordinate arrays
    MB_CHK_SET_ERR( read_coordinate_variables( startLatIdx, endLatIdx, startLonIdx, endLonIdx, arrays ), "Failed to read coordinate variables" );

    // Read data variables (htopo, landfract) and store on vertices
    MB_CHK_SET_ERR( read_data_variables( startLatIdx, endLatIdx, startLonIdx, endLonIdx, local_verts_range ), "Failed to read data variables" );

    // No faces created for point cloud
    faces.clear();

    return MB_SUCCESS;
}

ErrorCode NCHelperTOPO::read_coordinate_variables( int startLatIdx, int endLatIdx, int startLonIdx, int endLonIdx,
                                                   std::vector< double* >& arrays )
{
    DebugOutput& dbgOut = _readNC->dbgOut;
    int success;

    // Calculate local dimensions
    size_t nLocalLats = endLatIdx - startLatIdx + 1;
    size_t nLocalLons = endLonIdx - startLonIdx + 1;

    dbgOut.tprintf( 2, "Reading coordinates: lat[%d:%d] (%zu values), lon[%d:%d] (%zu values)\n",
                   startLatIdx, endLatIdx, nLocalLats, startLonIdx, endLonIdx, nLocalLons );

    // Read latitude values - only the local chunk
    int latVarId;
    success = NCFUNC( inq_varid )( _fileId, latDimName.c_str(), &latVarId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of " << latDimName );

    std::vector< float > latVals( nLocalLats );
    NCDF_SIZE lat_start = static_cast< NCDF_SIZE >( startLatIdx );
    NCDF_SIZE lat_count = static_cast< NCDF_SIZE >( nLocalLats );

#ifdef MOAB_HAVE_PNETCDF
    success = NCFUNCAG( _vara_float )( _fileId, latVarId, &lat_start, &lat_count, &latVals[0] );
#else
    success = NCFUNCAG( _vara_float )( _fileId, latVarId, &lat_start, &lat_count, &latVals[0] );
#endif
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read latitude values, error = " << success );

    // Read longitude values - only the local chunk
    int lonVarId;
    success = NCFUNC( inq_varid )( _fileId, lonDimName.c_str(), &lonVarId );
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to get variable id of " << lonDimName );

    std::vector< float > lonVals( nLocalLons );
    NCDF_SIZE lon_start = static_cast< NCDF_SIZE >( startLonIdx );
    NCDF_SIZE lon_count = static_cast< NCDF_SIZE >( nLocalLons );

#ifdef MOAB_HAVE_PNETCDF
    success = NCFUNCAG( _vara_float )( _fileId, lonVarId, &lon_start, &lon_count, &lonVals[0] );
#else
    success = NCFUNCAG( _vara_float )( _fileId, lonVarId, &lon_start, &lon_count, &lonVals[0] );
#endif
    if( success ) MB_SET_ERR( MB_FAILURE, "Failed to read longitude values, error = " << success );

    // Convert lat/lon to 3D Cartesian coordinates on unit sphere
    float factor = 1.0;
    if( degrees ) factor = pideg;

    // Convert lat/lon to 3D Cartesian coordinates on unit sphere
    size_t vert_idx = 0;
    for( size_t lat_idx = 0; lat_idx < nLocalLats; lat_idx++ )
    {
        float lat = latVals[lat_idx] * factor;
        double cosphi = cos( lat );
        double sinphi = sin( lat );

        for( size_t lon_idx = 0; lon_idx < nLocalLons; lon_idx++ )
        {
            double lon = lonVals[lon_idx] * factor;

            // Convert to 3D Cartesian coordinates on unit sphere
            arrays[0][vert_idx] = cosphi * cos( lon );  // x
            arrays[1][vert_idx] = cosphi * sin( lon );  // y
            arrays[2][vert_idx] = sinphi;               // z

            vert_idx++;
        }
    }

    dbgOut.tprintf( 2, "Successfully converted %zu coordinate pairs to 3D Cartesian\n", vert_idx );

    return MB_SUCCESS;
}

ErrorCode NCHelperTOPO::read_data_variables( int startLatIdx, int endLatIdx, int startLonIdx, int endLonIdx,
                                            const Range& local_verts )
{
    Interface*& mbImpl = _readNC->mbImpl;
    DebugOutput& dbgOut = _readNC->dbgOut;
    std::map< std::string, ReadNC::VarData >& varInfo = _readNC->varInfo;

    ErrorCode rval = MB_SUCCESS;
    int success;

    // Calculate local dimensions - use the same decomposition as mesh creation
    size_t nLocalLats = endLatIdx - startLatIdx + 1;
    size_t nLocalLons = endLonIdx - startLonIdx + 1;
    size_t nLocalVertices = nLocalLats * nLocalLons;

    // Verify we don't exceed PNetCDF limits (INT_MAX)
    if( nLocalVertices > INT_MAX )
    {
        MB_SET_ERR( MB_FAILURE, "Local vertex count " << nLocalVertices << " exceeds PNetCDF INT_MAX limit. Use more parallel processes." );
    }

    // Define chunk size for reading large datasets
    const size_t MAX_CHUNK_SIZE = 50000000; // ~50M elements to stay well under INT_MAX

    // Read htopo variable if it exists
    if( varInfo.find( "htopo" ) != varInfo.end() )
    {
        dbgOut.tprintf( 1, "Reading htopo variable...\n" );

        // Get variable ID
        int htopoVarId;
        success = NCFUNC( inq_varid )( _fileId, "htopo", &htopoVarId );
        if( success )
        {
            dbgOut.tprintf( 1, "Warning: Could not find htopo variable\n" );
        }
        else
        {
            // Create tag for htopo
            Tag htopoTag;
            int default_val = 0;
            rval = mbImpl->tag_get_handle( "htopo", 1, MB_TYPE_INTEGER, htopoTag,
                                         MB_TAG_DENSE | MB_TAG_CREAT, &default_val );
            MB_CHK_SET_ERR( rval, "Failed to create htopo tag" );

            // Read htopo data using chunked approach
            dbgOut.tprintf( 2, "Reading htopo data: start=(%d,%d), count=(%zu,%zu), total=%zu\n",
                           startLatIdx, startLonIdx, nLocalLats, nLocalLons, nLocalVertices );

            // Implement chunked reading to avoid PNetCDF limits
            std::vector< int > htopoVals( nLocalVertices );

            if( nLocalVertices > MAX_CHUNK_SIZE )
            {
                dbgOut.tprintf( 1, "Using chunked reading for htopo data (%zu vertices)\n", nLocalVertices );

                // Read in chunks along latitude dimension
                size_t chunk_lat_size = MAX_CHUNK_SIZE / nLocalLons;
                if( chunk_lat_size == 0 ) chunk_lat_size = 1;

                size_t data_offset = 0;
                for( size_t lat_start = 0; lat_start < nLocalLats; lat_start += chunk_lat_size )
                {
                    size_t lat_count = std::min( chunk_lat_size, nLocalLats - lat_start );
                    size_t chunk_size = lat_count * nLocalLons;

                    NCDF_SIZE start[2] = { static_cast<NCDF_SIZE>(startLatIdx + lat_start), static_cast<NCDF_SIZE>(startLonIdx) };
                    NCDF_SIZE count[2] = { static_cast<NCDF_SIZE>(lat_count), static_cast<NCDF_SIZE>(nLocalLons) };

                    dbgOut.tprintf( 3, "Reading htopo chunk: lat_start=%zu, lat_count=%zu, chunk_size=%zu\n",
                                   lat_start, lat_count, chunk_size );

#ifdef MOAB_HAVE_PNETCDF
                    success = NCFUNCAG( _vara_int )( _fileId, htopoVarId, start, count, &htopoVals[data_offset] );
#else
                    success = NCFUNCAG( _vara_int )( _fileId, htopoVarId, start, count, &htopoVals[data_offset] );
#endif
                    if( success )
                    {
                        dbgOut.tprintf( 1, "Warning: Failed to read htopo chunk, NetCDF error = %d\n", success );
                        break;
                    }
                    data_offset += chunk_size;
                }
            }
            else
            {
                // Read all at once if small enough
                NCDF_SIZE start[2] = { static_cast<NCDF_SIZE>(startLatIdx), static_cast<NCDF_SIZE>(startLonIdx) };
                NCDF_SIZE count[2] = { static_cast<NCDF_SIZE>(nLocalLats), static_cast<NCDF_SIZE>(nLocalLons) };

#ifdef MOAB_HAVE_PNETCDF
                success = NCFUNCAG( _vara_int )( _fileId, htopoVarId, start, count, &htopoVals[0] );
#else
                success = NCFUNCAG( _vara_int )( _fileId, htopoVarId, start, count, &htopoVals[0] );
#endif
            }

            if( !success )
            {
                // Set htopo data on vertices
                rval = mbImpl->tag_set_data( htopoTag, local_verts, &htopoVals[0] );
                MB_CHK_SET_ERR( rval, "Failed to set htopo data on vertices" );
                dbgOut.tprintf( 2, "Successfully read htopo data for %zu vertices\n", nLocalVertices );
            }
            else
            {
                dbgOut.tprintf( 1, "Warning: Failed to read htopo data, NetCDF error = %d\n", success );
            }
        }
    }

    // Read landfract variable if it exists
    if( varInfo.find( "landfract" ) != varInfo.end() )
    {
        dbgOut.tprintf( 1, "Reading landfract variable...\n" );

        // Get variable ID
        int landfractVarId;
        success = NCFUNC( inq_varid )( _fileId, "landfract", &landfractVarId );
        if( success )
        {
            dbgOut.tprintf( 1, "Warning: Failed to get landfract variable id\n" );
        }
        else
        {
            // Create tag for landfract
            Tag landfractTag;
            int default_val = 0;
            rval = mbImpl->tag_get_handle( "landfract", 1, MB_TYPE_INTEGER, landfractTag,
                                         MB_TAG_DENSE | MB_TAG_CREAT, &default_val );
            MB_CHK_SET_ERR( rval, "Failed to create landfract tag" );

            // Read landfract data using chunked approach
            dbgOut.tprintf( 2, "Reading landfract data: start=(%d,%d), count=(%zu,%zu), total=%zu\n",
                           startLatIdx, startLonIdx, nLocalLats, nLocalLons, nLocalVertices );

            // Implement chunked reading to avoid PNetCDF limits
            std::vector< int > landfractVals( nLocalVertices );

            if( nLocalVertices > MAX_CHUNK_SIZE )
            {
                dbgOut.tprintf( 1, "Using chunked reading for landfract data (%zu vertices)\n", nLocalVertices );

                // Read in chunks along latitude dimension
                size_t chunk_lat_size = MAX_CHUNK_SIZE / nLocalLons;
                if( chunk_lat_size == 0 ) chunk_lat_size = 1;

                size_t data_offset = 0;
                for( size_t lat_start = 0; lat_start < nLocalLats; lat_start += chunk_lat_size )
                {
                    size_t lat_count = std::min( chunk_lat_size, nLocalLats - lat_start );
                    size_t chunk_size = lat_count * nLocalLons;

                    NCDF_SIZE start[2] = { static_cast<NCDF_SIZE>(startLatIdx + lat_start), static_cast<NCDF_SIZE>(startLonIdx) };
                    NCDF_SIZE count[2] = { static_cast<NCDF_SIZE>(lat_count), static_cast<NCDF_SIZE>(nLocalLons) };

                    dbgOut.tprintf( 3, "Reading landfract chunk: lat_start=%zu, lat_count=%zu, chunk_size=%zu\n",
                                   lat_start, lat_count, chunk_size );

#ifdef MOAB_HAVE_PNETCDF
                    success = NCFUNCAG( _vara_int )( _fileId, landfractVarId, start, count, &landfractVals[data_offset] );
#else
                    success = NCFUNCAG( _vara_int )( _fileId, landfractVarId, start, count, &landfractVals[data_offset] );
#endif
                    if( success )
                    {
                        dbgOut.tprintf( 1, "Warning: Failed to read landfract chunk, NetCDF error = %d\n", success );
                        break;
                    }
                    data_offset += chunk_size;
                }
            }
            else
            {
                // Read all at once if small enough
                NCDF_SIZE start[2] = { static_cast<NCDF_SIZE>(startLatIdx), static_cast<NCDF_SIZE>(startLonIdx) };
                NCDF_SIZE count[2] = { static_cast<NCDF_SIZE>(nLocalLats), static_cast<NCDF_SIZE>(nLocalLons) };

#ifdef MOAB_HAVE_PNETCDF
                success = NCFUNCAG( _vara_int )( _fileId, landfractVarId, start, count, &landfractVals[0] );
#else
                success = NCFUNCAG( _vara_int )( _fileId, landfractVarId, start, count, &landfractVals[0] );
#endif
            }

            if( !success )
            {
                // Set landfract data on vertices
                rval = mbImpl->tag_set_data( landfractTag, local_verts, &landfractVals[0] );
                MB_CHK_SET_ERR( rval, "Failed to set landfract data on vertices" );
                dbgOut.tprintf( 2, "Successfully read landfract data for %zu vertices\n", nLocalVertices );
            }
            else
            {
                dbgOut.tprintf( 1, "Warning: Failed to read landfract data, NetCDF error = %d\n", success );
            }
        }
    }

    return MB_SUCCESS;
}

ErrorCode NCHelperTOPO::read_ucd_variables_to_nonset_allocate( std::vector< ReadNC::VarData >& vdatas,
                                                               std::vector< int >& /* tstep_nums */ )
{
    // Interface*& mbImpl = _readNC->mbImpl;
    DebugOutput& dbgOut = _readNC->dbgOut;

    ErrorCode rval = MB_SUCCESS;

    for( unsigned int i = 0; i < vdatas.size(); i++ )
    {
        // Skip coordinate variables (lat, lon) as they are handled separately
        if( vdatas[i].varName == latDimName || vdatas[i].varName == lonDimName )
            continue;

        // For USGS topo format, we expect 2D variables with dimensions (lat, lon)
        if( vdatas[i].varDims.size() != 2 )
        {
            dbgOut.tprintf( 1, "Warning: Variable %s has %d dimensions, expected 2 for USGS format\n",
                           vdatas[i].varName.c_str(), (int)vdatas[i].varDims.size() );
            continue;
        }

        // Check that dimensions match lat/lon
        bool isLatLonVar = false;
        for( unsigned int j = 0; j < vdatas[i].varDims.size(); j++ )
        {
            std::string dimName = _readNC->dimNames[vdatas[i].varDims[j]];
            if( dimName == latDimName || dimName == lonDimName )
            {
                isLatLonVar = true;
            }
        }

        if( !isLatLonVar )
        {
            dbgOut.tprintf( 1, "Warning: Variable %s does not have lat/lon dimensions\n",
                           vdatas[i].varName.c_str() );
            continue;
        }

        // Determine parallel reading parameters
        int startLatIdx = 0, endLatIdx = nLatVals - 1;
        int startLonIdx = 0, endLonIdx = nLonVals - 1;

#ifdef MOAB_HAVE_MPI
        bool& isParallel = _readNC->isParallel;
        if( isParallel )
        {
            ParallelComm* myPcomm = _readNC->myPcomm;
            int rank = myPcomm->proc_config().proc_rank();
            int procs = myPcomm->proc_config().proc_size();

            if( procs >= 2 )
            {
                // Use same 1D decomposition as in create_mesh
                int shifted_rank = rank;
                int& trivialPartitionShift = _readNC->trivialPartitionShift;
                if( trivialPartitionShift > 0 ) shifted_rank = ( rank + trivialPartitionShift ) % procs;

                int nLocalLats = int( std::floor( 1.0 * nLatVals / procs ) );
                startLatIdx = shifted_rank * nLocalLats;

                int iextra = nLatVals % procs;
                if( shifted_rank < iextra ) nLocalLats++;
                startLatIdx += std::min( shifted_rank, iextra );

                endLatIdx = startLatIdx + nLocalLats - 1;
            }
        }
#endif

        int nLocalLats = endLatIdx - startLatIdx + 1;
        int nLocalLons = endLonIdx - startLonIdx + 1;
        int nLocalVerts = nLocalLats * nLocalLons;

        dbgOut.tprintf( 1, "Reading variable %s with %d local values\n",
                       vdatas[i].varName.c_str(), nLocalVerts );

        // Allocate data for this variable
        switch( vdatas[i].varDataType )
        {
            case NC_BYTE:
            case NC_CHAR:
                vdatas[i].varDatas.resize( nLocalVerts );
                break;
            case NC_SHORT:
                vdatas[i].varDatas.resize( nLocalVerts * 2 );
                break;
            case NC_INT:
                vdatas[i].varDatas.resize( nLocalVerts * sizeof( int ) );
                break;
            case NC_FLOAT:
                vdatas[i].varDatas.resize( nLocalVerts * sizeof( float ) );
                break;
            case NC_DOUBLE:
                vdatas[i].varDatas.resize( nLocalVerts * sizeof( double ) );
                break;
            default:
                MB_SET_ERR( MB_FAILURE, "Unsupported data type for variable " << vdatas[i].varName );
        }

        // Store reading parameters for later use
        vdatas[i].readStarts.resize( 2 );
        vdatas[i].readCounts.resize( 2 );
        vdatas[i].readStarts[0] = startLatIdx;
        vdatas[i].readStarts[1] = startLonIdx;
        vdatas[i].readCounts[0] = nLocalLats;
        vdatas[i].readCounts[1] = nLocalLons;
    }

    return rval;
}

#ifdef MOAB_HAVE_PNETCDF
ErrorCode NCHelperTOPO::read_ucd_variables_to_nonset_async( std::vector< ReadNC::VarData >& vdatas,
                                                            std::vector< int >& tstep_nums )
{
    DebugOutput& dbgOut = _readNC->dbgOut;
    ErrorCode rval = MB_SUCCESS;

    std::vector< int > requests( vdatas.size() );
    std::vector< int > statuss( vdatas.size() );
    size_t idxReq = 0;

    for( unsigned int i = 0; i < vdatas.size(); i++ )
    {
        // Skip coordinate variables
        if( vdatas[i].varName == latDimName || vdatas[i].varName == lonDimName )
            continue;

        if( vdatas[i].varDatas.empty() )
            continue;

        // For time-independent variables, use first time step
        int tStep = 0;
        if( !tstep_nums.empty() )
            tStep = tstep_nums[0];

        NCDF_SIZE read_starts[3] = { static_cast< NCDF_SIZE >( tStep ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readStarts[0] ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readStarts[1] ) };
        NCDF_SIZE read_counts[3] = { 1,
                                     static_cast< NCDF_SIZE >( vdatas[i].readCounts[0] ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readCounts[1] ) };

        // Adjust for time-independent variables
        NCDF_SIZE* actual_starts = read_starts;
        NCDF_SIZE* actual_counts = read_counts;
        if( vdatas[i].varDims.size() == 2 )
        {
            actual_starts = &read_starts[1];
            actual_counts = &read_counts[1];
        }

        int success;
        switch( vdatas[i].varDataType )
        {
            case NC_BYTE:
            case NC_CHAR:
                success = NCFUNCREQG( _vara_uchar )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                     (unsigned char*)( &vdatas[i].varDatas[0] ), &requests[idxReq++] );
                break;
            case NC_SHORT:
                success = NCFUNCREQG( _vara_short )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                     (short*)( &vdatas[i].varDatas[0] ), &requests[idxReq++] );
                break;
            case NC_INT:
                success = NCFUNCREQG( _vara_int )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                   (int*)( &vdatas[i].varDatas[0] ), &requests[idxReq++] );
                break;
            case NC_FLOAT:
                success = NCFUNCREQG( _vara_float )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                     (float*)( &vdatas[i].varDatas[0] ), &requests[idxReq++] );
                break;
            case NC_DOUBLE:
                success = NCFUNCREQG( _vara_double )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                      (double*)( &vdatas[i].varDatas[0] ), &requests[idxReq++] );
                break;
            default:
                MB_SET_ERR( MB_FAILURE, "Unsupported data type for variable " << vdatas[i].varName );
        }

        if( success )
            MB_SET_ERR( MB_FAILURE, "Failed to read variable " << vdatas[i].varName );
    }

    // Wait for all requests to complete
    if( idxReq > 0 )
    {
        int success = NCFUNC( wait_all )( _fileId, idxReq, &requests[0], &statuss[0] );
        if( success )
            MB_SET_ERR( MB_FAILURE, "Failed on wait_all" );
    }

    dbgOut.tprintf( 1, "Successfully read %d variables using PNetCDF async\n", (int)idxReq );
    return rval;
}
#else
ErrorCode NCHelperTOPO::read_ucd_variables_to_nonset( std::vector< ReadNC::VarData >& vdatas,
                                                      std::vector< int >& tstep_nums )
{
    DebugOutput& dbgOut = _readNC->dbgOut;
    ErrorCode rval = MB_SUCCESS;

    for( unsigned int i = 0; i < vdatas.size(); i++ )
    {
        // Skip coordinate variables
        if( vdatas[i].varName == latDimName || vdatas[i].varName == lonDimName )
            continue;

        if( vdatas[i].varDatas.empty() )
            continue;

        // For time-independent variables, use first time step
        int tStep = 0;
        if( !tstep_nums.empty() )
            tStep = tstep_nums[0];

        NCDF_SIZE read_starts[3] = { static_cast< NCDF_SIZE >( tStep ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readStarts[0] ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readStarts[1] ) };
        NCDF_SIZE read_counts[3] = { 1,
                                     static_cast< NCDF_SIZE >( vdatas[i].readCounts[0] ),
                                     static_cast< NCDF_SIZE >( vdatas[i].readCounts[1] ) };

        // Adjust for time-independent variables
        NCDF_SIZE* actual_starts = read_starts;
        NCDF_SIZE* actual_counts = read_counts;
        if( vdatas[i].varDims.size() == 2 )
        {
            actual_starts = &read_starts[1];
            actual_counts = &read_counts[1];
        }

        int success;
        switch( vdatas[i].varDataType )
        {
            case NC_BYTE:
            case NC_CHAR:
                success = NCFUNCAG( _vara_uchar )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                   (unsigned char*)( &vdatas[i].varDatas[0] ) );
                break;
            case NC_SHORT:
                success = NCFUNCAG( _vara_short )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                   (short*)( &vdatas[i].varDatas[0] ) );
                break;
            case NC_INT:
                success = NCFUNCAG( _vara_int )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                 (int*)( &vdatas[i].varDatas[0] ) );
                break;
            case NC_FLOAT:
                success = NCFUNCAG( _vara_float )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                   (float*)( &vdatas[i].varDatas[0] ) );
                break;
            case NC_DOUBLE:
                success = NCFUNCAG( _vara_double )( _fileId, vdatas[i].varId, actual_starts, actual_counts,
                                                    (double*)( &vdatas[i].varDatas[0] ) );
                break;
            default:
                MB_SET_ERR( MB_FAILURE, "Unsupported data type for variable " << vdatas[i].varName );
        }

        if( success )
            MB_SET_ERR( MB_FAILURE, "Failed to read variable " << vdatas[i].varName );

        dbgOut.tprintf( 1, "Successfully read variable %s\n", vdatas[i].varName.c_str() );
    }

    return rval;
}
#endif

}  // namespace moab
