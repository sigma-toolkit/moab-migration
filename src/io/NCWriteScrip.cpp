//-------------------------------------------------------------------------
// Filename      : NCWriteScrip.cpp
//
// Purpose       : Implementation of the SCRIP grid file writer. See
//                 NCWriteScrip.hpp for design notes.
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#include "NCWriteScrip.hpp"
#include "MBTagConventions.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace moab
{

NCWriteScrip::~NCWriteScrip() {}

namespace
{

constexpr double kPi          = 3.14159265358979323846;
constexpr double kRadToDegree = 180.0 / kPi;

/// Convert Cartesian XYZ (on or near the unit sphere) to lat/lon in
/// degrees. Normalizes internally so callers don't have to. Longitude
/// is shifted to [0, 360).
inline void xyz_to_latlon_deg( double x, double y, double z, double& lat, double& lon )
{
    const double r = std::sqrt( x * x + y * y + z * z );
    if( r < 1.0e-30 )
    {
        lat = 0.0;
        lon = 0.0;
        return;
    }
    lat = std::asin( z / r ) * kRadToDegree;
    lon = std::atan2( y, x ) * kRadToDegree;
    if( lon < 0.0 ) lon += 360.0;
}

}  // namespace

// ============================================================================
// collect_mesh_info: gather local owned cells, compute per-cell center and
// corner lat/lon arrays, find the global max-corners value via Allreduce.
// ============================================================================
ErrorCode NCWriteScrip::collect_mesh_info()
{
    Interface*& mbImpl  = _writeNC->mbImpl;
    Tag& mGlobalIdTag   = _writeNC->mGlobalIdTag;
    DebugOutput& dbgOut = _writeNC->dbgOut;

    // Pick up all 2-D cells (polygonal: tris, quads, polygons).
    Range allCells;
    MB_CHK_SET_ERR( mbImpl->get_entities_by_dimension( _fileSet, 2, allCells ),
                    "Failed to gather 2-D cells for SCRIP write" );
    if( allCells.empty() ) MB_SET_ERR( MB_FAILURE, "No 2-D cells in file set; cannot write SCRIP grid" );

    // In parallel, restrict to owned cells so each rank contributes
    // exactly its share (no double-writing of ghosts).
    localCellsOwned = allCells;
#ifdef MOAB_HAVE_MPI
    bool& isParallel = _writeNC->isParallel;
    if( isParallel )
    {
        ParallelComm*& myPcomm = _writeNC->myPcomm;
        if( myPcomm && myPcomm->proc_config().proc_size() > 1 )
        {
            MB_CHK_SET_ERR( myPcomm->filter_pstatus( localCellsOwned, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                            "Failed to filter owned cells for SCRIP write" );
        }
    }
#endif

    mLocalCells = static_cast< long >( localCellsOwned.size() );

    // Allreduce to find the global maximum corner count so every rank's
    // arrays are sized identically — required for the final write because
    // grid_corner_lon/lat is a 2-D variable.
    int localMaxCorners = 0;
    for( Range::iterator cit = localCellsOwned.begin(); cit != localCellsOwned.end(); ++cit )
    {
        const EntityHandle* conn = nullptr;
        int numConn              = 0;
        ErrorCode rc             = mbImpl->get_connectivity( *cit, conn, numConn );
        if( MB_SUCCESS == rc && numConn > localMaxCorners ) localMaxCorners = numConn;
    }
#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm )
    {
        MPI_Allreduce( &localMaxCorners, &mMaxCornersGlobal, 1, MPI_INT, MPI_MAX,
                       _writeNC->myPcomm->proc_config().proc_comm() );
    }
    else
#endif
    {
        mMaxCornersGlobal = localMaxCorners;
    }
    if( mMaxCornersGlobal <= 0 )
        MB_SET_ERR( MB_FAILURE, "Mesh contains cells but all reported zero connectivity; cannot write SCRIP" );

    // Pre-allocate per-cell arrays.
    mLocalGids.resize( mLocalCells );
    mCenterLon.resize( mLocalCells );
    mCenterLat.resize( mLocalCells );
    mCornerLon.assign( static_cast< size_t >( mLocalCells ) * mMaxCornersGlobal, 0.0 );
    mCornerLat.assign( static_cast< size_t >( mLocalCells ) * mMaxCornersGlobal, 0.0 );
    mImask.assign( mLocalCells, 1 );

    // Optional: pull a GRID_IMASK tag if it exists (matches the reader).
    Tag maskTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "GRID_IMASK", 1, MB_TYPE_INTEGER, maskTag ) && maskTag )
    {
        // Read per-cell mask; if call fails (e.g. tag isn't densely set),
        // just keep the default 1's — the mask field is still well-formed.
        ErrorCode mrc = mbImpl->tag_get_data( maskTag, localCellsOwned, mImask.data() );
        if( MB_SUCCESS != mrc )
        {
            dbgOut.tprintf( 1, "  GRID_IMASK tag exists but read failed; writing all-1 mask\n" );
            std::fill( mImask.begin(), mImask.end(), 1 );
        }
    }

    // Optional: cell area tag.
    Tag areaTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "GRID_AREA", 1, MB_TYPE_DOUBLE, areaTag ) && areaTag )
    {
        mAreas.assign( mLocalCells, 0.0 );
        if( MB_SUCCESS == mbImpl->tag_get_data( areaTag, localCellsOwned, mAreas.data() ) )
            mHasAreas = true;
        else
            mAreas.clear();
    }

    // Fill the global-id, center, and corner arrays.
    long ci = 0;
    for( Range::iterator cit = localCellsOwned.begin(); cit != localCellsOwned.end(); ++cit, ++ci )
    {
        // Global ID for this cell — used to make the on-disk ordering
        // deterministic after a gather to rank 0.
        int gid = 0;
        if( mGlobalIdTag ) mbImpl->tag_get_data( mGlobalIdTag, &( *cit ), 1, &gid );
        mLocalGids[ci] = gid;

        const EntityHandle* conn = nullptr;
        int numConn              = 0;
        MB_CHK_SET_ERR( mbImpl->get_connectivity( *cit, conn, numConn ), "Cell connectivity lookup failed" );

        // Accumulate XYZ for the cell centroid as we go.
        double cx = 0.0, cy = 0.0, cz = 0.0;
        // First-corner fallback used to pad smaller polygons up to the
        // global max corner count (standard SCRIP convention).
        double firstLat = 0.0, firstLon = 0.0;

        for( int k = 0; k < numConn; ++k )
        {
            double vc[3] = { 0.0, 0.0, 0.0 };
            MB_CHK_SET_ERR( mbImpl->get_coords( &conn[k], 1, vc ), "Vertex coordinate lookup failed" );
            cx += vc[0];
            cy += vc[1];
            cz += vc[2];
            double vlat = 0.0, vlon = 0.0;
            xyz_to_latlon_deg( vc[0], vc[1], vc[2], vlat, vlon );
            mCornerLat[ci * mMaxCornersGlobal + k] = vlat;
            mCornerLon[ci * mMaxCornersGlobal + k] = vlon;
            if( k == 0 )
            {
                firstLat = vlat;
                firstLon = vlon;
            }
        }
        // Pad unused corners with the first corner (SCRIP convention).
        for( int k = numConn; k < mMaxCornersGlobal; ++k )
        {
            mCornerLat[ci * mMaxCornersGlobal + k] = firstLat;
            mCornerLon[ci * mMaxCornersGlobal + k] = firstLon;
        }

        // Center: average XYZ projected to the sphere.
        cx /= numConn;
        cy /= numConn;
        cz /= numConn;
        xyz_to_latlon_deg( cx, cy, cz, mCenterLat[ci], mCenterLon[ci] );
    }

    // Global cell count for the file-level grid_size dim.
    mGlobalCells = mLocalCells;
#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm )
    {
        long localN = mLocalCells;
        MPI_Allreduce( &localN, &mGlobalCells, 1, MPI_LONG, MPI_SUM,
                       _writeNC->myPcomm->proc_config().proc_comm() );
    }
#endif

    dbgOut.tprintf( 1, "  SCRIP write: local_cells=%ld global_cells=%ld max_corners=%d hasAreas=%d\n", mLocalCells,
                    mGlobalCells, mMaxCornersGlobal, (int)mHasAreas );

    return MB_SUCCESS;
}

// ============================================================================
// init_file: define dimensions and variables for the SCRIP schema.
// All ranks call collectively — the dispatch layer routes the def_dim /
// def_var calls to the right backend.
// ============================================================================
ErrorCode NCWriteScrip::init_file( std::vector< std::string >& /*var_names*/,
                                   std::vector< std::string >& /*desired_names*/,
                                   bool /*_append*/ )
{
    // Dimensions
    MB_CHK_SET_ERR( ( NCFUNC( def_dim )( _fileId, "grid_size", static_cast< size_t >( mGlobalCells ), &mDimGridSize ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_size dimension" );
    MB_CHK_SET_ERR(
        ( NCFUNC( def_dim )( _fileId, "grid_corners", static_cast< size_t >( mMaxCornersGlobal ), &mDimGridCorners ) )
                ? MB_FAILURE
                : MB_SUCCESS,
        "Failed to define grid_corners dimension" );
    MB_CHK_SET_ERR( ( NCFUNC( def_dim )( _fileId, "grid_rank", static_cast< size_t >( 1 ), &mDimGridRank ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_rank dimension" );

    int dim1[1]  = { mDimGridRank };
    int dim1c[1] = { mDimGridSize };
    int dim2[2]  = { mDimGridSize, mDimGridCorners };

    // Variables — match NCHelperScrip's reader expectations.
    MB_CHK_SET_ERR(
        ( NCFUNC( def_var )( _fileId, "grid_dims", NC_INT, 1, dim1, &mVarGridDims ) ) ? MB_FAILURE : MB_SUCCESS,
        "Failed to define grid_dims variable" );

    MB_CHK_SET_ERR( ( NCFUNC( def_var )( _fileId, "grid_center_lat", NC_DOUBLE, 1, dim1c, &mVarCenterLat ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_center_lat variable" );
    MB_CHK_SET_ERR( ( NCFUNC( def_var )( _fileId, "grid_center_lon", NC_DOUBLE, 1, dim1c, &mVarCenterLon ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_center_lon variable" );
    MB_CHK_SET_ERR( ( NCFUNC( def_var )( _fileId, "grid_corner_lat", NC_DOUBLE, 2, dim2, &mVarCornerLat ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_corner_lat variable" );
    MB_CHK_SET_ERR( ( NCFUNC( def_var )( _fileId, "grid_corner_lon", NC_DOUBLE, 2, dim2, &mVarCornerLon ) )
                            ? MB_FAILURE
                            : MB_SUCCESS,
                    "Failed to define grid_corner_lon variable" );
    MB_CHK_SET_ERR(
        ( NCFUNC( def_var )( _fileId, "grid_imask", NC_INT, 1, dim1c, &mVarImask ) ) ? MB_FAILURE : MB_SUCCESS,
        "Failed to define grid_imask variable" );

    if( mHasAreas )
    {
        MB_CHK_SET_ERR(
            ( NCFUNC( def_var )( _fileId, "grid_area", NC_DOUBLE, 1, dim1c, &mVarArea ) ) ? MB_FAILURE : MB_SUCCESS,
            "Failed to define grid_area variable" );
    }

    // Unit attributes — readers (including NCHelperScrip) key off these.
    const char* deg = "degrees";
    NCFUNC( put_att_text )( _fileId, mVarCenterLat, "units", 7, deg );
    NCFUNC( put_att_text )( _fileId, mVarCenterLon, "units", 7, deg );
    NCFUNC( put_att_text )( _fileId, mVarCornerLat, "units", 7, deg );
    NCFUNC( put_att_text )( _fileId, mVarCornerLon, "units", 7, deg );
    if( mHasAreas )
    {
        const char* rad2 = "radians^2";
        NCFUNC( put_att_text )( _fileId, mVarArea, "units", 9, rad2 );
    }

    // Global title attribute mirrors what the reader looks for in a few places.
    const char* title = "MOAB:NCWriteScrip generated SCRIP grid file";
    NCFUNC( put_att_text )( _fileId, NC_GLOBAL, "title", std::strlen( title ), title );

    // Close define mode so subsequent writes are data-mode operations.
    int endrc = NCFUNC( enddef )( _fileId );
    if( endrc ) MB_SET_ERR( MB_FAILURE, "enddef failed for SCRIP write" );

    return MB_SUCCESS;
}

// ============================================================================
// write_values: gather per-rank arrays to rank 0, sort by global id for
// deterministic ordering, then have rank 0 write the entire SCRIP grid
// via the dispatch layer.
//
// SCRIP grids are small enough (cells, not climate fields) that the
// gather-and-write pattern is fine; full parallel collective writes are
// a future enhancement if the file size ever grows large enough to matter.
// ============================================================================
ErrorCode NCWriteScrip::write_values( std::vector< std::string >& /*var_names*/,
                                      std::vector< int >& /*tstep_nums*/ )
{
    DebugOutput& dbgOut = _writeNC->dbgOut;

    // Helper to assemble the global-ordered arrays on rank 0.
    auto sortAndWrite = [&]( const std::vector< int >& allGids, const std::vector< double >& allCenterLat,
                             const std::vector< double >& allCenterLon, const std::vector< double >& allCornerLat,
                             const std::vector< double >& allCornerLon, const std::vector< int >& allImask,
                             const std::vector< double >& allAreas ) -> ErrorCode {
        const long N    = static_cast< long >( allGids.size() );
        const int  ncpc = mMaxCornersGlobal;

        // Permutation by global id for BfB output across rank counts.
        std::vector< long > perm( N );
        for( long i = 0; i < N; ++i )
            perm[i] = i;
        std::sort( perm.begin(), perm.end(),
                   [&]( long a, long b ) { return allGids[a] < allGids[b]; } );

        std::vector< double > sCenterLat( N ), sCenterLon( N ), sCornerLat( N * ncpc ), sCornerLon( N * ncpc );
        std::vector< int >    sImask( N );
        std::vector< double > sAreas;
        if( !allAreas.empty() ) sAreas.resize( N );

        for( long i = 0; i < N; ++i )
        {
            const long src = perm[i];
            sCenterLat[i]  = allCenterLat[src];
            sCenterLon[i]  = allCenterLon[src];
            sImask[i]      = allImask[src];
            if( !sAreas.empty() ) sAreas[i] = allAreas[src];
            for( int k = 0; k < ncpc; ++k )
            {
                sCornerLat[i * ncpc + k] = allCornerLat[src * ncpc + k];
                sCornerLon[i * ncpc + k] = allCornerLon[src * ncpc + k];
            }
        }

        // grid_dims = [grid_size]
        int gridDimsValue = static_cast< int >( N );
        size_t startD = 0, countD = 1;
        if( NCFUNCAP( _vara_int )( _fileId, mVarGridDims, &startD, &countD, &gridDimsValue ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_dims" );

        size_t s1 = 0, c1 = static_cast< size_t >( N );
        if( NCFUNCAP( _vara_double )( _fileId, mVarCenterLat, &s1, &c1, sCenterLat.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_center_lat" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarCenterLon, &s1, &c1, sCenterLon.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_center_lon" );
        if( NCFUNCAP( _vara_int )( _fileId, mVarImask, &s1, &c1, sImask.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_imask" );
        if( !sAreas.empty() )
        {
            if( NCFUNCAP( _vara_double )( _fileId, mVarArea, &s1, &c1, sAreas.data() ) )
                MB_SET_ERR( MB_FAILURE, "Failed to write grid_area" );
        }

        const size_t s2[2] = { 0, 0 };
        const size_t c2[2] = { static_cast< size_t >( N ), static_cast< size_t >( ncpc ) };
        if( NCFUNCAP( _vara_double )( _fileId, mVarCornerLat, s2, c2, sCornerLat.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_corner_lat" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarCornerLon, s2, c2, sCornerLon.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write grid_corner_lon" );

        return MB_SUCCESS;
    };

#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm && _writeNC->myPcomm->proc_config().proc_size() > 1 )
    {
        MPI_Comm comm = _writeNC->myPcomm->proc_config().proc_comm();
        const int rank = _writeNC->myPcomm->proc_config().proc_rank();
        const int size = _writeNC->myPcomm->proc_config().proc_size();
        const int ncpc = mMaxCornersGlobal;

        // Counts of cells per rank (Gatherv setup).
        int myN = static_cast< int >( mLocalCells );
        std::vector< int > counts( size, 0 ), displs( size, 0 );
        MPI_Gather( &myN, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm );
        if( rank == 0 )
        {
            int acc = 0;
            for( int i = 0; i < size; ++i )
            {
                displs[i] = acc;
                acc += counts[i];
            }
        }

        // Per-cell gather targets on rank 0.
        std::vector< int >    allGids;
        std::vector< double > allCenterLat, allCenterLon, allCornerLat, allCornerLon;
        std::vector< int >    allImask;
        std::vector< double > allAreas;
        if( rank == 0 )
        {
            allGids.resize( mGlobalCells );
            allCenterLat.resize( mGlobalCells );
            allCenterLon.resize( mGlobalCells );
            allCornerLat.resize( static_cast< size_t >( mGlobalCells ) * ncpc );
            allCornerLon.resize( static_cast< size_t >( mGlobalCells ) * ncpc );
            allImask.resize( mGlobalCells );
            if( mHasAreas ) allAreas.resize( mGlobalCells );
        }

        MPI_Gatherv( mLocalGids.data(), myN, MPI_INT, allGids.data(), counts.data(), displs.data(), MPI_INT, 0, comm );
        MPI_Gatherv( mCenterLat.data(), myN, MPI_DOUBLE, allCenterLat.data(), counts.data(), displs.data(), MPI_DOUBLE,
                     0, comm );
        MPI_Gatherv( mCenterLon.data(), myN, MPI_DOUBLE, allCenterLon.data(), counts.data(), displs.data(), MPI_DOUBLE,
                     0, comm );
        MPI_Gatherv( mImask.data(), myN, MPI_INT, allImask.data(), counts.data(), displs.data(), MPI_INT, 0, comm );

        // Corner arrays have a different per-rank count (myN * ncpc), so
        // build a second counts/displs.
        std::vector< int > countsC( size, 0 ), displsC( size, 0 );
        if( rank == 0 )
        {
            int acc = 0;
            for( int i = 0; i < size; ++i )
            {
                countsC[i] = counts[i] * ncpc;
                displsC[i] = acc;
                acc += countsC[i];
            }
        }
        MPI_Gatherv( mCornerLat.data(), myN * ncpc, MPI_DOUBLE, allCornerLat.data(), countsC.data(), displsC.data(),
                     MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mCornerLon.data(), myN * ncpc, MPI_DOUBLE, allCornerLon.data(), countsC.data(), displsC.data(),
                     MPI_DOUBLE, 0, comm );
        if( mHasAreas )
        {
            MPI_Gatherv( mAreas.data(), myN, MPI_DOUBLE, allAreas.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
                         comm );
        }

        // Rank 0 writes. Other ranks no-op (their file handle is still
        // open via the dispatch layer; close happens later in WriteNC).
        ErrorCode rc = MB_SUCCESS;
        if( rank == 0 )
        {
            rc = sortAndWrite( allGids, allCenterLat, allCenterLon, allCornerLat, allCornerLon, allImask, allAreas );
        }
        int rcInt = static_cast< int >( rc );
        MPI_Bcast( &rcInt, 1, MPI_INT, 0, comm );
        if( rcInt != MB_SUCCESS ) MB_SET_ERR( MB_FAILURE, "SCRIP write failed on rank 0" );

        dbgOut.tprintf( 1, "  SCRIP write: gathered+wrote %ld cells from %d ranks\n", mGlobalCells, size );
        return MB_SUCCESS;
    }
#endif

    // Serial path: rank-0-only arrays are already complete; sort + write.
    return sortAndWrite( mLocalGids, mCenterLat, mCenterLon, mCornerLat, mCornerLon, mImask, mAreas );
}

// SCRIP grid files don't have time-varying "nonset" variables; the abstract
// base requires us to implement this hook even though it has nothing to do.
ErrorCode NCWriteScrip::write_nonset_variables( std::vector< WriteNC::VarData >& /*vdatas*/,
                                                std::vector< int >& /*tstep_nums*/ )
{
    return MB_SUCCESS;
}

}  // namespace moab
