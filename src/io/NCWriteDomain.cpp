//-------------------------------------------------------------------------
// Filename      : NCWriteDomain.cpp
//
// Purpose       : CESM domain file writer. See NCWriteDomain.hpp for design.
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#include "NCWriteDomain.hpp"
#include "MBTagConventions.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace moab
{

NCWriteDomain::~NCWriteDomain() {}

namespace
{

constexpr double kPi          = 3.14159265358979323846;
constexpr double kRadToDegree = 180.0 / kPi;

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
// collect_mesh_info — gather owned cells, compute center + corner lat/lon
// per cell. (Domain is a per-cell layout; no vertex dedup needed, unlike
// NCWriteESMF.)
// ============================================================================
ErrorCode NCWriteDomain::collect_mesh_info()
{
    Interface*& mbImpl  = _writeNC->mbImpl;
    Tag& mGlobalIdTag   = _writeNC->mGlobalIdTag;
    DebugOutput& dbgOut = _writeNC->dbgOut;

    Range allCells;
    MB_CHK_SET_ERR( mbImpl->get_entities_by_dimension( _fileSet, 2, allCells ),
                    "Failed to gather 2-D cells for Domain write" );
    if( allCells.empty() ) MB_SET_ERR( MB_FAILURE, "No 2-D cells in file set; cannot write Domain grid" );

    localCellsOwned = allCells;
#ifdef MOAB_HAVE_MPI
    bool& isParallel = _writeNC->isParallel;
    if( isParallel )
    {
        ParallelComm*& myPcomm = _writeNC->myPcomm;
        if( myPcomm && myPcomm->proc_config().proc_size() > 1 )
        {
            MB_CHK_SET_ERR( myPcomm->filter_pstatus( localCellsOwned, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                            "Failed to filter owned cells for Domain write" );
        }
    }
#endif

    mLocalCells = static_cast< long >( localCellsOwned.size() );

    // Global max corners per cell.
    int localMaxCorners = 0;
    for( Range::iterator cit = localCellsOwned.begin(); cit != localCellsOwned.end(); ++cit )
    {
        const EntityHandle* conn = nullptr;
        int numConn              = 0;
        if( MB_SUCCESS == mbImpl->get_connectivity( *cit, conn, numConn ) && numConn > localMaxCorners )
            localMaxCorners = numConn;
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
    if( mMaxCornersGlobal <= 0 ) MB_SET_ERR( MB_FAILURE, "Cells reported zero connectivity; cannot write Domain" );

    mLocalGids.assign( mLocalCells, 0 );
    mXc.assign( mLocalCells, 0.0 );
    mYc.assign( mLocalCells, 0.0 );
    mXv.assign( static_cast< size_t >( mLocalCells ) * mMaxCornersGlobal, 0.0 );
    mYv.assign( static_cast< size_t >( mLocalCells ) * mMaxCornersGlobal, 0.0 );
    mMask.assign( mLocalCells, 1 );

    long ci = 0;
    for( Range::iterator cit = localCellsOwned.begin(); cit != localCellsOwned.end(); ++cit, ++ci )
    {
        int gid = 0;
        if( mGlobalIdTag ) mbImpl->tag_get_data( mGlobalIdTag, &( *cit ), 1, &gid );
        mLocalGids[ci] = gid;

        const EntityHandle* conn = nullptr;
        int numConn              = 0;
        MB_CHK_SET_ERR( mbImpl->get_connectivity( *cit, conn, numConn ), "Cell connectivity lookup failed" );

        double cx = 0.0, cy = 0.0, cz = 0.0;
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
            mYv[ci * mMaxCornersGlobal + k] = vlat;
            mXv[ci * mMaxCornersGlobal + k] = vlon;
            if( k == 0 )
            {
                firstLat = vlat;
                firstLon = vlon;
            }
        }
        for( int k = numConn; k < mMaxCornersGlobal; ++k )
        {
            mYv[ci * mMaxCornersGlobal + k] = firstLat;
            mXv[ci * mMaxCornersGlobal + k] = firstLon;
        }
        cx /= numConn;
        cy /= numConn;
        cz /= numConn;
        xyz_to_latlon_deg( cx, cy, cz, mYc[ci], mXc[ci] );
    }

    // Optional input tags
    Tag maskTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "DOMAIN_MASK", 1, MB_TYPE_INTEGER, maskTag ) && maskTag )
    {
        if( MB_SUCCESS == mbImpl->tag_get_data( maskTag, localCellsOwned, mMask.data() ) ) mHasMask = true;
    }
    // mask is required by the schema; if no tag found we still emit 1's.
    mHasMask = true;  // we always write mask

    Tag areaTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "GRID_AREA", 1, MB_TYPE_DOUBLE, areaTag ) && areaTag )
    {
        mAreas.assign( mLocalCells, 0.0 );
        if( MB_SUCCESS == mbImpl->tag_get_data( areaTag, localCellsOwned, mAreas.data() ) )
            mHasAreas = true;
        else
            mAreas.clear();
    }
    Tag fracTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "DOMAIN_FRAC", 1, MB_TYPE_DOUBLE, fracTag ) && fracTag )
    {
        mFrac.assign( mLocalCells, 1.0 );
        if( MB_SUCCESS == mbImpl->tag_get_data( fracTag, localCellsOwned, mFrac.data() ) )
            mHasFrac = true;
        else
            mFrac.clear();
    }
    // frac is also conventionally always present in CESM domain files; emit 1.0's by default.
    if( !mHasFrac )
    {
        mFrac.assign( mLocalCells, 1.0 );
        mHasFrac = true;
    }

    mGlobalCells = mLocalCells;
#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm )
    {
        long localN = mLocalCells;
        MPI_Allreduce( &localN, &mGlobalCells, 1, MPI_LONG, MPI_SUM,
                       _writeNC->myPcomm->proc_config().proc_comm() );
    }
#endif
    // Flat layout: nj=1, ni=total cells. Matches the CESM convention for
    // domain files derived from unstructured meshes (see e.g.
    // domain.ocn.ne4np4_oQU240 ships with ni=866, nj=1).
    mNj = 1;
    mNi = mGlobalCells;

    dbgOut.tprintf( 1, "  Domain write: local_cells=%ld global_cells=%ld ni=%ld nj=%ld max_corners=%d\n", mLocalCells,
                    mGlobalCells, mNi, mNj, mMaxCornersGlobal );

    return MB_SUCCESS;
}

ErrorCode NCWriteDomain::init_file( std::vector< std::string >& /*var_names*/,
                                    std::vector< std::string >& /*desired_names*/,
                                    bool /*_append*/ )
{
    if( NCFUNC( def_dim )( _fileId, "n", static_cast< size_t >( mGlobalCells ), &mDimN ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define n dim" );
    if( NCFUNC( def_dim )( _fileId, "ni", static_cast< size_t >( mNi ), &mDimNi ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define ni dim" );
    if( NCFUNC( def_dim )( _fileId, "nj", static_cast< size_t >( mNj ), &mDimNj ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define nj dim" );
    if( NCFUNC( def_dim )( _fileId, "nv", static_cast< size_t >( mMaxCornersGlobal ), &mDimNv ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define nv dim" );

    // Variable dims: arrays are (nj, ni) or (nj, ni, nv). j is the slow axis.
    int dimsCenter[2] = { mDimNj, mDimNi };
    int dimsVertex[3] = { mDimNj, mDimNi, mDimNv };

    if( NCFUNC( def_var )( _fileId, "xc", NC_DOUBLE, 2, dimsCenter, &mVarXc ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define xc var" );
    if( NCFUNC( def_var )( _fileId, "yc", NC_DOUBLE, 2, dimsCenter, &mVarYc ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define yc var" );
    if( NCFUNC( def_var )( _fileId, "xv", NC_DOUBLE, 3, dimsVertex, &mVarXv ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define xv var" );
    if( NCFUNC( def_var )( _fileId, "yv", NC_DOUBLE, 3, dimsVertex, &mVarYv ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define yv var" );
    if( NCFUNC( def_var )( _fileId, "mask", NC_INT, 2, dimsCenter, &mVarMask ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define mask var" );
    if( NCFUNC( def_var )( _fileId, "area", NC_DOUBLE, 2, dimsCenter, &mVarArea ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define area var" );
    if( NCFUNC( def_var )( _fileId, "frac", NC_DOUBLE, 2, dimsCenter, &mVarFrac ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define frac var" );

    // CF-style metadata matching the canonical CESM domain file schema.
    const char* deg_e   = "degrees_east";
    const char* deg_n   = "degrees_north";
    const char* rad2    = "radian2";
    const char* xv_name = "longitude of grid cell verticies";
    const char* yv_name = "latitude of grid cell verticies";
    const char* xc_name = "longitude of grid cell center";
    const char* yc_name = "latitude of grid cell center";
    const char* m_name  = "domain mask";
    const char* a_name  = "area of grid cell in radians squared";
    const char* f_name  = "fraction of grid cell that is active";

    NCFUNC( put_att_text )( _fileId, mVarXc, "long_name", std::strlen( xc_name ), xc_name );
    NCFUNC( put_att_text )( _fileId, mVarXc, "units", std::strlen( deg_e ), deg_e );
    NCFUNC( put_att_text )( _fileId, mVarXc, "bounds", 2, "xv" );

    NCFUNC( put_att_text )( _fileId, mVarYc, "long_name", std::strlen( yc_name ), yc_name );
    NCFUNC( put_att_text )( _fileId, mVarYc, "units", std::strlen( deg_n ), deg_n );
    NCFUNC( put_att_text )( _fileId, mVarYc, "bounds", 2, "yv" );

    NCFUNC( put_att_text )( _fileId, mVarXv, "long_name", std::strlen( xv_name ), xv_name );
    NCFUNC( put_att_text )( _fileId, mVarXv, "units", std::strlen( deg_e ), deg_e );

    NCFUNC( put_att_text )( _fileId, mVarYv, "long_name", std::strlen( yv_name ), yv_name );
    NCFUNC( put_att_text )( _fileId, mVarYv, "units", std::strlen( deg_n ), deg_n );

    NCFUNC( put_att_text )( _fileId, mVarMask, "long_name", std::strlen( m_name ), m_name );
    NCFUNC( put_att_text )( _fileId, mVarMask, "note", 8, "unitless" );
    NCFUNC( put_att_text )( _fileId, mVarMask, "coordinates", 5, "xc yc" );
    NCFUNC( put_att_text )( _fileId, mVarMask, "comment", 36, "0 value indicates cell is not active" );

    NCFUNC( put_att_text )( _fileId, mVarArea, "long_name", std::strlen( a_name ), a_name );
    NCFUNC( put_att_text )( _fileId, mVarArea, "units", std::strlen( rad2 ), rad2 );
    NCFUNC( put_att_text )( _fileId, mVarArea, "coordinates", 5, "xc yc" );

    NCFUNC( put_att_text )( _fileId, mVarFrac, "long_name", std::strlen( f_name ), f_name );
    NCFUNC( put_att_text )( _fileId, mVarFrac, "coordinates", 5, "xc yc" );
    NCFUNC( put_att_text )( _fileId, mVarFrac, "note", 8, "unitless" );

    const char* title = "MOAB:NCWriteDomain generated CESM-style domain file";
    NCFUNC( put_att_text )( _fileId, NC_GLOBAL, "title", std::strlen( title ), title );
    const char* conv = "CF-1.0";
    NCFUNC( put_att_text )( _fileId, NC_GLOBAL, "Conventions", std::strlen( conv ), conv );

    if( NCFUNC( enddef )( _fileId ) ) MB_SET_ERR( MB_FAILURE, "enddef failed for Domain write" );

    return MB_SUCCESS;
}

ErrorCode NCWriteDomain::write_values( std::vector< std::string >& /*var_names*/,
                                       std::vector< int >& /*tstep_nums*/ )
{
    DebugOutput& dbgOut = _writeNC->dbgOut;

    auto writeAll = [&]( const std::vector< int >& allGids, const std::vector< double >& allXc,
                         const std::vector< double >& allYc, const std::vector< double >& allXv,
                         const std::vector< double >& allYv, const std::vector< int >& allMask,
                         const std::vector< double >& allAreas, const std::vector< double >& allFrac ) -> ErrorCode {
        const long N    = static_cast< long >( allGids.size() );
        const int  ncpc = mMaxCornersGlobal;

        // Sort by global id for BfB output across rank counts.
        std::vector< long > perm( N );
        for( long i = 0; i < N; ++i )
            perm[i] = i;
        std::sort( perm.begin(), perm.end(), [&]( long a, long b ) { return allGids[a] < allGids[b]; } );

        std::vector< double > sXc( N ), sYc( N ), sXv( N * ncpc ), sYv( N * ncpc );
        std::vector< int >    sMask( N );
        std::vector< double > sArea( N, 0.0 ), sFrac( N, 1.0 );
        const bool            haveArea = !allAreas.empty();
        for( long i = 0; i < N; ++i )
        {
            const long src = perm[i];
            sXc[i]   = allXc[src];
            sYc[i]   = allYc[src];
            sMask[i] = allMask[src];
            if( haveArea ) sArea[i] = allAreas[src];
            if( !allFrac.empty() ) sFrac[i] = allFrac[src];
            for( int k = 0; k < ncpc; ++k )
            {
                sXv[i * ncpc + k] = allXv[src * ncpc + k];
                sYv[i * ncpc + k] = allYv[src * ncpc + k];
            }
        }

        // Writes are 2-D (nj=1, ni=N) for centers and mask/area/frac;
        // 3-D (nj=1, ni=N, nv=ncpc) for the vertex arrays.
        size_t s2[2] = { 0, 0 };
        size_t c2[2] = { 1, static_cast< size_t >( N ) };
        size_t s3[3] = { 0, 0, 0 };
        size_t c3[3] = { 1, static_cast< size_t >( N ), static_cast< size_t >( ncpc ) };

        if( NCFUNCAP( _vara_double )( _fileId, mVarXc, s2, c2, sXc.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write xc" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarYc, s2, c2, sYc.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write yc" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarXv, s3, c3, sXv.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write xv" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarYv, s3, c3, sYv.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write yv" );
        if( NCFUNCAP( _vara_int )( _fileId, mVarMask, s2, c2, sMask.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write mask" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarArea, s2, c2, sArea.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write area" );
        if( NCFUNCAP( _vara_double )( _fileId, mVarFrac, s2, c2, sFrac.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write frac" );
        return MB_SUCCESS;
    };

#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm && _writeNC->myPcomm->proc_config().proc_size() > 1 )
    {
        MPI_Comm comm = _writeNC->myPcomm->proc_config().proc_comm();
        const int rank = _writeNC->myPcomm->proc_config().proc_rank();
        const int size = _writeNC->myPcomm->proc_config().proc_size();
        const int ncpc = mMaxCornersGlobal;

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

        std::vector< int >    allGids;
        std::vector< double > allXc, allYc, allXv, allYv;
        std::vector< int >    allMask;
        std::vector< double > allAreas, allFrac;
        if( rank == 0 )
        {
            allGids.resize( mGlobalCells );
            allXc.resize( mGlobalCells );
            allYc.resize( mGlobalCells );
            allXv.resize( static_cast< size_t >( mGlobalCells ) * ncpc );
            allYv.resize( static_cast< size_t >( mGlobalCells ) * ncpc );
            allMask.resize( mGlobalCells );
            if( mHasAreas ) allAreas.resize( mGlobalCells );
            allFrac.resize( mGlobalCells );
        }

        MPI_Gatherv( mLocalGids.data(), myN, MPI_INT, allGids.data(), counts.data(), displs.data(), MPI_INT, 0, comm );
        MPI_Gatherv( mXc.data(), myN, MPI_DOUBLE, allXc.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mYc.data(), myN, MPI_DOUBLE, allYc.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mMask.data(), myN, MPI_INT, allMask.data(), counts.data(), displs.data(), MPI_INT, 0, comm );
        if( mHasAreas )
            MPI_Gatherv( mAreas.data(), myN, MPI_DOUBLE, allAreas.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
                         comm );
        MPI_Gatherv( mFrac.data(), myN, MPI_DOUBLE, allFrac.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, comm );

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
        MPI_Gatherv( mXv.data(), myN * ncpc, MPI_DOUBLE, allXv.data(), countsC.data(), displsC.data(), MPI_DOUBLE, 0,
                     comm );
        MPI_Gatherv( mYv.data(), myN * ncpc, MPI_DOUBLE, allYv.data(), countsC.data(), displsC.data(), MPI_DOUBLE, 0,
                     comm );

        ErrorCode rc = MB_SUCCESS;
        if( rank == 0 ) rc = writeAll( allGids, allXc, allYc, allXv, allYv, allMask, allAreas, allFrac );
        int rcInt = static_cast< int >( rc );
        MPI_Bcast( &rcInt, 1, MPI_INT, 0, comm );
        if( rcInt != MB_SUCCESS ) MB_SET_ERR( MB_FAILURE, "Domain write failed on rank 0" );

        dbgOut.tprintf( 1, "  Domain write: gathered+wrote %ld cells from %d ranks\n", mGlobalCells, size );
        return MB_SUCCESS;
    }
#endif

    return writeAll( mLocalGids, mXc, mYc, mXv, mYv, mMask, mAreas, mFrac );
}

ErrorCode NCWriteDomain::write_nonset_variables( std::vector< WriteNC::VarData >& /*vdatas*/,
                                                 std::vector< int >& /*tstep_nums*/ )
{
    return MB_SUCCESS;
}

}  // namespace moab
