//-------------------------------------------------------------------------
// Filename      : NCWriteESMF.cpp
//
// Purpose       : ESMF unstructured grid file writer. See NCWriteESMF.hpp
//                 for design notes.
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#include "NCWriteESMF.hpp"
#include "MBTagConventions.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <set>

namespace moab
{

NCWriteESMF::~NCWriteESMF() {}

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
// collect_mesh_info — gather owned cells; compute per-rank vertex lat/lon
// (de-duplicated within the rank), per-cell global-id connectivity, and
// center lat/lon.
// ============================================================================
ErrorCode NCWriteESMF::collect_mesh_info()
{
    Interface*& mbImpl  = _writeNC->mbImpl;
    Tag& mGlobalIdTag   = _writeNC->mGlobalIdTag;
    DebugOutput& dbgOut = _writeNC->dbgOut;

    Range allCells;
    MB_CHK_SET_ERR( mbImpl->get_entities_by_dimension( _fileSet, 2, allCells ),
                    "Failed to gather 2-D cells for ESMF write" );
    if( allCells.empty() ) MB_SET_ERR( MB_FAILURE, "No 2-D cells in file set; cannot write ESMF grid" );

    localCellsOwned = allCells;
#ifdef MOAB_HAVE_MPI
    bool& isParallel = _writeNC->isParallel;
    if( isParallel )
    {
        ParallelComm*& myPcomm = _writeNC->myPcomm;
        if( myPcomm && myPcomm->proc_config().proc_size() > 1 )
        {
            MB_CHK_SET_ERR( myPcomm->filter_pstatus( localCellsOwned, PSTATUS_NOT_OWNED, PSTATUS_NOT ),
                            "Failed to filter owned cells for ESMF write" );
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
    if( mMaxCornersGlobal <= 0 ) MB_SET_ERR( MB_FAILURE, "Cells reported zero connectivity; cannot write ESMF" );

    // First pass: collect the set of vertex handles referenced by owned
    // cells, plus per-cell (numNodes, connectivity-by-EntityHandle).
    std::set< EntityHandle > usedVerts;
    std::vector< std::vector< EntityHandle > > perCellConn( mLocalCells );
    mElementNumNodes.assign( mLocalCells, 0 );
    mLocalCellGids.assign( mLocalCells, 0 );
    mCenterLat.assign( mLocalCells, 0.0 );
    mCenterLon.assign( mLocalCells, 0.0 );

    long ci = 0;
    for( Range::iterator cit = localCellsOwned.begin(); cit != localCellsOwned.end(); ++cit, ++ci )
    {
        int gid = 0;
        if( mGlobalIdTag ) mbImpl->tag_get_data( mGlobalIdTag, &( *cit ), 1, &gid );
        mLocalCellGids[ci] = gid;

        const EntityHandle* conn = nullptr;
        int numConn              = 0;
        MB_CHK_SET_ERR( mbImpl->get_connectivity( *cit, conn, numConn ), "Cell connectivity lookup failed" );
        mElementNumNodes[ci] = numConn;
        perCellConn[ci].assign( conn, conn + numConn );

        // Compute the per-cell center directly from vertex coords.
        double cx = 0.0, cy = 0.0, cz = 0.0;
        for( int k = 0; k < numConn; ++k )
        {
            usedVerts.insert( conn[k] );
            double vc[3] = { 0.0, 0.0, 0.0 };
            MB_CHK_SET_ERR( mbImpl->get_coords( &conn[k], 1, vc ), "Vertex coordinate lookup failed" );
            cx += vc[0];
            cy += vc[1];
            cz += vc[2];
        }
        cx /= numConn;
        cy /= numConn;
        cz /= numConn;
        xyz_to_latlon_deg( cx, cy, cz, mCenterLat[ci], mCenterLon[ci] );
    }

    // Build per-rank node arrays (vertex GID + lat/lon).
    const size_t localNodeCount = usedVerts.size();
    mNodeGids.assign( localNodeCount, 0 );
    mNodeLon.assign( localNodeCount, 0.0 );
    mNodeLat.assign( localNodeCount, 0.0 );
    {
        size_t ni = 0;
        for( EntityHandle v : usedVerts )
        {
            int vgid = 0;
            if( mGlobalIdTag ) mbImpl->tag_get_data( mGlobalIdTag, &v, 1, &vgid );
            mNodeGids[ni] = vgid;
            double vc[3]  = { 0.0, 0.0, 0.0 };
            mbImpl->get_coords( &v, 1, vc );
            xyz_to_latlon_deg( vc[0], vc[1], vc[2], mNodeLat[ni], mNodeLon[ni] );
            ++ni;
        }
    }

    // Build per-rank connectivity as vertex GIDs (1-based index will be
    // assigned later at rank 0 after global dedup).
    mElementConn.assign( static_cast< size_t >( mLocalCells ) * mMaxCornersGlobal, 0 );
    for( long c = 0; c < mLocalCells; ++c )
    {
        const std::vector< EntityHandle >& cv = perCellConn[c];
        int firstGid                          = 0;
        for( int k = 0; k < (int)cv.size(); ++k )
        {
            int vgid = 0;
            if( mGlobalIdTag ) mbImpl->tag_get_data( mGlobalIdTag, &cv[k], 1, &vgid );
            mElementConn[c * mMaxCornersGlobal + k] = vgid;
            if( k == 0 ) firstGid = vgid;
        }
        // Pad unused slots with first GID (writers convention; the
        // numElementConn field tells the reader the real per-cell count).
        for( int k = (int)cv.size(); k < mMaxCornersGlobal; ++k )
            mElementConn[c * mMaxCornersGlobal + k] = firstGid;
    }

    // Optional masks / areas from familiar tags (matches what NCHelperESMF
    // would produce on read).
    Tag maskTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "ELEMENT_MASK", 1, MB_TYPE_INTEGER, maskTag ) && maskTag )
    {
        mMask.assign( mLocalCells, 1 );
        if( MB_SUCCESS == mbImpl->tag_get_data( maskTag, localCellsOwned, mMask.data() ) )
            mHasMask = true;
        else
            mMask.clear();
    }
    Tag areaTag = 0;
    if( MB_SUCCESS == mbImpl->tag_get_handle( "GRID_AREA", 1, MB_TYPE_DOUBLE, areaTag ) && areaTag )
    {
        mAreas.assign( mLocalCells, 0.0 );
        if( MB_SUCCESS == mbImpl->tag_get_data( areaTag, localCellsOwned, mAreas.data() ) )
            mHasAreas = true;
        else
            mAreas.clear();
    }

    // Global cell count.
    mGlobalCells = mLocalCells;
#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm )
    {
        long localN = mLocalCells;
        MPI_Allreduce( &localN, &mGlobalCells, 1, MPI_LONG, MPI_SUM,
                       _writeNC->myPcomm->proc_config().proc_comm() );
    }
#endif
    // Global node count is set after rank 0 does the dedup in write_values;
    // here we just stash a lower bound for the init_file dim size.
    mGlobalNodes = static_cast< long >( localNodeCount );
#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm )
    {
        // Worst-case (no sharing): sum local node counts. The reader
        // doesn't care if we declare a slightly larger nodeCount than
        // actually used, but for cleanliness we'll set the dim from the
        // actual deduplicated count once we know it (init_file uses
        // mGlobalNodes; we overwrite it before init_file runs).
        long localN = static_cast< long >( localNodeCount );
        MPI_Allreduce( &localN, &mGlobalNodes, 1, MPI_LONG, MPI_SUM,
                       _writeNC->myPcomm->proc_config().proc_comm() );
    }
#endif

    dbgOut.tprintf( 1, "  ESMF write: local_cells=%ld global_cells=%ld local_nodes=%zu max_corners=%d\n", mLocalCells,
                    mGlobalCells, localNodeCount, mMaxCornersGlobal );

    return MB_SUCCESS;
}

// ============================================================================
// init_file — define the ESMF schema. nodeCount is sized to a safe upper
// bound here (sum of per-rank node counts); write_values may reduce it
// after the rank-0 dedup. NetCDF allows defining a dim of a known size
// and writing fewer entries via the count[] argument.
// ============================================================================
ErrorCode NCWriteESMF::init_file( std::vector< std::string >& /*var_names*/,
                                  std::vector< std::string >& /*desired_names*/,
                                  bool /*_append*/ )
{
    // Dimensions
    if( NCFUNC( def_dim )( _fileId, "nodeCount", static_cast< size_t >( mGlobalNodes ), &mDimNodeCount ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define nodeCount dim" );
    if( NCFUNC( def_dim )( _fileId, "elementCount", static_cast< size_t >( mGlobalCells ), &mDimElementCount ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define elementCount dim" );
    if( NCFUNC( def_dim )( _fileId, "maxNodePElement", static_cast< size_t >( mMaxCornersGlobal ),
                           &mDimMaxNodePElement ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define maxNodePElement dim" );
    if( NCFUNC( def_dim )( _fileId, "coordDim", static_cast< size_t >( mCoordDim ), &mDimCoordDim ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define coordDim dim" );

    int dimsNode[2]   = { mDimNodeCount, mDimCoordDim };
    int dimsElem[2]   = { mDimElementCount, mDimMaxNodePElement };
    int dimsCenter[2] = { mDimElementCount, mDimCoordDim };
    int dimsElem1[1]  = { mDimElementCount };

    if( NCFUNC( def_var )( _fileId, "nodeCoords", NC_DOUBLE, 2, dimsNode, &mVarNodeCoords ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define nodeCoords var" );
    if( NCFUNC( def_var )( _fileId, "elementConn", NC_INT, 2, dimsElem, &mVarElementConn ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define elementConn var" );
    if( NCFUNC( def_var )( _fileId, "numElementConn", NC_INT, 1, dimsElem1, &mVarNumElementConn ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define numElementConn var" );
    if( NCFUNC( def_var )( _fileId, "centerCoords", NC_DOUBLE, 2, dimsCenter, &mVarCenterCoords ) )
        MB_SET_ERR( MB_FAILURE, "Failed to define centerCoords var" );

    if( mHasAreas )
    {
        if( NCFUNC( def_var )( _fileId, "elementArea", NC_DOUBLE, 1, dimsElem1, &mVarElementArea ) )
            MB_SET_ERR( MB_FAILURE, "Failed to define elementArea var" );
    }
    if( mHasMask )
    {
        if( NCFUNC( def_var )( _fileId, "elementMask", NC_INT, 1, dimsElem1, &mVarElementMask ) )
            MB_SET_ERR( MB_FAILURE, "Failed to define elementMask var" );
    }

    const char* deg = "degrees";
    NCFUNC( put_att_text )( _fileId, mVarNodeCoords, "units", 7, deg );
    NCFUNC( put_att_text )( _fileId, mVarCenterCoords, "units", 7, deg );
    if( mHasAreas )
    {
        const char* rad2 = "radians^2";
        NCFUNC( put_att_text )( _fileId, mVarElementArea, "units", 9, rad2 );
    }
    const char* gridtype = "unstructured";
    NCFUNC( put_att_text )( _fileId, NC_GLOBAL, "gridType", std::strlen( gridtype ), gridtype );
    const char* title = "MOAB:NCWriteESMF generated ESMF unstructured grid file";
    NCFUNC( put_att_text )( _fileId, NC_GLOBAL, "title", std::strlen( title ), title );

    if( NCFUNC( enddef )( _fileId ) ) MB_SET_ERR( MB_FAILURE, "enddef failed for ESMF write" );

    return MB_SUCCESS;
}

// ============================================================================
// write_values — gather to rank 0, dedup vertices by global id, write.
//
// The rank-0 dedup uses vertex GIDs (which are partition-invariant) as the
// canonical identity. After dedup we map each cell's vertex-GID array to
// 1-based indices into the unique-vertex list.
// ============================================================================
ErrorCode NCWriteESMF::write_values( std::vector< std::string >& /*var_names*/,
                                     std::vector< int >& /*tstep_nums*/ )
{
    DebugOutput& dbgOut = _writeNC->dbgOut;

    // Helper assembling the rank-0 unified arrays and writing the file.
    auto writeAll = [&]( const std::vector< int >& allCellGids, const std::vector< int >& allNumNodes,
                         const std::vector< int >& allConnGid, const std::vector< double >& allCenterLat,
                         const std::vector< double >& allCenterLon, const std::vector< int >& allNodeGids,
                         const std::vector< double >& allNodeLat, const std::vector< double >& allNodeLon,
                         const std::vector< double >& allAreas, const std::vector< int >& allMask ) -> ErrorCode {
        const int ncpc      = mMaxCornersGlobal;
        const long nCellsG  = static_cast< long >( allCellGids.size() );
        const long nNodesIn = static_cast< long >( allNodeGids.size() );

        // ---- Dedup nodes by GID ----------------------------------------
        // Build a sorted, unique vertex-GID list and a map gid -> 1-based
        // index for the elementConn rewrite.
        std::vector< int > nodeOrder( nNodesIn );
        for( long i = 0; i < nNodesIn; ++i )
            nodeOrder[i] = static_cast< int >( i );
        std::sort( nodeOrder.begin(), nodeOrder.end(),
                   [&]( int a, int b ) { return allNodeGids[a] < allNodeGids[b]; } );

        std::vector< int >    uniqNodeGids;
        std::vector< double > uniqNodeLat, uniqNodeLon;
        std::map< int, int >  gid2Idx;  // global vertex GID -> 1-based unique index
        uniqNodeGids.reserve( nNodesIn );
        uniqNodeLat.reserve( nNodesIn );
        uniqNodeLon.reserve( nNodesIn );
        int prevGid = -1;
        for( long i = 0; i < nNodesIn; ++i )
        {
            const int idx = nodeOrder[i];
            const int g   = allNodeGids[idx];
            if( g != prevGid )
            {
                uniqNodeGids.push_back( g );
                uniqNodeLat.push_back( allNodeLat[idx] );
                uniqNodeLon.push_back( allNodeLon[idx] );
                gid2Idx[g] = static_cast< int >( uniqNodeGids.size() );  // 1-based
                prevGid    = g;
            }
        }
        const long nNodesU = static_cast< long >( uniqNodeGids.size() );

        // ---- Sort cells by GID for deterministic on-disk ordering ------
        std::vector< long > cellOrder( nCellsG );
        for( long i = 0; i < nCellsG; ++i )
            cellOrder[i] = i;
        std::sort( cellOrder.begin(), cellOrder.end(),
                   [&]( long a, long b ) { return allCellGids[a] < allCellGids[b]; } );

        // ---- Build the final writable arrays in cell-sorted order ------
        std::vector< int >    numNodes( nCellsG );
        std::vector< int >    connIdx( static_cast< size_t >( nCellsG ) * ncpc );
        std::vector< double > centerLat( nCellsG ), centerLon( nCellsG );
        std::vector< double > areas;
        std::vector< int >    mask;
        if( !allAreas.empty() ) areas.resize( nCellsG );
        if( !allMask.empty() ) mask.resize( nCellsG );

        for( long i = 0; i < nCellsG; ++i )
        {
            const long src   = cellOrder[i];
            numNodes[i]      = allNumNodes[src];
            centerLat[i]     = allCenterLat[src];
            centerLon[i]     = allCenterLon[src];
            if( !areas.empty() ) areas[i] = allAreas[src];
            if( !mask.empty() ) mask[i] = allMask[src];

            for( int k = 0; k < ncpc; ++k )
            {
                const int g = allConnGid[src * ncpc + k];
                auto       it = gid2Idx.find( g );
                connIdx[i * ncpc + k] = ( it == gid2Idx.end() ) ? 0 : it->second;
            }
        }

        // ---- Interleave node and center coords as [lon, lat] -----------
        std::vector< double > nodeCoords( nNodesU * 2 ), centerCoords( nCellsG * 2 );
        for( long i = 0; i < nNodesU; ++i )
        {
            nodeCoords[i * 2 + 0] = uniqNodeLon[i];
            nodeCoords[i * 2 + 1] = uniqNodeLat[i];
        }
        for( long i = 0; i < nCellsG; ++i )
        {
            centerCoords[i * 2 + 0] = centerLon[i];
            centerCoords[i * 2 + 1] = centerLat[i];
        }

        // ---- Write to disk via the dispatch layer ----------------------
        size_t startN[2] = { 0, 0 };
        size_t countN[2] = { static_cast< size_t >( nNodesU ), static_cast< size_t >( mCoordDim ) };
        if( NCFUNCAP( _vara_double )( _fileId, mVarNodeCoords, startN, countN, nodeCoords.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write nodeCoords" );

        size_t startE[2] = { 0, 0 };
        size_t countE[2] = { static_cast< size_t >( nCellsG ), static_cast< size_t >( ncpc ) };
        if( NCFUNCAP( _vara_int )( _fileId, mVarElementConn, startE, countE, connIdx.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write elementConn" );

        size_t s1 = 0, c1 = static_cast< size_t >( nCellsG );
        if( NCFUNCAP( _vara_int )( _fileId, mVarNumElementConn, &s1, &c1, numNodes.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write numElementConn" );

        size_t startC[2] = { 0, 0 };
        size_t countC[2] = { static_cast< size_t >( nCellsG ), static_cast< size_t >( mCoordDim ) };
        if( NCFUNCAP( _vara_double )( _fileId, mVarCenterCoords, startC, countC, centerCoords.data() ) )
            MB_SET_ERR( MB_FAILURE, "Failed to write centerCoords" );

        if( !areas.empty() )
        {
            if( NCFUNCAP( _vara_double )( _fileId, mVarElementArea, &s1, &c1, areas.data() ) )
                MB_SET_ERR( MB_FAILURE, "Failed to write elementArea" );
        }
        if( !mask.empty() )
        {
            if( NCFUNCAP( _vara_int )( _fileId, mVarElementMask, &s1, &c1, mask.data() ) )
                MB_SET_ERR( MB_FAILURE, "Failed to write elementMask" );
        }

        return MB_SUCCESS;
    };

#ifdef MOAB_HAVE_MPI
    if( _writeNC->isParallel && _writeNC->myPcomm && _writeNC->myPcomm->proc_config().proc_size() > 1 )
    {
        MPI_Comm comm = _writeNC->myPcomm->proc_config().proc_comm();
        const int rank = _writeNC->myPcomm->proc_config().proc_rank();
        const int size = _writeNC->myPcomm->proc_config().proc_size();
        const int ncpc = mMaxCornersGlobal;

        // Per-rank cell counts + displacements for gather.
        int myCells = static_cast< int >( mLocalCells );
        std::vector< int > cellCounts( size, 0 ), cellDispls( size, 0 );
        MPI_Gather( &myCells, 1, MPI_INT, cellCounts.data(), 1, MPI_INT, 0, comm );

        // Per-rank node counts (these vary per rank — each rank reports
        // its own local node count).
        int myNodes = static_cast< int >( mNodeGids.size() );
        std::vector< int > nodeCounts( size, 0 ), nodeDispls( size, 0 );
        MPI_Gather( &myNodes, 1, MPI_INT, nodeCounts.data(), 1, MPI_INT, 0, comm );

        // Rank 0 sets up the displacement arrays.
        std::vector< int > cellCountsC( size, 0 ), cellDisplsC( size, 0 );
        if( rank == 0 )
        {
            int accC = 0, accCC = 0, accN = 0;
            for( int i = 0; i < size; ++i )
            {
                cellDispls[i]  = accC;
                accC += cellCounts[i];
                cellCountsC[i] = cellCounts[i] * ncpc;
                cellDisplsC[i] = accCC;
                accCC += cellCountsC[i];
                nodeDispls[i]  = accN;
                accN += nodeCounts[i];
            }
        }

        long nCellsGlobal = mGlobalCells;
        long nNodesGlobal = 0;
        if( rank == 0 ) for( int i = 0; i < size; ++i ) nNodesGlobal += nodeCounts[i];

        // Per-rank gather targets on rank 0.
        std::vector< int >    allCellGids, allNumNodes, allConnGid, allNodeGids;
        std::vector< double > allCenterLat, allCenterLon, allNodeLat, allNodeLon;
        std::vector< double > allAreas;
        std::vector< int >    allMask;
        if( rank == 0 )
        {
            allCellGids.resize( nCellsGlobal );
            allNumNodes.resize( nCellsGlobal );
            allConnGid.resize( static_cast< size_t >( nCellsGlobal ) * ncpc );
            allCenterLat.resize( nCellsGlobal );
            allCenterLon.resize( nCellsGlobal );
            allNodeGids.resize( nNodesGlobal );
            allNodeLat.resize( nNodesGlobal );
            allNodeLon.resize( nNodesGlobal );
            if( mHasAreas ) allAreas.resize( nCellsGlobal );
            if( mHasMask ) allMask.resize( nCellsGlobal );
        }

        MPI_Gatherv( mLocalCellGids.data(), myCells, MPI_INT, allCellGids.data(), cellCounts.data(),
                     cellDispls.data(), MPI_INT, 0, comm );
        MPI_Gatherv( mElementNumNodes.data(), myCells, MPI_INT, allNumNodes.data(), cellCounts.data(),
                     cellDispls.data(), MPI_INT, 0, comm );
        MPI_Gatherv( mCenterLat.data(), myCells, MPI_DOUBLE, allCenterLat.data(), cellCounts.data(), cellDispls.data(),
                     MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mCenterLon.data(), myCells, MPI_DOUBLE, allCenterLon.data(), cellCounts.data(), cellDispls.data(),
                     MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mElementConn.data(), myCells * ncpc, MPI_INT, allConnGid.data(), cellCountsC.data(),
                     cellDisplsC.data(), MPI_INT, 0, comm );

        MPI_Gatherv( mNodeGids.data(), myNodes, MPI_INT, allNodeGids.data(), nodeCounts.data(), nodeDispls.data(),
                     MPI_INT, 0, comm );
        MPI_Gatherv( mNodeLat.data(), myNodes, MPI_DOUBLE, allNodeLat.data(), nodeCounts.data(), nodeDispls.data(),
                     MPI_DOUBLE, 0, comm );
        MPI_Gatherv( mNodeLon.data(), myNodes, MPI_DOUBLE, allNodeLon.data(), nodeCounts.data(), nodeDispls.data(),
                     MPI_DOUBLE, 0, comm );

        if( mHasAreas )
            MPI_Gatherv( mAreas.data(), myCells, MPI_DOUBLE, allAreas.data(), cellCounts.data(), cellDispls.data(),
                         MPI_DOUBLE, 0, comm );
        if( mHasMask )
            MPI_Gatherv( mMask.data(), myCells, MPI_INT, allMask.data(), cellCounts.data(), cellDispls.data(), MPI_INT,
                         0, comm );

        ErrorCode rc = MB_SUCCESS;
        if( rank == 0 )
            rc = writeAll( allCellGids, allNumNodes, allConnGid, allCenterLat, allCenterLon, allNodeGids, allNodeLat,
                           allNodeLon, allAreas, allMask );
        int rcInt = static_cast< int >( rc );
        MPI_Bcast( &rcInt, 1, MPI_INT, 0, comm );
        if( rcInt != MB_SUCCESS ) MB_SET_ERR( MB_FAILURE, "ESMF write failed on rank 0" );

        dbgOut.tprintf( 1, "  ESMF write: gathered+wrote %ld cells from %d ranks\n", mGlobalCells, size );
        return MB_SUCCESS;
    }
#endif

    return writeAll( mLocalCellGids, mElementNumNodes, mElementConn, mCenterLat, mCenterLon, mNodeGids, mNodeLat,
                     mNodeLon, mAreas, mMask );
}

ErrorCode NCWriteESMF::write_nonset_variables( std::vector< WriteNC::VarData >& /*vdatas*/,
                                               std::vector< int >& /*tstep_nums*/ )
{
    return MB_SUCCESS;
}

}  // namespace moab
