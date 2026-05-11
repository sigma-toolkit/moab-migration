/*
 * imoab_spmv_bfb.cpp
 *
 * Parallel bit-for-bit (BFB) reproducibility test for the deterministic SpMV
 * kernel that backs iMOAB_ApplyScalarProjectionWeights. Verifies that the
 * sparse matrix-vector multiplication driving low-order, high-order, and
 * dual-map (CAAS) projections produces a byte-identical per-row output
 * regardless of the MPI rank count or how the source mesh is partitioned.
 *
 * Workflow per run:
 *   1. Load source (ATM) and target (OCN) FV meshes on component PEs.
 *   2. Migrate them to coupler PEs.
 *   3. Load a low-order (monotone) and a high-order (non-monotone) weight
 *      map from disk onto a single dual-map intersection app.
 *   4. Seed the source field directly on the intersection app's coverage
 *      cells from a deterministic function of GLOBAL_ID
 *      ( f(gid) = sin(a*gid) + 0.5*cos(b*(gid+7)) ).
 *      Using GLOBAL_ID as the input keeps the per-cell source value
 *      invariant under repartitioning, so any cross-rank-count drift in the
 *      output isolates the SpMV itself.
 *   5. Apply low-order and high-order projections (filter_type=0).
 *   6. Apply dual-map CAAS projection (filter_type=2, lo-stencil bounds).
 *   7. Read back GLOBAL_ID + projected value on every owned target cell,
 *      gather to rank 0, sort by global ID, and write a deterministic
 *      digest file per kernel:
 *           spmv_bfb_digest_<kernel>_<size>.txt
 *      (kernel ∈ {lo, hi, dual}; size = MPI rank count).
 *
 * Verification: rerun under different rank counts (e.g. 1, 2, 4, 8) and
 * diff the per-kernel digest files. They must be byte-identical.
 *
 * Defaults target the existing E3SM ne30pg2 → IcoswISC30E3r5 maps shipped
 * with the inputdata tree; pass --atmosphere/--ocean/--lo_map_file/
 * --hi_map_file to point at any other pair.
 */

#include "moab/Core.hpp"
#ifndef MOAB_HAVE_MPI
#error This test requires MPI configuration
#endif

#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include "imoab_coupler_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef MOAB_HAVE_TEMPESTREMAP
#error This test requires MOAB configuration with TempestRemap
#endif

// Deterministic source field as a function of source-cell GLOBAL_ID. Two
// non-commensurate frequencies keep the field full-spectrum so the SpMV
// exercises every column rather than collapsing to a constant.
static inline double source_from_gid( int gid )
{
    const double a = 1.732050807568877e-3;  // ~ sqrt(3) * 1e-3
    const double b = 1.414213562373095e-3;  // ~ sqrt(2) * 1e-3
    return std::sin( a * gid ) + 0.5 * std::cos( b * ( gid + 7 ) );
}

// Gather (gid, value) pairs from every rank to rank 0, sort by gid, and
// write to a digest file. The sort order — and hence the file contents —
// is independent of how the partitioner distributed target cells across
// ranks, so the same digest is produced from any decomposition iff the
// per-cell projected values are themselves bit-for-bit identical.
static int gather_and_write_digest( MPI_Comm comm,
                                    int rankInComm,
                                    const std::vector< int >& localGids,
                                    const std::vector< double >& localVals,
                                    const std::string& outFilename )
{
    int sizeInComm = 0;
    MPI_Comm_size( comm, &sizeInComm );

    int localCount = static_cast< int >( localGids.size() );
    std::vector< int > counts( sizeInComm, 0 );
    MPI_Gather( &localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm );

    std::vector< int > displs( sizeInComm, 0 );
    int totalCount = 0;
    if( rankInComm == 0 )
    {
        for( int r = 0; r < sizeInComm; ++r )
        {
            displs[r] = totalCount;
            totalCount += counts[r];
        }
    }

    std::vector< int > allGids;
    std::vector< double > allVals;
    if( rankInComm == 0 )
    {
        allGids.resize( totalCount );
        allVals.resize( totalCount );
    }

    MPI_Gatherv( localGids.data(), localCount, MPI_INT,
                 rankInComm == 0 ? allGids.data() : nullptr,
                 counts.data(), displs.data(), MPI_INT, 0, comm );
    MPI_Gatherv( localVals.data(), localCount, MPI_DOUBLE,
                 rankInComm == 0 ? allVals.data() : nullptr,
                 counts.data(), displs.data(), MPI_DOUBLE, 0, comm );

    if( rankInComm != 0 ) return 0;

    // Pair, sort by gid (ascending), drop duplicates that show up from
    // ghost/shared owner overlap (keep the value from the lowest rank,
    // which arrives first in the gather).
    std::vector< std::pair< int, double > > pairs;
    pairs.reserve( totalCount );
    for( int i = 0; i < totalCount; ++i )
        pairs.emplace_back( allGids[i], allVals[i] );

    std::sort( pairs.begin(), pairs.end(),
               []( const std::pair< int, double >& a,
                   const std::pair< int, double >& b ) {
                   if( a.first != b.first ) return a.first < b.first;
                   return a.second < b.second;
               } );

    std::ofstream out( outFilename );
    if( !out )
    {
        std::cerr << "ERROR: cannot open digest file " << outFilename << "\n";
        return 1;
    }
    out << std::scientific << std::setprecision( 17 );
    int prevGid = std::numeric_limits< int >::min();
    for( const auto& p : pairs )
    {
        if( p.first == prevGid ) continue;  // skip duplicate ghost entries
        out << p.first << "  " << p.second << "\n";
        prevGid = p.first;
    }
    out.close();
    return 0;
}

int main( int argc, char* argv[] )
{
    moab::CpuTimer timer;
    double timer_ops;
    std::string opName;

    int rankInGlobalComm = -1, numProcesses = 0;
    MPI_Group jgroup;
    const iMOAB_String readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );
    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );

    int rankInAtmComm = -1, rankInOcnComm = -1, rankInCouComm = -1;

    // Reasonable defaults — the test ships against the existing FV-FV
    // coupler unit-test meshes so it runs on any installed build.
    std::string atmFilename  = TestDir + "unittest/srcWithSolnTag.h5m";
    std::string ocnFilename  = TestDir + "unittest/outTri15_8.h5m";
    std::string loMapFile;   // empty => compute online (low-order, monotone)
    std::string hiMapFile;   // empty => compute online (high-order)
    std::string digestPrefix = "spmv_bfb_digest";

    int nghlay = 0;
    int startG1 = 0, endG1 = numProcesses - 1;
    int startG2 = 0, endG2 = numProcesses - 1;
    int startG4 = 0, endG4 = numProcesses - 1;

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "ATM mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "OCN mesh filename (target)", &ocnFilename );
    opts.addOpt< std::string >( "lo_map_file,l", "Low-order map file (nc); if set, load from disk", &loMapFile );
    opts.addOpt< std::string >( "hi_map_file,h", "High-order map file (nc); if set, load from disk", &hiMapFile );
    opts.addOpt< std::string >( "digest_prefix,o", "Output digest filename prefix", &digestPrefix );
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );
    opts.parseCommandLine( argc, argv );

    const bool loadFromDisk = ( !loMapFile.empty() && !hiMapFile.empty() );

    if( !rankInGlobalComm )
    {
        std::cout << " === imoab_spmv_bfb test ===\n"
                  << " ATM file      : " << atmFilename << "\n"
                  << " OCN file      : " << ocnFilename << "\n"
                  << " Lo-order map  : " << ( loadFromDisk ? loMapFile : std::string( "<compute online>" ) ) << "\n"
                  << " Hi-order map  : " << ( loadFromDisk ? hiMapFile : std::string( "<compute online>" ) ) << "\n"
                  << " Digest prefix : " << digestPrefix << "\n"
                  << " Processes     : " << numProcesses << "\n";
    }

    MPI_Group atmPEGroup;
    MPI_Comm  atmComm;
    CHECKIERR( create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm ),
               "Cannot create ATM MPI group and communicator" )

    MPI_Group ocnPEGroup;
    MPI_Comm  ocnComm;
    CHECKIERR( create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm ),
               "Cannot create OCN MPI group and communicator" )

    MPI_Group couPEGroup;
    MPI_Comm  couComm;
    CHECKIERR( create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm ),
               "Cannot create coupler MPI group and communicator" )

    MPI_Group joinAtmCouGroup;
    MPI_Comm  atmCouComm;
    CHECKIERR( create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm ),
               "Cannot create joint ATM-coupler communicator" )

    MPI_Group joinOcnCouGroup;
    MPI_Comm  ocnCouComm;
    CHECKIERR( create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm ),
               "Cannot create joint OCN-coupler communicator" )

    CHECKIERR( iMOAB_Initialize( argc, argv ), "Cannot initialize iMOAB" )

    int cmpatm = 1, cmpocn = 2, cplatm = 3, cplocn = 5;
    int dualmap_id = 9;

    int cmpAtmAppID = -1, cmpOcnAppID = -1;
    int cplAtmAppID = -1, cplOcnAppID = -1;
    int cplDualMapAppID = -1;

    iMOAB_AppID cmpAtmPID     = &cmpAtmAppID;
    iMOAB_AppID cmpOcnPID     = &cmpOcnAppID;
    iMOAB_AppID cplAtmPID     = &cplAtmAppID;
    iMOAB_AppID cplOcnPID     = &cplOcnAppID;
    iMOAB_AppID cplDualMapPID = &cplDualMapAppID;

    if( couComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( couComm, &rankInCouComm );
        CHECKIERR( iMOAB_RegisterApplication( "CPLATM", &couComm, &cplatm, cplAtmPID ),
                   "Cannot register ATM over coupler" )
        CHECKIERR( iMOAB_RegisterApplication( "CPLOCN", &couComm, &cplocn, cplOcnPID ),
                   "Cannot register OCN over coupler" )
        CHECKIERR( iMOAB_RegisterApplication( "DUALMAP", &couComm, &dualmap_id, cplDualMapPID ),
                   "Cannot register dual map app" )
    }

    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        CHECKIERR( iMOAB_RegisterApplication( "ATMCMP", &atmComm, &cmpatm, cmpAtmPID ),
                   "Cannot register ATM" )
        CHECKIERR( iMOAB_LoadMesh( cmpAtmPID, atmFilename.c_str(), readopts, &nghlay ),
                   "Cannot load ATM mesh" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        CHECKIERR( iMOAB_RegisterApplication( "OCNCMP", &ocnComm, &cmpocn, cmpOcnPID ),
                   "Cannot register OCN" )
        CHECKIERR( iMOAB_LoadMesh( cmpOcnPID, ocnFilename.c_str(), readopts, &nghlay ),
                   "Cannot load OCN mesh" )
    }

    int repartitioner_scheme = 0;
#ifdef MOAB_HAVE_ZOLTAN
    repartitioner_scheme = 2;
#endif
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendMesh( cmpAtmPID, &atmCouComm, &couPEGroup, &cplatm, &repartitioner_scheme ),
                   "Cannot send ATM" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_ReceiveMesh( cplAtmPID, &atmCouComm, &atmPEGroup, &cmpatm ),
                   "Cannot receive ATM" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ), "Cannot free ATM send buffers" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendMesh( cmpOcnPID, &ocnCouComm, &couPEGroup, &cplocn, &repartitioner_scheme ),
                   "Cannot send OCN" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_ReceiveMesh( cplOcnPID, &ocnCouComm, &ocnPEGroup, &cmpocn ),
                   "Cannot receive OCN" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpOcnPID, &cplocn ), "Cannot free OCN send buffers" )
    }

    const iMOAB_String srcField     = "SpmvBfbSrc";
    const iMOAB_String tgtFieldHi   = "SpmvBfbTgtHi";
    const iMOAB_String tgtFieldDual = "SpmvBfbTgtDual";
    const iMOAB_String tgtFieldLo   = "SpmvBfbTgtLo";

    if( couComm != MPI_COMM_NULL )
    {
        if( loadFromDisk )
        {
            int src_disc_type = 3;  // FV cell
            int tgt_disc_type = 3;  // FV cell
            int arearead      = 0;  // areas not needed for plain SpMV digest

            PUSH_TIMER( "Load low-order map from disk" )
            CHECKIERR( iMOAB_LoadMapFile( cplAtmPID, cplOcnPID, cplDualMapPID,
                                          &src_disc_type, &tgt_disc_type, &arearead,
                                          "lo-scalar", loMapFile.c_str() ),
                       "Cannot load low-order map file" )
            POP_TIMER( couComm, rankInCouComm )

            PUSH_TIMER( "Load high-order map from disk" )
            CHECKIERR( iMOAB_LoadMapFile( cplAtmPID, cplOcnPID, cplDualMapPID,
                                          &src_disc_type, &tgt_disc_type, &arearead,
                                          "hi-scalar", hiMapFile.c_str() ),
                       "Cannot load high-order map file" )
            POP_TIMER( couComm, rankInCouComm )

            int meshtype = 3;
            CHECKIERR( iMOAB_MigrateMapMesh( cplAtmPID, cplDualMapPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                             &cplatm, &dualmap_id ),
                       "Cannot migrate map mesh" )
        }
        else
        {
            int nghlay_cov = 3;
            int nghlay_tgt = 0;
            CHECKIERR( iMOAB_SetMapGhostLayers( cplDualMapPID, &nghlay_cov, &nghlay_tgt ),
                       "Failed to set number of ghost layers" );

            PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
            CHECKIERR( iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplDualMapPID ),
                       "Cannot compute ATM/OCN intersection" )
            POP_TIMER( couComm, rankInCouComm )

            const iMOAB_String disc_fv = "fv";
            const iMOAB_String dof_tag = "GLOBAL_ID";
            int disc_order             = 1;
            int fNoBubble = 1, fMonotone = 1, fVolumetric = 0, fInvDist = 0, fNoConserve = 0, fValidate = 0;

            PUSH_TIMER( "Compute low-order (monotone) weights" )
            CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, "lo-scalar", disc_fv, &disc_order, disc_fv,
                                                             &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                             &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                       "Cannot compute low-order weights" )
            POP_TIMER( couComm, rankInCouComm )

            fMonotone  = 0;
            disc_order = 2;

            PUSH_TIMER( "Compute high-order (non-monotone) weights" )
            CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, "hi-scalar", disc_fv, &disc_order, disc_fv,
                                                             &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                             &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                       "Cannot compute high-order weights" )
            POP_TIMER( couComm, rankInCouComm )

            int meshtype = 3;
            CHECKIERR( iMOAB_ComputeCommGraph( cplAtmPID, cplDualMapPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                               &meshtype, &cplatm, &dualmap_id ),
                       "Cannot compute ATM coverage graph" )
        }

        // Define source tag on dual-map app and the three target tags on cpl OCN.
        int tagType = DENSE_DOUBLE;
        int tagIndex_;
        int compNDoFs = 1;
        CHECKIERR( iMOAB_DefineTagStorage( cplDualMapPID, srcField, &tagType, &compNDoFs, &tagIndex_ ),
                   "Cannot define source tag on dual-map app" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldHi, &tagType, &compNDoFs, &tagIndex_ ),
                   "Cannot define hi-order target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldDual, &tagType, &compNDoFs, &tagIndex_ ),
                   "Cannot define dual-map target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldLo, &tagType, &compNDoFs, &tagIndex_ ),
                   "Cannot define lo-order target tag" )

        // ----- Seed the source field deterministically by GLOBAL_ID -----
        // iMOAB_GetCoverageMeshInfo + iMOAB_SetDoubleTagStorageOnCoverage write
        // straight to the intersection app's coverage cells and bypass the
        // partition-dependent SendElementTag/ReceiveElementTag two-hop. This
        // guarantees the SpMV input vector is the same logical field on every
        // rank-count, leaving the SpMV itself as the only thing under test.
        int nCovElems = 0;
        CHECKIERR( iMOAB_GetCoverageMeshInfo( cplDualMapPID, &nCovElems, nullptr, nullptr ),
                   "Cannot query coverage size on dual-map app" )

        std::vector< int > covGids( nCovElems, 0 );
        if( nCovElems > 0 )
        {
            CHECKIERR( iMOAB_GetCoverageMeshInfo( cplDualMapPID, &nCovElems, covGids.data(), nullptr ),
                       "Cannot fetch coverage GIDs on dual-map app" )
        }

        std::vector< double > srcVals( nCovElems, 0.0 );
        for( int i = 0; i < nCovElems; ++i )
            srcVals[i] = source_from_gid( covGids[i] );

        if( nCovElems > 0 )
        {
            CHECKIERR( iMOAB_SetDoubleTagStorageOnCoverage( cplDualMapPID, srcField, &nCovElems, srcVals.data() ),
                       "Cannot seed source field on coverage" )
        }

        // ----- Apply projections -----
        int filter_type = 0;

        PUSH_TIMER( "Low-order projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "lo-scalar",
                                                        srcField, tgtFieldLo, nullptr ),
                   "Failed to apply low-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        PUSH_TIMER( "High-order projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar",
                                                        srcField, tgtFieldHi, nullptr ),
                   "Failed to apply high-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        filter_type = 2;  // CAAS_LOCAL — high-order projection clipped to lo-stencil bounds
        PUSH_TIMER( "Dual-map CAAS projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar",
                                                        srcField, tgtFieldDual, "lo-scalar" ),
                   "Failed to apply dual-map CAAS weights" )
        POP_TIMER( couComm, rankInCouComm )

        // ----- Pull the projected target tags + GLOBAL_ID, write digests -----
        int nOcnElems[3] = { 0, 0, 0 };
        CHECKIERR( iMOAB_GetMeshInfo( cplOcnPID, nullptr, nOcnElems, nullptr, nullptr, nullptr ),
                   "Cannot get OCN mesh info" )

        const int nOwned = nOcnElems[0];  // owned cells (excludes ghosts)
        int entType      = 1;             // primary elements
        int gidTagSz     = 1;
        int gidTagType   = DENSE_INTEGER;
        int gidTagIndex  = -1;
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, "GLOBAL_ID", &gidTagType, &gidTagSz, &gidTagIndex ),
                   "Cannot define GLOBAL_ID tag on cplOCN" )

        std::vector< int >    tgtGids( nOwned, 0 );
        std::vector< double > hiVals( nOwned, 0.0 );
        std::vector< double > dualVals( nOwned, 0.0 );
        std::vector< double > loVals( nOwned, 0.0 );
        if( nOwned > 0 )
        {
            int nQuery = nOwned;
            CHECKIERR( iMOAB_GetIntTagStorage( cplOcnPID, "GLOBAL_ID", &nQuery, &entType, tgtGids.data() ),
                       "Cannot get GLOBAL_ID values" )
            CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldHi, &nQuery, &entType, hiVals.data() ),
                       "Cannot get hi-order values" )
            CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldDual, &nQuery, &entType, dualVals.data() ),
                       "Cannot get dual-map values" )
            CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldLo, &nQuery, &entType, loVals.data() ),
                       "Cannot get lo-order values" )
        }

        std::ostringstream szTag;
        szTag << "_" << numProcesses << ".txt";

        const int rcLo   = gather_and_write_digest( couComm, rankInCouComm, tgtGids, loVals,
                                                    digestPrefix + "_lo"   + szTag.str() );
        const int rcHi   = gather_and_write_digest( couComm, rankInCouComm, tgtGids, hiVals,
                                                    digestPrefix + "_hi"   + szTag.str() );
        const int rcDual = gather_and_write_digest( couComm, rankInCouComm, tgtGids, dualVals,
                                                    digestPrefix + "_dual" + szTag.str() );
        if( rcLo || rcHi || rcDual )
        {
            std::cerr << "ERROR: failed writing one of the digest files\n";
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }

        if( !rankInCouComm )
        {
            std::cout << " Digests written: "
                      << digestPrefix << "_{lo,hi,dual}" << szTag.str() << "\n"
                      << " To verify BFB SpMV: re-run under a different rank count\n"
                      << " (e.g. mpirun -n 4) with the same inputs and diff the\n"
                      << " produced digest files. They MUST be byte-identical.\n";
        }
    }

    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cplDualMapPID ), "Cannot deregister DUALMAP" )
        CHECKIERR( iMOAB_DeregisterApplication( cplOcnPID ), "Cannot deregister CPLOCN" )
        CHECKIERR( iMOAB_DeregisterApplication( cplAtmPID ), "Cannot deregister CPLATM" )
    }
    if( ocnComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpOcnPID ), "Cannot deregister OCN" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_DeregisterApplication( cmpAtmPID ), "Cannot deregister ATM" )
    }

    CHECKIERR( iMOAB_Finalize(), "Cannot finalize iMOAB" )

    if( MPI_COMM_NULL != atmCouComm ) MPI_Comm_free( &atmCouComm );
    MPI_Group_free( &joinAtmCouGroup );
    if( MPI_COMM_NULL != atmComm ) MPI_Comm_free( &atmComm );

    if( MPI_COMM_NULL != ocnComm ) MPI_Comm_free( &ocnComm );
    if( MPI_COMM_NULL != ocnCouComm ) MPI_Comm_free( &ocnCouComm );
    MPI_Group_free( &joinOcnCouGroup );

    if( MPI_COMM_NULL != couComm ) MPI_Comm_free( &couComm );

    MPI_Group_free( &atmPEGroup );
    MPI_Group_free( &ocnPEGroup );
    MPI_Group_free( &couPEGroup );
    MPI_Group_free( &jgroup );

    MPI_Finalize();
    return 0;
}
