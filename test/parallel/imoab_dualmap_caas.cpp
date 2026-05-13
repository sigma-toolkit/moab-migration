/*
 * imoab_dualmap_caas.cpp
 *
 * Test for dual-map nonlinear remapping (CAAS with low-order map bounds).
 *
 * The PRIMARY workflow loads pre-computed weight maps from disk via
 * --lo_map_file and --hi_map_file. This matches the standard E3SM coupler
 * workflow (maps are generated offline by TempestRemap and loaded at
 * runtime), and only this path produces strictly bit-for-bit reproducible
 * per-cell projection values across MPI rank counts. Online weight
 * generation via iMOAB_ComputeScalarProjectionWeights has a known
 * partition-dependent residual at the 1-3 ULP level on master and is
 * available as an alternative via --compute_online.
 *
 * Workflow:
 *   1. Load ATM (source) and OCN (target) meshes on all processes
 *   2. Migrate meshes to coupler communicator
 *   3. Set up two FV weight maps on the dual-map intersection app:
 *      - "lo-scalar" : low-order monotone map
 *      - "hi-scalar" : high-order non-monotone map
 *      Default: load both from disk (--lo_map_file, --hi_map_file).
 *      With --compute_online: compute both via iMOAB_ComputeScalarProjectionWeights.
 *   4. Define source field (degree-2 spherical harmonic on the ATM mesh)
 *   5. Apply lo, hi (no CAAS), and dual-map CAAS projections
 *   6. Verify: target dual-CAAS values are within per-row stencil bounds
 *   7. Optionally write per-cell BFB digest files (--digest_prefix)
 *
 * BFB digest workflow (cross-rank-count regression check):
 *   for n in 1 2 4 8; do mpirun -n $n ./imoab_dualmap_caas \
 *       -l <lo_map.nc> -h <hi_map.nc> -o digest ; done
 *   for k in lo hi dual; do diff -q digest_${k}_1.txt digest_${k}_4.txt; done
 *   All file pairs must be byte-identical when maps are loaded from disk.
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
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>
#include <cmath>

#ifndef MOAB_HAVE_TEMPESTREMAP
#error This test requires MOAB configuration with TempestRemap
#endif

// Gather (gid, value) pairs from every rank to rank 0, sort ascending by
// gid, drop duplicates that arise from ghost/shared owner overlap, and
// write to a digest file.  Same global IDs in any decomposition produce
// the same file iff the per-cell projected values are byte-identical
// across rank counts.  Diff the digest files from different mpirun -n
// invocations to verify cross-rank-count BFB reproducibility of the
// dual-map CAAS path.
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

    std::vector< int >    allGids;
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

    int rankInGlobalComm, numProcesses;
    MPI_Group jgroup;
    const iMOAB_String readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm_group( MPI_COMM_WORLD, &jgroup );

    int rankInAtmComm = -1, rankInOcnComm = -1, rankInCouComm = -1;

    // Use the same FV mesh files as imoab_read_compute_map.cpp so source field is FV.
    std::string atmFilename = TestDir + "unittest/srcWithSolnTag.h5m";
    std::string ocnFilename = TestDir + "unittest/recMeshOcn.h5m";
    std::string loMapFile;       // primary path: load from disk
    std::string hiMapFile;       // primary path: load from disk
    std::string digestPrefix;    // empty = skip digest dump
    std::string writeMapsPrefix; // if set, dump computed/loaded maps to <prefix>_{lo,hi}.nc
    bool compute_online = false; // alternative: compute weight maps online

    int nghlay = 0;

    // PE layout: default all tasks on all groups
    int startG1 = 0, endG1 = numProcesses - 1;  // ATM
    int startG2 = 0, endG2 = numProcesses - 1;  // OCN
    int startG4 = 0, endG4 = numProcesses - 1;  // Coupler

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "ATM mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "OCN mesh filename (target)", &ocnFilename );
    opts.addOpt< std::string >( "lo_map_file,l",
                                "Low-order map file (nc) — primary path: load from disk",
                                &loMapFile );
    opts.addOpt< std::string >( "hi_map_file,h",
                                "High-order map file (nc) — primary path: load from disk",
                                &hiMapFile );
    opts.addOpt< void >( "compute_online,n",
                         "Alternative: compute lo/hi weight maps online via iMOAB_ComputeScalarProjectionWeights "
                         "(non-BFB across rank counts at the ULP level — for testing only)",
                         &compute_online );
    opts.addOpt< std::string >( "write_maps,w",
                                "If set, write the active lo/hi weight maps to disk as <prefix>_lo.nc and "
                                "<prefix>_hi.nc after they are computed or loaded; intended for the BFB regression "
                                "workflow (serial compute → write → reload in parallel)",
                                &writeMapsPrefix );
    opts.addOpt< std::string >( "digest_prefix,o",
                                "If set, write per-cell BFB digest files <prefix>_{lo,hi,dual}_<np>.txt",
                                &digestPrefix );
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );
    opts.parseCommandLine( argc, argv );

    // Primary path: load lo/hi maps from disk (BFB across rank counts).
    // Fall back to online computation either when --compute_online is set
    // explicitly, or when no map files were supplied (so the test still
    // runs without arguments — but with a non-BFB warning).
    const bool haveBothMaps = ( !loMapFile.empty() && !hiMapFile.empty() );
    const bool fallbackOnline = ( !haveBothMaps && !compute_online );
    const bool loadFromDisk = ( haveBothMaps && !compute_online );

    if( !rankInGlobalComm )
    {
        std::cout << " === imoab_dualmap_caas test ===\n";
        std::cout << " ATM file: " << atmFilename << "\n";
        std::cout << " OCN file: " << ocnFilename << "\n";
        if( loadFromDisk )
        {
            std::cout << " Lo-order map: " << loMapFile << "\n";
            std::cout << " Hi-order map: " << hiMapFile << "\n";
            std::cout << " Mode: load maps from disk (BFB across rank counts)\n";
        }
        else if( compute_online )
        {
            std::cout << " Mode: compute maps online (--compute_online; non-BFB at ULP level)\n";
        }
        else
        {
            std::cout << " Mode: compute maps online (no map files supplied; fallback)\n";
            std::cout << " WARNING: online weight generation is NOT BFB across rank counts.\n"
                      << "          For BFB validation pass --lo_map_file and --hi_map_file\n"
                      << "          pointing at pre-computed netcdf weight files.\n";
            (void)fallbackOnline;
        }
        std::cout << " Processes: " << numProcesses << "\n";
        std::cout << " ATM tasks: " << startG1 << ":" << endG1
                  << ", OCN tasks: " << startG2 << ":" << endG2
                  << ", Coupler tasks: " << startG4 << ":" << endG4 << "\n";
    }

    // Create MPI communicators and groups using PE layout ranges
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    CHECKIERR( create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm ),
               "Cannot create ATM MPI group and communicator" )

    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    CHECKIERR( create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm ),
               "Cannot create OCN MPI group and communicator" )

    MPI_Group couPEGroup;
    MPI_Comm couComm;
    CHECKIERR( create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm ),
               "Cannot create coupler MPI group and communicator" )

    // Joint communicators for component-coupler data transfer
    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    CHECKIERR( create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm ),
               "Cannot create joint ATM-coupler communicator" )

    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    CHECKIERR( create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm ),
               "Cannot create joint OCN-coupler communicator" )

    CHECKIERR( iMOAB_Initialize( argc, argv ), "Cannot initialize iMOAB" )

    // Application IDs
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

    // Register applications on coupler
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

    // Register and load ATM on component PEs
    if( atmComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( atmComm, &rankInAtmComm );
        CHECKIERR( iMOAB_RegisterApplication( "ATMCMP", &atmComm, &cmpatm, cmpAtmPID ),
                   "Cannot register ATM" )
        CHECKIERR( iMOAB_LoadMesh( cmpAtmPID, atmFilename.c_str(), readopts, &nghlay ),
                   "Cannot load ATM mesh" )
    }

    // Register and load OCN on component PEs
    if( ocnComm != MPI_COMM_NULL )
    {
        MPI_Comm_rank( ocnComm, &rankInOcnComm );
        CHECKIERR( iMOAB_RegisterApplication( "OCNCMP", &ocnComm, &cmpocn, cmpOcnPID ),
                   "Cannot register OCN" )
        CHECKIERR( iMOAB_LoadMesh( cmpOcnPID, ocnFilename.c_str(), readopts, &nghlay ),
                   "Cannot load OCN mesh" )
    }

    // Migrate ATM mesh to coupler
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

    // Migrate OCN mesh to coupler
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

    // FV scalar field already present on srcWithSolnTag.h5m (matches imoab_read_compute_map.cpp).
    const iMOAB_String srcField     = "AnalyticalSolnSrcExact";
    const iMOAB_String tgtFieldHi   = "TargetHiOrder";
    const iMOAB_String tgtFieldDual = "TargetDualMap";
    const iMOAB_String tgtFieldLo   = "TargetLoOrder";
    if( couComm != MPI_COMM_NULL )
    {
        if( loadFromDisk )
        {
            // --- Load pre-computed weight maps from disk ---
            int src_disc_type = 3;  // FV cell
            int tgt_disc_type = 3;  // FV cell
            int arearead      = 3;  // read all areas

            PUSH_TIMER( "Load low-order map from disk" )
            CHECKIERR( iMOAB_LoadMapFile( cplAtmPID, cplOcnPID, cplDualMapPID,
                                          &src_disc_type, &tgt_disc_type, &arearead,
                                          "lo-scalar", loMapFile.c_str() ),
                       "Cannot load low-order map file" )
            POP_TIMER( couComm, rankInCouComm )

            arearead = 0;  // do not read areas
            PUSH_TIMER( "Load high-order map from disk" )
            CHECKIERR( iMOAB_LoadMapFile( cplAtmPID, cplOcnPID, cplDualMapPID,
                                          &src_disc_type, &tgt_disc_type, &arearead,
                                          "hi-scalar", hiMapFile.c_str() ),
                       "Cannot load high-order map file" )
            POP_TIMER( couComm, rankInCouComm )

            // Migrate the coverage mesh so source tag data can be transferred
            int meshtype = 3;
            CHECKIERR( iMOAB_MigrateMapMesh( cplAtmPID, cplDualMapPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                             &cplatm, &dualmap_id ),
                       "Cannot migrate map mesh for lo-scalar" )
        }
        else
        {
            // --- Compute weight maps online ---

            // Set the ghost layers on the coupler for the ATM mesh
            int nghlay = 3;
            int nghlay_tgt = 0;
            CHECKIERR( iMOAB_SetMapGhostLayers( cplDualMapPID, &nghlay, &nghlay_tgt ),
                       "Failed to set number of ghost layers on ATM mesh" );

            // Compute mesh intersection between ATM and OCN on coupler
            PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
            CHECKIERR( iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplDualMapPID ),
                       "Cannot compute ATM/OCN intersection" )
            POP_TIMER( couComm, rankInCouComm )

            // Compute LOW-ORDER weights (monotone FV, 1st order)
            const iMOAB_String disc_fv  = "fv";
            const iMOAB_String dof_tag  = "GLOBAL_ID";
            int disc_order              = 1;
            int fNoBubble = 1, fMonotone = 1, fVolumetric = 0, fInvDist = 0, fNoConserve = 0, fValidate = 0;

            PUSH_TIMER( "Compute low-order (monotone) weights" )
            CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, "lo-scalar", disc_fv, &disc_order, disc_fv,
                                                             &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                             &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                       "Cannot compute low-order weights" )
            POP_TIMER( couComm, rankInCouComm )

            // Compute HIGH-ORDER weights (non-monotone FV)
            fMonotone  = 0;
            disc_order = 2;

            PUSH_TIMER( "Compute high-order (non-monotone) weights" )
            CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, "hi-scalar", disc_fv, &disc_order, disc_fv,
                                                             &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                             &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                       "Cannot compute high-order weights" )
            POP_TIMER( couComm, rankInCouComm )

            // Compute coverage comm graph for tag migration
            int meshtype = 3;
            CHECKIERR( iMOAB_ComputeCommGraph( cplAtmPID, cplDualMapPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                               &meshtype, &cplatm, &dualmap_id ),
                       "Cannot compute ATM coverage graph" )
        }

        // Optional: persist the active maps to disk for the BFB regression
        // workflow. Called collectively. Two separate netcdf files are
        // emitted, one per map identifier, so a subsequent run can reload
        // them via --lo_map_file / --hi_map_file and verify cross-rank-count
        // BFB through diff'ing per-cell digests.
        if( !writeMapsPrefix.empty() )
        {
            const std::string loOut = writeMapsPrefix + "_lo.nc";
            const std::string hiOut = writeMapsPrefix + "_hi.nc";
            CHECKIERR( iMOAB_WriteMapFile( cplDualMapPID, "lo-scalar", loOut.c_str() ),
                       "Cannot write low-order map file" )
            CHECKIERR( iMOAB_WriteMapFile( cplDualMapPID, "hi-scalar", hiOut.c_str() ),
                       "Cannot write high-order map file" )
            if( !rankInCouComm )
                std::cout << " Wrote weight maps to " << loOut << " and " << hiOut << "\n";
        }

        // Define source and target tags
        int tagType                     = 1;  // DENSE_DOUBLE
        int tagIndex;
        int atmCompNDoFs = 1, ocnCompNDoFs = 1;

        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, srcField, &tagType, &atmCompNDoFs, &tagIndex ),
                   "Cannot define source tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplDualMapPID, srcField, &tagType, &atmCompNDoFs, &tagIndex ),
                   "Cannot define source tag on dual-map app" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldHi, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define high-order target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldDual, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define dual-map target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldLo, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define low-order target tag" )
    }

    // Set source field values directly on the dual-map coverage mesh:
    // The source tag (srcField = "a2oTbot") is already present on the on-disk
    // ATM mesh (wholeATM_T.h5m). We follow the exact tag-migration pattern
    // used by imoab_read_compute_map.cpp: define-on-cpl-side, then
    // SendElementTag/ReceiveElementTag in two hops:
    //   cmpAtm -> cplAtm  (via atmCouComm, context cplatm)
    //   cplAtm -> cplDualMap (via couComm, context dualmap_id)
    // Both maps (lo + hi) live on cplDualMapPID and consume the same coverage,
    // so we migrate the source tag once.
    if( couComm != MPI_COMM_NULL )
    {
        int tagType = DENSE_DOUBLE, atmCompNDoFs = 1, tagIndex_ = -1;
        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, srcField, &tagType, &atmCompNDoFs, &tagIndex_ ),
                   "Cannot define src tag on cplAtm" )
    }

    // First hop: cmpAtm -> cplAtm
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cmpAtmPID, srcField, &atmCouComm, &cplatm ),
                   "Cannot send src tag from cmpAtm to cplAtm" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_ReceiveElementTag( cplAtmPID, srcField, &atmCouComm, &cmpatm ),
                   "Cannot receive src tag on cplAtm" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ), "Cannot free sender buffers (cmpAtm)" )
    }

    // Second hop: cplAtm -> cplDualMap (the intersection app's coverage mesh).
    // Mirrors the COMPUTE_FILE_MAP / COMPUTE_ONLINE_MAP send-tag block in
    // imoab_read_compute_map.cpp.
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cplAtmPID, srcField, &couComm, &dualmap_id ),
                   "Cannot send src tag from cplAtm to cplDualMap" )
        CHECKIERR( iMOAB_ReceiveElementTag( cplDualMapPID, srcField, &couComm, &cplatm ),
                   "Cannot receive src tag on cplDualMap coverage" )
        CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &dualmap_id ), "Cannot free sender buffers (cplAtm)" )
    }

    // === Apply projections and test ===
    if( couComm != MPI_COMM_NULL )
    {
        int filter_type = 0;  // no CAAS for plain high-order

        // 1) Apply high-order projection WITHOUT CAAS (baseline)
        PUSH_TIMER( "Apply high-order projection (no CAAS)" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar",
                                                        srcField, tgtFieldHi, nullptr ),
                   "Failed to apply high-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 2) Apply low-order projection (reference, should already be bounded)
        filter_type = 0;
        PUSH_TIMER( "Apply low-order projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "lo-scalar",
                                                        srcField, tgtFieldLo, nullptr ),
                   "Failed to apply low-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 3) Apply high-order projection WITH dual-map CAAS bounds from low-order map
        filter_type = 2;  // CAAS_LOCAL
        PUSH_TIMER( "Apply dual-map CAAS projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar", srcField, tgtFieldDual,
                                                       "lo-scalar" ),
                   "Failed to apply dual-map CAAS weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 4) Verify bounds preservation: the dual-map CAAS result must stay within
        //    the per-row stencil bounds computed from the low-order weight map.
        //    Note: we check stencil bounds (not global source range) because the CAAS
        //    algorithm guarantees per-row bounds, not global bounds.
        int nOcnElems[3];
        CHECKIERR( iMOAB_GetMeshInfo( cplOcnPID, nullptr, nOcnElems, nullptr, nullptr, nullptr ),
                   "Cannot get OCN mesh info" )
        int nDualElems[3];
        CHECKIERR( iMOAB_GetMeshInfo( cplDualMapPID, nullptr, nDualElems, nullptr, nullptr, nullptr ),
                   "Cannot get DualMap mesh info" )

        std::vector< double > hiVals( nOcnElems[2] ), dualVals( nOcnElems[2] ), loVals( nOcnElems[2] );
        int entity_type = 1;

        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldHi, &nOcnElems[2], &entity_type, hiVals.data() ),
                   "Cannot get hi-order values" )
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldDual, &nOcnElems[2], &entity_type, dualVals.data() ),
                   "Cannot get dual-map values" )
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, tgtFieldLo, &nOcnElems[2], &entity_type, loVals.data() ),
                   "Cannot get lo-order values" )

        // Get per-row stencil bounds from the low-order map (stored by ComputeRowBounds)
        // These tags were created by iMOAB_ApplyScalarProjectionWeights on the intersection
        // app's target entities; register them on cplOcnPID so we can read them.
        std::string loBoundName = std::string(srcField) + "_DualMapLoBound";
        std::string hiBoundName = std::string(srcField) + "_DualMapHiBound";
        {
            int tagType_dbl = 1, tagIndex, nDoFs = 1;
            CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, loBoundName.c_str(), &tagType_dbl, &nDoFs, &tagIndex ),
                       "Cannot define stencil lo bound tag" )
            CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, hiBoundName.c_str(), &tagType_dbl, &nDoFs, &tagIndex ),
                       "Cannot define stencil hi bound tag" )
        }
        std::vector< double > stencilLo( nOcnElems[2] ), stencilHi( nOcnElems[2] );
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, loBoundName.c_str(), &nOcnElems[2], &entity_type, stencilLo.data() ),
                   "Cannot get stencil lower bounds" )
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, hiBoundName.c_str(), &nOcnElems[2], &entity_type, stencilHi.data() ),
                   "Cannot get stencil upper bounds" )

        // Count violations against per-row stencil bounds
        int localHiViolations   = 0;
        int localDualViolations = 0;
        double maxHiExceedance  = 0.0;
        double maxDualExceedance = 0.0;
        const double tol_bounds = 1e-10;

        for( int i = 0; i < nOcnElems[2]; i++ )
        {
            // Check high-order against stencil bounds (may violate — that's expected)
            if( hiVals[i] < stencilLo[i] - tol_bounds || hiVals[i] > stencilHi[i] + tol_bounds )
            {
                localHiViolations++;
                double exc = std::max( stencilLo[i] - hiVals[i], hiVals[i] - stencilHi[i] );
                maxHiExceedance = std::max( maxHiExceedance, exc );
            }

            // Check dual-map against stencil bounds (should NOT violate)
            if( dualVals[i] < stencilLo[i] - tol_bounds || dualVals[i] > stencilHi[i] + tol_bounds )
            {
                localDualViolations++;
                double exc = std::max( stencilLo[i] - dualVals[i], dualVals[i] - stencilHi[i] );
                maxDualExceedance = std::max( maxDualExceedance, exc );
                if( localDualViolations <= 5 )
                    printf( "  [rank %d] dual violation #%d: i=%d dualVal=%.15e bounds=[%.15e,%.15e]\n",
                            rankInCouComm, localDualViolations, i, dualVals[i], stencilLo[i], stencilHi[i] );
            }
        }

        int globalHiViolations = 0, globalDualViolations = 0;
        double globalMaxHiExc = 0.0, globalMaxDualExc = 0.0;
        MPI_Reduce( &localHiViolations, &globalHiViolations, 1, MPI_INT, MPI_SUM, 0, couComm );
        MPI_Reduce( &localDualViolations, &globalDualViolations, 1, MPI_INT, MPI_SUM, 0, couComm );
        MPI_Reduce( &maxHiExceedance, &globalMaxHiExc, 1, MPI_DOUBLE, MPI_MAX, 0, couComm );
        MPI_Reduce( &maxDualExceedance, &globalMaxDualExc, 1, MPI_DOUBLE, MPI_MAX, 0, couComm );

        if( !rankInCouComm )
        {
            std::cout << "\n === Dual-Map CAAS Bounds Test Results ===\n";
            std::cout << " High-order violations:  " << globalHiViolations
                      << " (max exceedance: " << globalMaxHiExc << ")\n";
            std::cout << " Dual-map violations:    " << globalDualViolations
                      << " (max exceedance: " << globalMaxDualExc << ")\n";

            if( globalDualViolations == 0 )
                std::cout << " RESULT: PASS - Dual-map CAAS preserved bounds\n\n";
            else
                std::cout << " RESULT: FAIL - Dual-map CAAS violated bounds\n\n";
        }

        // Return failure if dual-map violated bounds
        if( globalDualViolations > 0 )
        {
            MPI_Abort( MPI_COMM_WORLD, 1 );
            return 1;
        }

        // 5) BFB digest dump (optional). Gather (target_gid, value) per
        //    OWNED OCN cell to root, sort by gid, write digest_{lo,hi,dual}_<np>.txt.
        //    Source field comes from srcWithSolnTag.h5m so per-cell input
        //    values are partition-independent by construction. Running this
        //    binary under different mpirun -n values must produce byte-identical
        //    digest files for each kernel.
        //
        //    We restrict to OWNED cells (nOcnElems[0], not nOcnElems[2]). MOAB
        //    orders entities owned-first in the visible range, so the first
        //    nOcnElems[0] tag entries belong to cells this rank owns. Ghost
        //    cells (the next nOcnElems[1] entries) are written by
        //    ApplyWeightsWithDualMap with stale lcl_lo/lcl_hi=0 bounds and
        //    therefore receive a dM_total/cap_g-scaled garbage value, so
        //    including them would make the dual digest non-BFB across rank
        //    counts even though the underlying SpMV result on owned cells
        //    is bit-identical.
        if( !digestPrefix.empty() )
        {
            int gidTagType = DENSE_INTEGER;
            int gidNDoFs   = 1;
            int gidIndex   = -1;
            int entType    = 1;
            CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, "GLOBAL_ID", &gidTagType, &gidNDoFs, &gidIndex ),
                       "Cannot define GLOBAL_ID tag on cplOcn" )
            const int nOwned = nOcnElems[0];
            std::vector< int > tgtGids( nOwned );
            int nQuery = nOwned;
            CHECKIERR( iMOAB_GetIntTagStorage( cplOcnPID, "GLOBAL_ID", &nQuery, &entType, tgtGids.data() ),
                       "Cannot get GLOBAL_ID values on cplOcn" )

            // Truncate the projected-value vectors to OWNED only too, so the
            // digest gathers bit-identical (gid, value) pairs across rank counts.
            std::vector< double > loValsOwned( loVals.begin(), loVals.begin() + nOwned );
            std::vector< double > hiValsOwned( hiVals.begin(), hiVals.begin() + nOwned );
            std::vector< double > dualValsOwned( dualVals.begin(), dualVals.begin() + nOwned );

            std::ostringstream szTag;
            szTag << "_" << numProcesses << ".txt";

            const int rcLo   = gather_and_write_digest( couComm, rankInCouComm, tgtGids, loValsOwned,
                                                        digestPrefix + "_lo"   + szTag.str() );
            const int rcHi   = gather_and_write_digest( couComm, rankInCouComm, tgtGids, hiValsOwned,
                                                        digestPrefix + "_hi"   + szTag.str() );
            const int rcDual = gather_and_write_digest( couComm, rankInCouComm, tgtGids, dualValsOwned,
                                                        digestPrefix + "_dual" + szTag.str() );
            if( rcLo || rcHi || rcDual )
            {
                std::cerr << "ERROR: failed writing one of the digest files\n";
                MPI_Abort( MPI_COMM_WORLD, 1 );
            }
            if( !rankInCouComm )
            {
                std::cout << " Digests written: " << digestPrefix
                          << "_{lo,hi,dual}" << szTag.str() << "\n"
                          << " Verify cross-rank-count BFB by re-running with a\n"
                          << " different mpirun -n and diff'ing the digest files.\n";
            }
        }
    }

    if( couComm != MPI_COMM_NULL )
    {
        // Write the OCN coupler-side mesh, which carries the projected target
        // tags (TargetHiOrder, TargetLoOrder, TargetDualMap) plus the per-row
        // bound diagnostics (SourceAnalytical_DualMapLoBound/HiBound).
        char outputFileCpl[]    = "cplOcnProjFile.h5m";
        char fileWriteOptions[] = "PARALLEL=WRITE_PART";
        CHECKIERR( iMOAB_WriteMesh( cplOcnPID, outputFileCpl, fileWriteOptions ),
                   "could not write cplOcnProjFile.h5m to disk" )
    }

    if( ocnComm != MPI_COMM_NULL )
    {
        // write only for n==1 case
        char outputFileRecvd[]  = "cmpOcnProjFile.h5m";
        char fileWriteOptions[] = "PARALLEL=WRITE_PART";
        CHECKIERR( iMOAB_WriteMesh( cmpOcnPID, outputFileRecvd, fileWriteOptions ),
                   "could not write cmpOcnProjFile.h5m to disk" )
    }

    // Cleanup
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

    // Free MPI communicators and groups (matching imoab_read_compute_map.cpp pattern)
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
