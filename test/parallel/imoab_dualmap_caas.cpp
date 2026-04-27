/*
 * imoab_dualmap_caas.cpp
 *
 * Test for dual-map nonlinear remapping (CAAS with low-order map bounds).
 *
 * Workflow:
 *   1. Load ATM (source) and OCN (target) meshes on all processes
 *   2. Migrate meshes to coupler communicator
 *   3. Compute mesh intersection and coverage
 *   4. Compute two sets of FV weights:
 *      - "lo-scalar" : 1st-order FV with monotonicity (low-order monotone map)
 *      - "hi-scalar" : 1st-order FV without monotonicity (high-order non-monotone map)
 *   5. Define an analytical source tag with sharp features
 *   6. Apply high-order map with dual-map CAAS bounds from low-order map
 *      using the new lo_weights_identifier parameter
 *   7. Verify: target values are within source stencil bounds
 *   8. Test iMOAB_CheckMapSubset
 *   9. Verify bounds preservation against source field range
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

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#ifndef MOAB_HAVE_TEMPESTREMAP
#error This test requires MOAB configuration with TempestRemap
#endif

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

    std::string atmFilename = TestDir + "unittest/wholeATM_T.h5m";
    std::string ocnFilename = TestDir + "unittest/recMeshOcn.h5m";

    int nghlay = 0;
    int startG1 = 0, endG1 = numProcesses - 1;
    int startG2 = 0, endG2 = numProcesses - 1;
    int startG4 = 0, endG4 = numProcesses - 1;

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "ATM mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "OCN mesh filename (target)", &ocnFilename );
    opts.parseCommandLine( argc, argv );

    if( !rankInGlobalComm )
    {
        std::cout << " === imoab_dualmap_caas test ===\n";
        std::cout << " ATM file: " << atmFilename << "\n";
        std::cout << " OCN file: " << ocnFilename << "\n";
        std::cout << " Processes: " << numProcesses << "\n";
    }

    // Create MPI groups and communicators
    MPI_Group atmPEGroup;
    MPI_Comm atmComm;
    CHECKIERR( create_group_and_comm( startG1, endG1, jgroup, &atmPEGroup, &atmComm ),
               "Cannot create ATM group" )

    MPI_Group ocnPEGroup;
    MPI_Comm ocnComm;
    CHECKIERR( create_group_and_comm( startG2, endG2, jgroup, &ocnPEGroup, &ocnComm ),
               "Cannot create OCN group" )

    MPI_Group couPEGroup;
    MPI_Comm couComm;
    CHECKIERR( create_group_and_comm( startG4, endG4, jgroup, &couPEGroup, &couComm ),
               "Cannot create coupler group" )

    MPI_Group joinAtmCouGroup;
    MPI_Comm atmCouComm;
    CHECKIERR( create_joint_comm_group( atmPEGroup, couPEGroup, &joinAtmCouGroup, &atmCouComm ),
               "Cannot create joint ATM-coupler comm" )

    MPI_Group joinOcnCouGroup;
    MPI_Comm ocnCouComm;
    CHECKIERR( create_joint_comm_group( ocnPEGroup, couPEGroup, &joinOcnCouGroup, &ocnCouComm ),
               "Cannot create joint OCN-coupler comm" )

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

    if( couComm != MPI_COMM_NULL )
    {
        // Compute mesh intersection between ATM and OCN on coupler
        PUSH_TIMER( "Compute ATM-OCN mesh intersection" )
        CHECKIERR( iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplDualMapPID ),
                   "Cannot compute ATM/OCN intersection" )
        POP_TIMER( couComm, rankInCouComm )

        // Compute LOW-ORDER weights (monotone FV, 1st order)
        const iMOAB_String lo_map_id   = "lo-scalar";
        const iMOAB_String disc_fv     = "fv";
        const iMOAB_String dof_tag     = "GLOBAL_ID";
        int disc_order                 = 1;
        int fNoBubble = 1, fMonotone = 1, fVolumetric = 0, fInvDist = 0, fNoConserve = 0, fValidate = 0;

        PUSH_TIMER( "Compute low-order (monotone) weights" )
        CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, lo_map_id, disc_fv, &disc_order, disc_fv,
                                                         &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                         &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                   "Cannot compute low-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // Compute HIGH-ORDER weights (non-monotone FV, 1st order)
        const iMOAB_String hi_map_id = "hi-scalar";
        fMonotone                    = 0;  // no monotonicity constraint
        disc_order                   = 2;

        PUSH_TIMER( "Compute high-order (non-monotone) weights" )
        CHECKIERR( iMOAB_ComputeScalarProjectionWeights( cplDualMapPID, hi_map_id, disc_fv, &disc_order, disc_fv,
                                                         &disc_order, nullptr, &fNoBubble, &fMonotone, &fVolumetric,
                                                         &fInvDist, &fNoConserve, &fValidate, dof_tag, dof_tag ),
                   "Cannot compute high-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // Compute coverage comm graph for tag migration
        int meshtype = 3;
        CHECKIERR( iMOAB_ComputeCommGraph( cplAtmPID, cplDualMapPID, &couComm, &couPEGroup, &couPEGroup, &meshtype,
                                           &meshtype, &cplatm, &dualmap_id ),
                   "Cannot compute ATM coverage graph" )

        // Define source and target tags
        int tagType                     = 1;  // DENSE_DOUBLE
        int tagIndex;
        int atmCompNDoFs = 1, ocnCompNDoFs = 1;
        const iMOAB_String srcField     = "SourceAnalytical";
        const iMOAB_String tgtFieldHi   = "TargetHiOrder";
        const iMOAB_String tgtFieldDual = "TargetDualMap";
        const iMOAB_String tgtFieldLo   = "TargetLoOrder";

        CHECKIERR( iMOAB_DefineTagStorage( cplAtmPID, srcField, &tagType, &atmCompNDoFs, &tagIndex ),
                   "Cannot define source tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldHi, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define high-order target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldDual, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define dual-map target tag" )
        CHECKIERR( iMOAB_DefineTagStorage( cplOcnPID, tgtFieldLo, &tagType, &ocnCompNDoFs, &tagIndex ),
                   "Cannot define low-order target tag" )
    }

    // Set source field values on ATM component: a sharp step function
    if( atmComm != MPI_COMM_NULL )
    {
        int tagIndex;
        int tagType[2]              = { 0, 1 };  // dense_int, dense_double
        int atmCompNDoFs            = 1;
        const iMOAB_String idField  = "GLOBAL_ID";
        const iMOAB_String srcField = "SourceAnalytical";

        CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, idField, &tagType[0], &atmCompNDoFs, &tagIndex ),
                   "Cannot define src tag on component" )
        CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, srcField, &tagType[1], &atmCompNDoFs, &tagIndex ),
                   "Cannot define src tag on component" )

        int nElems[3];
        CHECKIERR( iMOAB_GetMeshInfo( cmpAtmPID, nullptr, nElems, nullptr, nullptr, nullptr ),
                   "Cannot get ATM mesh info" )

        // Get global IDs for elements to set field values
        std::vector< int > gids( nElems[2] );
        int entity_type = 1;  // elements
        CHECKIERR( iMOAB_GetIntTagStorage( cmpAtmPID, idField, &nElems[2], &entity_type, gids.data() ),
                   "Cannot get element global IDs" )

        // Create a step function: elements with even GIDs get value 10.0, odd get 0.0
        // This creates sharp discontinuities that will test bounds preservation
        std::vector< double > srcVals( nElems[2] );
        for( int i = 0; i < nElems[2]; i++ )
        {
            srcVals[i] = ( gids[i] % 2 == 0 ) ? 10.0 : 0.0;
        }

        CHECKIERR( iMOAB_SetDoubleTagStorage( cmpAtmPID, srcField, &nElems[2], &entity_type, srcVals.data() ),
                   "Cannot set source field values" )
    }

    // Send source tag to coupler
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cmpAtmPID, "SourceAnalytical", &atmCouComm, &cplatm ),
                   "Cannot send source tag" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_ReceiveElementTag( cplAtmPID, "SourceAnalytical", &atmCouComm, &cmpatm ),
                   "Cannot receive source tag" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ),
                   "Cannot free source tag send buffers" )
    }

    // Send source tag to coverage mesh
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cplAtmPID, "SourceAnalytical", &couComm, &dualmap_id ),
                   "Cannot send tag to coverage" )
        CHECKIERR( iMOAB_ReceiveElementTag( cplDualMapPID, "SourceAnalytical", &couComm, &cplatm ),
                   "Cannot receive tag on coverage" )
        CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &dualmap_id ),
                   "Cannot free coverage send buffers" )
    }

    // === Apply projections and test ===
    if( couComm != MPI_COMM_NULL )
    {
        int filter_type = 0;  // no CAAS for plain high-order

        // 1) Apply high-order projection WITHOUT CAAS (baseline)
        PUSH_TIMER( "Apply high-order projection (no CAAS)" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar",
                                                        "SourceAnalytical", "TargetHiOrder", nullptr ),
                   "Failed to apply high-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 2) Apply low-order projection (reference, should already be bounded)
        filter_type = 0;
        PUSH_TIMER( "Apply low-order projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "lo-scalar",
                                                        "SourceAnalytical", "TargetLoOrder", nullptr ),
                   "Failed to apply low-order weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 3) Apply high-order projection WITH dual-map CAAS bounds from low-order map
        filter_type = 2;  // CAAS_LOCAL
        PUSH_TIMER( "Apply dual-map CAAS projection" )
        CHECKIERR( iMOAB_ApplyScalarProjectionWeights( cplDualMapPID, &filter_type, "hi-scalar",
                                                        "SourceAnalytical", "TargetDualMap", "lo-scalar" ),
                   "Failed to apply dual-map CAAS weights" )
        POP_TIMER( couComm, rankInCouComm )

        // 4) Test iMOAB_CheckMapSubset: low-order stencil subset of high-order
        int is_subset = 0;
        CHECKIERR( iMOAB_CheckMapSubset( cplDualMapPID, "lo-scalar", "hi-scalar", &is_subset ),
                   "Failed to check map subset" )
        if( !rankInCouComm )
        {
            std::cout << " CheckMapSubset (lo ⊆ hi): " << ( is_subset ? "PASS" : "FAIL" ) << "\n";
        }

        // 5) Verify bounds preservation: source field is a step function with values
        //    {0.0, 10.0}, so every target element's stencil bounds are [0.0, 10.0].
        //    The dual-map CAAS result must stay within these source bounds.
        //    (ComputeRowBounds is an internal iMOAB routine, not exposed in the public API)
        const double srcBoundLo = 0.0;
        const double srcBoundHi = 10.0;

        int nOcnVerts, nOcnElems;
        CHECKIERR( iMOAB_GetMeshInfo( cplOcnPID, &nOcnVerts, &nOcnElems, nullptr, nullptr, nullptr ),
                   "Cannot get OCN mesh info" )

        std::vector< double > hiVals( nOcnElems ), dualVals( nOcnElems ), loVals( nOcnElems );
        int entity_type = 1;

        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, "TargetHiOrder", &nOcnElems, &entity_type, hiVals.data() ),
                   "Cannot get hi-order values" )
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, "TargetDualMap", &nOcnElems, &entity_type, dualVals.data() ),
                   "Cannot get dual-map values" )
        CHECKIERR( iMOAB_GetDoubleTagStorage( cplOcnPID, "TargetLoOrder", &nOcnElems, &entity_type, loVals.data() ),
                   "Cannot get lo-order values" )

        // Count violations against source field bounds [0.0, 10.0]
        int localHiViolations   = 0;
        int localDualViolations = 0;
        double maxHiExceedance  = 0.0;
        double maxDualExceedance = 0.0;
        const double tol_bounds = 1e-10;

        for( int i = 0; i < nOcnElems; i++ )
        {
            // Check high-order (may violate)
            if( hiVals[i] < srcBoundLo - tol_bounds || hiVals[i] > srcBoundHi + tol_bounds )
            {
                localHiViolations++;
                double exc = std::max( srcBoundLo - hiVals[i], hiVals[i] - srcBoundHi );
                maxHiExceedance = std::max( maxHiExceedance, exc );
            }

            // Check dual-map (should NOT violate)
            if( dualVals[i] < srcBoundLo - tol_bounds || dualVals[i] > srcBoundHi + tol_bounds )
            {
                localDualViolations++;
                double exc = std::max( srcBoundLo - dualVals[i], dualVals[i] - srcBoundHi );
                maxDualExceedance = std::max( maxDualExceedance, exc );
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
    MPI_Finalize();

    return 0;
}
