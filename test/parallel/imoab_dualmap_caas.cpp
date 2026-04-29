/*
 * imoab_dualmap_caas.cpp
 *
 * Test for dual-map nonlinear remapping (CAAS with low-order map bounds).
 *
 * Workflow:
 *   1. Load ATM (source) and OCN (target) meshes on all processes
 *   2. Migrate meshes to coupler communicator
 *   3. Either compute mesh intersection and weight maps online, or
 *      load pre-computed weight maps from disk (--lo_map_file, --hi_map_file)
 *   4. Two sets of FV weights are used:
 *      - "lo-scalar" : low-order monotone map
 *      - "hi-scalar" : high-order non-monotone map
 *   5. Define an analytical source field (degree-2 spherical harmonic)
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
    std::string loMapFile;  // empty = compute online
    std::string hiMapFile;  // empty = compute online

    int nghlay = 0;

    // PE layout: default all tasks on all groups
    int startG1 = 0, endG1 = numProcesses - 1;  // ATM
    int startG2 = 0, endG2 = numProcesses - 1;  // OCN
    int startG4 = 0, endG4 = numProcesses - 1;  // Coupler

    ProgOptions opts;
    opts.addOpt< std::string >( "atmosphere,t", "ATM mesh filename (source)", &atmFilename );
    opts.addOpt< std::string >( "ocean,m", "OCN mesh filename (target)", &ocnFilename );
    opts.addOpt< std::string >( "lo_map_file,l", "Low-order map file (nc); if set, load from disk", &loMapFile );
    opts.addOpt< std::string >( "hi_map_file,h", "High-order map file (nc); if set, load from disk", &hiMapFile );
    opts.addOpt< int >( "startAtm,a", "start task for atmosphere layout", &startG1 );
    opts.addOpt< int >( "endAtm,b", "end task for atmosphere layout", &endG1 );
    opts.addOpt< int >( "startOcn,c", "start task for ocean layout", &startG2 );
    opts.addOpt< int >( "endOcn,d", "end task for ocean layout", &endG2 );
    opts.addOpt< int >( "startCoupler,g", "start task for coupler layout", &startG4 );
    opts.addOpt< int >( "endCoupler,j", "end task for coupler layout", &endG4 );
    opts.parseCommandLine( argc, argv );

    bool loadFromDisk = ( !loMapFile.empty() && !hiMapFile.empty() );

    if( !rankInGlobalComm )
    {
        std::cout << " === imoab_dualmap_caas test ===\n";
        std::cout << " ATM file: " << atmFilename << "\n";
        std::cout << " OCN file: " << ocnFilename << "\n";
        if( loadFromDisk )
        {
            std::cout << " Lo-order map: " << loMapFile << "\n";
            std::cout << " Hi-order map: " << hiMapFile << "\n";
            std::cout << " Mode: load maps from disk\n";
        }
        else
        {
            std::cout << " Mode: compute maps online\n";
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

    const iMOAB_String srcField     = "SourceAnalytical";
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
            int arearead      = 0;  // do not read areas

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
            int nghlay = 0;
            int nghlay_tgt = 0;
            CHECKIERR( iMOAB_SetMapGhostLayers( cplAtmPID, &nghlay, &nghlay_tgt ),
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

    // Set source field values on ATM component: spherical harmonic evaluated at element centroids.
    // f(x,y,z) = 1.0 + 0.5*(3z^2 - 1) + 0.8*(x^2 - y^2)
    //          = 1.0 + P_2(z) + 0.8*Y_2^2(x,y)     [unnormalized]
    // This is a smooth degree-2 polynomial on the unit sphere that exercises the
    // remapping well and produces deterministic results independent of mesh partitioning.
    double localSrcMin = 1e308, localSrcMax = -1e308;

    if( atmComm != MPI_COMM_NULL )
    {
        int tagIndex;
        int tagType_dbl = 1;  // DENSE_DOUBLE
        int atmCompNDoFs = 1;

        CHECKIERR( iMOAB_DefineTagStorage( cmpAtmPID, srcField, &tagType_dbl, &atmCompNDoFs, &tagIndex ),
                   "Cannot define src tag on component" )

        int nVerts[3], nElems[3], nBlocks[3];
        CHECKIERR( iMOAB_GetMeshInfo( cmpAtmPID, nVerts, nElems, nBlocks, nullptr, nullptr ),
                   "Cannot get ATM mesh info" )

        // Get vertex coordinates (interleaved x,y,z)
        int coordsLen = nVerts[2] * 3;
        std::vector< double > coords( coordsLen );
        CHECKIERR( iMOAB_GetVisibleVerticesCoordinates( cmpAtmPID, &coordsLen, coords.data() ),
                   "Cannot get ATM vertex coordinates" )

        // Get block IDs
        std::vector< int > blockIDs( nBlocks[2] );
        CHECKIERR( iMOAB_GetBlockID( cmpAtmPID, &nBlocks[2], blockIDs.data() ),
                   "Cannot get block IDs" )

        // Compute element centroids and evaluate spherical harmonic
        std::vector< double > srcVals( nElems[2] );
        int elemOffset = 0;

        for( int b = 0; b < nBlocks[2]; b++ )
        {
            int vertsPerElem, numElemsInBlock;
            CHECKIERR( iMOAB_GetBlockInfo( cmpAtmPID, &blockIDs[b], &vertsPerElem, &numElemsInBlock ),
                       "Cannot get block info" )

            int connLen = vertsPerElem * numElemsInBlock;
            std::vector< int > conn( connLen );
            CHECKIERR( iMOAB_GetBlockElementConnectivities( cmpAtmPID, &blockIDs[b], &connLen, conn.data() ),
                       "Cannot get element connectivity" )

            for( int e = 0; e < numElemsInBlock; e++ )
            {
                // Compute centroid of this element
                double cx = 0.0, cy = 0.0, cz = 0.0;
                for( int v = 0; v < vertsPerElem; v++ )
                {
                    int vidx = conn[e * vertsPerElem + v] - 1;  // 1-based to 0-based
                    cx += coords[3 * vidx + 0];
                    cy += coords[3 * vidx + 1];
                    cz += coords[3 * vidx + 2];
                }
                cx /= vertsPerElem;
                cy /= vertsPerElem;
                cz /= vertsPerElem;

                // Normalize to unit sphere (in case mesh radius != 1)
                double r = std::sqrt( cx * cx + cy * cy + cz * cz );
                if( r > 1e-14 )
                {
                    cx /= r;
                    cy /= r;
                    cz /= r;
                }

                // Spherical harmonic: f = 1.0 + 0.5*(3z^2 - 1) + 0.8*(x^2 - y^2)
                double val = 1.0 + 0.5 * ( 3.0 * cz * cz - 1.0 ) + 0.8 * ( cx * cx - cy * cy );
                srcVals[elemOffset + e] = val;

                localSrcMin = std::min( localSrcMin, val );
                localSrcMax = std::max( localSrcMax, val );
            }
            elemOffset += numElemsInBlock;
        }

        int entity_type = 1;  // elements
        CHECKIERR( iMOAB_SetDoubleTagStorage( cmpAtmPID, srcField, &nElems[2], &entity_type, srcVals.data() ),
                   "Cannot set source field values" )
    }

    // Get global min/max of source field for bounds checking later
    double globalSrcMin, globalSrcMax;
    MPI_Allreduce( &localSrcMin, &globalSrcMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( &localSrcMax, &globalSrcMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

    if( !rankInGlobalComm )
    {
        std::cout << " Source field range: [" << globalSrcMin << ", " << globalSrcMax << "]\n";
    }

    // Send source tag to coupler
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cmpAtmPID, srcField, &atmCouComm, &cplatm ),
                   "Cannot send source tag" )
    }
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_ReceiveElementTag( cplAtmPID, srcField, &atmCouComm, &cmpatm ),
                   "Cannot receive source tag" )
    }
    if( atmComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_FreeSenderBuffers( cmpAtmPID, &cplatm ),
                   "Cannot free source tag send buffers" )
    }

    // Send source tag from coupler-atm coverage to dual-map coverage so that
    // the high/low order maps can apply weights against actual source values
    // (otherwise the dual-map app's source tag is the default fill value).
    if( couComm != MPI_COMM_NULL )
    {
        CHECKIERR( iMOAB_SendElementTag( cplAtmPID, srcField, &couComm, &dualmap_id ),
                   "Cannot send tag to coverage" )
        CHECKIERR( iMOAB_ReceiveElementTag( cplDualMapPID, srcField, &couComm, &cplatm ),
                   "Cannot receive tag on coverage" )
        CHECKIERR( iMOAB_FreeSenderBuffers( cplAtmPID, &dualmap_id ),
                   "Cannot free coverage send buffers" )
    }

    if( couComm != MPI_COMM_NULL  )
    {
        // Pre-projection snapshot: dual-map coverage with the migrated source field.
        // Projected target fields live on cplOcnPID and are written below after
        // projection; writing this here only captures the source-field state.
        char outputFileRecvd[] = "cplAtmFile.h5m";
        char fileWriteOptions[] = "PARALLEL=WRITE_PART";
        CHECKIERR( iMOAB_WriteMesh( cplDualMapPID, outputFileRecvd, fileWriteOptions ),
                   "could not write cplAtmFile.h5m to disk" )
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

        // 4) Test iMOAB_CheckMapSubset: low-order stencil subset of high-order
        int is_subset = 0;
        CHECKIERR( iMOAB_CheckMapSubset( cplDualMapPID, "lo-scalar", "hi-scalar", &is_subset ),
                   "Failed to check map subset" )
        if( !rankInCouComm )
        {
            std::cout << " CheckMapSubset (lo ⊆ hi): " << ( is_subset ? "PASS" : "FAIL" ) << "\n";
        }

        // 5) Verify bounds preservation: the dual-map CAAS result must stay within
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
