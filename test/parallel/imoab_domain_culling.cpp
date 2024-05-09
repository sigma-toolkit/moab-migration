/*
 * imoab_domain_culling.cpp
 *  Created on: Jan. 16, 2024
 *
 *   This imoab_domain_culling test will test reading of a domain with or without culling , and see what is the
 * effect on intersection and map generation
 * basically, 1 moab file will read the domain with or without culling, and see what is the effect on intersection
 * and projection of a field
 * in a data model for ocean, the ocean domain file is read directly from an nc domain file;
 */

#include "moab/Core.hpp"

// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"

#include "moab/iMOAB.h"
#include "TestUtil.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ProgOptions.hpp"
#include <iostream>
#include <sstream>
// CHECKIERR
#include "imoab_coupler_utils.hpp"

using namespace moab;

int main( int argc, char* argv[] )
{
    int ierr;
    int rankInGlobalComm, numProcesses;
    std::string readoptsDomain( "PARALLEL=READ_PART;PARTITION_METHOD=SQIJ;NO_CULLING" );
    // options for a regular atmosphere file
    std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

    // Timer data
    moab::CpuTimer timer;
    std::string opName;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rankInGlobalComm );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcesses );

    MPI_Comm dup_comm_world;
    MPI_Comm_dup( MPI_COMM_WORLD, &dup_comm_world );

    MPI_Group mpigrp_CPLID;
    MPI_Comm_group( dup_comm_world, &mpigrp_CPLID );

    std::string domainFile = TestDir + "unittest/io/domain.ocn.ne4np4_oQU240.160614.nc";
    std::string atmFile    = TestDir + "unittest/atm_c2x.h5m";
    std::string tagname( "Sa_pbot:Sa_dens" );
    int nghlay = 0;  // no ghosts

    ierr = iMOAB_Initialize( argc, argv );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKIERR( ierr, "Cannot initialize iMOAB" )

    ProgOptions opts;
    opts.addOpt< std::string >( "source,s", "atm mesh filename (source)", &atmFile );
    opts.addOpt< std::string >( "target,t", "domain filename (target)", &domainFile );
    opts.addOpt< std::string >( "tags,g", "list of tags, colon separated", &tagname );

    opts.parseCommandLine( argc, argv );

    int cplAtmAppID       = -1;            // -1 means it is not initialized
    iMOAB_AppID cplAtmPID = &cplAtmAppID;  // atm on coupler PEs, moab dist

    int cplOcnAppID = -1, cplAtmOcnAppID = -1;   // -1 means it is not initialized
    iMOAB_AppID cplOcnPID    = &cplOcnAppID;     // ocn on coupler PEs
    iMOAB_AppID cplAtmOcnPID = &cplAtmOcnAppID;  // intx atm -ocn on coupler PEs

    int cplatm = 6;
    ierr       = iMOAB_RegisterApplication( "ATMX", &dup_comm_world, &cplatm,
                                            cplAtmPID );  // mesh on coupler pes
    CHECKIERR( ierr, "Cannot register ATMX over coupler PEs" )

    ierr = iMOAB_LoadMesh( cplAtmPID, atmFile.c_str(), readopts.c_str(), &nghlay );  // moab mesh can be cells
    CHECKIERR( ierr, "Cannot load moab mesh on coupler pes" )

    int cplocn = 18;
    ierr = iMOAB_RegisterApplication( "OCNX", &dup_comm_world, &cplocn,
                                      cplOcnPID );  // ocn on coupler pes
    CHECKIERR( ierr, "Cannot register OCN over coupler PEs" )
    // use domain file options
    ierr = iMOAB_LoadMesh( cplOcnPID, domainFile.c_str(), readoptsDomain.c_str(), &nghlay );  // moab mesh can be cells
    CHECKIERR( ierr, "Cannot load ocean domain mesh on coupler pes" )

    char fileWriteOptions[] = "PARALLEL=WRITE_PART";
    char outputFile[]       = "OcnDomMesh.h5m";
    ierr                    = iMOAB_WriteMesh( cplOcnPID, outputFile, fileWriteOptions );
    CHECKIERR( ierr, "Cannot write ocean domain mesh from coupler pes" )

    int sizeTag    = 1;
    int tagIndex   = -1;
    int tagType = DENSE_DOUBLE;
    ierr = iMOAB_DefineTagStorage( cplAtmPID, tagname.c_str(), &tagType, &sizeTag, &tagIndex );
    CHECKIERR( ierr, "Cannot define source tags tag" )

    ierr = iMOAB_DefineTagStorage( cplOcnPID, tagname.c_str(), &tagType, &sizeTag, &tagIndex );
    CHECKIERR( ierr, "Cannot define target tags tag" )

    int atmocnid = 618;
    // now compute intersection between OCNx and ATMx on coupler PEs
    ierr = iMOAB_RegisterApplication( "ATMOCN", &dup_comm_world, &atmocnid, cplAtmOcnPID );
    CHECKIERR( ierr, "Cannot register ocn_atm intx over coupler pes " )

    ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID, cplAtmOcnPID );
    // coverage mesh was computed here, for cplAtmPID, atm on coupler pes
    // basically, atm was redistributed according to target (ocean) partition, to "cover" the
    // ocean partitions check if intx valid, write some h5m intx file
    CHECKIERR( ierr, "cannot compute intersection" )

#ifdef VERBOSE
    char prefix[] = "intx_atmocn";
    ierr          = iMOAB_WriteLocalMesh( cplAtmOcnPID, prefix, strlen( prefix ) );
    CHECKIERR( ierr, "failed to write local intx mesh" );
#endif

    int type1 = 3, type2 = 3;
    ierr = iMOAB_ComputeCommGraph( cplAtmPID, cplAtmOcnPID, &dup_comm_world, &mpigrp_CPLID, &mpigrp_CPLID, &type1,
                                   &type2, &cplatm, &atmocnid );
    CHECKIERR( ierr, "failed to compute comm graph" );
    const std::string weights_identifiers[2] = { "scalar", "scalar-pc" };
    const std::string disc_methods[3]        = { "cgll", "fv", "pcloud" };
    const std::string dof_tag_names[3]       = { "GLOBAL_DOFS", "GLOBAL_ID", "GLOBAL_ID" };

    int fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 0, fNoConserve = 0, fNoBubble = 1, fInverseDistanceMap = 0;
    int filter_type = 0;
    int disc_orders[3] = { 4, 1, 1 };
    ierr = iMOAB_ComputeScalarProjectionWeights( cplAtmOcnPID, weights_identifiers[0].c_str(), disc_methods[1].c_str(),
                                                 &disc_orders[1], disc_methods[1].c_str(), &disc_orders[1], nullptr,
                                                 &fNoBubble, &fMonotoneTypeID, &fVolumetric, &fInverseDistanceMap,
                                                 &fNoConserve, &fValidate, dof_tag_names[1].c_str(),
                                                 dof_tag_names[1].c_str() );

    CHECKIERR( ierr, "cannot compute scalar projection weights" )

    const std::string atmocn_map_file_name = "atmDOCN_map.nc";
    ierr =
        iMOAB_WriteMappingWeightsToFile( cplAtmOcnPID, weights_identifiers[0].c_str(), atmocn_map_file_name.c_str() );
    CHECKIERR( ierr, "failed to write map file to disk" );

    // as always, use nonblocking sends
    // this is for projection to ocean:
    ierr = iMOAB_SendElementTag( cplAtmPID, tagname.c_str(), &dup_comm_world, &atmocnid );
    CHECKIERR( ierr, "cannot send tag values" )

    // receive on atm on coupler pes, that was redistributed according to coverage
    // receive in the coverage mesh, basically
    ierr = iMOAB_ReceiveElementTag( cplAtmOcnPID, tagname.c_str(), &dup_comm_world, &cplatm );
    CHECKIERR( ierr, "cannot receive tag values" )

    ierr = iMOAB_FreeSenderBuffers( cplAtmPID, &atmocnid );  // context is intx external id
    CHECKIERR( ierr, "cannot free buffers used to resend atm tag towards the coverage mesh" )

    ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, &filter_type, weights_identifiers[0].c_str(), tagname.c_str(),
            tagname.c_str() );
    CHECKIERR( ierr, "failed to compute projection weight application" );

    char outputFile2[]       = "OcnDomMeshProj.h5m";
    ierr                    = iMOAB_WriteMesh( cplOcnPID, outputFile2, fileWriteOptions );
    CHECKIERR( ierr, "Cannot write ocean domain mesh from coupler pes" )

    ierr = iMOAB_DeregisterApplication( cplAtmOcnPID );
    CHECKIERR( ierr, "cannot deregister app LNDX2" )
    ierr = iMOAB_DeregisterApplication( cplOcnPID );
    CHECKIERR( ierr, "cannot deregister app LNDX2" )
    ierr = iMOAB_DeregisterApplication( cplAtmPID );
    CHECKIERR( ierr, "cannot deregister app LNDX" )

    ierr = iMOAB_Finalize();
    CHECKIERR( ierr, "did not finalize iMOAB" )

    MPI_Comm_free( &dup_comm_world );
    MPI_Finalize();
    return 0;
}
