/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning( disable : 4786 )
#endif
#endif

#include <iostream>
#include "moab/Core.hpp"
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include "Internals.hpp"
#include "TestRunner.hpp"

using namespace moab;

static const double radius               = 1.0;
const double MOAB_PI                     = 3.1415926535897932384626433832795028841971693993751058209749445923;
static const double surface_area         = 4.0 * MOAB_PI * radius * radius;
static const std::string outFilenames[5] = { "outTempestCS.g", "outTempestRLL.g", "outTempestICO.g", "outTempestICOD.g",
                                             "outTempestOV.g" };

void test_tempest_cs_create();
void test_tempest_rll_create();
void test_tempest_ico_create();
void test_tempest_mpas_create();
void test_tempest_overlap_combinations();
void test_tempest_to_moab_convert();
#ifdef MOAB_HAVE_PNETCDF
void test_pnetcdf_map_roundtrip();
void test_pnetcdf_parallel_read();
#endif

int main( int argc, char** argv )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    // The mesh-generation tests below write TempestRemap meshes in NetCDF formats
    // (GenerateCSMesh/RLLMesh/... with "NetCDF4"), so they only run when NetCDF is enabled.
#ifdef MOAB_HAVE_NETCDF
    REGISTER_TEST( test_tempest_cs_create );
    REGISTER_TEST( test_tempest_rll_create );
    REGISTER_TEST( test_tempest_ico_create );
    REGISTER_TEST( test_tempest_mpas_create );
    REGISTER_TEST( test_tempest_overlap_combinations );
    REGISTER_TEST( test_tempest_to_moab_convert );
#endif

    // Full SCRIP map I/O round-trip via PnetCDF (no NetCDF required): generate meshes
    // in-memory, compute weights, write + read the map with PnetCDF, and transfer a field.
#ifdef MOAB_HAVE_PNETCDF
    REGISTER_TEST( test_pnetcdf_map_roundtrip );
    REGISTER_TEST( test_pnetcdf_parallel_read );
#endif

    int result = RUN_TESTS( argc, argv );
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return result;
}

// All tests below generate/write/load TempestRemap meshes in NetCDF formats and are only
// compiled (and run) when NetCDF is available; otherwise main() above skips the suite.
#ifdef MOAB_HAVE_NETCDF

void test_tempest_cs_create()
{
    NcError error( NcError::verbose_nonfatal );
    const int blockSize           = 30;
    const std::string outFilename = outFilenames[0];

    std::cout << "Creating TempestRemap Cubed-Sphere Mesh ...\n";
    Mesh tempest_mesh;
    int ierr = GenerateCSMesh( tempest_mesh, blockSize, outFilename, "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of CS mesh
    const double sphere_area = tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( sphere_area, surface_area, 1e-10 );
}

void test_tempest_rll_create()
{
#ifdef MOAB_HAVE_NETCDF
    NcError error( NcError::verbose_nonfatal );
#endif
    const int blockSize           = 30;
    const std::string outFilename = outFilenames[1];

    std::cout << "Creating TempestRemap Latitude-Longitude Mesh ...\n";
    Mesh tempest_mesh;
    int ierr =
        GenerateRLLMesh( tempest_mesh, blockSize * 2, blockSize, 0.0, 360.0, -90.0, 90.0, false, false, true, "", "",
                         "",  // std::string strInputFile, std::string strInputFileLonName, std::string
                              // strInputFileLatName,
                         outFilename, "NetCDF4",  // std::string strOutputFile, std::string strOutputFormat
                         false );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of RLL mesh
    const double sphere_area = tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( sphere_area, surface_area, 1e-10 );
}

void test_tempest_ico_create()
{
#ifdef MOAB_HAVE_NETCDF
    NcError error( NcError::verbose_nonfatal );
#endif
    const int blockSize           = 30;
    const bool computeDual        = false;
    const std::string outFilename = outFilenames[2];

    std::cout << "Creating TempestRemap Icosahedral Mesh ...\n";
    Mesh tempest_mesh;
    int ierr = GenerateICOMesh( tempest_mesh, blockSize, computeDual, outFilename, "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of ICO mesh
    const double sphere_area = tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( sphere_area, surface_area, 1e-10 );
}

void test_tempest_mpas_create()
{
#ifdef MOAB_HAVE_NETCDF
    NcError error( NcError::verbose_nonfatal );
#endif
    const int blockSize           = 30;
    const bool computeDual        = true;
    const std::string outFilename = outFilenames[3];

    std::cout << "Creating TempestRemap MPAS Mesh (dual of the Icosahedral) ...\n";
    Mesh tempest_mesh;
    int ierr = GenerateICOMesh( tempest_mesh, blockSize, computeDual, outFilename, "NetCDF4" );
    CHECK_EQUAL( ierr, 0 );

    // Compute the surface area of MPAS mesh
    const double sphere_area = tempest_mesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( sphere_area, surface_area, 1e-10 );
}

void test_tempest_overlap_combinations()
{
#ifdef MOAB_HAVE_NETCDF
    NcError error( NcError::verbose_nonfatal );
#endif
    const std::string outFilename = outFilenames[4];

    Mesh inpMesh( outFilenames[0] );
    // verify input mesh area first
    const double inpArea = inpMesh.CalculateFaceAreas( false );
    CHECK_REAL_EQUAL( inpArea, surface_area, 1e-10 );

    for( int isrc = 0; isrc < 4; ++isrc )
    {
        for( int jsrc = 0; jsrc < 4; ++jsrc )
        {
            std::cout << "Computing Overlap between " << outFilenames[isrc] << " and " << outFilenames[jsrc]
                      << " ...\n";
            Mesh tempest_mesh;
            int ierr = GenerateOverlapMesh( outFilenames[isrc], outFilenames[jsrc], tempest_mesh, outFilename,
                                            "NetCDF4", "exact", false, false, false, false, false );
            CHECK_EQUAL( ierr, 0 );
            // verify overlap mesh area
            const double ovArea = tempest_mesh.CalculateFaceAreas( false );
            CHECK_REAL_EQUAL( ovArea, surface_area, 1e-10 );
        }
    }
}

void test_tempest_to_moab_convert()
{
#ifdef MOAB_HAVE_NETCDF
    NcError error( NcError::verbose_nonfatal );
#endif

    // Allocate and create MOAB Remapper object
    moab::ErrorCode rval;
    moab::Interface* mbCore = new( std::nothrow ) moab::Core;
    CHECK( NULL != mbCore );

#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm       = new moab::ParallelComm( mbCore, MPI_COMM_WORLD, 0 );
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mbCore, pcomm );
#else
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mbCore );
#endif
    remapper->meshValidate     = true;
    remapper->constructEdgeMap = true;
    remapper->initialize();

#ifdef MOAB_HAVE_MPI
    rval = pcomm->check_all_shared_handles();CHECK_ERR( rval );
#endif

    rval = remapper->LoadMesh( moab::Remapper::SourceMesh, outFilenames[0], moab::TempestRemapper::CS );CHECK_ERR( rval );

    // Load the meshes and validate
    rval = remapper->ConvertTempestMesh( moab::Remapper::SourceMesh );CHECK_ERR( rval );

    Mesh* srcTempest = remapper->GetMesh( moab::Remapper::SourceMesh );

    moab::EntityHandle srcset = remapper->GetMeshSet( moab::Remapper::SourceMesh );

    moab::EntityHandle& tgtset = remapper->GetMeshSet( moab::Remapper::TargetMesh );

    tgtset = srcset;

    // Load the meshes and validate
    rval = remapper->ConvertMeshToTempest( moab::Remapper::TargetMesh );CHECK_ERR( rval );

    Mesh* tgtTempest = remapper->GetMesh( moab::Remapper::TargetMesh );

    const size_t tempest_nodes_src = srcTempest->nodes.size(), tempest_elems_src = srcTempest->faces.size();
    const size_t tempest_nodes_tgt = tgtTempest->nodes.size(), tempest_elems_tgt = tgtTempest->faces.size();
    CHECK_EQUAL( tempest_nodes_src, tempest_nodes_tgt );
    CHECK_EQUAL( tempest_elems_src, tempest_elems_tgt );

    delete remapper;
#ifdef MOAB_HAVE_MPI
    delete pcomm;
#endif
    delete mbCore;
}

#endif  // MOAB_HAVE_NETCDF

#ifdef MOAB_HAVE_PNETCDF
// End-to-end PnetCDF round-trip that verifies the whole toolchain without NetCDF:
//   (1) generate CS5 (source) + ICOD5 (target) meshes in memory (empty filename => no
//       file write, so no NetCDF dependency),
//   (2) compute FV-FV remapping weights with TempestOnlineMap,
//   (3) write the SCRIP (.nc) weight map to disk using the PnetCDF writer,
//   (4) read it back with the PnetCDF reader and confirm the dimensions round-trip,
//   (5) transfer a unit source field through the map and confirm it maps to a unit
//       target field (consistency of a conservative FV-FV map).
// The remapper is built on MPI_COMM_SELF and each rank writes a rank-unique file, so the
// test is correct for any number of MPI processes (each performs an independent round-trip).
void test_pnetcdf_map_roundtrip()
{
    moab::ErrorCode rval;
    int wrank = 0;
#ifdef MOAB_HAVE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
#endif

    // (1) In-memory meshes. The remapper takes ownership (its destructor deletes them),
    //     so heap-allocate and do not delete here.
    Mesh* srcMesh = new Mesh();
    Mesh* tgtMesh = new Mesh();
    CHECK_EQUAL( 0, GenerateCSMesh( *srcMesh, 5, "", "" ) );          // CS5   source
    CHECK_EQUAL( 0, GenerateICOMesh( *tgtMesh, 5, true, "", "" ) );   // ICOD5 target (dual)

    // (2) Remapper on MPI_COMM_SELF (independent per-rank serial round-trip)
    moab::Interface* mb = new( std::nothrow ) moab::Core;
    CHECK( NULL != mb );
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm       = new moab::ParallelComm( mb, MPI_COMM_SELF );
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mb, pcomm );
#else
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mb );
#endif
    remapper->meshValidate     = true;
    remapper->constructEdgeMap = false;
    rval                       = remapper->initialize();CHECK_ERR( rval );

    remapper->SetMesh( moab::Remapper::SourceMesh, srcMesh );
    rval = remapper->ConvertTempestMesh( moab::Remapper::SourceMesh );CHECK_ERR( rval );
    remapper->SetMesh( moab::Remapper::TargetMesh, tgtMesh );
    rval = remapper->ConvertTempestMesh( moab::Remapper::TargetMesh );CHECK_ERR( rval );

    rval = remapper->ConstructCoveringSet( 1e-8, 1.0, 1.0, 0.1, false, true, 0 );CHECK_ERR( rval );
    rval = remapper->ComputeOverlapMesh( true, false );CHECK_ERR( rval );

    // (3) Compute FV-FV weights
    moab::TempestOnlineMap* onlinemap = new moab::TempestOnlineMap( remapper );
    GenerateOfflineMapAlgorithmOptions mapOptions;
    mapOptions.nPin            = 1;
    mapOptions.nPout           = 1;
    mapOptions.fNoConservation = false;
    mapOptions.fMonotone       = false;
    mapOptions.fNoCheck        = true;
    rval = onlinemap->GenerateRemappingWeights( "fv", "fv", mapOptions );CHECK_ERR( rval );

    const int nSrc = onlinemap->GetSourceGlobalNDofs();
    const int nTgt = onlinemap->GetDestinationGlobalNDofs();
    CHECK( nSrc > 0 );
    CHECK( nTgt > 0 );
    std::cout << "[roundtrip rank " << wrank << "] computed map: nA=" << nSrc << " nB=" << nTgt << "\n";

    // (4) Write the SCRIP map via PnetCDF
    std::stringstream fname;
    fname << "pnetcdf_roundtrip_map_" << wrank << ".nc";
    std::map< std::string, std::string > attrMap;
    attrMap["Title"]   = "MOAB PnetCDF round-trip test map";
    attrMap["Creator"] = "test_remapping::test_pnetcdf_map_roundtrip";
    rval               = onlinemap->WriteParallelMap( fname.str(), attrMap );CHECK_ERR( rval );

    // (5) Transfer a unit source field through the map: consistent map => unit target field
    moab::EntityHandle srcSet = remapper->GetMeshSet( moab::Remapper::SourceMesh );
    moab::EntityHandle tgtSet = remapper->GetMeshSet( moab::Remapper::TargetMesh );
    moab::Range srcCells, tgtCells;
    rval = mb->get_entities_by_dimension( srcSet, 2, srcCells );CHECK_ERR( rval );
    rval = mb->get_entities_by_dimension( tgtSet, 2, tgtCells );CHECK_ERR( rval );

    moab::Tag srcTag, tgtTag;
    double defval = 0.0;
    rval          = mb->tag_get_handle( "src_field", 1, moab::MB_TYPE_DOUBLE, srcTag,
                                        moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defval );CHECK_ERR( rval );
    rval = mb->tag_get_handle( "tgt_field", 1, moab::MB_TYPE_DOUBLE, tgtTag,
                               moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, &defval );CHECK_ERR( rval );
    std::vector< double > ones( srcCells.size(), 1.0 );
    rval = mb->tag_set_data( srcTag, srcCells, ones.data() );CHECK_ERR( rval );

    rval = onlinemap->ApplyWeights( srcTag, tgtTag );CHECK_ERR( rval );

    std::vector< double > tvals( tgtCells.size(), 0.0 );
    rval          = mb->tag_get_data( tgtTag, tgtCells, tvals.data() );CHECK_ERR( rval );
    double tmin = tvals.empty() ? 0.0 : tvals[0], tmax = tvals.empty() ? 0.0 : tvals[0];
    for( size_t i = 0; i < tvals.size(); ++i )
    {
        if( tvals[i] < tmin ) tmin = tvals[i];
        if( tvals[i] > tmax ) tmax = tvals[i];
    }
    std::cout << "[roundtrip rank " << wrank << "] field transfer target min/max = " << tmin << " / " << tmax << "\n";
    CHECK_REAL_EQUAL( 1.0, tmin, 1e-6 );
    CHECK_REAL_EQUAL( 1.0, tmax, 1e-6 );

    // (6) Read the map back via PnetCDF and confirm dimensions round-trip
    moab::Interface* mb2 = new( std::nothrow ) moab::Core;
    CHECK( NULL != mb2 );
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm2       = new moab::ParallelComm( mb2, MPI_COMM_SELF );
    moab::TempestRemapper* remapper2 = new moab::TempestRemapper( mb2, pcomm2 );
#else
    moab::TempestRemapper* remapper2 = new moab::TempestRemapper( mb2 );
#endif
    rval = remapper2->initialize( false );CHECK_ERR( rval );
    moab::TempestOnlineMap* readmap = new moab::TempestOnlineMap( remapper2 );
    std::vector< int > tgt_owned_ids;
    std::vector< double > areaA, areaB;
    int rnA = 0, rnB = 0;
    rval    = readmap->ReadParallelMap( fname.str().c_str(), tgt_owned_ids, 0, areaA, rnA, areaB, rnB );CHECK_ERR( rval );
    std::cout << "[roundtrip rank " << wrank << "] read-back map: nA=" << rnA << " nB=" << rnB << "\n";
    CHECK_EQUAL( nSrc, rnA );
    CHECK_EQUAL( nTgt, rnB );

    delete readmap;
    delete remapper2;
#ifdef MOAB_HAVE_MPI
    delete pcomm2;
#endif
    delete mb2;

    delete onlinemap;
    delete remapper;  // deletes srcMesh and tgtMesh (owned)
#ifdef MOAB_HAVE_MPI
    delete pcomm;
#endif
    delete mb;
}

// Multi-rank parallel read through the PnetCDF reader. Rank 0 computes an FV-FV map
// (CS5->ICOD5) in-memory and writes it as a shared SCRIP file with the PnetCDF writer
// (CDF-5, classic family). All ranks then read it back collectively with ReadParallelMap.
// With n_s = 8976 (< buffered-read threshold) this exercises the rank-0 PnetCDF dimension
// probe + buffered read + MPI scatter across all ranks. (The stock reference map
// outCS5ICOD5_map.nc is HDF5-based NetCDF-4, which PnetCDF cannot read, so we generate a
// PnetCDF/CDF map here instead.)
void test_pnetcdf_parallel_read()
{
    moab::ErrorCode rval;
    int wrank = 0, wsize = 1;
#ifdef MOAB_HAVE_MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &wrank );
    MPI_Comm_size( MPI_COMM_WORLD, &wsize );
#endif
    const std::string mapfile = "pnetcdf_shared_map.nc";

    // --- Rank 0: compute + write a CDF-5 SCRIP map via the PnetCDF writer ---
    if( wrank == 0 )
    {
        Mesh* srcMesh = new Mesh();
        Mesh* tgtMesh = new Mesh();
        CHECK_EQUAL( 0, GenerateCSMesh( *srcMesh, 5, "", "" ) );
        CHECK_EQUAL( 0, GenerateICOMesh( *tgtMesh, 5, true, "", "" ) );

        moab::Interface* mb = new( std::nothrow ) moab::Core;
        CHECK( NULL != mb );
#ifdef MOAB_HAVE_MPI
        moab::ParallelComm* pcomm       = new moab::ParallelComm( mb, MPI_COMM_SELF );
        moab::TempestRemapper* remapper = new moab::TempestRemapper( mb, pcomm );
#else
        moab::TempestRemapper* remapper = new moab::TempestRemapper( mb );
#endif
        remapper->meshValidate     = true;
        remapper->constructEdgeMap = false;
        rval                       = remapper->initialize();CHECK_ERR( rval );
        remapper->SetMesh( moab::Remapper::SourceMesh, srcMesh );
        rval = remapper->ConvertTempestMesh( moab::Remapper::SourceMesh );CHECK_ERR( rval );
        remapper->SetMesh( moab::Remapper::TargetMesh, tgtMesh );
        rval = remapper->ConvertTempestMesh( moab::Remapper::TargetMesh );CHECK_ERR( rval );
        rval = remapper->ConstructCoveringSet( 1e-8, 1.0, 1.0, 0.1, false, true, 0 );CHECK_ERR( rval );
        rval = remapper->ComputeOverlapMesh( true, false );CHECK_ERR( rval );

        moab::TempestOnlineMap* onlinemap = new moab::TempestOnlineMap( remapper );
        GenerateOfflineMapAlgorithmOptions mapOptions;
        mapOptions.nPin     = 1;
        mapOptions.nPout    = 1;
        mapOptions.fNoCheck = true;
        rval                = onlinemap->GenerateRemappingWeights( "fv", "fv", mapOptions );CHECK_ERR( rval );
        std::map< std::string, std::string > attrMap;
        attrMap["Title"] = "MOAB PnetCDF parallel-read test map";
        rval             = onlinemap->WriteParallelMap( mapfile, attrMap );CHECK_ERR( rval );

        delete onlinemap;
        delete remapper;
#ifdef MOAB_HAVE_MPI
        delete pcomm;
#endif
        delete mb;
    }

#ifdef MOAB_HAVE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    // --- All ranks: read the shared map back in parallel ---
    moab::Interface* mb = new( std::nothrow ) moab::Core;
    CHECK( NULL != mb );
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm       = new moab::ParallelComm( mb, MPI_COMM_WORLD );
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mb, pcomm );
#else
    moab::TempestRemapper* remapper = new moab::TempestRemapper( mb );
#endif
    rval = remapper->initialize( false );CHECK_ERR( rval );
    moab::TempestOnlineMap* onlinemap = new moab::TempestOnlineMap( remapper );
    std::vector< int > tgt_owned_ids;
    std::vector< double > areaA, areaB;
    int nA = 0, nB = 0;
    rval = onlinemap->ReadParallelMap( mapfile.c_str(), tgt_owned_ids, 0, areaA, nA, areaB, nB );CHECK_ERR( rval );
    std::cout << "[parallel-read rank " << wrank << "/" << wsize << "] read map: nA=" << nA << " nB=" << nB << "\n";
    CHECK_EQUAL( 150, nA );
    CHECK_EQUAL( 252, nB );

    delete onlinemap;
    delete remapper;
#ifdef MOAB_HAVE_MPI
    delete pcomm;
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    delete mb;
    if( wrank == 0 ) remove( mapfile.c_str() );
}
#endif  // MOAB_HAVE_PNETCDF
