/**
 * @file generate_meshes.cpp
 * @brief Generate MOAB-format mesh files for mbtempest regression tests.
 *
 * Usage:
 *   generate_meshes <type> <resolution> <output_file>
 *
 * Types: cs, rll, ico, icod (icosahedral dual)
 *
 * Example:
 *   generate_meshes cs 5 cs5.h5m
 *   generate_meshes rll 10 rll10.h5m
 *   generate_meshes ico 5 ico5.h5m
 *   generate_meshes icod 5 icod5.h5m
 */

#include "moab/Core.hpp"
#include "moab/earthsystem/remapping/TempestRemapper.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

#include <iostream>
#include <string>

int main( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    if( argc < 4 )
    {
        std::cerr << "Usage: " << argv[0] << " <type> <resolution> <output_file>\n";
        std::cerr << "  type: cs, rll, ico, icod\n";
        std::cerr << "  resolution: integer (e.g. 5, 10, 25)\n";
        std::cerr << "  output_file: path to write (e.g. cs5.h5m)\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    std::string meshType  = argv[1];
    int resolution        = std::atoi( argv[2] );
    std::string outputFile = argv[3];

    moab::Core mbCore;
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm pcomm( &mbCore, MPI_COMM_WORLD );
#endif

    moab::TempestRemapper remapper( &mbCore,
#ifdef MOAB_HAVE_MPI
                                    &pcomm
#else
                                    nullptr
#endif
    );

    moab::ErrorCode rval;

    // Determine mesh type
    moab::TempestRemapper::TempestMeshType trType;
    if( meshType == "cs" )
        trType = moab::TempestRemapper::CS;
    else if( meshType == "rll" )
        trType = moab::TempestRemapper::RLL;
    else if( meshType == "ico" )
        trType = moab::TempestRemapper::ICO;
    else if( meshType == "icod" )
        trType = moab::TempestRemapper::ICOD;
    else
    {
        std::cerr << "Unknown mesh type: " << meshType << "\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    // Initialize remapper
    rval = remapper.Initialize();
    if( rval != moab::MB_SUCCESS )
    {
        std::cerr << "Failed to initialize TempestRemapper\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    // Generate the mesh using the source mesh slot
    rval = remapper.GenerateMesh( moab::Remapper::SourceMesh, trType, resolution, outputFile );
    if( rval != moab::MB_SUCCESS )
    {
        std::cerr << "Failed to generate " << meshType << " mesh at resolution " << resolution << "\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    // Write the MOAB mesh to h5m file
    moab::EntityHandle meshSet = remapper.GetMeshSet( moab::Remapper::SourceMesh );
    rval                       = mbCore.write_file( outputFile.c_str(), nullptr, nullptr, &meshSet, 1 );
    if( rval != moab::MB_SUCCESS )
    {
        std::cerr << "Failed to write mesh to " << outputFile << "\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    std::cout << "Generated " << meshType << " mesh (res " << resolution << ") -> " << outputFile << "\n";

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
