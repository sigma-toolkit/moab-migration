#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/ParallelMergeMesh.hpp"
#include "TestUtil.hpp"
#include <iostream>

using namespace moab;

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    int nproc, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    std::string filename0 = TestDir + "unittest/brick1.vtk";
    std::string filename1 = TestDir + "unittest/brick2.vtk";

    moab::Core* mb         = new moab::Core();
    moab::ParallelComm* pc = new moab::ParallelComm( mb, MPI_COMM_WORLD );
    ErrorCode rval         = MB_SUCCESS;

    if( 0 == rank%2 )
        rval = mb->load_file( filename0.c_str() );
    else
        rval = mb->load_file( filename1.c_str() );

    if( rval != MB_SUCCESS )
    {
        std::cout << "fail to load file\n";
        delete pc;
        delete mb;
        MPI_Finalize();
        return 1;
    }
    // if rank > 2, translate in z direction, with the distance ( (rank)/2 ) * 10
    if (rank >= 2)
    {
        Range verts;
        rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_ERR( rval );
        int num_verts=(int)verts.size();
        std::vector<double> coords;
        coords.resize(num_verts*3);
        rval = mb->get_coords(verts, &coords[0]);MB_CHK_ERR( rval );
        int steps=rank/2;
        double z_translate = steps*10.;// the 2 bricks are size 10
        for (int i=0; i<num_verts; i++)
            coords[3*i+2] += z_translate;
        rval = mb->set_coords(verts, &coords[0]);MB_CHK_ERR( rval );
    }


    ParallelMergeMesh pm( pc, 0.001 );
    rval = pm.merge();
    if( rval != MB_SUCCESS )
    {
        std::cout << "fail to merge in parallel \n";
        delete pc;
        delete mb;
        MPI_Finalize();
        return 1;
    }
    if (nproc == 2)
    {
        // check number of shared entities
        Range shared_ents;
        // Get entities shared with all other processors
        rval = pc->get_shared_entities( -1, shared_ents );
        if( rval != MB_SUCCESS )
        {
            delete pc;
            delete mb;
            MPI_Finalize();
            return 1;
        }
        // there should be exactly 9 vertices, 12 edges, 4 faces
        unsigned numV = shared_ents.num_of_dimension( 0 );
        unsigned numE = shared_ents.num_of_dimension( 1 );
        unsigned numF = shared_ents.num_of_dimension( 2 );

        if( numV != 9 || numE != 12 || numF != 4 )
        {
            std::cout << " wrong number of shared entities on proc " << rank << " v:" << numV << " e:" << numE
                      << " f:" << numF << "\n";
            delete pc;
            delete mb;
            MPI_Finalize();
            return 1;
        }
    }
    rval = mb->write_file( "testpm.h5m", 0, "PARALLEL=WRITE_PART" );
    if( rval != MB_SUCCESS )
    {
        std::cout << "fail to write output file \n";
        delete pc;
        delete mb;
        MPI_Finalize();
        return 1;
    }

    delete pc;
    delete mb;

    MPI_Finalize();
    return 0;
}
