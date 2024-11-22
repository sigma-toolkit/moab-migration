#include <cmath>
#include <iostream>
#include <numeric>

#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/MergeMesh.hpp"

using namespace moab;

#define MBERRORR( rval, STR )                  \
    {                                          \
        if( MB_SUCCESS != ( rval ) )           \
        {                                      \
            std::cout << ( STR ) << std::endl; \
            std::cout.flush();                 \
            std::exit( 1 );                    \
        }                                      \
    }

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );

    // let us start the computation
    {
        moab::Core mb;
        // moab::Interface& mb = moab;

        MPI_Comm global_comm = MPI_COMM_WORLD;
        moab::ParallelComm pcomm( &mb, global_comm );

        ErrorCode rval;
        std::string output_file;

        auto rank = pcomm.rank();
        auto size = pcomm.size();

        std::string file_name = std::string( MESH_DIR ) + std::string( "/64bricks_1khex.h5m" ) ;
        std::string options   = "PARALLEL=READ_PART;"
                                "PARALLEL_RESOLVE_SHARED_ENTS;"
                                "PARTITION=PARALLEL_PARTITION;";

        // let us load the user file
        rval = mb.load_file( file_name.c_str(), 0, options.c_str() );
        MBERRORR( rval, "load file" )

        mb.write_file( "current_mesh.h5m", 0, "PARALLEL=WRITE_PART" );
        output_file = "original_local_mesh_" + std::to_string( rank ) + ".h5m";
        mb.write_file( output_file.c_str() );

        std::vector< int > from_procs( size - 1 );  // size-1 sending processes
        std::iota( from_procs.begin(), from_procs.end(), 1 );
        std::vector< int > to_procs( 1, 0 );  // just 1-receiving process, the global root

        auto send_to_root = [&]( moab::Range& ents ) -> moab::Range {

            int sendrecv_tag = 100;

            pcomm.set_debug_verbosity( 5 );

            // two requests per communication typically: 1 - message size/acknowledge, 2 - buffer
            MPI_Request recv_remoteh_reqs = MPI_REQUEST_NULL;
            MPI_Status status;

            std::vector< int > buffer_sizes( size, 0 );
            int local_send_buffer_send = 0;

            // if the current rank is in the sending end...
            ParallelComm::Buffer* send_buffer = nullptr;
            if( std::find( from_procs.begin(), from_procs.end(), rank ) != from_procs.end() )
            {
                send_buffer = new ParallelComm::Buffer( ParallelComm::INITIAL_BUFF_SIZE );
                send_buffer->reset_ptr( sizeof( int ) );
                MBERRORR( pcomm.pack_buffer( ents, false, false, false, to_procs[0], send_buffer ),
                          " can't pack buffer for entities to send" );
                local_send_buffer_send = send_buffer->get_current_size();

                MBERRORR( MPI_Isend( send_buffer->mem_ptr, local_send_buffer_send, MPI_UNSIGNED_CHAR, to_procs[0],
                                     sendrecv_tag, global_comm, &recv_remoteh_reqs ),
                          "sending buffer failed" );  // we have to use global communicator
            }

            // if the current rank is in the receiving end...
            if( std::find( to_procs.begin(), to_procs.end(), rank ) != to_procs.end() )
            {
                Range entities;
                MPI_Status recvstatus;
                for( size_t irank = 0; irank < size; ++irank )
                {
                    if( irank == rank ) continue;  // no self messages!

                    MBERRORR( MPI_Probe( irank, sendrecv_tag, global_comm, &recvstatus ),
                              " MPI_Probe failure in MPI_Probe " );

                    // get the count of data received from the MPI_Status structure
                    int local_recv_buffer_send;
                    MBERRORR( MPI_Get_count( &recvstatus, MPI_CHAR, &local_recv_buffer_send ),
                              " MPI_Get_count failure" );

                    // now resize the buffer, then receive it
                    ParallelComm::Buffer* buffer = new ParallelComm::Buffer( local_recv_buffer_send );

                    MBERRORR( MPI_Recv( buffer->mem_ptr, local_recv_buffer_send, MPI_UNSIGNED_CHAR, irank, sendrecv_tag,
                                        global_comm, &recvstatus ),
                              " MPI_Recv failure" );

                    // now unpack the buffer we just received
                    std::vector< std::vector< EntityHandle > > L1hloc, L1hrem;
                    std::vector< std::vector< int > > L1p;
                    std::vector< EntityHandle > L2hloc, L2hrem;
                    std::vector< unsigned int > L2p;

                    buffer->reset_ptr( sizeof( int ) );
                    std::vector< EntityHandle > entities_vec;
                    MBERRORR( pcomm.unpack_buffer( buffer->buff_ptr, false, -1, -1, L1hloc, L1hrem, L1p, L2hloc, L2hrem,
                                                   L2p, entities_vec ),
                              "unable to unpack buffer" );
                    delete buffer;

                    // accumulate the received entities into the local range of entities
                    std::copy( entities_vec.begin(), entities_vec.end(), range_inserter( entities ) );
                }

                moab::EntityHandle root_set = 0;
                MergeMesh merger( &mb );
                // merge all the new handles with any existing shared entities that are already locally available
                MBERRORR( merger.merge_all( root_set, 1e-8 ), "merging local mesh failed" );
            }

            // wait for all messages to be done.
            MPI_Wait( &recv_remoteh_reqs, &status );

            if( send_buffer ) delete send_buffer;

            Range new_ents;
            MBERRORR( mb.get_entities_by_handle( 0, new_ents, true ), "getting entities by handle failed" );

            return subtract( new_ents, ents );
        };

        Range ents;
        MBERRORR( mb.get_entities_by_handle( 0, ents, true ), "failed get_entities_by_handle " );

        Range new_ents = send_to_root( ents );

        if( std::find( to_procs.begin(), to_procs.end(), rank ) != to_procs.end() )
        {
            std::cout << rank << ": Number of new entities obtained = " << new_ents.size() << "\n";
        }

        // write out the updated mesh to disk for verification
        output_file = "updated_local_mesh_" + std::to_string( rank ) + ".h5m";
        mb.write_file( output_file.c_str() );
    }

    MPI_Finalize();
    return 0;
}
