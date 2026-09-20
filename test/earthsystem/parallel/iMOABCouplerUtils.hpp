/*
 * iMOABCouplerUtils.hpp
 *
 *  Created on: Aug. 22, 2020
 *  \brief Utility methods shared by iMOAB coupler tests to avoid repetitive tasks
 *  \ even migrate tests can use some of these utilities
 *  1) create_comm_group(int start, int end, int tag, MPI_Group& group, MPI_Comm& comm)
 *
 */

#ifndef TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_
#define TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>

#define CHECKIERR( rc, message )                                              \
    if( 0 != ( rc ) )                                                         \
    {                                                                         \
        printf( "Error at line number %d in file %s\n", __LINE__, __FILE__ ); \
        printf( "%s. ErrorCode = %d\n", message, rc );                        \
        return 1;                                                             \
    }

#define PUSH_TIMER( operation )               \
    {                                         \
        timer_ops = timer.time_since_birth(); \
        opName    = operation;                \
    }
#define POP_TIMER( localcomm, localrank )                                                         \
    {                                                                                             \
        double locElapsed = timer.time_since_birth() - timer_ops, minElapsed = 0, maxElapsed = 0; \
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, localcomm );             \
        MPI_Reduce( &locElapsed, &minElapsed, 1, MPI_DOUBLE, MPI_MIN, 0, localcomm );             \
        if( !( localrank ) )                                                                      \
            std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed     \
                      << ", avg = " << ( maxElapsed + minElapsed ) / 2 << "\n";                   \
        opName.clear();                                                                           \
    }

/*
 *  \brief create an MPI group and an MPI communicator for the group, in the global communicator
 */
inline int create_group_and_comm( int start, int end, MPI_Group worldGroup, MPI_Group* group, MPI_Comm* comm )
{
    std::vector< int > groupTasks;
    groupTasks.resize( end - start + 1, 0 );
    for( int i = start; i <= end; i++ )
        groupTasks[i - start] = i;

    int ierr = MPI_Group_incl( worldGroup, end - start + 1, &groupTasks[0], group );
    CHECKIERR( ierr, "Cannot create group" )

    ierr = MPI_Comm_create( MPI_COMM_WORLD, *group, comm );
    CHECKIERR( ierr, "Cannot create comm" )

    return 0;
}

inline int create_joint_comm_group( MPI_Group agroup, MPI_Group bgroup, MPI_Group* abgroup, MPI_Comm* abcomm )
{
    int ierr = MPI_Group_union( agroup, bgroup, abgroup );
    CHECKIERR( ierr, "Cannot create joint union group" )

    ierr = MPI_Comm_create( MPI_COMM_WORLD, *abgroup, abcomm );
    CHECKIERR( ierr, "Cannot create joint communicator from union group" )

    return 0;
}

inline int setup_component_coupler_meshes( iMOAB_AppID cmpId,
                                    int cmpTag,
                                    iMOAB_AppID cplCmpId,
                                    int cmpcouTag,
                                    MPI_Comm* cmpcomm,
                                    MPI_Group* cmpPEGroup,
                                    MPI_Comm* coucomm,
                                    MPI_Group* cplPEGroup,
                                    MPI_Comm* cmpcoucomm,
                                     const std::string& filename,
                                     const std::string& readopts,
                                    int nghlay,
                                    int repartitioner_scheme )
{
    int ierr = 0;
    if( *cmpcomm != MPI_COMM_NULL )
    {
        // load first mesh
        ierr = iMOAB_LoadMesh( cmpId, filename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load component mesh" )

        // then send mesh to coupler pes
        ierr = iMOAB_SendMesh( cmpId, cmpcoucomm, cplPEGroup, &cmpcouTag,
                               &repartitioner_scheme );  // send to  coupler pes
        CHECKIERR( ierr, "cannot send elements" )
    }
    // now, receive mesh, on coupler communicator; first mesh 1, atm
    if( *coucomm != MPI_COMM_NULL )
    {

        ierr = iMOAB_ReceiveMesh( cplCmpId, cmpcoucomm, cmpPEGroup,
                                  &cmpTag );  // receive from component
        CHECKIERR( ierr, "cannot receive elements on coupler app" )
    }

    // we can now free the sender buffers
    if( *cmpcomm != MPI_COMM_NULL )
    {
        int context_id = cmpcouTag;
        ierr           = iMOAB_FreeSenderBuffers( cmpId, &context_id );
        CHECKIERR( ierr, "cannot free buffers used to send atm mesh" )
    }
    return 0;
}

// Gather (GLOBAL_ID, <tagName>) pairs from every rank in comm to
// rank 0, sort by GID, and write to a digest file.  The sort-order is
// decomposition-independent so the digest is byte-identical iff the
// per-cell projected values are bit-for-bit identical across rank counts.
inline int gather_and_write_proj_tag( MPI_Comm comm,
                               int rankInComm,
                               iMOAB_AppID pid,
                               const std::string& tagName,
                               const std::string& outFilename )
{
    int nverts[3], nelem[3];
    int ierr = iMOAB_GetMeshInfo( pid, nverts, nelem, 0, 0, 0 );
    if( ierr ) return 1;

    int tag_type = DENSE_INTEGER, ncomp = 1, tagInd = 0;
    ierr = iMOAB_DefineTagStorage( pid, "GLOBAL_ID", &tag_type, &ncomp, &tagInd );
    if( ierr ) return 1;

    int ent_type = 1;  // elements
    int sz       = nelem[2];
    std::vector< int > gids( sz, 0 );
    std::vector< double > vals( sz, 0.0 );
    ierr = iMOAB_GetIntTagStorage( pid, "GLOBAL_ID", &sz, &ent_type, gids.data() );
    if( ierr ) return 1;
    ierr = iMOAB_GetDoubleTagStorage( pid, tagName.c_str(), &sz, &ent_type, vals.data() );
    if( ierr ) return 1;

    int sizeInComm = 0;
    MPI_Comm_size( comm, &sizeInComm );

    std::vector< int > counts( sizeInComm, 0 );
    MPI_Gather( &sz, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm );

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

    MPI_Gatherv( gids.data(), sz, MPI_INT, rankInComm == 0 ? allGids.data() : nullptr, counts.data(), displs.data(),
                 MPI_INT, 0, comm );
    MPI_Gatherv( vals.data(), sz, MPI_DOUBLE, rankInComm == 0 ? allVals.data() : nullptr, counts.data(), displs.data(),
                 MPI_DOUBLE, 0, comm );

    if( rankInComm != 0 ) return 0;

    // Pair, sort by GID (ascending), dedup
    std::vector< std::pair< int, double > > pairs;
    pairs.reserve( totalCount );
    for( int i = 0; i < totalCount; ++i )
        pairs.emplace_back( allGids[i], allVals[i] );

    std::sort( pairs.begin(), pairs.end(), []( const std::pair< int, double >& a, const std::pair< int, double >& b ) {
        return a.first < b.first;
    } );

    auto last = std::unique( pairs.begin(), pairs.end(),
                             []( const std::pair< int, double >& a, const std::pair< int, double >& b ) {
                                 return a.first == b.first;
                             } );
    pairs.erase( last, pairs.end() );

    std::ofstream fs( outFilename );
    if( !fs.is_open() ) return 1;
    fs << std::fixed << std::setprecision( 16 );
    for( auto& p : pairs )
        fs << p.first << " " << p.second << "\n";
    return 0;
}

#endif /* TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_ */
