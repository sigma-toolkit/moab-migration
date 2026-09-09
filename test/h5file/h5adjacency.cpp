/**\file h5adjacency.cpp
 *
 * Regression test for ReadHDF5::read_adjacencies() chunked reads.
 *
 * The adjacency table is read in chunks of the internal data buffer.  When a
 * record straddles a chunk boundary the trailing partial record is carried
 * over to the next iteration.  That path used to be broken in three ways:
 *
 *  1. \c count was computed as \c min(chunk_size,remaining) and then had
 *     \c left_over subtracted from it, which underflows \c size_t when
 *     \c left_over is larger.  mhdf takes a signed count, so the wrapped value
 *     arrived negative ("Invalid input for read: offset = 0, count = -6947677").
 *  2. The file \c offset was never advanced, so every chunk re-read from 0.
 *  3. The memmove of the carried-over record passed a handle *count* where a
 *     *byte* count was required, moving only 1/8 of the data and silently
 *     corrupting the record spanning the boundary.
 *
 * Reproducing this with the default 128 MB buffer would need an adjacency
 * table of more than 16.7M entries.  Instead we use the (test-only)
 * BUFFER_SIZE read option to shrink the buffer so a small mesh still produces
 * many chunks.  Record lengths are deliberately varied so that records
 * straddle the chunk boundaries at differing offsets.
 *
 * Note that only *stored* adjacencies are written to the file.  Adjacencies
 * MOAB can recompute on demand are not serialized, so the test must set them
 * explicitly via Interface::add_adjacencies() -- the same store that
 * "mbpart -j" fills, which is how the file that triggered the original
 * failure was produced.
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "TestUtil.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <algorithm>
#include <ctime>
#include <sstream>
#include <vector>

using namespace moab;

static const char filename[] = "h5adjacency_tmp.h5m";

/* A quad strip: (NQUAD+1) x 2 vertices, NQUAD quads, plus a pool of edges to
 * point at.  Sized so the adjacency table spans many chunks at the reduced
 * buffer sizes used below. */
static const int NQUAD = 400;

void test_adjacency_chunked_read();

int main( int argc, char* argv[] )
{
#ifdef MOAB_HAVE_MPI
    int fail = MPI_Init( &argc, &argv );
    if( fail ) return fail;
#else
    argv[0] = argv[argc - argc];  // silence unused-parameter warning in serial
#endif

    int exitval = RUN_TEST( test_adjacency_chunked_read );

#ifdef MOAB_HAVE_MPI
    fail = MPI_Finalize();
    if( fail ) return fail;
#endif

    return exitval;
}

void test_adjacency_chunked_read()
{
    ErrorCode rval;

    // ---- build ----------------------------------------------------------
    Core mbcore;
    Interface& mb = mbcore;

    std::vector< EntityHandle > verts( 2 * ( NQUAD + 1 ) );
    for( int i = 0; i <= NQUAD; ++i )
    {
        double c0[3] = { static_cast< double >( i ), 0.0, 0.0 };
        double c1[3] = { static_cast< double >( i ), 1.0, 0.0 };
        rval         = mb.create_vertex( c0, verts[2 * i] );CHECK_ERR( rval );
        rval = mb.create_vertex( c1, verts[2 * i + 1] );CHECK_ERR( rval );
    }

    // Edge pool for the adjacency lists to reference.
    std::vector< EntityHandle > edges( NQUAD );
    for( int i = 0; i < NQUAD; ++i )
    {
        EntityHandle conn[2] = { verts[2 * i], verts[2 * i + 2] };
        rval                 = mb.create_element( MBEDGE, conn, 2, edges[i] );CHECK_ERR( rval );
    }

    Tag gid = mb.globalId_tag();

    // Quads, each with an explicitly stored adjacency list whose length varies
    // from 1 to 7.  Varying the length keeps record boundaries from lining up
    // with the power-of-two chunk sizes, so records straddle chunk ends.
    std::vector< EntityHandle > quads( NQUAD );
    std::vector< int > expected_counts( NQUAD );  // lengths written to the file
    for( int q = 0; q < NQUAD; ++q )
    {
        EntityHandle conn[4] = { verts[2 * q], verts[2 * q + 2], verts[2 * q + 3], verts[2 * q + 1] };
        rval                 = mb.create_element( MBQUAD, conn, 4, quads[q] );CHECK_ERR( rval );

        int qid = 1000000 + q;
        rval    = mb.tag_set_data( gid, &quads[q], 1, &qid );CHECK_ERR( rval );

        const int nadj = 1 + ( q % 7 );
        std::vector< EntityHandle > adj;
        for( int a = 0; a < nadj; ++a )
        {
            EntityHandle e = edges[( q + a * 13 ) % NQUAD];
            if( std::find( adj.begin(), adj.end(), e ) == adj.end() ) adj.push_back( e );
        }
        rval = mb.add_adjacencies( quads[q], &adj[0], adj.size(), false );CHECK_ERR( rval );
        expected_counts[q] = (int)adj.size();
    }

    rval = mb.write_file( filename, "MOAB" );CHECK_ERR( rval );

    // ---- read back at several buffer sizes ------------------------------
    //
    // The invariant under test: the mesh that comes back must not depend on
    // the read buffer size.  The default buffer (128 MB) reads this table in a
    // single chunk, so the carry-over path is never taken -- that is the
    // reference.  The small sizes below force many chunks with records
    // straddling the boundaries, which is the path that was broken.
    //
    // Comparing against the reference read (rather than against what was
    // stored) keeps the test focused on the chunking bug and independent of
    // how get_adjacencies() blends stored and derived adjacencies.
    //
    // 1024 bytes = 128 handles per chunk, far smaller than the 2397-entry
    // table.  The other sizes place chunk boundaries at different offsets
    // within the records.
    const int buffer_sizes[] = { 1024, 2048, 4099, 8192 };
    const int num_sizes      = (int)( sizeof( buffer_sizes ) / sizeof( buffer_sizes[0] ) );

    // Reference: default buffer size, single chunk.
    std::vector< int > reference_counts;
    {
        Core refcore;
        Interface& ref = refcore;
        rval           = ref.load_file( filename );CHECK_ERR( rval );

        Range rquads;
        rval = ref.get_entities_by_type( 0, MBQUAD, rquads );CHECK_ERR( rval );
        CHECK_EQUAL( (size_t)NQUAD, rquads.size() );

        Tag rgid = ref.globalId_tag();
        reference_counts.resize( NQUAD, -1 );
        for( Range::iterator it = rquads.begin(); it != rquads.end(); ++it )
        {
            int qid = 0;
            rval    = ref.tag_get_data( rgid, &( *it ), 1, &qid );CHECK_ERR( rval );
            const int q = qid - 1000000;
            CHECK( q >= 0 && q < NQUAD );

            std::vector< EntityHandle > adj;
            rval = ref.get_adjacencies( &( *it ), 1, 1, false, adj );CHECK_ERR( rval );
            reference_counts[q] = (int)adj.size();
        }

        for( int q = 0; q < NQUAD; ++q )
            CHECK( reference_counts[q] >= 0 );
    }

    for( int s = 0; s < num_sizes; ++s )
    {
        Core readcore;
        Interface& rmb = readcore;

        std::ostringstream opts;
        opts << "BUFFER_SIZE=" << buffer_sizes[s];

        // Pre-fix this call does not return.  With the file offset never
        // advancing, the loop reaches left_over == count == chunk_size, so
        // count becomes 0, `remaining` stops decreasing and it spins forever.
        // (Observed directly: `remaining` frozen at 2140 after 4 iterations.)
        // On a larger table the same arithmetic instead underflows and mhdf
        // rejects the negative count, which is how this first surfaced.
        // Either way the CTest timeout below turns it into a failure.
        const std::clock_t t0 = std::clock();
        rval                  = rmb.load_file( filename, 0, opts.str().c_str() );CHECK_ERR( rval );
        // Guard in case this is ever run outside CTest: a correct read of this
        // tiny file is essentially instant, so anything near a wall second
        // means the chunk loop is not making progress.
        CHECK( ( std::clock() - t0 ) < 30 * CLOCKS_PER_SEC );

        Range rquads;
        rval = rmb.get_entities_by_type( 0, MBQUAD, rquads );CHECK_ERR( rval );
        CHECK_EQUAL( (size_t)NQUAD, rquads.size() );

        Tag rgid = rmb.globalId_tag();

        // Every quad must match the reference exactly.  A truncated memmove
        // shows up here as a wrong count on the record that straddled a
        // chunk boundary.
        for( Range::iterator it = rquads.begin(); it != rquads.end(); ++it )
        {
            int qid = 0;
            rval    = rmb.tag_get_data( rgid, &( *it ), 1, &qid );CHECK_ERR( rval );
            const int q = qid - 1000000;
            CHECK( q >= 0 && q < NQUAD );

            std::vector< EntityHandle > adj;
            rval = rmb.get_adjacencies( &( *it ), 1, 1, false, adj );CHECK_ERR( rval );
            CHECK_EQUAL( reference_counts[q], (int)adj.size() );
        }
    }

    remove( filename );
}
