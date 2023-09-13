#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ProgOptions.hpp"
#include "MBParallelConventions.h"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#include <sstream>

using namespace moab;

std::string example = TestDir + "unittest/io/ne4np4-esmf.nc";

std::string read_options;
const double eps = 1e-20;


void read_mesh_parallel( bool rcbzoltan, bool no_mixed_elements )
{
    Core moab;
    Interface& mb = moab;

    read_options = "PARALLEL=READ_PART;PARTITION_METHOD=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";
    if( rcbzoltan )
        read_options = "PARALLEL=READ_PART;PARTITION_METHOD=RCBZOLTAN;PARALLEL_RESOLVE_SHARED_ENTS;VARIABLE=";

    if( no_mixed_elements ) read_options += ";NO_MIXED_ELEMENTS";

    ErrorCode rval = mb.load_file( example.c_str(), NULL, read_options.c_str() );CHECK_ERR( rval );

    ParallelComm* pcomm = ParallelComm::get_pcomm( &mb, 0 );
    int procs           = pcomm->proc_config().proc_size();
    int rank            = pcomm->proc_config().proc_rank();

    rval = pcomm->check_all_shared_handles();CHECK_ERR( rval );

    // Get local vertices
    Range local_verts;
    rval = mb.get_entities_by_type( 0, MBVERTEX, local_verts );CHECK_ERR( rval );

    int verts_num = local_verts.size();
  /*  if( 2 == procs )
    {
        if( rcbzoltan )
        {
            if( 0 == rank )
                CHECK_EQUAL( 685, verts_num );
            else if( 1 == rank )
                CHECK_EQUAL( 685, verts_num );  // Not owned vertices included
        }
        else
        {
            if( 0 == rank )
                CHECK_EQUAL( 1120, verts_num );
            else if( 1 == rank )
                CHECK_EQUAL( 1122, verts_num );  // Not owned vertices included
        }
    }

    rval = pcomm->filter_pstatus( local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT );CHECK_ERR( rval );

    verts_num = local_verts.size();
    if( 2 == procs )
    {
        if( rcbzoltan )
        {
            if( 0 == rank )
                CHECK_EQUAL( 685, verts_num );
            else if( 1 == rank )
                CHECK_EQUAL( 595, verts_num );  // Not owned vertices excluded
        }
        else
        {
            if( 0 == rank )
                CHECK_EQUAL( 1120, verts_num );
            else if( 1 == rank )
                CHECK_EQUAL( 160, verts_num );  // Not owned vertices excluded
        }
    }

    // Get local edges
    Range local_edges;
    rval = mb.get_entities_by_type( 0, MBEDGE, local_edges );CHECK_ERR( rval );

    int edges_num = local_edges.size();
    if( 2 == procs )
    {
        if( rcbzoltan )
        {
            if( 0 == rank )
                CHECK_EQUAL( 1005, edges_num );
            else if( 1 == rank )
                CHECK_EQUAL( 1005, edges_num );  // Not owned edges included
        }
        else
        {
            if( 0 == rank )
                CHECK_EQUAL( 1438, edges_num );
            else if( 1 == rank )
                CHECK_EQUAL( 1444, edges_num );  // Not owned edges included
        }
    }

    rval = pcomm->filter_pstatus( local_edges, PSTATUS_NOT_OWNED, PSTATUS_NOT );CHECK_ERR( rval );

    edges_num = local_edges.size();
    if( 2 == procs )
    {
        if( rcbzoltan )
        {
            if( 0 == rank )
                CHECK_EQUAL( 1005, edges_num );
            else if( 1 == rank )
                CHECK_EQUAL( 915, edges_num );  // Not owned edges excluded
        }
        else
        {
            if( 0 == rank )
                CHECK_EQUAL( 1438, edges_num );
            else if( 1 == rank )
                CHECK_EQUAL( 482, edges_num );  // Not owned edges excluded
        }
    }

    // Get local cells
    Range local_cells;
    rval = mb.get_entities_by_type( 0, MBPOLYGON, local_cells );CHECK_ERR( rval );

    int cells_num = local_cells.size();
    if( 2 == procs )
    {
        CHECK_EQUAL( 321, cells_num );
        CHECK_EQUAL( (size_t)1, local_cells.psize() );
    }

    rval = pcomm->filter_pstatus( local_cells, PSTATUS_NOT_OWNED, PSTATUS_NOT );CHECK_ERR( rval );

    cells_num = local_cells.size();
    if( 2 == procs )
    {
        CHECK_EQUAL( 321, cells_num );
        CHECK_EQUAL( (size_t)1, local_cells.psize() );
    }

    std::cout << "proc: " << rank << " verts:" << verts_num << "\n";

    int total_verts_num;
    MPI_Reduce( &verts_num, &total_verts_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm() );
    if( 0 == rank )
    {
        std::cout << "total vertices: " << total_verts_num << "\n";
        CHECK_EQUAL( 1280, total_verts_num );
    }

    std::cout << "proc: " << rank << " edges:" << edges_num << "\n";

    int total_edges_num;
    MPI_Reduce( &edges_num, &total_edges_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm() );
    if( 0 == rank )
    {
        std::cout << "total edges: " << total_edges_num << "\n";
        CHECK_EQUAL( 1920, total_edges_num );
    }

    std::cout << "proc: " << rank << " cells:" << cells_num << "\n";

    int total_cells_num;
    MPI_Reduce( &cells_num, &total_cells_num, 1, MPI_INT, MPI_SUM, 0, pcomm->proc_config().proc_comm() );
    if( 0 == rank )
    {
        std::cout << "total cells: " << total_cells_num << "\n";
        CHECK_EQUAL( 642, total_cells_num );
    }*/

#ifdef MOAB_HAVE_HDF5_PARALLEL
    std::string write_options( "PARALLEL=WRITE_PART;" );

    std::string output_file = "test_esmf";
    if( rcbzoltan ) output_file += "_rcbzoltan";
    if( no_mixed_elements ) output_file += "_no_mixed_elements";
    output_file += ".h5m";

    mb.write_file( output_file.c_str(), NULL, write_options.c_str() );
#endif
}


void test_read_mesh_parallel_trivial()
{
    read_mesh_parallel( false, false );
}

void test_read_mesh_parallel_trivial_no_mixed_elements()
{
    read_mesh_parallel( false, true );
}

void test_read_mesh_parallel_rcbzoltan()
{
    read_mesh_parallel( true, false );
}

void test_read_mesh_parallel_rcbzoltan_no_mixed_elements()
{
    read_mesh_parallel( true, true );
}


int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    int result = 0;

    result += RUN_TEST( test_read_mesh_parallel_trivial );
    result += RUN_TEST( test_read_mesh_parallel_trivial_no_mixed_elements );
#if defined( MOAB_HAVE_MPI ) && defined( MOAB_HAVE_ZOLTAN )
    result += RUN_TEST( test_read_mesh_parallel_rcbzoltan );
    result += RUN_TEST( test_read_mesh_parallel_rcbzoltan_no_mixed_elements );
#endif

    MPI_Finalize();
    return result;
}
