/** @example CoupleMGen.cpp
 * Example demonstrating parallel mesh coupling and interpolation
 *
 * \details This example demonstrates:
 * - Parallel mesh generation with different partitioning schemes
 * - Online coupling between two meshes covering the same domain
 * - Field interpolation from source to target mesh
 * - Performance timing of mesh generation and coupling operations
 * - L-infinity norm calculation for verification
 *
 * The program creates two different meshes with configurable partitioning
 * and demonstrates coupling between them using the MBCoupler interface.
 *
 * \note Requires MOAB to be built with MBCoupler support
 *
 * \par Algorithm:
 * 1. Initialize MPI and create two different meshes
 * 2. Set up interpolation field on source mesh
 * 3. Create coupler instance
 * 4. Interpolate field from source to target mesh
 * 5. Calculate L-infinity norm for verification
 * 6. Report timing information
 *
 * \par Usage:
 * \code
 * mpiexec -np 16 CoupleMGen -K 4 -N 4 -print
 * \endcode
 *
 * \param[in] -blockSize Block size of mesh (default: 4)
 * \param[in] -xproc Number of processors in x direction (default: 1)
 * \param[in] -yproc Number of processors in y direction (default: 1)
 * \param[in] -zproc Number of processors in z direction (default: 1)
 * \param[in] -xblocks Number of blocks on a task in x direction (default: 2)
 * \param[in] -yblocks Number of blocks on a task in y direction (default: 2)
 * \param[in] -zblocks Number of blocks on a task in z direction (default: 2)
 * \param[in] -xsize Total size in x direction (default: 1.0)
 * \param[in] -ysize Total size in y direction (default: 1.0)
 * \param[in] -zsize Total size in z direction (default: 1.0)
 * \param[in] -eps Tolerance for coupling (default: 1e-6)
 * \param[in] -print Write meshes to files
 *
 * \return 0 on success, 1 on failure
 *
 * \par Example:
 * \code
 * mpiexec -np 8 CoupleMGen -M 2 -N 2 -K 2 -blockSize 8
 * \endcode
 *
 * \see Core, ParallelComm, MeshGeneration, Coupler
 * \author MOAB Development Team
 * \date Last Updated 2025
 *
 */
// MOAB includes
#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#ifdef MOAB_HAVE_MBCOUPLER
#include "mbcoupler/Coupler.hpp"
#include "mbcoupler/ElemUtil.hpp"
#else
#error Requires MOAB to be built with MBCoupler
#endif
#include "moab/MeshGeneration.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;
using std::string;

double physField( double x, double y, double z )
{

    double out = sin( M_PI * x ) * cos( M_PI * y ) * sin( M_PI * z );

    return out;
}

int main( int argc, char* argv[] )
{
    int proc_id = 0, size = 1;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    Core mcore;
    Interface* mb = &mcore;
    EntityHandle fileset1, fileset2;  // for 2 different meshes
    MeshGeneration::BrickOpts opts;
    // default options
    opts.A = opts.B = opts.C = 1;
    opts.M = opts.N = opts.K = 1;
    opts.blockSize           = 4;
    opts.xsize = opts.ysize = opts.zsize = 1.;
    opts.ui                              = CartVect( 1., 0, 0. );
    opts.uj                              = CartVect( 0., 1., 0. );
    opts.uk                              = CartVect( 0., 0., 1. );
    opts.newMergeMethod = opts.quadratic = opts.keep_skins = opts.tetra = false;
    opts.adjEnts = opts.parmerge = false;
    opts.GL                      = 0;

    ProgOptions popts;

    popts.addOpt< int >( string( "blockSize,b" ), string( "Block size of mesh (default=4)" ), &opts.blockSize );
    popts.addOpt< int >( string( "xproc,M" ), string( "Number of processors in x dir (default=1)" ), &opts.M );
    popts.addOpt< int >( string( "yproc,N" ), string( "Number of processors in y dir (default=1)" ), &opts.N );
    popts.addOpt< int >( string( "zproc,K" ), string( "Number of processors in z dir (default=1)" ), &opts.K );

    popts.addOpt< int >( string( "xblocks,A" ), string( "Number of blocks on a task in x dir (default=2)" ), &opts.A );
    popts.addOpt< int >( string( "yblocks,B" ), string( "Number of blocks on a task in y dir (default=2)" ), &opts.B );
    popts.addOpt< int >( string( "zblocks,C" ), string( "Number of blocks on a task in x dir (default=2)" ), &opts.C );

    popts.addOpt< double >( string( "xsize,x" ), string( "Total size in x direction (default=1.)" ), &opts.xsize );
    popts.addOpt< double >( string( "ysize,y" ), string( "Total size in y direction (default=1.)" ), &opts.ysize );
    popts.addOpt< double >( string( "zsize,z" ), string( "Total size in z direction (default=1.)" ), &opts.zsize );

    popts.addOpt< void >( "newMerge,w", "use new merging method", &opts.newMergeMethod );

    popts.addOpt< void >( "quadratic,q", "use hex 27 elements", &opts.quadratic );

    popts.addOpt< void >( "keep_skins,k", "keep skins with shared entities", &opts.keep_skins );

    popts.addOpt< void >( "tetrahedrons,t", "generate tetrahedrons", &opts.tetra );

    popts.addOpt< void >( "faces_edges,f", "create all faces and edges", &opts.adjEnts );

    popts.addOpt< int >( string( "ghost_layers,g" ), string( "Number of ghost layers (default=0)" ), &opts.GL );

    popts.addOpt< void >( "parallel_merge,p", "use parallel mesh merge, not vertex ID based merge", &opts.parmerge );

    Coupler::Method method = Coupler::LINEAR_FE;

    double toler = 1.e-6;
    popts.addOpt< double >( string( "eps,e" ), string( "tolerance for coupling, used in locating points" ), &toler );

    bool writeMeshes = false;
    popts.addOpt< void >( "print,p", "write meshes", &writeMeshes );

    popts.parseCommandLine( argc, argv );

    double start_time = MPI_Wtime();

    MB_CHK_ERR( mb->create_meshset( MESHSET_SET, fileset1 ) );
    MB_CHK_ERR( mb->create_meshset( MESHSET_SET, fileset2 ) );

    ParallelComm* pc1     = new ParallelComm( mb, MPI_COMM_WORLD );
    MeshGeneration* mgen1 = new MeshGeneration( mb, pc1, fileset1 );

    MB_CHK_ERR( mgen1->BrickInstance( opts ) );  // this will generate first mesh on fileset1

    double instance_time = MPI_Wtime();
    double current       = instance_time;
    if( !proc_id ) std::cout << " instantiate first mesh " << instance_time - start_time << "\n";
    // set an interpolation tag on source mesh, from phys field
    std::string interpTag( "interp_tag" );
    Tag tag;
    MB_CHK_ERR( mb->tag_get_handle( interpTag.c_str(), 1, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT | MB_TAG_DENSE ) );

    Range src_elems;
    MB_CHK_ERR( pc1->get_part_entities( src_elems, 3 ) );
    Range src_verts;
    MB_CHK_ERR( mb->get_connectivity( src_elems, src_verts ) );
    for( Range::iterator vit = src_verts.begin(); vit != src_verts.end(); ++vit )
    {
        EntityHandle vert = *vit;  //?

        double vertPos[3];
        mb->get_coords( &vert, 1, vertPos );

        double fieldValue = physField( vertPos[0], vertPos[1], vertPos[2] );

        MB_CHK_ERR( mb->tag_set_data( tag, &vert, 1, &fieldValue ) );
    }

    double setTag_time = MPI_Wtime();
    if( !proc_id ) std::cout << " set tag " << setTag_time - current;
    current = instance_time;
    // change some options, so it is a different mesh
    int tmp1   = opts.K;
    opts.K     = opts.M;
    opts.M     = tmp1;  // swap (opts.K, opts.M)
    opts.tetra = !opts.tetra;
    opts.blockSize++;

    ParallelComm* pc2     = new ParallelComm( mb, MPI_COMM_WORLD );
    MeshGeneration* mgen2 = new MeshGeneration( mb, pc2, fileset2 );

    MB_CHK_ERR( mgen2->BrickInstance( opts ) );  // this will generate second mesh on fileset2

    double instance_second = MPI_Wtime();
    if( !proc_id ) std::cout << " instance second mesh" << instance_second - current << "\n";
    current = instance_second;

    // test the sets are fine
    if( writeMeshes )
    {
        MB_CHK_SET_ERR( mb->write_file( "mesh1.h5m", 0, ";;PARALLEL=WRITE_PART;CPUTIME;PARALLEL_COMM=0;", &fileset1,
                                        1 ),
                        "Can't write in parallel mesh 1" );
        MB_CHK_SET_ERR( mb->write_file( "mesh2.h5m", 0, ";;PARALLEL=WRITE_PART;CPUTIME;PARALLEL_COMM=1;", &fileset2,
                                        1 ),
                        "Can't write in parallel mesh 1" );
        double write_files = MPI_Wtime();
        if( !proc_id ) std::cout << " write files " << write_files - current << "\n";
        current = write_files;
    }

    // Instantiate a coupler, which also initializes the tree
    Coupler mbc( mb, pc1, src_elems, 0 );

    double instancecoupler = MPI_Wtime();
    if( !proc_id ) std::cout << " instance coupler " << instancecoupler - current << "\n";
    current = instancecoupler;

    // Get points from the target mesh to interpolate
    // We have to treat differently the case when the target is a spectral mesh
    // In that case, the points of interest are the GL points, not the vertex nodes
    std::vector< double > vpos;  // This will have the positions we are interested in
    int numPointsOfInterest = 0;

    Range targ_elems;
    Range targ_verts;

    // First get all vertices adj to partition entities in target mesh
    MB_CHK_ERR( pc2->get_part_entities( targ_elems, 3 ) );

    MB_CHK_ERR( mb->get_adjacencies( targ_elems, 0, false, targ_verts, Interface::UNION ) );
    Range tmp_verts;
    // Then get non-owned verts and subtract
    MB_CHK_ERR( pc2->get_pstatus_entities( 0, PSTATUS_NOT_OWNED, tmp_verts ) );
    targ_verts = subtract( targ_verts, tmp_verts );
    // get position of these entities; these are the target points
    numPointsOfInterest = (int)targ_verts.size();
    vpos.resize( 3 * targ_verts.size() );
    MB_CHK_ERR( mb->get_coords( targ_verts, &vpos[0] ) );
    // Locate those points in the source mesh
    // std::cout<<"rank "<< proc_id<< " points of interest: " << numPointsOfInterest << "\n";
    MB_CHK_ERR( mbc.locate_points( &vpos[0], numPointsOfInterest, 0, toler ) );

    double locatetime = MPI_Wtime();
    if( !proc_id ) std::cout << " locate points: " << locatetime - current << "\n";
    current = locatetime;

    // Now interpolate tag onto target points
    std::vector< double > field( numPointsOfInterest );

    MB_CHK_ERR( mbc.interpolate( method, interpTag, &field[0] ) );

    // compare with the actual phys field
    double err_max = 0;
    for( int i = 0; i < numPointsOfInterest; i++ )
    {
        double trval = physField( vpos[3 * i], vpos[3 * i + 1], vpos[3 * i + 2] );
        double err2  = fabs( trval - field[i] );
        if( err2 > err_max ) err_max = err2;
    }

    double interpolateTime = MPI_Wtime();
    if( !proc_id ) std::cout << " interpolate points: " << interpolateTime - current << "\n";
    current = interpolateTime;

    double gerr;
    MPI_Allreduce( &err_max, &gerr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    if( 0 == proc_id ) std::cout << "max err  " << gerr << "\n";

    delete mgen1;
    delete mgen2;

    MPI_Finalize();

    return 0;
}
