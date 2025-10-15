/*
 * visuMap.cpp
 * this tool will take a source file, target file (h5m) and a map file in nc format, and will visualize weights
 *
 * example of usage:
 * ./mbvisumap -s source.h5m -t target.h5m -m map.nc -b startSourceID \
 *          -e endSourceID  -c startTargetID -f endTargetID -o 1
 *  will associate row i, corresponding to target DOF i, in the map, with a partial mesh with entities from source mesh that
 *      target the DOF i; i.e. the weights w(i,j)!=0 , j=1,n_b, will be displayed on source cells with global DOF j
 *  will associate column j corresponding to source DOF j, in the map, with a partial mesh with entities from the target mesh
 *      that are affected by the source DOF j; i.e., the weights w(i,j)!=0, i=1,n_a, will be displayed on target cells with
 *      global DOF i
 *
 *      The option -o controls if the row and columns files are output in vtk or in h5m format
 *
 * can be built only if netcdf and hdf5 and eigen3 are available
 *
 * default option is now -o 2, which will create an h5m edge mesh file, with the each edge corresponding to  w(i,j)!=0 in the
 *  map file, connecting source center i with target center j. The map will be displayed on a sphere of radius 1,
 *  with the source centers highly elevated from the surfaces, to differentiate them from the target vertices, which
 *  stay on the sphere of radius 1; the elevation is controlled by a new option, -r, with a default value of .05
 *  which means that the source vertices will be put on a sphere of radius 1.05, creating an umbrella for each source center
 *
 *  only the map file is needed, positions for source and target centers are taken from the map file itself
 *  example of usage:
 *   ./mbvisumap  -m map.nc  -r 0.01
 *
 */
#include "moab/MOABConfig.h"

#ifndef MOAB_HAVE_EIGEN3
#error mbvisumap tool requires eigen3 configuration
#endif

#ifndef MOAB_HAVE_HDF5
#error mbvisumap tool requires hdf5 configuration
#endif

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/ReadUtilIface.hpp"

#include "netcdf.h"
#include <cmath>
#include <sstream>
#include <map>
#include <Eigen/Sparse>

#define ERR_NC( e )                                \
    {                                              \
        printf( "Error: %s\n", nc_strerror( e ) ); \
        exit( 2 );                                 \
    }

// copy from ReadNCDF.cpp some useful macros for reading from a netcdf file
// ncFile1 is an integer initialized when opening the nc file in read mode

int ncFile1;

#define GET_DIM1( ncdim, name, val )                            \
    {                                                           \
        int gdfail = nc_inq_dimid( ncFile1, name, &( ncdim ) ); \
        if( NC_NOERR == gdfail )                                \
        {                                                       \
            size_t tmp_val;                                     \
            gdfail = nc_inq_dimlen( ncFile1, ncdim, &tmp_val ); \
            if( NC_NOERR != gdfail )                            \
            {                                                   \
                ERR_NC( gdfail )                                \
            }                                                   \
            else                                                \
                ( val ) = tmp_val;                              \
        }                                                       \
        else                                                    \
            ( val ) = 0;                                        \
    }

#define GET_VAR1( name, id, dims )                                     \
    {                                                                  \
        ( id )     = -1;                                               \
        int gvfail = nc_inq_varid( ncFile1, name, &( id ) );           \
        if( NC_NOERR == gvfail )                                       \
        {                                                              \
            int ndims;                                                 \
            gvfail = nc_inq_varndims( ncFile1, id, &ndims );           \
            if( NC_NOERR == gvfail )                                   \
            {                                                          \
                ( dims ).resize( ndims );                              \
                gvfail = nc_inq_vardimid( ncFile1, id, &( dims )[0] ); \
                if( NC_NOERR != gvfail )                               \
                {                                                      \
                    ERR_NC( gvfail )                                   \
                }                                                      \
            }                                                          \
        }                                                              \
    }

#define GET_1D_INT_VAR1( name, id, vals )                                               \
    {                                                                                   \
        GET_VAR1( name, id, vals );                                                     \
        if( -1 != ( id ) )                                                              \
        {                                                                               \
            size_t ntmp;                                                                \
            int ivfail = nc_inq_dimlen( ncFile1, ( vals )[0], &ntmp );                  \
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
            ( vals ).resize( ntmp );                                                    \
            size_t ntmp1 = 0;                                                           \
            ivfail       = nc_get_vara_int( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
        }                                                                               \
    }

#define GET_1D_DBL_VAR1( name, id, vals )                                                  \
    {                                                                                      \
        std::vector< int > dum_dims;                                                       \
        GET_VAR1( name, id, dum_dims );                                                    \
        if( -1 != ( id ) )                                                                 \
        {                                                                                  \
            size_t ntmp;                                                                   \
            int dvfail = nc_inq_dimlen( ncFile1, dum_dims[0], &ntmp );                     \
            if( NC_NOERR != dvfail )                                                       \
            {                                                                              \
                ERR_NC( dvfail )                                                           \
            }                                                                              \
            ( vals ).resize( ntmp );                                                       \
            size_t ntmp1 = 0;                                                              \
            dvfail       = nc_get_vara_double( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != dvfail )                                                       \
            {                                                                              \
                ERR_NC( dvfail )                                                           \
            }                                                                              \
        }                                                                                  \
    }

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{

    ProgOptions opts;
    int dimSource   = 2;  // for FV meshes is 2; for SE meshes, use fine mesh, dim will be 0
    int dimTarget   = 2;  //
    int otype       = 2;
    double fraction = 0.05;
    std::string inputfile1, inputSource, inputTarget;
    opts.addOpt< std::string >( "map,m", "input map ", &inputfile1 );
    opts.addOpt< std::string >( "source,s", "source mesh", &inputSource );
    opts.addOpt< std::string >( "target,t", "target mesh", &inputTarget );
    opts.addOpt< int >( "dimSource,d", "dimension of source  ", &dimSource );
    opts.addOpt< int >( "dimTarget,g", "dimension of target  ", &dimTarget );
    opts.addOpt< int >( "typeOutput,o", " output type vtk(0), h5m(1), view(default = 2) ", &otype );

    int startSourceID = -1, endSourceID = -1, startTargetID = -1, endTargetID = -1;
    opts.addOpt< int >( "startSourceID,b", "start source id ", &startSourceID );
    opts.addOpt< int >( "endSourceID,e", "end source id ", &endSourceID );
    opts.addOpt< int >( "startTargetID,c", "start target id ", &startTargetID );
    opts.addOpt< int >( "endTargetID,f", "end target id ", &endTargetID );
    opts.addOpt< double >( "raiseFraction,r", "fraction for raising source points height (default 0.05)", &fraction );
    //  -b startSourceID -e endSourceID  -c startTargetID -f endTargetID

    opts.parseCommandLine( argc, argv );

    std::string extension = ".vtk";
    if( 1 <= otype ) extension = ".h5m";

    // Open netcdf map file
    int fail = nc_open( inputfile1.c_str(), 0, &ncFile1 );
    if( NC_NOWRITE != fail )
    {
        ERR_NC( fail )
    }

    std::cout << " opened " << inputfile1 << " for map 1 \n";

    int temp_dim;
    int na1, nb1, ns1;
    GET_DIM1( temp_dim, "n_a", na1 );
    GET_DIM1( temp_dim, "n_b", nb1 );
    GET_DIM1( temp_dim, "n_s", ns1 );
    std::cout << " n_a, n_b, n_s : " << na1 << ", " << nb1 << ", " << ns1 << " for map 1 \n";
    std::vector< int > col1( ns1 ), row1( ns1 );
    std::vector< double > val1( ns1 );
    int idrow1, idcol1, ids1;
    GET_1D_INT_VAR1( "row", idrow1, row1 );
    GET_1D_INT_VAR1( "col", idcol1, col1 );
    GET_1D_DBL_VAR1( "S", ids1, val1 );

    // we read the matrix; now read moab source and target
    Core core;
    Interface* mb = &core;
    ErrorCode rval;
    Tag gtag = mb->globalId_tag();

    // a dense tag for weights
    Tag wtag;
    double defVal = 0;

    std::string name_map = inputfile1;
    // strip last 3 chars (.nc extension)
    name_map.erase( name_map.begin() + name_map.length() - 3, name_map.end() );
    // if path , remove from name
    size_t pos = name_map.rfind( '/', name_map.length() );
    if( pos != std::string::npos ) name_map = name_map.erase( 0, pos + 1 );

    MB_CHK_SET_ERR( mb->tag_get_handle( "weight", 1, MB_TYPE_DOUBLE, wtag, MB_TAG_CREAT | MB_TAG_DENSE, &defVal ),
                    "Failed to create weight" );

    if( 2 == otype )
    {
        // create a view of the full map, in which each weight is shown on an edge that starts at the source cell center
        // and ends at the target cell center
        // source cell centers are raised a little, let's say a fraction 0.05 * radius, which is 1
        // each edge gets the associated weight as a tag
        // first read the cell centers for source and target meshes, directly from the map file

        std::vector< double > xc_a( na1 ), yc_a( na1 );
        std::vector< double > xc_b( nb1 ), yc_b( nb1 );
        int idxc_a, idxc_b, idyc_a, idyc_b;
        GET_1D_DBL_VAR1( "xc_a", idxc_a, xc_a );
        GET_1D_DBL_VAR1( "xc_b", idxc_b, xc_b );
        GET_1D_DBL_VAR1( "yc_a", idyc_a, yc_a );
        GET_1D_DBL_VAR1( "yc_b", idyc_b, yc_b );
        // create source vertices, and target vertices
        std::vector< double > vertex_coords_src( 3 * na1 );
        // xc_a and yc_a are in degrees, usually
        for( int i = 0; i < na1; i++ )
        {
            IntxUtils::SphereCoords sph;
            sph.R                        = 1 + fraction;  // slightly higher
            sph.lon                      = xc_a[i] * M_PI / 180;
            sph.lat                      = yc_a[i] * M_PI / 180;
            CartVect pos                 = IntxUtils::spherical_to_cart( sph );
            vertex_coords_src[3 * i]     = pos[0];
            vertex_coords_src[3 * i + 1] = pos[1];
            vertex_coords_src[3 * i + 2] = pos[2];
        }
        Range source_verts;
        MB_CHK_SET_ERR( mb->create_vertices( &vertex_coords_src[0], na1, source_verts ),
                        "can't create source vertices" );
        // create a set with source vertices
        EntityHandle srcSet;
        MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, srcSet ), "can't create source set for vertices" );
        MB_CHK_SET_ERR( mb->add_entities( srcSet, source_verts ), "can't add vertices" );
        std::vector< int > vgid( na1 );
        for( int i = 0; i < na1; i++ )
            vgid[i] = i + 1;
        MB_CHK_SET_ERR( mb->tag_set_data( gtag, source_verts, &vgid[0] ), "can't set global id on source verts" );

        std::vector< double > vertex_coords_tgt( 3 * nb1 );
        // xc_a and yc_a are in degrees, usually
        for( int i = 0; i < nb1; i++ )
        {
            IntxUtils::SphereCoords sph;
            sph.R                        = 1;  //
            sph.lon                      = xc_b[i] * M_PI / 180;
            sph.lat                      = yc_b[i] * M_PI / 180;
            CartVect pos                 = IntxUtils::spherical_to_cart( sph );
            vertex_coords_tgt[3 * i]     = pos[0];
            vertex_coords_tgt[3 * i + 1] = pos[1];
            vertex_coords_tgt[3 * i + 2] = pos[2];
        }
        Range target_verts;
        MB_CHK_SET_ERR( mb->create_vertices( &vertex_coords_tgt[0], nb1, target_verts ),
                        "can't create target vertices" );
        EntityHandle tgtSet;
        MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, tgtSet ), "can't create target set for vertices" );
        MB_CHK_SET_ERR( mb->add_entities( tgtSet, target_verts ), "can't add vertices" );
        vgid.resize( nb1 );
        for( int i = 0; i < nb1; i++ )
            vgid[i] = i + 1;
        MB_CHK_SET_ERR( mb->tag_set_data( gtag, target_verts, &vgid[0] ), "can't set global id on target verts" );
        // create ns1 edges

        ReadUtilIface* read_iface;
        MB_CHK_ERR( mb->query_interface( read_iface ) );

        EntityHandle actual_start_handle;
        EntityHandle* array = nullptr;
        MB_CHK_ERR( read_iface->get_element_connect( ns1, 2, MBEDGE, 1, actual_start_handle, array ) );

        for( int i = 0; i < ns1; i++ )
        {
            array[2 * i]     = source_verts[col1[i] - 1];  // 1 based to 0 based index
            array[2 * i + 1] = target_verts[row1[i] - 1];  // 1 based to 0 based index
        }
        Range edges( actual_start_handle, actual_start_handle + ns1 - 1 );

        MB_CHK_SET_ERR( mb->tag_set_data( wtag, edges, &val1[0] ), "can't set tag on edges" );

        vgid.resize( ns1 );
        for( int i = 0; i < ns1; i++ )
            vgid[i] = i + 1;
        MB_CHK_SET_ERR( mb->tag_set_data( gtag, edges, &vgid[0] ), "can't set global id on edges" );

        std::string name_file = name_map + extension;
        MB_CHK_ERR( mb->write_mesh( name_file.c_str() ) );
        std::cout << " wrote view map file " << name_file << " with source fraction height: " << fraction << "\n";

        return 0;  // do not bother with other files created, just one file with weights on edges
    }
    // first matrix
    typedef Eigen::Triplet< double > Triplet;
    std::vector< Triplet > tripletList;
    tripletList.reserve( ns1 );
    for( int iv = 0; iv < ns1; iv++ )
    {
        // all row and column indices are 1-based in the map file.
        // inside Eigen3, we will use 0-based; then we will have to add back 1 when
        //     we dump out the files
        tripletList.push_back( Triplet( row1[iv] - 1, col1[iv] - 1, val1[iv] ) );
    }
    Eigen::SparseMatrix< double > weight1( nb1, na1 );

    weight1.setFromTriplets( tripletList.begin(), tripletList.end() );
    weight1.makeCompressed();
    EntityHandle sourceSet, targetSet;
    // those are the maps from global ids to the moab entity handles corresponding to those global ids
    //   which are corresponding to the global DOFs
    map< int, EntityHandle > sourceHandles;
    map< int, EntityHandle > targetHandles;
    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, sourceSet ), "can't create source mesh set" );
    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, targetSet ), "can't create target mesh set" );
    const char* readopts = "";
    MB_CHK_SET_ERR( mb->load_file( inputSource.c_str(), &sourceSet, readopts ), "Failed to read" );
    MB_CHK_SET_ERR( mb->load_file( inputTarget.c_str(), &targetSet, readopts ), "Failed to read" );
    Range sRange;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( sourceSet, dimSource, sRange ), "Failed to get sRange" );
    vector< int > sids;
    sids.resize( sRange.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( gtag, sRange, &sids[0] ), "Failed to get ids for srange" );
    // all global ids are 1 based in the source file, and they correspond to the dofs in the map file
    for( size_t i = 0; i < sids.size(); i++ )
    {
        EntityHandle eh    = sRange[i];
        int gid            = sids[i];
        sourceHandles[gid] = eh;
    }
    Range tRange;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( targetSet, dimTarget, tRange ), "Failed to get tRange" );
    vector< int > tids;
    tids.resize( tRange.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( gtag, tRange, &tids[0] ), "Failed to get ids for trange" );
    // all global ids are 1 based in the target file, and they correspond to the dofs in the map file
    for( size_t i = 0; i < tids.size(); i++ )
    {
        EntityHandle eh    = tRange[i];
        int gid            = tids[i];
        targetHandles[gid] = eh;
    }
    EntityHandle partialSet;
    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, partialSet ), "can't create partial set" );
    // how to get a complete row in sparse matrix? Or a complete column ?
    for( int col = startSourceID - 1; col <= endSourceID - 1; col++ )
    {
        Range targetEnts;  // will find entries for column col-1 in sparse matrix weight1, and its entries
        // will assign a dense tag with values, and write out the file
        if( col < 0 ) continue;
        Eigen::SparseVector< double > colVect = weight1.col( col );
        // the row indices correspond to target cells
        for( Eigen::SparseVector< double >::InnerIterator it( colVect ); it; ++it )
        {
            double weight = it.value();  // == vec[ it.index() ]
                                         // we add back the 1 that we subtract
            int globalIdRow = it.index() + 1;
            EntityHandle th = targetHandles[globalIdRow];
            targetEnts.insert( th );
            MB_CHK_SET_ERR( mb->tag_set_data( wtag, &th, 1, &weight ), "Failed to set weight tag on target" );
        }

        if( dimTarget == 0 )
        {
            Range adjCells;
            MB_CHK_SET_ERR( mb->get_adjacencies( targetEnts, 2, false, adjCells, Interface::UNION ),
                            " can't get adj cells " );
            targetEnts.merge( adjCells );
        }

        MB_CHK_SET_ERR( mb->add_entities( partialSet, targetEnts ), "Failed to add target entities to partial set" );
        // write now the set in a numbered file
        std::stringstream fff;
        fff << name_map << "_column" << col + 1 << extension;
        MB_CHK_ERR( mb->write_mesh( fff.str().c_str(), &partialSet, 1 ) );
        // remove from partial set the entities it has
        MB_CHK_SET_ERR( mb->clear_meshset( &partialSet, 1 ), "Failed to empty partial set" );
    }

    // how to get a complete row in sparse matrix?
    for( int row = startTargetID - 1; row <= endTargetID - 1; row++ )
    {
        Range sourceEnts;  // will find entries for row in sparse matrix weight1, and its entries
        // will assign a dense tag with values, and write out the file
        if( row < 0 ) continue;
        Eigen::SparseVector< double > rowVect = weight1.row( row );
        // the row indices correspond to target cells
        for( Eigen::SparseVector< double >::InnerIterator it( rowVect ); it; ++it )
        {
            double weight   = it.value();  // == vec[ it.index() ]
            int globalIdCol = it.index() + 1;
            EntityHandle sh = sourceHandles[globalIdCol];
            sourceEnts.insert( sh );
            MB_CHK_SET_ERR( mb->tag_set_data( wtag, &sh, 1, &weight ), "Failed to set weight tag on source" );
        }
        if( dimSource == 0 )
        {
            Range adjCells;
            MB_CHK_SET_ERR( mb->get_adjacencies( sourceEnts, 2, false, adjCells, Interface::UNION ),
                            " can't get adj cells " );
            sourceEnts.merge( adjCells );
        }
        MB_CHK_SET_ERR( mb->add_entities( partialSet, sourceEnts ), "Failed to add source entities" );
        // write now the set in a numbered file
        std::stringstream fff;
        fff << name_map << "_row" << row + 1 << extension;
        MB_CHK_ERR( mb->write_mesh( fff.str().c_str(), &partialSet, 1 ) );
        MB_CHK_SET_ERR( mb->clear_meshset( &partialSet, 1 ), "Failed to empty partial set" );
    }
    return 0;
}
