#ifndef __compute_delaunay_hpp__
#define __compute_delaunay_hpp__

#include "RemapMPASROMS.hpp"
#include "tetgen.h"  // Defined tetgenio, tetrahedralize().
#include "ComputeNN.hpp"

// constexpr double rescale_factor = 6.37122E4;
constexpr double rescale_factor = 2E4;

typedef Eigen::Matrix< double, 4, 1 > Vec4d;
/*============================================================================*/
void crossProduct( double* ans, double* v1, double* v2 )
/*==============================================================================
 This function gives the cross product of two vectors.

 Input  Type[len]  Description
 -----  ---------  -----------
 *ans   double[3]  Pointer to an answer array. Values will be overwritten
 *v1    double[3]  Pointer to array of vector 1 components
 *v2    double[3]  Pointer to array of vector 2 components

 No outputs other than *ans

*=============================================================================*/
{
    ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return;
}

/*============================================================================*/
double dotProduct( double* v1, double* v2 )
/*==============================================================================
 This function gives the dot product of two vectors.

 Input  Type[len]  Description
 -----  ---------  -----------
 *v1    double[3]  Pointer to array of vector 1 components
 *v2    double[3]  Pointer to array of vector 2 components

 Output  Type      Description
 ------ ------     -----------
 result double     the dot product

*=============================================================================*/
{
    double result = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    return ( result );
}

/*============================================================================*/
bool tetrahedron_barycentric2( const double* a,
                               const double* b,
                               const double* c,
                               const double* d,
                               const double p[3],
                               Vec4d& bcoords )
/*==============================================================================
 This function gives the barycentric coordinates for a tetrahedron in 3 dimensions.
 It works by calculating the volume of the subtetrahedron using the scalar triple
 product. { V_tet = 1/6 * v1 * ( v2 x v2 ) }

 This works for points outside the tetrahedron as well.

 Input  Type[len]  Description
 -----  ---------  -----------
 *ans    double[4]  Pointer to an answer array. Values 0-3 will be overwritten
 *p      double[3]  Pointer to array of (x,y,z) coordinates for the test point
 *a      double[3]  Pointer to array of (x,y,z) coordinates for tet node 1
 *b      double[3]  Pointer to array of (x,y,z) coordinates for tet node 2
 *c      double[3]  Pointer to array of (x,y,z) coordinates for tet node 3
 *d      double[3]  Pointer to array of (x,y,z) coordinates for tet node 4

 No outputs other than *ans

 The node order for a tetrahedron is conterclockwise around the base, then up.

         4               3
          /|\             |.
         / | \            | . .
        /  |  \           |  .   .  2
       /   |   \          |    .  /\
      /    |    \         |     ./  \
     /     |     \        |     / .  \
    /      |      \       |    /   .  \
 0 /_______|_______\ 2    |   /     .  \
   \       |       /      |  /        . \
    \      |      /       | /          . \
     \     |     /        |/_____________.\
      \    |    /        0                  1
       \   |   /
        \  |  /
         \ | /
          \|/
            1


*=============================================================================*/
{
    // const double* a = &vcoords[0];
    // const double* b = &vcoords[3];
    // const double* c = &vcoords[6];
    // const double* d = &vcoords[9];
    double vap[3];
    double vbp[3];
    double vcp[3];
    double vdp[3];
    double vab[3];
    double vac[3];
    double vad[3];
    double vbc[3];
    double vbd[3];
    double va;
    double vb;
    double vc;
    double vd;
    double v;
    double temp[3];

    for( int i = 0; i < 3; i++ )
    {
        vap[i] = p[i] - a[i];
        vbp[i] = p[i] - b[i];
        vcp[i] = p[i] - c[i];
        vdp[i] = p[i] - d[i];
        vab[i] = b[i] - a[i];
        vac[i] = c[i] - a[i];
        vad[i] = d[i] - a[i];
        vbc[i] = c[i] - b[i];
        vbd[i] = d[i] - b[i];
    }
    crossProduct( temp, vbd, vbc );
    va = dotProduct( vbp, temp ) / 6.0;
    crossProduct( temp, vac, vad );
    vb = dotProduct( vap, temp ) / 6.0;
    crossProduct( temp, vad, vab );
    vc = dotProduct( vap, temp ) / 6.0;
    crossProduct( temp, vab, vac );
    vd = dotProduct( vap, temp ) / 6.0;
    crossProduct( temp, vac, vad );
    v = dotProduct( vab, temp ) / 6.0;

    bcoords( 0 ) = va / v;
    bcoords( 1 ) = vb / v;
    bcoords( 2 ) = vc / v;
    bcoords( 3 ) = vd / v;

    return !( ( bcoords.array() < 0 ).any() );
}

moab::ErrorCode SetupDelaunayInterpolant( RuntimeContext& context, std::vector< double >& xyzd )
{
    tetgenio in, out;

    // All indices start from 1.
    in.firstnumber = 0;

    in.numberofpoints = xyzd.size() / 3;
    in.pointlist      = new REAL[in.numberofpoints * 3];
    // Set node coordinates: memcpy the data
    std::copy( xyzd.begin(), xyzd.end(), in.pointlist );
    in.numberofpointattributes = 0;

    for( auto ix = 0; ix < in.numberofpoints; ix++ )
        in.pointlist[ix * 3 + 2] /= rescale_factor;

    // Output the PLC to files 'mpas3d.node'
    // in.save_nodes( "mpas3d" );

    //////////////////////////////////////////////////////////////
    // Tetrahedralize the PLC with switches are chosen as follows:
    //
    // -p: Tetrahedralizes a piecewise linear complex (PLC).
    // -Y: Preserves the input surface mesh (does not modify it).
    // -M: No merge of coplanar facets or very close vertices.
    // -c: Retains the convex hull of the PLC.
    // -z: Numbers all output items starting from zero.
    // -f: Outputs all faces to .face file.
    // -e: Outputs all edges to .edge file.
    // -n: Outputs tetrahedra neighbors to .neigh file.
    // -k: Outputs mesh to .vtk file for viewing by Paraview.
    // -J: No jettison of unused vertices from output .node file.
    tetgenbehavior tetgen_be;
    std::vector< char > options = { 'p', 'Y', 'c', 'z', 'J' };  // 'f', 'e', 'k'
    tetgen_be.parse_commandline( options.data() );
    tetgen_be.opt_max_flip_level = 0;
    tetgen_be.opt_iterations     = 1;
    tetrahedralize( &tetgen_be, &in, &out );

    assert( out.numberofcorners == 4 );

    int nvertices      = out.numberofpoints;
    int ntrianglefaces = out.numberoftrifaces;
    int ntetrahedrons  = out.numberoftetrahedra;
    int ncorners       = out.numberofcorners;

    printf( "Number of points: %d, triangles: %d, corners: %d, tetrahedra: %d\n", nvertices, ntrianglefaces, ncorners,
            ntetrahedrons );

    // moab::ReadUtilIface* read_iface;
    // runchk( context.moab_interface->query_interface( read_iface ), "Error in query_interface" );

    // // Create quads
    // moab::EntityHandle start_elem;
    // moab::EntityHandle* connect;
    // runchk( read_iface->get_element_connect( ntetrahedrons, 4, moab::MBTET, 0, start_elem, connect ),
    //         "Error in get_element_connect" );
    // moab::Range& tetrahedrons = context.mpas3d_dual_elems;
    // tetrahedrons = moab::Range( start_elem, start_elem + ntetrahedrons );
    // std::copy( out.tetrahedronlist, out.tetrahedronlist + ntetrahedrons * 4, connect );

    // context.dual_mpas_tetrahedron =
    //     Eigen::Map< Eigen::Matrix< int, Eigen::Dynamic, 4 > >( out.tetrahedronlist, ntetrahedrons, 4 );

    assert( ncorners == 4 );
    context.tetrahedraconn.resize( ntetrahedrons * ncorners );
    std::copy( out.tetrahedronlist, out.tetrahedronlist + ntetrahedrons * ncorners, context.tetrahedraconn.begin() );

#ifdef USE_DUAL_TETS
    std::vector< double > tetcentroids( ntetrahedrons * 3, 0.0 );
    // #pragma omp parallel for shared( tetcentroids, xyzd )
    for( auto index = 0; index < ntetrahedrons; index++ )
    {
        const int offset = index * 3;
        for( auto eindex = 0; eindex < ncorners; eindex++ )
        {
            for( auto cindex = 0; cindex < 3; cindex++ )
                tetcentroids[offset + cindex] +=
                    in.pointlist[context.tetrahedraconn[index * ncorners + eindex] * 3 + cindex];
        }
        for( auto cindex = 0; cindex < 3; cindex++ )
            tetcentroids[offset + cindex] /= ncorners;  // centroid of the tetrahedron
        // tetcentroids[offset + 2] /= rescale_factor;
        // std::cout << "Tetrahedron: " << index << " Coordinates: " << tetcentroids[offset] << ", "
        //           << tetcentroids[offset + 1] << ", " << tetcentroids[offset + 2] << std::endl;
    }
    context.dual_mpas_tetrahedron_centroids = tetcentroids;
#else
    context.nvertcache        = std::vector< int >( nvertices, 0 );
    context.vertex_to_element = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic >( nvertices, N_MAX_V2EADJACENCIES );
    context.vertex_to_element.setZero();
    // context.vertex_to_element += 1;
    // Store vertex to element adjacency list
    for( auto ie = 0; ie < ntetrahedrons; ++ie )
    {
        for( auto ic = 0; ic < ncorners; ic++ )
        {
            const int iv = context.tetrahedraconn[ie * ncorners + ic];
            // const int size                            = context.vertex_to_element( iv, 0 );
            context.vertex_to_element( iv, context.nvertcache[iv] ) = ie;
            // context.vertex_to_element( iv, 0 )        = size + 1;
            context.nvertcache[iv]++;
        }
    }
#endif

    // for( auto ie = 0; ie < nvertices; ++ie )
    //     std::cout << "Vertex: " << ie << ",  nAdjacentelements: " << context.vertex_to_element( ie, 0 ) << std::endl;

    // std::cin.get();
    // context.dual_mpas_tetrahedron_centroids = std::vector< double >( ntetrahedrons, 0.0 );
    // for( auto index = 0; index < ntetrahedrons; index++ )
    // {
    //     const int offset = index * 3;
    //     for( auto eindex = 0; eindex < 4; eindex++ )
    //     {
    //         for( auto cindex = 0; cindex < 3; cindex++ )
    //             context.dual_mpas_tetrahedron_centroids[offset + cindex] +=
    //                 xyzd[context.dual_mpas_tetrahedron( index, eindex ) * 3 + cindex];
    //     }
    //     for( auto cindex = 0; cindex < 3; cindex++ )
    //         context.dual_mpas_tetrahedron_centroids[offset + cindex] /= 4;  // centroid of the tetrahedron
    // }

    // Output mesh to files 'mpasdel3d.node', 'mpasdel3d.ele' and 'mpasdel3d.face'.
    // out.save_nodes( "mpasdel3d" );
    // out.save_faces( "mpasdel3d" );
    // out.save_elements( "mpasdel3d" );

    return moab::MB_SUCCESS;
}

typedef Eigen::Matrix< double, 3, 1 > Vec3d;
typedef Eigen::Matrix< double, 4, 1 > Vec4d;
typedef Eigen::Matrix< int, 4, 1 > Vec4i;
typedef Eigen::Matrix< int, Eigen::Dynamic, 1 > VecXi;
typedef Eigen::Matrix< double, 3, 3 > Mat3d;

/// @brief Compute the barycentric coordinates of a tetrahedron in physical space
/// @param vcoords Vertex coordinates in 3D; size(12):= 4 vertices x 3D
/// @param point Point for which barycentric coordinates are to be evaluated
/// @param bcoords Barycentric coordinates to be returned
/// @return True if point inside tetrahedron; False otherwise
bool tetrahedron_barycentric( const double* ar,
                              const double* br,
                              const double* cr,
                              const double* dr,
                              const double* point,
                              Vec4d& bcoords )
{
    Eigen::Map< const Vec3d > a( ar ), b( br ), c( cr ), d( dr ), p( point );

    // form the jacobian matrix for the tetrahedron
    Mat3d LHS;
    LHS << b - a, c - a, d - a;
    Vec3d RHS = p - a;

    bcoords.setZero();

    // Barycentric coordinates
    // Eigen::ColPivHouseholderQR< Mat3d > dec( LHS );
    // bcoords.head( 3 ) = dec.solve( RHS );
    bcoords.head( 3 ) = LHS.inverse() * RHS;
    bcoords( 3 )      = 1.0 - bcoords( 0 ) - bcoords( 1 ) - bcoords( 2 );

    // ensure that all barycentric coordinates are positive to
    // indicate point \p is inside the tetrahedron
    // std::cout << "Query: " << "point: " << p.transpose() << ", Barycentric coords: " << bcoords << std::endl;
    return !( ( bcoords.array() < 0 ).any() );
}

moab::ErrorCode ComputeDelaunayInterpolant( RuntimeContext& context,
                                            const std::vector< double >& xyzd,
                                            std::vector< double >& fd,
                                            const std::vector< double >& xyzi,
                                            std::vector< double >& fi )
{
    // runchk( SetupDelaunayInterpolant( context, xyzd ), "Computing delaunay 3D triangulation failed" );

    size_t nd = fd.size();
    size_t ni = fi.size();

    std::vector< double > xyzld( xyzd );
    for( size_t in = 0; in < nd; in++ )
        xyzld[3 * in + 2] /= rescale_factor;
#ifdef USE_DUAL_TETS
    PC3D< double > cloud( context.dual_mpas_tetrahedron_centroids );
    const size_t num_results = 64;
#else
    PC3D< double > cloud( xyzld );
    const size_t num_results = 1;
#endif

    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    KdTree::BoundingBox bbox_src;
    tree.computeBoundingBox( bbox_src );
    printf( "Source bounding boxes: (%f, %f), (%f, %f), (%3.10e, %3.10e)\n", bbox_src[0].low, bbox_src[0].high,
            bbox_src[1].low, bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

    nanoflann::SearchParameters sparams;
    nanoflann::KNNResultSet< double > resultSet( num_results );
    int num_not_found = 0;
    for( size_t index = 0; index < ni; index++ )
    {
        // double* query_pt = &xyzi[index * 3];
        double query_pt[3] = { xyzi[index * 3], xyzi[index * 3 + 1], xyzi[index * 3 + 2] / rescale_factor };

        // Do a KNN search
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        resultSet.init( srcindx.data(), srcdist.data() );

        tree.findNeighbors( resultSet, query_pt, sparams );
        // bool treefound = tree.closest( resultSet, query_pt, sparams );
        assert( treefound );

        bool found      = false;
        double mindist  = srcdist[0];
        size_t minindex = srcindx[0];
        for( size_t jindex = 0; jindex < num_results; ++jindex )
        {
            // printf( "srcindex: %zu --  jindex=%zu, eindex=%lu, edistance=%f\n", index, jindex, srcindx[jindex], srcdist[jindex] );
            // Perform the natural interpolation on the element
            size_t element_index = srcindx[jindex];

#ifdef USE_DUAL_TETS
            const int* connectivity = context.tetrahedraconn.data() + element_index * 4;
            const double* a         = xyzld.data() + connectivity[0] * 3;
            const double* b         = xyzld.data() + connectivity[1] * 3;
            const double* c         = xyzld.data() + connectivity[2] * 3;
            const double* d         = xyzld.data() + connectivity[3] * 3;

            Vec4d bcoords;
            if( tetrahedron_barycentric2( a, b, c, d, query_pt,
                                          bcoords ) )  // tetrahedron_barycentric, tetrahedron_barycentric2
            {
                // std::cout << "Query: " << "Connectivity: " << connectivity << ", Barycentric coords: " << bcoords << std::endl;
                fi[index] = 0.0;
                for( auto ic = 0; ic < 4; ++ic )
                {
                    fi[index] += fd[connectivity[ic]] * bcoords( ic );
                }
                // if( index < 100000 )
                // printf( "%d: Found point at result %d with distance=%3.8e\n", index, jindex, srcdist[jindex] );
                found = true;
                break;
            }
            else
            {
                if( srcdist[jindex] < mindist )
                {
                    mindist  = srcdist[jindex];
                    minindex = srcindx[jindex];
                }
            }
#else
            VecXi v2e =
                context.vertex_to_element( element_index,
                                           Eigen::all );  //.head( context.vertex_to_element( element_index, 0 ) + 1 );

            // printf( "%d: Found nearest vertex with element adjacencies: %d\n ", element_index, v2e( 0 ) );

            for( int it = 0; it < context.nvertcache[element_index]; ++it )
            {
                const int* connectivity = context.tetrahedraconn.data() + v2e( it ) * 4;
                const double* a         = xyzld.data() + connectivity[0] * 3;
                const double* b         = xyzld.data() + connectivity[1] * 3;
                const double* c         = xyzld.data() + connectivity[2] * 3;
                const double* d         = xyzld.data() + connectivity[3] * 3;

                // printf( "Vertex{%d}: element=%d, connectivity=[%d, %d, %d, %d]\n", index, v2e( it ), connectivity[0],
                //         connectivity[1], connectivity[2], connectivity[3] );

                Vec4d bcoords;
                if( tetrahedron_barycentric2( a, b, c, d, query_pt,
                                              bcoords ) )  // tetrahedron_barycentric, tetrahedron_barycentric2
                {
                    Eigen::Map< const Eigen::Matrix< int, 1, 4 > > conn( connectivity, 1, 4 );
                    fi[index] = 0.0;
                    // std::cout << "Query point: " << query_pt[0] << ", " << query_pt[1] << ", " << query_pt[2]
                    //   << std::endl;
                    for( auto ic = 0; ic < 4; ++ic )
                    {
                        // std::cout << "\t ic = " << connectivity[ic] << ", bcoords( ic ) = " << bcoords( ic )
                        //           << ", pvalue = " << fd[connectivity[ic]] << std::endl;
                        fi[index] += fd[connectivity[ic]] * bcoords( ic );
                    }
                    found = true;
                    // std::cout << " Interpolated value = " << fi[index] << std::endl << std::endl;
                    break;

                    // // std::cout << "Query: " << "Connectivity: " << connectivity << ", Barycentric coords: " << bcoords << std::endl;
                    // fi[index] = 0.0;
                    // for( auto ic = 0; ic < 4; ++ic )
                    // {
                    //     fi[index] += fd[connectivity( ic )] * bcoords( ic );
                    // }
                    // found = true;
                    // break;
                }
                else
                {
                    if( srcdist[jindex] < mindist )
                    {
                        mindist  = srcdist[jindex];
                        minindex = v2e( it );
                    }
                }
            }

            if( found ) break;
#endif
        }
        // assert( found );
        if( !found )
        {
            // std::cout << "Query point: (" << query_pt[0] << ", " << query_pt[1] << ", " << query_pt[2]
            //           << ") not found in the domain; mindist = " << mindist << ", minindex = " << minindex << std::endl;

            // std::cin.get();
            const int* connectivity = context.tetrahedraconn.data() + minindex * 4;
            fi[index]               = 1.0 / 4 *
                        ( fd[connectivity[0]] + fd[connectivity[1]] + fd[connectivity[2]] +
                          fd[connectivity[3]] );  // assign value of first vertex (extrapolant); this is arbitrary
            // fi[index] = fd[connectivity[0]];
            num_not_found++;
        }
    }

    printf( "Number of query points not found = %d/%zu\n", num_not_found, ni );

    return moab::MB_SUCCESS;
}

#endif  // __compute_delaunay_hpp__
