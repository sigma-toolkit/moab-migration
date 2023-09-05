
#include "RemapMPASROMS.hpp"
#include "tetgen.h"  // Defined tetgenio, tetrahedralize().
#include "ComputeNN.hpp"
#include <Eigen/Dense>

constexpr double mpas_radius = 1.0;

typedef Eigen::Matrix< double, 3, 1 > Vec3d;
typedef Eigen::Matrix< double, 4, 1 > Vec4d;
typedef Eigen::Matrix< int, 4, 1 > Vec4i;
typedef Eigen::Matrix< int, Eigen::Dynamic, 1 > VecXi;
typedef Eigen::Matrix< double, 3, 3 > Mat3d;

std::vector< int > tetrahedraconn;
Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > vertex_to_element;

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

void Setup( std::vector< double >& xyzd )
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
        in.pointlist[ix * 3 + 2] /= mpas_radius;

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
    std::vector< char > options = { 'p', 'Y', 'c', 'z', 'J', 'k' };  // 'f', 'e',
    tetgen_be.parse_commandline( options.data() );
    tetrahedralize( &tetgen_be, &in, &out );

    assert( out.numberofcorners == 4 );

    int nvertices      = out.numberofpoints;
    int ntrianglefaces = out.numberoftrifaces;
    int ntetrahedrons  = out.numberoftetrahedra;

    printf( "Number of points: %d, triangles: %d, tetrahedra: %d\n", nvertices, ntrianglefaces, ntetrahedrons );

    tetrahedraconn.resize( ntetrahedrons * 4 );
    std::copy( out.tetrahedronlist, out.tetrahedronlist + ntetrahedrons * 4, tetrahedraconn.begin() );

    vertex_to_element = Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic >( nvertices, 128 );
    vertex_to_element.setZero();
    // Store vertex to element adjacency list
    for( auto ie = 0; ie < ntetrahedrons; ++ie )
    {
        for( auto ic = 0; ic < 4; ic++ )
        {
            const int iv                              = tetrahedraconn[ ie * 4 + ic ];
            const int size                            = vertex_to_element( iv, 0 );
            vertex_to_element( iv, size + 1 ) = ie;
            vertex_to_element( iv, 0 )        = size + 1;
        }
    }

    for( auto ie = 0; ie < nvertices; ++ie )
        if( ie < 10 )
            std::cout << "Vertex: " << ie << ",  nAdjacentelements: " << vertex_to_element( ie, 0 )
                      << std::endl;

    // Output mesh to files 'mpasdel3d.node', 'mpasdel3d.ele' and 'mpasdel3d.face'.
    // out.save_nodes( "mpasdel3d" );
    // out.save_faces( "mpasdel3d" );
    // out.save_elements( "mpasdel3d" );
}

/*============================================================================*/
static bool bary_tet( const double vcoords[12], const double p[3], Vec4d& bcoords )
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
    const double* a = &vcoords[0];
    const double* b = &vcoords[3];
    const double* c = &vcoords[6];
    const double* d = &vcoords[9];
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

    int i;

    for( i = 0; i < 3; i++ )
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

    bcoords(0) = va / v;
    bcoords(1) = vb / v;
    bcoords(2) = vc / v;
    bcoords(3) = vd / v;

    return !( ( bcoords.array() < 0 ).any() );
}

/// @brief Compute the barycentric coordinates of a tetrahedron in physical space
/// @param vcoords Vertex coordinates in 3D; size(12):= 4 vertices x 3D
/// @param point Point for which barycentric coordinates are to be evaluated
/// @param bcoords Barycentric coordinates to be returned
/// @return True if point inside tetrahedron; False otherwise
static bool tetrahedron_barycentric( const double vcoords[12], const double point[3], Vec4d& bcoords )
{
    Eigen::Map< const Vec3d > a( vcoords ), b( vcoords + 3 ), c( vcoords + 6 ), d( vcoords + 9 ), p( point );

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

void ComputeDelaunayInterpolant( std::vector< double >& xyzd,
                                 std::vector< double >& fd,
                                 std::vector< double >& xyzi,
                                 std::vector< double >& fi )
{
    size_t ni = fi.size();

    std::vector< double > xyzld( xyzd );
    // for( size_t in = 0; in < ni; in++ )
    //     xyzld[3 * in + 2] /= mpas_radius;
    // PC3D< double > cloud( context.dual_mpas_tetrahedron_centroids );
    PC3D< double > cloud( xyzld );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    KdTree::BoundingBox bbox_src;
    tree.computeBoundingBox( bbox_src );
    printf( "Source bounding boxes: (%f, %f), (%f, %f), (%3.10e, %3.10e)\n", bbox_src[0].low, bbox_src[0].high, bbox_src[1].low,
            bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

    const size_t num_results = 1;
    nanoflann::KNNResultSet< double > resultSet( num_results );
    int num_not_found = 0;
    for( size_t index = 0; index < ni; index++ )
    {
        // double* query_pt = &xyzi[index * 3];
        double query_pt[3] = { xyzi[index * 3], xyzi[index * 3 + 1], xyzi[index * 3 + 2] / mpas_radius };

        // Do a KNN search
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        resultSet.init( srcindx.data(), srcdist.data() );
        tree.findNeighbors( resultSet, query_pt );

        bool found      = false;
        double mindist  = 1E6;
        size_t minindex = num_results;
        for( size_t jindex = 0; jindex < num_results; ++jindex )
        {
            // Perform the natural interpolation on the element
            size_t element_index = srcindx[jindex];

            VecXi v2e =
                vertex_to_element( element_index, Eigen::all );//.head( context.vertex_to_element( element_index, 0 ) + 1 );

            printf( "%d: Found nearest vertex with element adjacencies: %d\n ", element_index, v2e( 0 ) );

            for (int it = 1; it < v2e(0)+1; ++it)
            {
                const int* connectivity = &tetrahedraconn[ v2e( it ) * 4 ];
                // int nnodes;
                // runchk( context.moab_interface->get_connectivity( element, connectivity, nnodes, true ) );
                // assert( nnodes == 4 ); // we are expecting only tetrahedrons

                double vcoords[12];
                std::copy( xyzld.data() + connectivity[0] * 3, xyzld.data() + connectivity[0] * 3 + 3, vcoords );
                std::copy( xyzld.data() + connectivity[1] * 3, xyzld.data() + connectivity[1] * 3 + 3, vcoords + 3 );
                std::copy( xyzld.data() + connectivity[2] * 3, xyzld.data() + connectivity[2] * 3 + 3, vcoords + 6 );
                std::copy( xyzld.data() + connectivity[3] * 3, xyzld.data() + connectivity[3] * 3 + 3, vcoords + 9 );
                // runchk( context.moab_interface->get_coords( connectivity, nnodes, vcoords ) );

                Vec4d bcoords;
                if( tetrahedron_barycentric( vcoords, query_pt, bcoords ) )
                {
                    Eigen::Map< const Eigen::Matrix< int, 1, 4 > > conn( connectivity, 1, 4 );
                    fi[index] = 0.0;
                    std::cout << "Query point: " << query_pt[0] << ", " << query_pt[1] << ", " << query_pt[2]
                              << std::endl;
                    for( auto ic = 0; ic < 4; ++ic )
                    {
                        std::cout << "\t ic = " << connectivity[ic] << ", bcoords( ic ) = " << bcoords( ic )
                                  << ", pvalue = " << fd[connectivity[ic]] << std::endl;
                        fi[index] += fd[connectivity[ic]] * bcoords( ic );
                    }
                    found = true;
                    std::cout << " Interpolated value = " << fi[index] << std::endl << std::endl;
                    break;
                }
            }

            if( found ) break;
        }
        // assert( found );
        if( !found )
        {
            // std::cout << "Query point: (" << query_pt[0] << ", " << query_pt[1] << ", " << query_pt[2]
            //           << ") not found in the domain; mindist = " << mindist << ", minindex = " << minindex << std::endl;

            // std::cin.get();
            const int* connectivity = &tetrahedraconn[minindex * 4];
            fi[index]          = 1.0 / 4 *
                        ( fd[connectivity[0]] + fd[connectivity[1]] + fd[connectivity[2]] +
                          fd[connectivity[3]] );  // assign value of first vertex (extrapolant); this is arbitrary
            num_not_found++;
        }
    }

    printf( "Number of query points not found = %d/%d\n", num_not_found,  ni );
}

int main()
{
    // Tensor product hexahedra made up of two hexagons (top/bottom):
    constexpr int Nd           = 12;
    std::vector< double > xyzd = { -2, 0, 0, -1, 1, 0, 1, 1, 0, 2, 0, 0, 1, -1, 0, -1, -1, 0,
                                   -2, 0, 1, -1, 1, 1, 1, 1, 1, 2, 0, 1, 1, -1, 1, -1, -1, 1 };
    // std::vector< double > fd   = { 10, 12, 8, 3, 0.5, 0, 20, 22, 15, 6, 1, 0.1 };
    std::vector< double > fd   = { 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2 };
    constexpr int Ni           = 2;
    std::vector< double > xyzi = { -1.1, 1, 0, -1, 0, 0.5 };
    std::vector< double > fi( Ni );
    std::vector< double > fi_expected = { 1, 1.5 };

    Setup( xyzd );

    ComputeDelaunayInterpolant( xyzd, fd, xyzi, fi );

    printf( "Value at point (%f, %f, %f) = %f; Expected = %f\n", xyzi[0], xyzi[1], xyzi[2], fi[0], fi_expected[0] );
    printf( "Value at point (%f, %f, %f) = %f; Expected = %f\n", xyzi[3], xyzi[4], xyzi[5], fi[1], fi_expected[1] );

    return 0;
}
