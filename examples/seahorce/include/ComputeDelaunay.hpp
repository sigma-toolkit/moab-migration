#ifndef __compute_delaunay_hpp__
#define __compute_delaunay_hpp__

#include "RemapMPASROMS.hpp"
#include "tetgen.h"  // Defined tetgenio, tetrahedralize().
#include "ComputeNN.hpp"
#include <Eigen/Dense>

// moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& xyzd,
//                                             std::vector< double >& fd,
//                                             std::vector< double >& xyzi,
//                                             std::vector< double >& fi );

constexpr double mpas_radius = 6371220.0;

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

    context.dual_mpas_tetrahedron =
        Eigen::Map< Eigen::Matrix< int, Eigen::Dynamic, 4 > >( out.tetrahedronlist, ntetrahedrons, 4 );
    context.dual_mpas_tetrahedron_centroids = std::vector< double >( ntetrahedrons, 0.0 );
    for( auto index = 0; index < ntetrahedrons; index++ )
    {
        const int offset = index * 3;
        for( auto eindex = 0; eindex < 4; eindex++ )
        {
            for( auto cindex = 0; cindex < 3; cindex++ )
                context.dual_mpas_tetrahedron_centroids[offset + cindex] +=
                    xyzd[context.dual_mpas_tetrahedron( index, eindex ) * 3 + cindex];
        }
        for( auto cindex = 0; cindex < 3; cindex++ )
            context.dual_mpas_tetrahedron_centroids[offset + cindex] /= 4;  // centroid of the tetrahedron
    }
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
typedef Eigen::Matrix< double, 3, 3 > Mat3d;

/// @brief Compute the barycentric coordinates of a tetrahedron in physical space
/// @param vcoords Vertex coordinates in 3D; size(12):= 4 vertices x 3D
/// @param point Point for which barycentric coordinates are to be evaluated
/// @param bcoords Barycentric coordinates to be returned
/// @return True if point inside tetrahedron; False otherwise
bool tetrahedron_barycentric( const double vcoords[12], const double point[3], Vec4d& bcoords )
{
    Eigen::Map< const Vec3d > a( vcoords ), b( vcoords + 3 ), c( vcoords + 6 ), d( vcoords + 9 ), p( point );

    // form the jacobian matrix for the tetrahedron
    Mat3d LHS;
    LHS << b - a, c - a, d - a;
    Vec3d RHS = p - a;

    bcoords.setZero();

    // Barycentric coordinates
    Eigen::ColPivHouseholderQR< Mat3d > dec( LHS );
    bcoords.head( 3 ) = dec.solve( RHS );
    bcoords( 3 )      = 1.0 - bcoords( 0 ) - bcoords( 1 ) - bcoords( 2 );

    // ensure that all barycentric coordinates are positive to
    // indicate point \p is inside the tetrahedron
    return !( ( bcoords.array() < 0 ).any() );
}

moab::ErrorCode ComputeDelaunayInterpolant( RuntimeContext& context,
                                            std::vector< double >& xyzd,
                                            std::vector< double >& fd,
                                            std::vector< double >& xyzi,
                                            std::vector< double >& fi )
{
    assert( context.mpas3d_dual_elems.size() );

    // moab::Range tetrahedrons;
    // runchk( SetupDelaunayInterpolant( context, xyzd ), "Computing delaunay 3D triangulation failed" );

    size_t ni = fi.size();

    // construct a kd-tree index:
    // moab::Range& tetrahedrons = context.mpas3d_dual_elems;
    // std::vector< double > tet_xyz( tetrahedrons.size() * 3 );
    // runchk( context.moab_interface->get_coords( tetrahedrons, tet_xyz.data() ) );
    PC3D< double > cloud( context.dual_mpas_tetrahedron_centroids );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    const size_t num_results = 3;
    nanoflann::KNNResultSet< double > resultSet( num_results );
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

            Vec4i connectivity = context.dual_mpas_tetrahedron( element_index, Eigen::all );
            // int nnodes;
            // runchk( context.moab_interface->get_connectivity( element, connectivity, nnodes, true ) );
            // assert( nnodes == 4 ); // we are expecting only tetrahedrons

            double vcoords[12];
            std::copy( xyzd.data() + connectivity( 0 ) * 3, xyzd.data() + connectivity( 0 ) * 3 + 3, vcoords );
            std::copy( xyzd.data() + connectivity( 1 ) * 3, xyzd.data() + connectivity( 1 ) * 3 + 3, vcoords + 3 );
            std::copy( xyzd.data() + connectivity( 2 ) * 3, xyzd.data() + connectivity( 2 ) * 3 + 3, vcoords + 6 );
            std::copy( xyzd.data() + connectivity( 3 ) * 3, xyzd.data() + connectivity( 3 ) * 3 + 3, vcoords + 9 );
            // runchk( context.moab_interface->get_coords( connectivity, nnodes, vcoords ) );
            for( auto ic = 0; ic < 4; ++ic )
                vcoords[ic * 3 + 2] /= mpas_radius;

            Vec4d bcoords;
            if( tetrahedron_barycentric( vcoords, query_pt, bcoords ) )
            {
                // std::cout << "Query: " << "Connectivity: " << connectivity << ", Barycentric coords: " << bcoords << std::endl;
                fi[index] = 0.0;
                for( auto ic = 0; ic < 4; ++ic )
                {
                    fi[index] += fd[connectivity( ic )] * bcoords( ic );
                }
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
        }
        // assert( found );
        if( !found )
        {
            std::cout << "Query point: (" << query_pt[0] << ", " << query_pt[1] << ", " << query_pt[2]
                      << ") not found in the domain" << std::endl;
            Vec4i connectivity = context.dual_mpas_tetrahedron( minindex, Eigen::all );
            fi[index]          = 1.0 / 4 *
                        ( fd[connectivity( 0 )] + fd[connectivity( 1 )] + fd[connectivity( 2 )] +
                          fd[connectivity( 3 )] );  // assign value of first vertex (extrapolant); this is arbitrary
        }
    }

    return moab::MB_SUCCESS;
}

#endif  // __compute_delaunay_hpp__