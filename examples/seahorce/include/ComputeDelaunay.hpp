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

moab::ErrorCode SetupDelaunayInterpolant( RuntimeContext& context,
                                          std::vector< double >& xyzd )
{
    tetgenio in, out;

    // All indices start from 1.
    in.firstnumber = 1;

    in.numberofpoints = xyzd.size() / 3;
    in.pointlist      = new REAL[in.numberofpoints * 3];
    // Set node coordinates: memcpy the data
    std::copy( xyzd.begin(), xyzd.end(), in.pointlist );
    in.numberofpointattributes = 0;

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
    std::vector< char > options = { 'p', 'Y', 'c', 'z', 'J', 'f', 'e', 'k' };
    tetgen_be.parse_commandline( options.data() );
    tetrahedralize( &tetgen_be, &in, &out );

    assert( out.numberofcorners == 4 );

    int nvertices      = out.numberofpoints;
    int ntrianglefaces = out.numberoftrifaces;
    int ntetrahedrons  = out.numberoftetrahedra;

    printf( "Number of points: %d, triangles: %d, tetrahedra: %d\n", nvertices, ntrianglefaces, ntetrahedrons );

    moab::ReadUtilIface* read_iface;
    runchk( context.moab_interface->query_interface( read_iface ), "Error in query_interface" );

    // Create quads
    moab::EntityHandle start_elem;
    moab::EntityHandle* connect;
    runchk( read_iface->get_element_connect( ntetrahedrons, 4, moab::MBTET, 0, start_elem, connect ),
            "Error in get_element_connect" );
    moab::Range& tetrahedrons = context.mpas3d_dual_elems;
    tetrahedrons = moab::Range( start_elem, start_elem + ntetrahedrons );
    std::copy( out.tetrahedronlist, out.tetrahedronlist + ntetrahedrons * 4, connect );

    // Output mesh to files 'mpasdel3d.node', 'mpasdel3d.ele' and 'mpasdel3d.face'.
    // out.save_nodes( "mpasdel3d" );
    // out.save_faces( "mpasdel3d" );
    // out.save_elements( "mpasdel3d" );

    return moab::MB_SUCCESS;
}

typedef Eigen::Matrix< double, 3, 1 > Vec3d;
typedef Eigen::Matrix< double, 4, 1 > Vec4d;
typedef Eigen::Matrix< double, 3, 3 > Mat3d;

/// @brief Compute the barycentric coordinates of a tetrahedron in physical space
/// @param vcoords Vertex coordinates in 3D; size(12):= 4 vertices x 3D
/// @param point Point for which barycentric coordinates are to be evaluated
/// @param bcoords Barycentric coordinates to be returned
/// @return True if point inside tetrahedron; False otherwise
bool tetrahedron_barycentric( const double vcoords[12], const double point[3], Vec4d& bcoords )
{
    Eigen::Map< const Vec3d > a( vcoords ), b( vcoords + 3 ), c( vcoords + 6 ), d( vcoords + 9 ), p(point);

    // form the jacobian matrix for the tetrahedron
    Mat3d LHS;
    LHS << b-a, c-a, d-a;
    Vec3d RHS = p - a;

    bcoords.setZero();

    // Barycentric coordinates
    Eigen::ColPivHouseholderQR< Mat3d > dec( LHS );
    bcoords.head(3) = dec.solve( RHS );
    bcoords( 3 )    = 1.0 - bcoords( 0 ) - bcoords( 1 ) - bcoords( 2 );

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

    size_t nd = fd.size();
    size_t ni = fi.size();

    // construct a kd-tree index:
    moab::Range& tetrahedrons = context.mpas3d_dual_elems;
    std::vector< double > tet_xyz( tetrahedrons.size() * 3 );
    runchk( context.moab_interface->get_coords( tetrahedrons, tet_xyz.data() ) );
    PC3D< double > cloud( tet_xyz );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    const size_t num_results = 3;
    std::vector< size_t > srcindx( num_results );
    std::vector< double > srcdist( num_results );
    nanoflann::KNNResultSet< double > resultSet( num_results );
    for( size_t index = 0; index < ni; index++ )
    {
        double* query_pt = &xyzi[index * 3];

        // Do a KNN search
        resultSet.init( srcindx.data(), srcdist.data() );
        tree.findNeighbors( resultSet, query_pt );

        bool found = false;
        for( size_t jindex = 0; jindex < num_results; ++jindex )
        {
            // Perform the natural interpolation on the element
            moab::EntityHandle element = context.mpas3d_dual_elems[srcindx[jindex]];

            const moab::EntityHandle* connectivity;
            int nnodes;
            runchk( context.moab_interface->get_connectivity( element, connectivity, nnodes, true ) );
            assert( nnodes == 4 ); // we are expecting only tetrahedrons

            double vcoords[12];
            runchk( context.moab_interface->get_coords( connectivity, nnodes, vcoords ) );

            Vec4d bcoords;
            if( tetrahedron_barycentric( vcoords, query_pt, bcoords ) )
            {
                fi[index] = 0.0;
                for (auto ic = 0; ic < 4; ++ic)
                    fi[index] += fd[connectivity[ic]] * bcoords(ic);
                found = true;
                break;
            }
        }
        assert( found );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ComputeDelaunayInterpolantTetgen( std::vector< double >& xyzd,
                                            std::vector< double >& fd,
                                            std::vector< double >& /*xyzi*/,
                                            std::vector< double >& /*fi*/ )
{
    // int nd = fd.size();
    // int ni = fi.size();

    tetgenio in, out;
    tetgenio::facet* f;
    tetgenio::polygon* p;

    // All indices start from 1.
    in.firstnumber = 0;

    in.numberofpoints = fd.size();
    in.pointlist      = new REAL[in.numberofpoints * 3];
    // Set node coordinates: memcpy the data
    std::copy( xyzd.begin(), xyzd.end(), in.pointlist );
    in.numberofpointattributes = 0;

    // Output the PLC to files 'barin.node' and 'barin.poly'.
    in.save_nodes( "mpas3d" );

    // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //   do quality mesh generation (q) with a specified quality bound
    //   (1.414), and apply a maximum volume constraint (a0.1).

    tetgenbehavior tetgen_be;
    char* options = "pYczJ";
    tetgen_be.parse_commandline( options );
    tetrahedralize( &tetgen_be, &in, &out );

    printf( "Number of points: %d, triangles: %d, tetrahedra: %d\n", out.numberofpoints, out.numberoftrifaces, out.numberoftetrahedra );

    // Output mesh to files 'mpasdel3d.node', 'mpasdel3d.ele' and 'mpasdel3d.face'.
    out.save_nodes( "mpasdel3d" );
    out.save_elements( "mpasdel3d" );
    out.save_faces( "mpasdel3d" );

    return moab::MB_SUCCESS;
}


#endif  // __compute_delaunay_hpp__