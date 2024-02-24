#ifndef __compute_mba_hpp__
#define __compute_mba_hpp__

#include "RemapMPASROMS.hpp"
#include "moab/Remapping/MBA.hpp"
#include "ComputeNN.hpp"

/**
 * Computes the MBA interpolant for a given set of data points.
 *
 * @param context The runtime context.
 * @param xyzd The input coordinates of the data points.
 * @param fd The input function values at the data points.
 * @param xyzi The output interpolated coordinates.
 * @param fi The output interpolated function values.
 * @param is_threed Flag indicating whether the data points are in 3D.
 * @param order The order of the interpolant.
 * @param grid The size of the grid used for interpolation.
 * @param bbox_user The user-defined bounding box for the data points.
 * @param nlevels The number of levels used for multilevel interpolation. Default is 7.
 *
 * @return The error code indicating the success or failure of the computation.
 */
moab::ErrorCode ComputeMBAInterpolant( RuntimeContext& context,
                                       std::vector< double >& xyzd,
                                       std::vector< double >& fd,
                                       std::vector< double >& xyzi,
                                       std::vector< double >& fi,
                                       bool is_threed,
                                       int order,
                                       std::array< size_t, 3 >& grid,
                                       std::array< double, 6 >& bbox_user,
                                       int nlevels = 7 )
{
    const size_t nd = fd.size();
    const size_t ni = fi.size();

    if( false )
    {
        if( is_threed )
        {
            std::vector< double > xd( nd ), yd( nd ), zd( nd ), xi( ni ), yi( ni ), zi( ni );
#pragma omp parallel for
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                xd[k]               = xyzd[offset];
                yd[k]               = xyzd[offset + 1];
                zd[k]               = xyzd[offset + 2];
            }
#pragma omp parallel for
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                xi[k]               = xyzi[offset];
                yi[k]               = xyzi[offset + 1];
                zi[k]               = xyzi[offset + 2];
            }
            // mlinterp::interp( &nd, ni, fd, fi, xd, xi, yd, yi, zd, zi );

            return moab::MB_SUCCESS;
        }
        else
        {
            std::vector< double > xd( nd ), yd( nd ), xi( ni ), yi( ni );
#pragma omp parallel for
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                xd[k]               = xyzd[offset];
                yd[k]               = xyzd[offset + 1];
            }
#pragma omp parallel for
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                xi[k]               = xyzi[offset];
                yi[k]               = xyzi[offset + 1];
            }
            // mlinterp::interp( &nd, ni, fd, fi, xd, xi, yd, yi );

            return moab::MB_SUCCESS;
        }
    }

    constexpr bool use_recursive = false;
    std::vector< double > filocal;
    if( use_recursive ) filocal.resize( fi.size() );

    // Algorithm setup.
    if( order == 1 )
    {
        mba::linear_approximation< 3 >* interp;
        if( use_recursive )
        {
            moab::ErrorCode err = ComputeNNInterpolant( context, xyzd, fd, xyzi, filocal );MB_CHK_ERR( err );

            std::vector< mba::point< 3 > > coords( ni );
#pragma omp parallel for shared( coords )
            for( size_t k = 0; k < ni; k++ )
            {
                const size_t offset = k * 3;
                coords[k]           = mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] };
            }

            interp = new mba::linear_approximation< 3 >( coords.begin(), coords.end(), filocal.begin() );
        }
        else
        {
            std::vector< mba::point< 3 > > coords( nd );
#pragma omp parallel for shared( coords )
            for( size_t k = 0; k < nd; k++ )
            {
                const size_t offset = k * 3;
                coords[k]           = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };
                // coords[k] = mba::point< 2 >{ xyzd[offset], xyzd[offset + 1] };
            }

            interp = new mba::linear_approximation< 3 >( coords.begin(), coords.end(), fd.begin() );
        }

        // Get interpolated value at arbitrary location.
#pragma omp parallel for shared( fi, interp )
        for( size_t k = 0; k < ni; k++ )
        {
            const size_t offset = k * 3;
            fi[k]               = ( *interp )( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );
        }
        // fi[k] = interp( mba::point< 2 >{ xyzi[offset], xyzi[offset + 1] } );

        delete interp;
    }
    else
    {
        std::vector< mba::point< 3 > > coords( nd );
#pragma omp parallel for shared( coords, xyzd )
        for( size_t k = 0; k < nd; k++ )
        {
            const size_t offset = k * 3;
            coords[k]           = mba::point< 3 >{ xyzd[offset], xyzd[offset + 1], xyzd[offset + 2] };
        }

        // construct a kd-tree index:
        if( false )
        {
            PC3D< double > cloud_src( xyzd );
            KdTree tree_src( 3 /*dim*/, cloud_src, { 5 /* max leaf */ } );
            // construct a kd-tree index:
            PC3D< double > cloud_tgt( xyzi );
            KdTree tree_tgt( 3 /*dim*/, cloud_tgt, { 5 /* max leaf */ } );

            KdTree::BoundingBox bbox_src, bbox_tgt;
            tree_src.computeBoundingBox( bbox_src );
            tree_tgt.computeBoundingBox( bbox_tgt );

            printf( "Source bounding boxes: (%f, %f), (%f, %f), (%f, %f)\n", bbox_src[0].low, bbox_src[0].high,
                    bbox_src[1].low, bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

            printf( "Target bounding boxes: (%f, %f), (%f, %f), (%f, %f)\n", bbox_tgt[0].low, bbox_tgt[0].high,
                    bbox_tgt[1].low, bbox_tgt[1].high, bbox_tgt[2].low, bbox_tgt[2].high );
        }

        // Bounding box containing the data points.
        mba::point< 3 > lo    = { bbox_user[0], bbox_user[1], bbox_user[2] };
        mba::point< 3 > hi    = { bbox_user[3], bbox_user[4], bbox_user[5] };
        mba::MBA< 3 >* interp = nullptr;
        if( use_recursive )
        {
            moab::ErrorCode err = ComputeNNInterpolant( context, xyzd, fd, xyzi, filocal );MB_CHK_ERR( err );

            PC3D< double > cloud_src( xyzd );
            KdTree tree_src( 3 /*dim*/, cloud_src, { 5 /* max leaf */ } );
            const auto& tmpTree = tree_src;

            std::function< double( mba::point< 3 > ) > initFn = [&tmpTree, fd]( mba::point< 3 > query_pt ) {
                size_t srcindx;
                double srcdist;
                nanoflann::KNNResultSet< double > resultSet( 1 );

                // Do a KNN search
                resultSet.init( &srcindx, &srcdist );
                // printf( "Querying point: %f %f %f\n", query_pt[0], query_pt[1], query_pt[2] );
                tmpTree.findNeighbors( resultSet, query_pt.data() );

                return fd[srcindx];
            };

            interp = new mba::MBA< 3 >( lo, hi, grid, coords, fd, nlevels /*levels*/, 1e-12 /*tolerance*/,
                                        0.5 /*min_fill*/, initFn );
        }
        else
        {
            interp = new mba::MBA< 3 >( lo, hi, grid, coords, fd, nlevels /*levels*/, 1e-12 /*tolerance*/,
                                        0.5 /*min_fill*/ );
        }

        // Get interpolated value at arbitrary location.
        std::cout << "\nEvaluating the interpolant now...\n";
#pragma omp parallel for shared( fi, interp, xyzi ) schedule( guided, 64 )
        for( size_t k = 0; k < ni; k++ )
        {
            auto offset = k * 3;
            fi[k]       = ( *interp )( mba::point< 3 >{ xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] } );
        }

        delete interp;
    }

    return moab::MB_SUCCESS;
}

#endif  // __compute_mba_hpp__
