#ifndef __computeNN_hpp__
#define __computeNN_hpp__

#include "RemapMPASROMS.hpp"
#include "moab/nanoflann.hpp"

template < typename T >
struct PC3D
{
    using coord_t = T;  //!< The type of each coordinate

    const static int dimension = 3;
    const std::vector< T >& xyz;
    const size_t count;

    PC3D( const std::vector< T >& pxyz ) : xyz( pxyz ), count( pxyz.size() / dimension ) {}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return count;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt( const size_t idx, const size_t dim ) const
    {
        return xyz[idx * dimension + dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template < class BBOX >
    bool kdtree_get_bbox( BBOX& /* bb */ ) const
    {
        return false;
    }
};

using KdTree = nanoflann::
    KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >, PC3D< double >, 3 /* dim */
                              >;

/**
 * Computes the nearest neighbor interpolant for the given source and destination data.
 *
 * @param context The runtime context.
 * @param src_xyz The source coordinates.
 * @param src_tdata The source data values.
 * @param dst_xyz The destination coordinates.
 * @param dst_tdata The interpolated destination data values.
 * @return The error code indicating the success or failure of the computation.
 */
moab::ErrorCode ComputeNNInterpolant( RuntimeContext&,
                                      const std::vector< double >& src_xyz,
                                      const std::vector< double >& src_tdata,
                                      const std::vector< double >& dst_xyz,
                                      std::vector< double >& dst_tdata )
{
    constexpr double power       = 2.0;
    constexpr size_t num_results = 1;

    // construct a kd-tree index:
    PC3D< double > cloud( src_xyz );
    KdTree tree( 3 /*dim*/, cloud, { 15 /* max leaf */ } );
    KdTree::BoundingBox bbox_src;
    tree.computeBoundingBox( bbox_src );
    printf( "Source bounding boxes: (%f, %f), (%f, %f), (%3.10e, %3.10e)\n", bbox_src[0].low, bbox_src[0].high,
            bbox_src[1].low, bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

#pragma omp parallel for shared( tree, dst_xyz, src_tdata, dst_tdata )
    for( size_t i = 0; i < dst_tdata.size(); i++ )
    {
        const size_t offset = i * 3;
        std::vector< size_t > srcindx( num_results );
        std::vector< double > srcdist( num_results );
        nanoflann::KNNResultSet< double > resultSet( num_results );

        const double* query_pt = dst_xyz.data() + offset;

        // Do a KNN search
        resultSet.init( srcindx.data(), srcdist.data() );
        tree.findNeighbors( resultSet, query_pt );

        // check if the point is outside the bounding box
        if( ( query_pt[0] < bbox_src[0].low && query_pt[1] < bbox_src[1].low && query_pt[2] < bbox_src[2].low ) ||
            ( query_pt[0] > bbox_src[0].high && query_pt[1] > bbox_src[1].high && query_pt[2] > bbox_src[2].high ) )
        {
            // data needs to be extrapolated
            if( query_pt[2] < bbox_src[2].low )  // point is below the MPAS sea bed
                dst_tdata[i] = src_tdata[srcindx[0]];
            else  // point is above the MPAS sea surface
                dst_tdata[i] = src_tdata[srcindx[0]];
        }
        else
        {
            double value = 0.0, weights = 0.0;
            for( size_t j = 0; j < num_results; ++j )
            {
                const double idw = std::max( 1e-8, std::pow( srcdist[j], power ) );
                value += src_tdata[srcindx[j]] / idw;
                weights += 1.0 / idw;
            }

            // dst_tdata[i] = src_tdata[srcindx[0]];
            dst_tdata[i] = value / weights;
        }
    }

    return moab::MB_SUCCESS;
}

#endif  // __computeNN_hpp__