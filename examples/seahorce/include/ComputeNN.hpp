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
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );
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

        double value = 0.0, weights = 0.0;
        for( size_t j = 0; j < num_results; ++j )
        {
            value += src_tdata[srcindx[j]] / std::pow( srcdist[j], power );
            weights += 1.0 / std::pow( srcdist[j], power );
        }

        // dst_tdata[i] = src_tdata[srcindx[0]];
        dst_tdata[i] = value / weights;
    }

    return moab::MB_SUCCESS;
}

#endif  // __computeNN_hpp__