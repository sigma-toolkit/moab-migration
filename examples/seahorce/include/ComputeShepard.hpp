#ifndef __compute_shepard_hpp__
#define __compute_shepard_hpp__

#include "RemapMPASROMS.hpp"

moab::ErrorCode ShepardInterpolator( int dimension,
                                     std::vector< double >& xyzd,
                                     std::vector< double >& fd,
                                     double power,
                                     std::vector< double >& xyzi,
                                     std::vector< double >& fi )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    size_t nd  = xyzd.size() / dimension;
    size_t ni  = xyzi.size() / dimension;
    double ind = 1.0 / static_cast< double >( nd );

    for( size_t i = 0; i < ni; i++ )
    {
        std::vector< double > w( nd, 0.0 );
        int z;
        const int ioffset = i * dimension;
        if( power < 1 )
        {
            for( size_t j = 0; j < nd; j++ )
                w[j] = ind;
        }
        else
        {
            z = -1;
            for( size_t j = 0; j < nd; j++ )
            {
                double t          = 0.0;
                const int joffset = j * dimension;
                for( int i2 = 0; i2 < dimension; i2++ )
                {
                    t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                }
                w[j] = std::sqrt( t );
                if( w[j] < 1e-12 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < nd; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
            }
            else
            {
                double s = 0.0;
                for( size_t j = 0; j < nd; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], power );
                    s += w[j];
                }

                for( size_t j = 0; j < nd; j++ )
                    w[j] /= s;
            }
        }

        fi[i] = 0.0;
        for( size_t k = 0; k < nd; k++ )
        {
            fi[i] += w[k] * fd[k];
        }
        w.clear();
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode ModifiedShepardInterpolator( int dimension,
                                             std::vector< double >& xyzd,
                                             std::vector< double >& fd,
                                             std::vector< double >& xyzi,
                                             std::vector< double >& fi,
                                             int order )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7302837/
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    const double shepard_power = order + 1;
    // First call the regular inverse-distance weighting method
    ShepardInterpolator( dimension, xyzd, fd, shepard_power, xyzi, fi );

    size_t nd = xyzd.size() / dimension;
    size_t ni = xyzi.size() / dimension;
    std::vector< double > w( nd, 0.0 );

    double fisecnum = 0.0;
    for( size_t k = 0; k < nd; k++ )
        fisecnum += fd[k];

    assert( shepard_power >= 1.0 );

    for( size_t i = 0; i < ni; i++ )
    {
        const int ioffset = i * dimension;
        int z;
        double s = 0.0, is = 0.0;
        {
            z = -1;
            for( size_t j = 0; j < nd; j++ )
            {
                double t          = 0.0;
                const int joffset = j * dimension;
                for( int i2 = 0; i2 < dimension; i2++ )
                {
                    t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                }
                w[j] = std::sqrt( t );
                if( w[j] < 1e-6 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < nd; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
                s = is = 1.0;
            }
            else
            {
                for( size_t j = 0; j < nd; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], shepard_power );
                    s += w[j];
                    is += 1.0 / w[j];
                }

                for( size_t j = 0; j < nd; j++ )
                    w[j] /= s;
            }
        }

        fi[i] += nd * ( fisecnum - nd * fi[i] ) / ( nd * nd - s * is );
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode modified_shepard_interpolate2( int dimension,
                                               std::vector< double >& xyzd,
                                               std::vector< double >& fd,
                                               std::vector< double >& xyzi,
                                               std::vector< double >& fi )
//****************************************************************************80
//  Original source from:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//    https://people.sc.fsu.edu/~jburkardt/f_src/shepard_interp_nd/shepard_interp_nd.html
//
//  Reference:
//
//    Donald Shepard,
//    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7302837/
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
{
    const double shepard_power = 2;
    // size_t nd  = xyzd.size() / dimension;
    size_t ni = xyzi.size() / dimension;

    // construct a kd-tree index:
    using KdTree = nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PC3D< double > >,
                                                        PC3D< double >, 3 /* dim */
                                                        >;

    PC3D< double > cloud( xyzd );
    KdTree tree( 3 /*dim*/, cloud, { 10 /* max leaf */ } );

    double query_pt[3];  // dimension

    const size_t num_results =
        static_cast< size_t >( ( 2 * shepard_power + 1 ) * ( 2 * shepard_power + 1 ) - shepard_power );
    double ind = ( 1.0 / num_results );
    std::vector< size_t > srcindx( num_results );
    std::vector< double > srcdist( num_results );
    for( size_t i = 0; i < ni; i++ )
    {
        const int ioffset = i * dimension;
        {
            for( int dd = 0; dd < dimension; dd++ )
                query_pt[dd] = xyzi[ioffset + dd];
            // IntxUtils::transform_coordinates( query_pt, PointCloud< double >::projection );

            // Do a KNN search
            nanoflann::KNNResultSet< double > resultSet( num_results );
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            // tree.knnSearch( &query_pt[0], num_results, &srcindx[0], &srcdist[0] );

            // std::cout << "knnSearch(nn=" << num_results << "): \n";
            for( size_t ll = 0; ll < num_results; ++ll )
            {
                double dist = 0.0;
                for( int dd = 0; dd < dimension; dd++ )
                    dist += ( xyzi[ioffset + dd] - xyzd[srcindx[ll] * dimension + dd] ) *
                            ( xyzi[ioffset + dd] - xyzd[srcindx[ll] * dimension + dd] );
                dist = std::sqrt( dist );
                std::cout << i << "\tret_index=" << srcindx[ll] << " out_dist_sqr=" << srcdist[ll] << ", " << dist
                          << std::endl;
                srcdist[ll] = dist;
            }
        }

        std::vector< double > w( num_results, 0.0 );
        int z;
        if( shepard_power < 1 )
        {
            for( size_t j = 0; j < num_results; j++ )
                w[j] = ind;
        }
        else
        {
            z = -1;
            for( size_t j = 0; j < num_results; j++ )
            {
                // double t          = 0.0;
                // const int joffset = j * dimension;
                // for( int i2 = 0; i2 < dimension; i2++ )
                // {
                //     t += std::pow( xyzi[i2 + ioffset] - xyzd[i2 + joffset], 2.0 );
                // }
                w[j] = srcdist[j];  // std::sqrt( t );
                if( w[j] < 1e-12 )
                {
                    z = j;
                    break;
                }
            }

            if( z != -1 )
            {
                for( size_t j = 0; j < num_results; j++ )
                    w[j] = 0.0;
                w[z] = 1.0;
            }
            else
            {
                double s = 0.0;
                for( size_t j = 0; j < num_results; j++ )
                {
                    w[j] = 1.0 / std::pow( w[j], shepard_power );
                    s += w[j];
                }

                for( size_t j = 0; j < num_results; j++ )
                    w[j] /= s;
            }
        }

        fi[i] = 0.0;
        for( size_t k = 0; k < num_results; k++ )
        {
            fi[i] += w[k] * fd[k];
        }
        w.clear();
    }

    return moab::MB_SUCCESS;
}

#endif  // __compute_shepard_hpp__
