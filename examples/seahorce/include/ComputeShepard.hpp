#ifndef __compute_shepard_hpp__
#define __compute_shepard_hpp__

#include "RemapMPASROMS.hpp"

#ifdef MOAB_HAVE_ALGLIB

// alglib includes
#include "interpolation.h"

moab::ErrorCode ShepardInterpolatorAlgLib( int dimension,
                                           std::vector< double >& xyzd,
                                           std::vector< double >& fd,
                                           std::vector< double >& xyzi,
                                           std::vector< double >& fi )
{
    using namespace alglib;
    size_t nd  = xyzd.size() / dimension;
    size_t ni  = xyzi.size() / dimension;

    // NOTE: we can work with N-dimensional models and vector-valued functions too :)
    //
    // Typical sequence of steps is given below:
    // 1. we create IDW builder object
    // 2. we attach our dataset to the IDW builder and tune algorithm settings
    // 3. we generate IDW model
    // 4. we use IDW model instance (evaluate, serialize, etc.)
    //

    //
    // Step 1: IDW builder creation.
    //
    // We have to specify dimensionality of the space (2 or 3) and
    // dimensionality of the function (scalar or vector).
    //
    // New builder object is empty - it has not dataset and uses
    // default model construction settings
    //
    idwbuilder builder;
    idwbuildercreate( dimension, 1, builder );

    alglib::real_2d_array indata;
    // ae_matrix_wrapper indata;
    indata.setlength( nd, dimension + 1 );
    for( size_t i = 0; i < nd; i++ )
    {
        const auto offset = i * 3;
        for( int d = 0; d < dimension; d++ )
            indata( i, d ) = xyzd[offset + d];
        indata( i, dimension ) = fd[i];
    }

    idwbuildersetpoints( builder, indata );

    //
    // Step 3: choose IDW algorithm and generate model
    //
    // We use modified stabilized IDW algorithm with following parameters:
    // * SRad - set to 5.0 (search radius must be large enough)
    //
    // IDW-MSTAB algorithm is a state-of-the-art implementation of IDW which
    // is competitive with RBFs and bicubic splines. See comments on the
    // idwbuildersetalgomstab() function for more information.
    //
    idwmodel model;
    idwreport rep;
    idwbuildersetalgomstab( builder, 50.0 );
    idwfit( builder, model, rep );

    //
    // Step 4: model was built, evaluate its value
    //
    for( size_t i = 0; i < ni; i++ )
    {
        const auto offset = i * 3;
        fi[i]             = idwcalc3( model, xyzi[offset], xyzi[offset + 1], xyzi[offset + 2] );
    }

    return moab::MB_SUCCESS;
}

#endif

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
                // std::cout << i << "\tret_index=" << srcindx[ll] << " out_dist_sqr=" << srcdist[ll] << ", " << dist
                //           << std::endl;
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

moab::ErrorCode ComputeHierarchicalShepardInterpolant( RuntimeContext&,
                                                       const std::vector< double >& src_xyz,
                                                       const std::vector< double >& src_tdata,
                                                       const std::vector< double >& dst_xyz,
                                                       std::vector< double >& dst_tdata,
                                                       int order )
{
    const double power = 1.0 * order;

    // construct a kd-tree index:
    PC3D< double > cloud( src_xyz );
    KdTree tree( 3 /*dim*/, cloud, { 15 /* max leaf */ } );
    KdTree::BoundingBox bbox_src;
    tree.computeBoundingBox( bbox_src );
    printf( "Source bounding boxes: (%f, %f), (%f, %f), (%3.10e, %3.10e)\n", bbox_src[0].low, bbox_src[0].high,
            bbox_src[1].low, bbox_src[1].high, bbox_src[2].low, bbox_src[2].high );

    // constexpr double Radius = 1.0E-6;
#pragma omp parallel for shared( tree, dst_xyz, src_tdata, dst_tdata )
    for( size_t i = 0; i < dst_tdata.size(); i++ )
    {
        const size_t offset = i * 3;
        dst_tdata[i]        = 0;

        // constexpr size_t num_results[] = {1, 8, 13, 27};
        // constexpr double num_res_weights[] = { 0.125, 0.33, 0.33, 0.215 };
        // for( size_t ires = 0; ires < 4; ++ires )
        // {
        //     std::vector< size_t > srcindx( num_results[ires] );
        //     std::vector< double > srcdist( num_results[ires] );
        //     nanoflann::KNNResultSet< double > resultSet( num_results[ires] );

        //     const double* query_pt = dst_xyz.data() + offset;

        //     // Do a KNN search
        //     resultSet.init( srcindx.data(), srcdist.data() );
        //     tree.findNeighbors( resultSet, query_pt );

        //     // check if the point is outside the bounding box
        //     if( ( query_pt[0] < bbox_src[0].low || query_pt[1] < bbox_src[1].low || query_pt[2] < bbox_src[2].low ) ||
        //         ( query_pt[0] > bbox_src[0].high || query_pt[1] > bbox_src[1].high || query_pt[2] > bbox_src[2].high ) )
        //     {
        //         // data needs to be extrapolated
        //         if( query_pt[2] < bbox_src[2].low )  // point is below the MPAS sea bed
        //             dst_tdata[i] += num_res_weights[ires] * src_tdata[srcindx[0]];
        //         else  // point is above the MPAS sea surface
        //             dst_tdata[i] += num_res_weights[ires] * src_tdata[srcindx[0]];
        //     }
        //     else
        //     {
        //         double value = 0.0, weights = 0.0;
        //         for( size_t j = 0; j < num_results[ires]; ++j )
        //         {
        //             // const double idw = 1.0 / std::pow( srcdist[j], power );
        //             const double idw = std::pow( fmax( 0.0, Radius - srcdist[j] ) / Radius / srcdist[j], power );
        //             value += src_tdata[srcindx[j]] * idw;
        //             weights += idw;
        //         }

        //         // dst_tdata[i] = src_tdata[srcindx[0]];
        //         if( weights > Radius ) dst_tdata[i] += num_res_weights[ires] * value / weights;
        //         else
        //             dst_tdata[i] += num_res_weights[ires] * src_tdata[srcindx[0]];
        //     }
        // }
        // constexpr int numlevels = 4;
        // constexpr size_t num_results[]     = { 1, 8, 13, 27 };
        // constexpr double num_res_weights[] = { 1.0/49, 8.0/49, 13.0/49, 27.0/49 };
        // constexpr int numlevels            = 2;
        // constexpr size_t num_results[]     = { 1, 8 };
        // constexpr double num_res_weights[] = { 1.0 / 9, 8.0 / 9 };
        constexpr int numlevels            = 1;
        constexpr size_t num_results[]     = { 8 };
        constexpr double num_res_weights[] = { 1.0 };
        {
            std::vector< size_t > srcindx( num_results[numlevels - 1] );
            std::vector< double > srcdist( num_results[numlevels - 1] );
            nanoflann::KNNResultSet< double > resultSet( num_results[numlevels - 1] );

            const double* query_pt = dst_xyz.data() + offset;

            // Do a KNN search
            resultSet.init( srcindx.data(), srcdist.data() );
            tree.findNeighbors( resultSet, query_pt );

            // check if the point is outside the bounding box
            if( ( query_pt[0] < bbox_src[0].low || query_pt[1] < bbox_src[1].low || query_pt[2] < bbox_src[2].low ) ||
                ( query_pt[0] > bbox_src[0].high || query_pt[1] > bbox_src[1].high || query_pt[2] > bbox_src[2].high ) )
            {
                // data needs to be extrapolated
                if( query_pt[2] < bbox_src[2].low )  // point is below the MPAS sea bed
                    dst_tdata[i] = src_tdata[srcindx[0]];
                else  // point is above the MPAS sea surface
                    dst_tdata[i] = src_tdata[srcindx[0]];
            }
            else
            {
                // constexpr double alpha = 0.666;//2.0 / ( num_results[3] + 1 );
                // double value = 0.0, weights = 1.0;
                // for( size_t j = 0; j < num_results[3]; ++j )
                // // for( int j = num_results[3] - 1; j >= 0; j-- )
                // {
                //     // const double idw = 1.0 / std::pow( srcdist[j], power );
                //     // const double idw = std::pow( fmax( 0.0, Radius - srcdist[j] ) / Radius / srcdist[j], power );

                //     value = alpha * src_tdata[srcindx[j]] + (1.0 - alpha) * value;
                // }
                // dst_tdata[i] = value;

                double finalvalue = 0.0;
                for( size_t k = 0; k < numlevels; ++k )
                {
                    double value = 0.0, weights = 0.0;
                    for( size_t j = 0; j < num_results[k]; ++j )
                    {
                        const double idw = std::max( 1e-12, std::pow( srcdist[j], power ) );
                        value += src_tdata[srcindx[j]] / idw;
                        weights += 1.0 / idw;
                    }
                    finalvalue += num_res_weights[k] * value / weights;
                }

                // dst_tdata[i] = src_tdata[srcindx[0]];
                dst_tdata[i] = finalvalue;
            }
        }
    }

    return moab::MB_SUCCESS;
}

#endif  // __compute_shepard_hpp__
