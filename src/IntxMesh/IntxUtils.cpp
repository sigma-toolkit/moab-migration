/*
 * IntxUtils.cpp
 *
 *  Created on: Oct 3, 2012
 */
#if defined( _MSC_VER ) || defined( WIN32 ) /* windows */
#define _USE_MATH_DEFINES                   // For M_PI
#endif

#include <cmath>
#include <cassert>
#include <iostream>

#include "moab/IntxMesh/IntxUtils.hpp"
// this is from mbcoupler; maybe it should be moved somewhere in moab src
// right now, add a dependency to mbcoupler
// #include "ElemUtil.hpp"
#include "moab/MergeMesh.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"

#define CHECKNEGATIVEAREA
#ifdef CHECKNEGATIVEAREA
#include <iomanip>
#endif
#include <queue>
#include <map>

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "GridElements.h"
#endif

#ifdef MOAB_HAVE_EIGEN3
#define EIGEN_NO_DEBUG
#define EIGEN_MAX_CPP_VER 11
#include "Eigen/Dense"
#endif

namespace moab
{
/**
 * This code defines several utility functions for computing edge intersections and performing geometric operations.
 *
 * - `borderPointsOfXinY2`: Computes the border points of a set of points `X` inside another set of points `Y`.
 * - `SortAndRemoveDoubles2`: Sorts a set of points `P` according to their angles and removes duplicate points.
 * - `EdgeIntersections2`: Computes the intersections between the edges of two sets of points `blue` and `red`.
 * - `EdgeIntxRllCs`: Computes the intersections between the edges of a set of points `blue` and a set of points `red` on a specific plane.
 *
 * The code also defines some helper structs and functions used by these utility functions.
 */

#define CORRTAGNAME "__correspondent"
#define MAXEDGES    10

/**
 * Computes the border points of X in Y2.
 *
 * @param X The array of points representing X.
 * @param nX The number of points in X.
 * @param Y The array of points representing Y.
 * @param nY The number of points in Y.
 * @param P The array to store the border points of X in Y2.
 * @param side The array to store the side information for each point in X.
 * @param epsilon_area The epsilon value for area comparison.
 * @return The number of extra points found.
 */
int IntxUtils::borderPointsOfXinY2( double* X, int nX, double* Y, int nY, double* P, int* side, double epsilon_area )
{
    // 2 triangles, 3 corners, is the corner of X in Y?
    // Y must have a positive area
    /*
     */
    int extraPoint = 0;
    for( int i = 0; i < nX; i++ )
    {
        // compute double the area of all nY triangles formed by a side of Y and a corner of X; if
        // one is negative, stop (negative means it is outside; X and Y are all oriented such that
        // they are positive oriented;
        //  if one area is negative, it means it is outside the convex region, for sure)
        double* A = X + 2 * i;

        int inside = 1;
        for( int j = 0; j < nY; j++ )
        {
            double* B = Y + 2 * j;
            int j1    = ( j + 1 ) % nY;
            double* C = Y + 2 * j1;  // no copy of data

            double area2 = ( B[0] - A[0] ) * ( C[1] - A[1] ) - ( C[0] - A[0] ) * ( B[1] - A[1] );
            if( area2 < -epsilon_area )
            {
                inside = 0;
                break;
            }
        }
        if( inside )
        {
            side[i] = 1;  // so vertex i of X is inside the convex region formed by Y
            // so side has nX dimension (first array)
            P[extraPoint * 2]     = A[0];
            P[extraPoint * 2 + 1] = A[1];
            extraPoint++;
        }
    }
    return extraPoint;
}

// used to order according to angle, so it can remove double easily
struct angleAndIndex
{
    double angle;
    int index;
};

bool angleCompare( angleAndIndex lhs, angleAndIndex rhs )
{
    return lhs.angle < rhs.angle;
}

/**
 * Sorts and removes duplicate points in the given array.
 *
 * Note: nP might be modified too, we will remove duplicates if found
 *
 * @param P The array of points to be sorted and checked for duplicates.
 * @param nP The number of points in P.
 * @param epsilon_1 The epsilon value for distance comparison.
 * @return 0 if successful.
 */
int IntxUtils::SortAndRemoveDoubles2( double* P, int& nP, double epsilon_1 )
{
    if( nP < 2 ) return 0;  // nothing to do

    // center of gravity for the points
    double c[2] = { 0., 0. };
    int k       = 0;
    for( k = 0; k < nP; k++ )
    {
        c[0] += P[2 * k];
        c[1] += P[2 * k + 1];
    }
    c[0] /= nP;
    c[1] /= nP;

    // how many? we dimensioned P at MAXEDGES*10; so we imply we could have at most 5*MAXEDGES
    // intersection points
    struct angleAndIndex pairAngleIndex[5 * MAXEDGES];

    for( k = 0; k < nP; k++ )
    {
        double x = P[2 * k] - c[0], y = P[2 * k + 1] - c[1];
        if( x != 0. || y != 0. )
        {
            pairAngleIndex[k].angle = atan2( y, x );
        }
        else
        {
            pairAngleIndex[k].angle = 0;
            // it would mean that the cells are touching at a vertex
        }
        pairAngleIndex[k].index = k;
    }

    // this should be faster than the bubble sort we had before
    std::sort( pairAngleIndex, pairAngleIndex + nP, angleCompare );
    // copy now to a new double array
    double PCopy[10 * MAXEDGES];  // the same dimension as P; very conservative, but faster than
                                  // reallocate for a vector
    for( k = 0; k < nP; k++ )     // index will show where it should go now;
    {
        int ck           = pairAngleIndex[k].index;
        PCopy[2 * k]     = P[2 * ck];
        PCopy[2 * k + 1] = P[2 * ck + 1];
    }
    // now copy from PCopy over original P location
    std::copy( PCopy, PCopy + 2 * nP, P );

    // eliminate duplicates, finally

    int i = 0, j = 1;  // the next one; j may advance faster than i
    // check the unit
    // double epsilon_1 = 1.e-5; // these are cm; 2 points are the same if the distance is less
    // than 1.e-5 cm
    while( j < nP )
    {
        double d2 = dist2( &P[2 * i], &P[2 * j] );
        if( d2 > epsilon_1 )
        {
            i++;
            P[2 * i]     = P[2 * j];
            P[2 * i + 1] = P[2 * j + 1];
        }
        j++;
    }
    // test also the last point with the first one (index 0)
    // the first one could be at -PI; last one could be at +PI, according to atan2 span

    double d2 = dist2( P, &P[2 * i] );  // check the first and last points (ordered from -pi to +pi)
    if( d2 > epsilon_1 )
    {
        nP = i + 1;
    }
    else
        nP = i;            // effectively delete the last point (that would have been the same with first)
    if( nP == 0 ) nP = 1;  // we should be left with at least one point we already tested if nP is 0 originally
    return 0;
}

/**
 * Computes the edge intersections of two elements.
 *
 * @param blue The array of points representing the blue element.
 * @param nsBlue The number of points in the blue element.
 * @param red The array of points representing the red element.
 * @param nsRed The number of points in the red element.
 * @param markb The array to mark the intersecting edges of the blue element.
 * @param markr The array to mark the intersecting edges of the red element.
 * @param points The array to store the intersection points.
 * @param nPoints The number of intersection points found.
 * @return The error code.
 */
ErrorCode IntxUtils::EdgeIntersections2( double* blue,
                                         int nsBlue,
                                         double* red,
                                         int nsRed,
                                         int* markb,
                                         int* markr,
                                         double* points,
                                         int& nPoints )
{
    /* EDGEINTERSECTIONS computes edge intersections of two elements
     [P,n]=EdgeIntersections(X,Y) computes for the two given elements  * red
     and blue ( stored column wise )
     (point coordinates are stored column-wise, in counter clock
     order) the points P where their edges intersect. In addition,
     in n the indices of which neighbors of red  are also intersecting
     with blue are given.
     */

    // points is an array with enough slots   (24 * 2 doubles)
    nPoints = 0;
    for( int i = 0; i < MAXEDGES; i++ )
    {
        markb[i] = markr[i] = 0;
    }

    for( int i = 0; i < nsBlue; i++ )
    {
        for( int j = 0; j < nsRed; j++ )
        {
            double b[2];
            double a[2][2];  // 2*2
            int iPlus1 = ( i + 1 ) % nsBlue;
            int jPlus1 = ( j + 1 ) % nsRed;
            for( int k = 0; k < 2; k++ )
            {
                b[k] = red[2 * j + k] - blue[2 * i + k];
                // row k of a: a(k, 0), a(k, 1)
                a[k][0] = blue[2 * iPlus1 + k] - blue[2 * i + k];
                a[k][1] = red[2 * j + k] - red[2 * jPlus1 + k];
            }
            double delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
            if( fabs( delta ) > 1.e-14 )  // this is close to machine epsilon
            {
                // not parallel
                double alfa = ( b[0] * a[1][1] - a[0][1] * b[1] ) / delta;
                double beta = ( -b[0] * a[1][0] + b[1] * a[0][0] ) / delta;
                if( 0 <= alfa && alfa <= 1. && 0 <= beta && beta <= 1. )
                {
                    // the intersection is good
                    for( int k = 0; k < 2; k++ )
                    {
                        points[2 * nPoints + k] = blue[2 * i + k] + alfa * ( blue[2 * iPlus1 + k] - blue[2 * i + k] );
                    }
                    markb[i] = 1;  // so neighbor number i of blue will be considered too.
                    markr[j] = 1;  // this will be used in advancing red around blue quad
                    nPoints++;
                }
            }
            // the case delta ~ 0. will be considered by the interior points logic
        }
    }
    return MB_SUCCESS;
}

/**
 * Computes the edge intersections between a RLL and CS quad.
 *
 * Note: Special function.
 *
 * @param blue The array of points representing the blue element.
 * @param bluec The array of Cartesian coordinates for the blue element.
 * @param blueEdgeType The array of edge types for the blue element.
 * @param nsBlue The number of points in the blue element.
 * @param red The array of points representing the red element.
 * @param redc The array of Cartesian coordinates for the red element.
 * @param nsRed The number of points in the red element.
 * @param markb The array to mark the intersecting edges of the blue element.
 * @param markr The array to mark the intersecting edges of the red element.
 * @param plane The plane of intersection.
 * @param R The radius of the sphere.
 * @param points The array to store the intersection points.
 * @param nPoints The number of intersection points found.
 * @return The error code.
 */
ErrorCode IntxUtils::EdgeIntxRllCs( double* blue,
                                    CartVect* bluec,
                                    int* blueEdgeType,
                                    int nsBlue,
                                    double* red,
                                    CartVect* redc,
                                    int nsRed,
                                    int* markb,
                                    int* markr,
                                    int plane,
                                    double R,
                                    double* points,
                                    int& nPoints )
{
    // if blue edge type is 1, intersect in 3d then project to 2d by gnomonic projection
    // everything else the same (except if there are 2 points resulting, which is rare)
    for( int i = 0; i < 4; i++ )
    {  // always at most 4 , so maybe don't bother
        markb[i] = markr[i] = 0;
    }

    for( int i = 0; i < nsBlue; i++ )
    {
        int iPlus1 = ( i + 1 ) % nsBlue;
        if( blueEdgeType[i] == 0 )  // old style, just 2d
        {
            for( int j = 0; j < nsRed; j++ )
            {
                double b[2];
                double a[2][2];  // 2*2

                int jPlus1 = ( j + 1 ) % nsRed;
                for( int k = 0; k < 2; k++ )
                {
                    b[k] = red[2 * j + k] - blue[2 * i + k];
                    // row k of a: a(k, 0), a(k, 1)
                    a[k][0] = blue[2 * iPlus1 + k] - blue[2 * i + k];
                    a[k][1] = red[2 * j + k] - red[2 * jPlus1 + k];
                }
                double delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
                if( fabs( delta ) > 1.e-14 )
                {
                    // not parallel
                    double alfa = ( b[0] * a[1][1] - a[0][1] * b[1] ) / delta;
                    double beta = ( -b[0] * a[1][0] + b[1] * a[0][0] ) / delta;
                    if( 0 <= alfa && alfa <= 1. && 0 <= beta && beta <= 1. )
                    {
                        // the intersection is good
                        for( int k = 0; k < 2; k++ )
                        {
                            points[2 * nPoints + k] =
                                blue[2 * i + k] + alfa * ( blue[2 * iPlus1 + k] - blue[2 * i + k] );
                        }
                        markb[i] = 1;  // so neighbor number i of blue will be considered too.
                        markr[j] = 1;  // this will be used in advancing red around blue quad
                        nPoints++;
                    }
                }  // if the edges are too "parallel", skip them
            }
        }
        else  // edge type is 1, so use 3d intersection
        {
            CartVect& C = bluec[i];
            CartVect& D = bluec[iPlus1];
            for( int j = 0; j < nsRed; j++ )
            {
                int jPlus1  = ( j + 1 ) % nsRed;  // nsRed is just 4, forget about it, usually
                CartVect& A = redc[j];
                CartVect& B = redc[jPlus1];
                int np      = 0;
                double E[9];
                intersect_great_circle_arc_with_clat_arc( A.array(), B.array(), C.array(), D.array(), R, E, np );
                if( np == 0 ) continue;
                if( np >= 2 )
                {
                    std::cout << "intersection with 2 points :" << A << B << C << D << "\n";
                }
                for( int k = 0; k < np; k++ )
                {
                    gnomonic_projection( CartVect( E + k * 3 ), R, plane, points[2 * nPoints],
                                         points[2 * nPoints + 1] );
                    nPoints++;
                }
                markb[i] = 1;  // so neighbor number i of blue will be considered too.
                markr[j] = 1;  // this will be used in advancing red around blue quad
            }
        }
    }
    return MB_SUCCESS;
}

// vec utils related to gnomonic projection on a sphere

// vec utils

/*
 *
 * position on a sphere of radius R
 * if plane specified, use it; if not, return the plane, and the point in the plane
 * there are 6 planes, numbered 1 to 6
 * plane 1: x=R, plane 2: y=R, 3: x=-R, 4: y=-R, 5: z=-R, 6: z=R
 *
 * projection on the plane will preserve the orientation, such that a triangle, quad pointing
 * outside the sphere will have a positive orientation in the projection plane
 */
void IntxUtils::decide_gnomonic_plane( const CartVect& pos, int& plane )
{
    // decide plane, based on max x, y, z
    if( fabs( pos[0] ) < fabs( pos[1] ) )
    {
        if( fabs( pos[2] ) < fabs( pos[1] ) )
        {
            // pos[1] is biggest
            if( pos[1] > 0 )
                plane = 2;
            else
                plane = 4;
        }
        else
        {
            // pos[2] is biggest
            if( pos[2] < 0 )
                plane = 5;
            else
                plane = 6;
        }
    }
    else
    {
        if( fabs( pos[2] ) < fabs( pos[0] ) )
        {
            // pos[0] is the greatest
            if( pos[0] > 0 )
                plane = 1;
            else
                plane = 3;
        }
        else
        {
            // pos[2] is biggest
            if( pos[2] < 0 )
                plane = 5;
            else
                plane = 6;
        }
    }
    return;
}

// point on a sphere is projected on one of six planes, decided earlier
ErrorCode IntxUtils::gnomonic_projection( const CartVect& pos, double R, int plane, double& c1, double& c2 )
{
    double alfa = 1.;  // the new point will be on line alfa*pos

    switch( plane )
    {
        case 1: {
            // the plane with x = R; x>y, x>z
            // c1->y, c2->z
            alfa = R / pos[0];
            c1   = alfa * pos[1];
            c2   = alfa * pos[2];
            break;
        }
        case 2: {
            // y = R -> zx
            alfa = R / pos[1];
            c1   = alfa * pos[2];
            c2   = alfa * pos[0];
            break;
        }
        case 3: {
            // x=-R, -> yz
            alfa = -R / pos[0];
            c1   = -alfa * pos[1];  // the sign is to preserve orientation
            c2   = alfa * pos[2];
            break;
        }
        case 4: {
            // y = -R
            alfa = -R / pos[1];
            c1   = -alfa * pos[2];  // the sign is to preserve orientation
            c2   = alfa * pos[0];
            break;
        }
        case 5: {
            // z = -R
            alfa = -R / pos[2];
            c1   = -alfa * pos[0];  // the sign is to preserve orientation
            c2   = alfa * pos[1];
            break;
        }
        case 6: {
            alfa = R / pos[2];
            c1   = alfa * pos[0];
            c2   = alfa * pos[1];
            break;
        }
        default:
            return MB_FAILURE;  // error
    }

    return MB_SUCCESS;  // no error
}

// given the position on plane (one out of 6), find out the position on sphere
ErrorCode IntxUtils::reverse_gnomonic_projection( const double& c1,
                                                  const double& c2,
                                                  double R,
                                                  int plane,
                                                  CartVect& pos )
{

    // the new point will be on line beta*pos
    double len  = sqrt( c1 * c1 + c2 * c2 + R * R );
    double beta = R / len;  // it is less than 1, in general

    switch( plane )
    {
        case 1: {
            // the plane with x = R; x>y, x>z
            // c1->y, c2->z
            pos[0] = beta * R;
            pos[1] = c1 * beta;
            pos[2] = c2 * beta;
            break;
        }
        case 2: {
            // y = R -> zx
            pos[1] = R * beta;
            pos[2] = c1 * beta;
            pos[0] = c2 * beta;
            break;
        }
        case 3: {
            // x=-R, -> yz
            pos[0] = -R * beta;
            pos[1] = -c1 * beta;  // the sign is to preserve orientation
            pos[2] = c2 * beta;
            break;
        }
        case 4: {
            // y = -R
            pos[1] = -R * beta;
            pos[2] = -c1 * beta;  // the sign is to preserve orientation
            pos[0] = c2 * beta;
            break;
        }
        case 5: {
            // z = -R
            pos[2] = -R * beta;
            pos[0] = -c1 * beta;  // the sign is to preserve orientation
            pos[1] = c2 * beta;
            break;
        }
        case 6: {
            pos[2] = R * beta;
            pos[0] = c1 * beta;
            pos[1] = c2 * beta;
            break;
        }
        default:
            return MB_FAILURE;  // error
    }

    return MB_SUCCESS;  // no error
}

void IntxUtils::gnomonic_unroll( double& c1, double& c2, double R, int plane )
{
    double tmp;
    switch( plane )
    {
        case 1:
            break;  // nothing
        case 2:     // rotate + 90
            tmp = c1;
            c1  = -c2;
            c2  = tmp;
            c1  = c1 + 2 * R;
            break;
        case 3:
            c1 = c1 + 4 * R;
            break;
        case 4:  // rotate with -90 x-> -y; y -> x

            tmp = c1;
            c1  = c2;
            c2  = -tmp;
            c1  = c1 - 2 * R;
            break;
        case 5:  // South Pole
            // rotate 180 then move to (-2, -2)
            c1 = -c1 - 2. * R;
            c2 = -c2 - 2. * R;
            break;
            ;
        case 6:  // North Pole
            c1 = c1 - 2. * R;
            c2 = c2 + 2. * R;
            break;
    }
    return;
}
// given a mesh on the sphere, project all centers in 6 gnomonic planes, or project mesh too
ErrorCode IntxUtils::global_gnomonic_projection( Interface* mb,
                                                 EntityHandle inSet,
                                                 double R,
                                                 bool centers_only,
                                                 EntityHandle& outSet )
{
    std::string parTagName( "PARALLEL_PARTITION" );
    Tag part_tag;
    Range partSets;
    ErrorCode rval = mb->tag_get_handle( parTagName.c_str(), part_tag );
    if( MB_SUCCESS == rval && part_tag != 0 )
    {
        rval = mb->get_entities_by_type_and_tag( inSet, MBENTITYSET, &part_tag, NULL, 1, partSets, Interface::UNION );MB_CHK_ERR( rval );
    }
    rval = ScaleToRadius( mb, inSet, 1.0 );MB_CHK_ERR( rval );
    // Get all entities of dimension 2
    Range inputRange;  // get
    rval = mb->get_entities_by_dimension( inSet, 1, inputRange );MB_CHK_ERR( rval );
    rval = mb->get_entities_by_dimension( inSet, 2, inputRange );MB_CHK_ERR( rval );

    std::map< EntityHandle, int > partsAssign;
    std::map< int, EntityHandle > newPartSets;
    if( !partSets.empty() )
    {
        // get all cells, and assign parts
        for( Range::iterator setIt = partSets.begin(); setIt != partSets.end(); ++setIt )
        {
            EntityHandle pSet = *setIt;
            Range ents;
            rval = mb->get_entities_by_handle( pSet, ents );MB_CHK_ERR( rval );
            int val;
            rval = mb->tag_get_data( part_tag, &pSet, 1, &val );MB_CHK_ERR( rval );
            // create a new set with the same part id tag, in the outSet
            EntityHandle newPartSet;
            rval = mb->create_meshset( MESHSET_SET, newPartSet );MB_CHK_ERR( rval );
            rval = mb->tag_set_data( part_tag, &newPartSet, 1, &val );MB_CHK_ERR( rval );
            newPartSets[val] = newPartSet;
            rval             = mb->add_entities( outSet, &newPartSet, 1 );MB_CHK_ERR( rval );
            for( Range::iterator it = ents.begin(); it != ents.end(); ++it )
            {
                partsAssign[*it] = val;
            }
        }
    }

    if( centers_only )
    {
        for( Range::iterator it = inputRange.begin(); it != inputRange.end(); ++it )
        {
            CartVect center;
            EntityHandle cell = *it;
            rval              = mb->get_coords( &cell, 1, center.array() );MB_CHK_ERR( rval );
            int plane;
            decide_gnomonic_plane( center, plane );
            double c[3];
            c[2] = 0.;
            gnomonic_projection( center, R, plane, c[0], c[1] );

            gnomonic_unroll( c[0], c[1], R, plane );

            EntityHandle vertex;
            rval = mb->create_vertex( c, vertex );MB_CHK_ERR( rval );
            rval = mb->add_entities( outSet, &vertex, 1 );MB_CHK_ERR( rval );
        }
    }
    else
    {
        // distribute the cells to 6 planes, based on the center
        Range subranges[6];
        for( Range::iterator it = inputRange.begin(); it != inputRange.end(); ++it )
        {
            CartVect center;
            EntityHandle cell = *it;
            rval              = mb->get_coords( &cell, 1, center.array() );MB_CHK_ERR( rval );
            int plane;
            decide_gnomonic_plane( center, plane );
            subranges[plane - 1].insert( cell );
        }
        for( int i = 1; i <= 6; i++ )
        {
            Range verts;
            rval = mb->get_connectivity( subranges[i - 1], verts );MB_CHK_ERR( rval );
            std::map< EntityHandle, EntityHandle > corr;
            for( Range::iterator vt = verts.begin(); vt != verts.end(); ++vt )
            {
                CartVect vect;
                EntityHandle v = *vt;
                rval           = mb->get_coords( &v, 1, vect.array() );MB_CHK_ERR( rval );
                double c[3];
                c[2] = 0.;
                gnomonic_projection( vect, R, i, c[0], c[1] );
                gnomonic_unroll( c[0], c[1], R, i );
                EntityHandle vertex;
                rval = mb->create_vertex( c, vertex );MB_CHK_ERR( rval );
                corr[v] = vertex;  // for new connectivity
            }
            EntityHandle new_conn[20];  // max edges in 2d ?
            for( Range::iterator eit = subranges[i - 1].begin(); eit != subranges[i - 1].end(); ++eit )
            {
                EntityHandle eh          = *eit;
                const EntityHandle* conn = NULL;
                int num_nodes;
                rval = mb->get_connectivity( eh, conn, num_nodes );MB_CHK_ERR( rval );
                // build a new vertex array
                for( int j = 0; j < num_nodes; j++ )
                    new_conn[j] = corr[conn[j]];
                EntityType type = mb->type_from_handle( eh );
                EntityHandle newCell;
                rval = mb->create_element( type, new_conn, num_nodes, newCell );MB_CHK_ERR( rval );
                rval = mb->add_entities( outSet, &newCell, 1 );MB_CHK_ERR( rval );
                std::map< EntityHandle, int >::iterator mit = partsAssign.find( eh );
                if( mit != partsAssign.end() )
                {
                    int val = mit->second;
                    rval    = mb->add_entities( newPartSets[val], &newCell, 1 );MB_CHK_ERR( rval );
                }
            }
        }
    }

    return MB_SUCCESS;
}

ErrorCode IntxUtils::EdgeMap(Interface* mb, EntityHandle inputSet, EntityHandle intx_set, bool sourceMap)
{
    Tag parentTag;
    ErrorCode rval;
    if (sourceMap)
    {
        rval = mb->tag_get_handle("SourceParent", parentTag);MB_CHK_SET_ERR( rval, "can't get parent source tag in edge map" );
    }
    else
    {
        rval = mb->tag_get_handle("TargetParent", parentTag);MB_CHK_SET_ERR( rval, "can't get parent target tag in edge map" );
    }
    Tag fractionTag;
    rval = mb->tag_get_handle( "EdgeRecoveryFraction", fractionTag);
    Tag subTag;
    rval = mb->tag_get_handle( "NumSubEdges", subTag);
    // get all polygons in the intx set
    Range cells;
    rval = mb->get_entities_by_dimension(intx_set, 2, cells);MB_CHK_SET_ERR( rval, "can't get intersection cells" );
    // create all edges adjacent to the cells in the intx set
    // some might be original edges from edge or target meshes
    Range intxEdges;
    rval = mb->get_adjacencies(cells, 1, true, intxEdges, Interface::UNION);MB_CHK_SET_ERR( rval, "can't get intersection edges" );
    std::cout << " number of intx edges:" << intxEdges.size() << "\n";
    Range parentCells;
    rval = mb->get_entities_by_dimension(inputSet, 2, parentCells);MB_CHK_SET_ERR( rval, "can't get intersection cells" );
    Range initialEdges;
    rval = mb->get_adjacencies(parentCells, 1, true, initialEdges, Interface::UNION);MB_CHK_SET_ERR( rval, "can't get intersection edges" );
    std::cout << " number of original input edges:" << initialEdges.size() << "\n";

    // first, identify input polygons that are recovered fully
    // get global ids of the initial cells
    std::vector<int> parentGids(parentCells.size());
    Tag gidTag = mb->globalId_tag();
    rval = mb->tag_get_data(gidTag, parentCells, &parentGids[0]);MB_CHK_SET_ERR( rval, "can't get parent global ids" );
    std::map< int, double > initAreas; // get areas of those initial cells
    int i=0;
    std::vector< double > coords( 3 * 20 );// enough vertices, at most 30, good for intx polygons too
    int num_nodes = 0; // num nodes in cells
    const EntityHandle* verts;
    IntxAreaUtils areaAdaptor( IntxAreaUtils::lHuiller ); // GaussQuadrature , lHuiller
    std::map< int, double > recoveredAreas; // get areas of from intersection cells
    std::map<int, EntityHandle> mapFromGIDToParent;
    std::map<int, std::vector<EntityHandle>> mapFromParentGIDToIntxCells;
    for (Range::iterator it=parentCells.begin(); it!=parentCells.end(); it++, i++)
    {
        EntityHandle parentCell=*it;
        rval = mb->get_connectivity( parentCell, verts, num_nodes );MB_CHK_SET_ERR( rval, "can't get connectivity of parent cell" );
        // get coordinates
        rval = mb->get_coords( verts, num_nodes, &coords[0] );MB_CHK_SET_ERR( rval, "can't get coords of parent cell" );

        double area = areaAdaptor.area_spherical_polygon( &coords[0], num_nodes, 1. );
        int parentID = parentGids[i];
        initAreas[parentID] = area;
        recoveredAreas[parentID] = 0.;
        mapFromGIDToParent[parentID] = parentCell;
        mapFromParentGIDToIntxCells[parentID]; // just initialize it with empty vector
    }


    for (Range::iterator it=cells.begin(); it!=cells.end(); it++)
    {
        EntityHandle cell=*it;
        rval = mb->get_connectivity( cell, verts, num_nodes );MB_CHK_SET_ERR( rval, "can't get connectivity of intx cell" );
        // get coordinates
        rval = mb->get_coords( verts, num_nodes, &coords[0] );MB_CHK_SET_ERR( rval, "can't get coods of intx cell" );
        double intx_area = areaAdaptor.area_spherical_polygon( &coords[0], num_nodes, 1. );
        int parentID = 0;
        rval = mb->tag_get_data(parentTag, &cell, 1, &parentID);MB_CHK_SET_ERR( rval, "can't get parent Tag" );
        recoveredAreas[parentID] += intx_area;
        mapFromParentGIDToIntxCells[parentID].push_back(cell);
    }
    int recovered = 0, notRecovered =0;
    Range recoveredCells;
    for (size_t j=0; j<parentGids.size(); j++)
    {
        int parentID = parentGids[j];
        double areaDiff = fabs(initAreas[parentID] - recoveredAreas[parentID]);
        if (areaDiff < 1.e-12)
        {
            recovered++;
            recoveredCells.insert(parentCells[j]); // should we use a std::vector, that will be ordered already ?
        }
        else
            notRecovered++;
    }
    std::cout << "recovered initial cells: " << recovered << " vs:" << notRecovered << " not recovered \n";
    // initial edges that should be decomposable from intx edges
    Range recoverableEdges;
    rval = mb->get_adjacencies(recoveredCells, 1, false, recoverableEdges,  Interface::UNION);MB_CHK_SET_ERR( rval, "can't get recoverable edges" );

    // now recover each edge, looking at intx polygons forming one of the adjacent cells, that is in initial recoverable set
    int recoveredEdges = 0;
    int unrecovered = 0;
    int identity_edges = 0; // edges that are formed by one intx edge, itself, the original
    std::map<EntityHandle, std::vector<EntityHandle>> mapEdges;
    for (Range::iterator eit=recoverableEdges.begin(); eit!=recoverableEdges.end(); eit++)
    {
        EntityHandle initialEdge = *eit;
        // if this edge is among intxEdges, we are done
        int index = intxEdges.index(initialEdge);
        if (index>=0)
        {
            mapEdges[initialEdge].push_back(initialEdge);
            identity_edges++;
            recoveredEdges++;
            continue; // no need to sweat it anymore
        }
        // get adjacent polygons; if there is an adj polygon in intx cells, we are done; if not, get one in the initial recoveredCells range
        Range initialAdjCells;
        rval = mb-> get_adjacencies(&initialEdge, 1, 2, false, initialAdjCells, Interface::UNION);MB_CHK_SET_ERR( rval, "can't get adjacent initial cells" );
        if (initialAdjCells.size()==0) MB_CHK_SET_ERR( MB_FAILURE, "no adjacent initial cells" );
        // intersect with recoveredCells
        Range adjRecoveredInitialCell = intersect(initialAdjCells, recoveredCells);
        if(adjRecoveredInitialCell.size() == 0) MB_CHK_SET_ERR( MB_FAILURE, "no adjacent initial cells that are recovered" );
        // get all edges adjacent to intersection polygons that form one of the recovered initial cell
        EntityHandle recoveredCell = adjRecoveredInitialCell[0]; // first one
        int nv = 0;
        const EntityHandle * connCell;
        rval = mb->get_connectivity(recoveredCell, connCell, nv);MB_CHK_SET_ERR( rval, "can't get connectivity of parent cell" );
        int cellGlobalId = 0;
        rval = mb->tag_get_data(gidTag, &recoveredCell, 1, &cellGlobalId);MB_CHK_SET_ERR( rval, "can't get id of parent cell" );
        // list of intx polygons in it:
        std::vector<EntityHandle>  & listIntxCells = mapFromParentGIDToIntxCells[cellGlobalId];
        std::vector<EntityHandle> candidateEdges;
        rval = mb->get_adjacencies(&listIntxCells[0], listIntxCells.size(), 1, false, candidateEdges,  Interface::UNION);MB_CHK_SET_ERR( rval, "can't get adj edges" );
        // add all edges that have both vertices on the original edge
        CartVect verticesOriginal[2];
        const EntityHandle *connEdge;
        int numVerts = 0;
        rval = mb->get_connectivity(initialEdge, connEdge, numVerts);MB_CHK_SET_ERR( rval, "can't get connectivity of initial edge" );
        // need to find out if edge oriented positive
        bool found = false;
        bool reverted = false;
        for (int k=0; k<nv && !found; k++)
        {
            EntityHandle v1=connCell[k];
            EntityHandle v2=connCell[ (k+1)%nv ];
            if (connEdge[0] == v1 && connEdge[1] == v2)
            {
                found = true; // positive
            }
            else if(connEdge[0] == v2 && connEdge[1] == v1)
            {
                found = true; // positive
                reverted = true;
            }
        }
        if (!found) MB_CHK_SET_ERR( MB_FAILURE, "can't get orientation of adjacent edge" );
        rval = mb->get_coords( connEdge, numVerts,  verticesOriginal[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of initial edge" );
        if (reverted)
        {
            // change vertices cart coord , so we will always have positive areas when computing
            CartVect tmp  = verticesOriginal[0];
            verticesOriginal[0] = verticesOriginal[1];
            verticesOriginal[1] = tmp;
        }
        double edgeLength = angle( verticesOriginal[0], verticesOriginal[1] ); // distance on sphere of radius 1 is the angle
        // loop now over candidateEdges
        double recoveredLength = 0.;
        for (size_t i=0; i< candidateEdges.size(); i++)
        {
            EntityHandle subedge = candidateEdges[i];
            // get its vertex coordinates:
            rval = mb->get_connectivity(subedge, connEdge, numVerts);MB_CHK_SET_ERR( rval, "can't get connectivity of candidate edge" );
            CartVect verticesSubEdge[2];
            rval = mb->get_coords( connEdge, numVerts,  verticesSubEdge[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of initial edge" );
            // decide if they are on the initial edge (form an angle)
            bool onEdge = true;
            for (int j=0; j<2 && onEdge; j++)
            {
                double dist0 = angle_robust(verticesOriginal[0], verticesSubEdge[j]);
                double dist1 = angle_robust(verticesOriginal[1], verticesSubEdge[j]);
                if (dist0 < 1.e-11 || dist1 < 1.e-11) continue; // end  points
                if (dist1+dist0 - edgeLength > 1.e-11) // triangle inequality
                {
                    onEdge = false;
                    continue;
                }

                double areaTriangle = areaAdaptor.area_spherical_triangle( verticesOriginal[0].array(), verticesOriginal[1].array(), verticesSubEdge[j].array() , 1.0 );
                //double angle = IntxUtils::oriented_spherical_angle( verticesOriginal[0].array(), verticesOriginal[1].array(), verticesSubEdge[j].array() );
                if ( fabs(areaTriangle) > 1.e-12 ) onEdge = false;
            }
            if (onEdge)
            {
                double subEdgeLen = angle_robust(verticesSubEdge[0], verticesSubEdge[1]);
                recoveredLength += subEdgeLen;
                mapEdges[initialEdge].push_back(subedge);
            }
        }
        if ( fabs(edgeLength-recoveredLength) < 1.e-10)
            recoveredEdges++;
        else
            unrecovered++;
        double fraction = recoveredLength/edgeLength;
        rval = mb->tag_set_data(fractionTag, &initialEdge, 1, &fraction);MB_CHK_SET_ERR( rval, "can't set fraction on initial edge" );
        fraction = double(mapEdges[initialEdge].size());
        rval = mb->tag_set_data(subTag, &initialEdge, 1, &fraction);MB_CHK_SET_ERR( rval, "can't set number of subedges on initial edge" );
    }

    std::cout << " recoveredEdges:"  << recoveredEdges << " identity edges:" << identity_edges <<  " unrecovered edges:" <<
            unrecovered << "\n";
    return MB_SUCCESS;
}

void IntxUtils::transform_coordinates( double* avg_position, int projection_type )
{
    if( projection_type == 1 )
    {
        double R =
            avg_position[0] * avg_position[0] + avg_position[1] * avg_position[1] + avg_position[2] * avg_position[2];
        R               = sqrt( R );
        double lat      = asin( avg_position[2] / R );
        double lon      = atan2( avg_position[1], avg_position[0] );
        avg_position[0] = lon;
        avg_position[1] = lat;
        avg_position[2] = R;
    }
    else if( projection_type == 2 )  // gnomonic projection
    {
        CartVect pos( avg_position );
        int gplane;
        IntxUtils::decide_gnomonic_plane( pos, gplane );

        IntxUtils::gnomonic_projection( pos, 1.0, gplane, avg_position[0], avg_position[1] );
        avg_position[2] = 0;
        IntxUtils::gnomonic_unroll( avg_position[0], avg_position[1], 1.0, gplane );
    }
}
/*
 *
 use physical_constants, only : dd_pi
 type(cartesian3D_t), intent(in) :: cart
 type(spherical_polar_t)         :: sphere

 sphere%r=distance(cart)
 sphere%lat=ASIN(cart%z/sphere%r)
 sphere%lon=0

 ! ==========================================================
 ! enforce three facts:
 !
 ! 1) lon at poles is defined to be zero
 !
 ! 2) Grid points must be separated by about .01 Meter (on earth)
 !    from pole to be considered "not the pole".
 !
 ! 3) range of lon is { 0<= lon < 2*pi }
 !
 ! ==========================================================

 if (distance(cart) >= DIST_THRESHOLD) then
 sphere%lon=ATAN2(cart%y,cart%x)
 if (sphere%lon<0) then
 sphere%lon=sphere%lon+2*DD_PI
 end if
 end if

 end function cart_to_spherical
 */
IntxUtils::SphereCoords IntxUtils::cart_to_spherical( CartVect& cart3d )
{
    SphereCoords res;
    res.R = cart3d.length();
    if( res.R < 0 )
    {
        res.lon = res.lat = 0.;
        return res;
    }
    res.lat = asin( cart3d[2] / res.R );
    res.lon = atan2( cart3d[1], cart3d[0] );
    if( res.lon < 0 ) res.lon += 2 * M_PI;  // M_PI is defined in math.h? it seems to be true, although
    // there are some defines it depends on :(
    // #if defined __USE_BSD || defined __USE_XOPEN ???

    return res;
}

/*
 * ! ===================================================================
 ! spherical_to_cart:
 ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
 ! on unit sphere
 ! ===================================================================

 function spherical_to_cart(sphere) result (cart)

 type(spherical_polar_t), intent(in) :: sphere
 type(cartesian3D_t)                 :: cart

 cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
 cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
 cart%z=sphere%r*SIN(sphere%lat)

 end function spherical_to_cart
 */
CartVect IntxUtils::spherical_to_cart( IntxUtils::SphereCoords& sc )
{
    CartVect res;
    res[0] = sc.R * cos( sc.lat ) * cos( sc.lon );  // x coordinate
    res[1] = sc.R * cos( sc.lat ) * sin( sc.lon );  // y
    res[2] = sc.R * sin( sc.lat );                  // z
    return res;
}

ErrorCode IntxUtils::ScaleToRadius( Interface* mb, EntityHandle set, double R )
{
    Range nodes;
    ErrorCode rval = mb->get_entities_by_type( set, MBVERTEX, nodes, true );  // recursive
    if( rval != moab::MB_SUCCESS ) return rval;

    // one by one, get the node and project it on the sphere, with a radius given
    // the center of the sphere is at 0,0,0
    for( Range::iterator nit = nodes.begin(); nit != nodes.end(); ++nit )
    {
        EntityHandle nd = *nit;
        CartVect pos;
        rval = mb->get_coords( &nd, 1, (double*)&( pos[0] ) );
        if( rval != moab::MB_SUCCESS ) return rval;
        double len = pos.length();
        if( len == 0. ) return MB_FAILURE;
        pos  = R / len * pos;
        rval = mb->set_coords( &nd, 1, (double*)&( pos[0] ) );
        if( rval != moab::MB_SUCCESS ) return rval;
    }
    return MB_SUCCESS;
}

// assume they are one the same sphere
double IntxAreaUtils::spherical_angle( const double* A, const double* B, const double* C, double Radius )
{
    // the angle by definition is between the planes OAB and OBC
    CartVect a( A );
    CartVect b( B );
    CartVect c( C );
    double err1 = a.length_squared() - Radius * Radius;
    if( fabs( err1 ) > 0.0001 )
    {
        std::cout << " error in input " << a << " radius: " << Radius << " error:" << err1 << "\n";
    }
    CartVect normalOAB = a * b;
    CartVect normalOCB = c * b;
    return angle( normalOAB, normalOCB );
}

// could be bigger than M_PI;
// angle at B could be bigger than M_PI, if the orientation is such that ABC points toward the
// interior
double IntxUtils::oriented_spherical_angle( const double* A, const double* B, const double* C )
{
    // assume the same radius, sphere at origin
    CartVect a( A ), b( B ), c( C );
    CartVect normalOAB = a * b;
    CartVect normalOCB = c * b;
    CartVect orient    = ( c - b ) * ( a - b );
    double ang         = angle( normalOAB, normalOCB );  // this is between 0 and M_PI
    if( ang != ang )
    {
        // signal of a nan
        std::cout << a << " " << b << " " << c << "\n";
        std::cout << ang << "\n";
    }
    if( orient % b < 0 ) return ( 2 * M_PI - ang );  // the other angle, supplement

    return ang;
}

double IntxAreaUtils::area_spherical_triangle( const double* A, const double* B, const double* C, double Radius )
{
    switch( m_eAreaMethod )
    {
        case Girard:
            return area_spherical_triangle_girard( A, B, C, Radius );
#ifdef MOAB_HAVE_TEMPESTREMAP
        case GaussQuadrature:
            return area_spherical_triangle_GQ( A, B, C );
#endif
        case lHuiller:
        default:
            return area_spherical_triangle_lHuiller( A, B, C, Radius );
    }
}

double IntxAreaUtils::area_spherical_polygon( const double* A, int N, double Radius, int* sign )
{
    switch( m_eAreaMethod )
    {
        case Girard:
            return area_spherical_polygon_girard( A, N, Radius );
#ifdef MOAB_HAVE_TEMPESTREMAP
        case GaussQuadrature:
            return area_spherical_polygon_GQ( A, N ) * Radius * Radius;  //area_spherical_polygon_GQ normalizes
#endif
        case lHuiller:
        default:
            return area_spherical_polygon_lHuiller( A, N, Radius, sign );
    }
}

double IntxAreaUtils::area_spherical_triangle_girard( const double* A, const double* B, const double* C, double Radius )
{
    double correction = spherical_angle( A, B, C, Radius ) + spherical_angle( B, C, A, Radius ) +
                        spherical_angle( C, A, B, Radius ) - M_PI;
    double area = Radius * Radius * correction;
    // now, is it negative or positive? is it pointing toward the center or outward?
    CartVect a( A ), b( B ), c( C );
    CartVect abc = ( b - a ) * ( c - a );
    if( abc % a > 0 )  // dot product positive, means ABC points out
        return area;
    else
        return -area;
}

double IntxAreaUtils::area_spherical_polygon_girard( const double* A, int N, double Radius )
{
    // this should work for non-convex polygons too
    // assume that the A, A+3, ..., A+3*(N-1) are the coordinates
    //
    if( N <= 2 ) return 0.;
    double sum_angles = 0.;
    for( int i = 0; i < N; i++ )
    {
        int i1 = ( i + 1 ) % N;
        int i2 = ( i + 2 ) % N;
        sum_angles += IntxUtils::oriented_spherical_angle( A + 3 * i, A + 3 * i1, A + 3 * i2 );
    }
    double correction = sum_angles - ( N - 2 ) * M_PI;
    return Radius * Radius * correction;
}

double IntxAreaUtils::area_spherical_polygon_lHuiller( const double* A, int N, double Radius, int* sign )
{
    // This should work for non-convex polygons too
    // In the input vector A, assume that the A, A+3, ..., A+3*(N-1) are the coordinates
    // We also assume that the orientation is positive;
    // If negative orientation, the area will be negative
    if( N <= 2 ) return 0.;

    int lsign   = 1;  // assume positive orientain
    double area = 0.;
    for( int i = 1; i < N - 1; i++ )
    {
        int i1              = i + 1;
        double areaTriangle = area_spherical_triangle_lHuiller( A, A + 3 * i, A + 3 * i1, Radius );
        if( areaTriangle < 0 )
            lsign = -1;  // signal that we have at least one triangle with negative orientation ;
                         // possible nonconvex polygon
        area += areaTriangle;
    }
    if( sign ) *sign = lsign;

    return area;
}

#ifdef MOAB_HAVE_TEMPESTREMAP

double IntxAreaUtils::area_spherical_polygon_GQ( const double* A, int N )
{
    // this should work for non-convex polygons too
    // In the input vector A, assume that the A, A+3, ..., A+3*(N-1) are the coordinates
    // We also assume that the orientation is positive;
    // If negative orientation, the area can be negative
    if( N <= 2 ) return 0.;

    // assume positive orientation
    double area = 0.;
    for( int i = 1; i < N - 1; i++ )
    {
        area += area_spherical_triangle_GQ( A, A + 3 * i, A + 3 * ( i + 1 ) );
    }
    return area;
}

template < typename Derived >
Eigen::Array< typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime > shift(
    const Eigen::ArrayBase< Derived >& array,
    int positions )
{
    Eigen::Array< typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime > result = array;
    if( positions > 0 )
    {
        result.segment( positions, array.size() - positions ) = array.head( array.size() - positions );
        result.head( positions ).setZero();
    }
    else if( positions < 0 )
    {
        result.head( array.size() + positions ) = array.tail( array.size() + positions );
        result.tail( -positions ).setZero();
    }
    return result;
}

double IntxAreaUtils::area_spherical_triangle_GQ( const double* inode1, const double* inode2, const double* inode3 )
{
#if defined( MOAB_HAVE_EIGEN3 )
    typedef Eigen::Map< const Eigen::Vector3d > V3d;
    const V3d node1( inode1 );
    const V3d node2( inode2 );
    const V3d node3( inode3 );
    const int nOrder = 6;

    // If we change the quadrature order, use the call: GaussQuadrature::GetPoints(nOrder, 0.0, 1.0, dG, dW);
    const double dG[6] = { 0.03376524289842397, 0.1693953067668678, 0.3806904069584016,
                           0.6193095930415985,  0.8306046932331322, 0.966234757101576 };
    const double dW[6] = { 0.08566224618958521, 0.1803807865240693, 0.2339569672863455,
                           0.2339569672863455,  0.1803807865240693, 0.08566224618958521 };

    double dFaceArea = 0.0;
    Eigen::Vector3d dF, dF2, dDaF, dDbF, dDaG, dDbG;
    double nodeCross[3];

    // Calculate area at quadrature node and sum it up
    for( int p = 0; p < nOrder; p++ )
    {
        for( int q = 0; q < nOrder; q++ )
        {

            const double dA = dG[p];
            const double dB = dG[q];

            // V3d dF = (1.0 - dB) * (((1.0 - dA) * node1) + (dA * node2)) + (dB * node3);
            dF  = ( ( ( 1.0 - dB ) * ( 1.0 - dA ) ) * node1 ) + ( ( ( 1.0 - dB ) * dA ) * node2 ) + ( dB * node3 );
            dF2 = dF.array().square();

            dDaF = ( node1 - node2 );

            // dDbF = -( 1.0 - dA ) * node1 - dA * node2 + node3;
            dDbF = ( node3 - node1 ) + dA * dDaF;
            dDaF *= ( dB - 1.0 );

            const double dDenomTerm = std::pow( dF.norm(), -3.0 );

            // Eigen::Vector3d temp1 = dF2;
            // temp1( 2 ) += dF2( 0 );                           // temp1 = [dF2(0), dF2(1), dF2(0)+dF2(2)]
            // Eigen::Vector3d temp2 = dDaF.cwiseProduct( dF );  // temp2 = [dDaF(0)*dF(0), dDaF(1)*dF(1), dDaF(2)*dF(2)]
            // dDaG                  = dDaF.cwiseProduct( temp1.segment( 1, 2 ) ) - temp2.segment( 1, 2 );
            // Eigen::Vector3d temp3 = dDbF.cwiseProduct( dF );  // temp3 = [dDbF(0)*dF(0), dDbF(1)*dF(1), dDbF(2)*dF(2)]
            // dDbG                  = dDbF.cwiseProduct( temp1.segment( 1, 2 ) ) - temp3.segment( 1, 2 );

            // Eigen::Vector3d dF2_shifted = dF2 + Eigen::Vector3d( dF2( 1 ), dF2( 2 ), dF2( 0 ) );
            // Eigen::Vector3d dDaF_shifted = Eigen::Vector3d( dDaF( 1 ), dDaF( 2 ), dDaF( 0 ) ).cwiseProduct( dF );

            // dDaG = dDaF.cwiseProduct( dF2_shifted ) - dF.cwiseProduct( dDaF_shifted );

            // Eigen::Vector3d dDbF_shifted = Eigen::Vector3d( dDbF( 1 ), dDbF( 2 ), dDbF( 0 ) ).cwiseProduct( dF );

            // dDbG = dDbF.cwiseProduct( dF2_shifted ) - dF.cwiseProduct( dDbF_shifted );

            dDaG( 0 ) = dDaF( 0 ) * ( dF2( 1 ) + dF2( 2 ) ) - dF( 0 ) * ( dDaF( 1 ) * dF( 1 ) + dDaF( 2 ) * dF( 2 ) );
            dDaG( 1 ) = dDaF( 1 ) * ( dF2( 0 ) + dF2( 2 ) ) - dF( 1 ) * ( dDaF( 0 ) * dF( 0 ) + dDaF( 2 ) * dF( 2 ) );
            dDaG( 2 ) = dDaF( 2 ) * ( dF2( 0 ) + dF2( 1 ) ) - dF( 2 ) * ( dDaF( 0 ) * dF( 0 ) + dDaF( 1 ) * dF( 1 ) );

            // dDbG      = dDbF.cwiseProduct( ( dF2.array().shift( -1 ) + dF2.array().shift( -2 ) ) ) + dF.cwiseProduct( dDbF );
            dDbG( 0 ) = dDbF( 0 ) * ( dF2( 1 ) + dF2( 2 ) ) - dF( 0 ) * ( dDbF( 1 ) * dF( 1 ) + dDbF( 2 ) * dF( 2 ) );
            dDbG( 1 ) = dDbF( 1 ) * ( dF2( 0 ) + dF2( 2 ) ) - dF( 1 ) * ( dDbF( 0 ) * dF( 0 ) + dDbF( 2 ) * dF( 2 ) );
            dDbG( 2 ) = dDbF( 2 ) * ( dF2( 0 ) + dF2( 1 ) ) - dF( 2 ) * ( dDbF( 0 ) * dF( 0 ) + dDbF( 1 ) * dF( 1 ) );

            // Scale the vectors by the denominator
            dDaG *= dDenomTerm;
            dDbG *= dDenomTerm;

            // Cross product gives local Jacobian: dGaG x dDbG
            nodeCross[0] = dDaG( 1 ) * dDbG( 2 ) - dDaG( 2 ) * dDbG( 1 );
            nodeCross[1] = dDaG( 2 ) * dDbG( 0 ) - dDaG( 0 ) * dDbG( 2 );
            nodeCross[2] = dDaG( 0 ) * dDbG( 1 ) - dDaG( 1 ) * dDbG( 0 );

            const double dJacobian =
                std::sqrt( nodeCross[0] * nodeCross[0] + nodeCross[1] * nodeCross[1] + nodeCross[2] * nodeCross[2] );

            // dFaceArea += 2.0 * dW[p] * dW[q] * (1.0 - dG[q]) * dJacobian;
            dFaceArea += dW[p] * dW[q] * dJacobian;
        }
    }

    return dFaceArea;
#else
    /* compute the area by using Gauss-Quadratures; use TR interfaces directly */
    Face face( 3 );
    NodeVector nodes( 3 );
    nodes[0] = Node( inode1[0], inode1[1], inode1[2] );
    nodes[1] = Node( inode2[0], inode2[1], inode2[2] );
    nodes[2] = Node( inode3[0], inode3[1], inode3[2] );
    face.SetNode( 0, 0 );
    face.SetNode( 1, 1 );
    face.SetNode( 2, 2 );
    return CalculateFaceArea( face, nodes );
#endif
}

#endif

/*
 *  l'Huiller's formula for spherical triangle
 *  http://williams.best.vwh.net/avform.htm
 *  a, b, c are arc measures in radians, too
 *  A, B, C are angles on the sphere, for which we already have formula
 *               c
 *         A -------B
 *          \       |
 *           \      |
 *            \b    |a
 *             \    |
 *              \   |
 *               \  |
 *                \C|
 *                 \|
 *
 *  (The angle at B is not necessarily a right angle)
 *
 *    sin(a)  sin(b)   sin(c)
 *    ----- = ------ = ------
 *    sin(A)  sin(B)   sin(C)
 *
 * In terms of the sides (this is excess, as before, but numerically stable)
 *
 *  E = 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)))
 */
double IntxAreaUtils::area_spherical_triangle_lHuiller( const double* ptA,
                                                        const double* ptB,
                                                        const double* ptC,
                                                        double Radius )
{

    // now, a is angle BOC, O is origin
    CartVect vA( ptA ), vB( ptB ), vC( ptC );
    double a = angle_robust( vB, vC );
    double b = angle_robust( vC, vA );
    double c = angle_robust( vA, vB );
    int sign = 1;
    // if( fabs( ( vA * vB ) % vC ) < 1e-17 ) sign = -1;
    if( ( vA * vB ) % vC < 0 ) sign = -1;
    double s  = ( a + b + c ) / 2;
    double a1 = ( s - a ) / 2;
    double b1 = ( s - b ) / 2;
    double c1 = ( s - c ) / 2;
#ifdef MOAB_HAVE_TEMPESTREMAP
    if( fabs( a1 ) < 1.e-14 || fabs( b1 ) < 1.e-14 || fabs( c1 ) < 1.e-14 )
    {
        double area = area_spherical_triangle_GQ( ptA, ptB, ptC ) * sign;
#ifdef VERBOSE
        std::cout << " very obtuse angle, use TR to compute area " << " a1:" << a1 << " b1:" << b1 << " c1:" << c1
                  << "\n";
        std::cout << " area with TR: " << area << "\n";
#endif
        return area;
    }
#endif
    double tmp = tan( s / 2 ) * tan( a1 ) * tan( b1 ) * tan( c1 );
    if( tmp < 0. ) tmp = 0.;

    double E = 4 * atan( sqrt( tmp ) );
    if( E != E ) std::cout << " NaN at spherical triangle area \n";

    double area = sign * E * Radius * Radius;

#ifdef CHECKNEGATIVEAREA
    if( area < 0 )
    {
        std::cout << "negative area: " << area << "\n";
        std::cout << std::setprecision( 15 );
        std::cout << "vA: " << vA << "\n";
        std::cout << "vB: " << vB << "\n";
        std::cout << "vC: " << vC << "\n";
        std::cout << "sign: " << sign << "\n";
        std::cout << " a: " << a << "\n";
        std::cout << " b: " << b << "\n";
        std::cout << " c: " << c << "\n";
    }
#endif

    return area;
}

double IntxAreaUtils::area_on_sphere( Interface* mb, EntityHandle set, double R )
{
    // Get all entities of dimension 2
    Range inputRange;
    ErrorCode rval = mb->get_entities_by_dimension( set, 2, inputRange );MB_CHK_ERR_RET_VAL( rval, -1.0 );

    // Filter by elements that are owned by current process
    std::vector< int > ownerinfo( inputRange.size(), -1 );
    Tag intxOwnerTag;
    rval = mb->tag_get_handle( "ORIG_PROC", intxOwnerTag );
    if( MB_SUCCESS == rval )
    {
        rval = mb->tag_get_data( intxOwnerTag, inputRange, &ownerinfo[0] );MB_CHK_ERR_RET_VAL( rval, -1.0 );
    }

    // compare total area with 4*M_PI * R^2
    int ie            = 0;
    double total_area = 0.;
    for( Range::iterator eit = inputRange.begin(); eit != inputRange.end(); ++eit )
    {

        // All zero/positive owner data represents ghosted elems
        if( ownerinfo[ie++] >= 0 ) continue;

        EntityHandle eh        = *eit;
        const double elem_area = this->area_spherical_element( mb, eh, R );

        // check whether the area of the spherical element is positive.
        if( elem_area <= 0 )
        {
            std::cout << "Area of element " << mb->id_from_handle( eh ) << " is = " << elem_area << "\n";
            mb->list_entity( eh );
        }
        assert( elem_area > 0 );

        // sum up the contribution
        total_area += elem_area;
    }

    // return total mesh area
    return total_area;
}

double IntxAreaUtils::area_spherical_element( Interface* mb, EntityHandle elem, double R )
{
    // get the nodes, then the coordinates
    const EntityHandle* verts;
    int nsides;
    ErrorCode rval = mb->get_connectivity( elem, verts, nsides );MB_CHK_ERR_RET_VAL( rval, -1.0 );

    // account for possible padded polygons
    while( verts[nsides - 2] == verts[nsides - 1] && nsides > 3 )
        nsides--;

    // get coordinates
    std::vector< double > coords( 3 * nsides );
    rval = mb->get_coords( verts, nsides, &coords[0] );MB_CHK_ERR_RET_VAL( rval, -1.0 );

    // compute and return the area of the polygonal element
    return area_spherical_polygon( &coords[0], nsides, R );
}

double IntxUtils::distance_on_great_circle( CartVect& p1, CartVect& p2 )
{
    SphereCoords sph1 = cart_to_spherical( p1 );
    SphereCoords sph2 = cart_to_spherical( p2 );
    // radius should be the same
    return sph1.R *
           acos( sin( sph1.lon ) * sin( sph2.lon ) + cos( sph1.lat ) * cos( sph2.lat ) * cos( sph2.lon - sph2.lon ) );
}

//
/**
 * @brief Enforces convexity for a given set of polygons.
 *
 * This function checks each polygon in the input set and computes the angles of each vertex.
 * If a reflex angle is found, the polygon is broken into triangles and added back to the set.
 * This process continues until all polygons in the set are convex.
 *
 * @param mb The interface to the MOAB instance.
 * @param lset The handle of the input set containing the polygons.
 * @param my_rank The rank of the local process.
 * @return The error code indicating the success or failure of the operation.
 */
ErrorCode IntxUtils::enforce_convexity( Interface* mb, EntityHandle lset, int my_rank )
{
    Range inputRange;
    ErrorCode rval = mb->get_entities_by_dimension( lset, 2, inputRange );MB_CHK_ERR( rval );

    Tag corrTag       = 0;
    EntityHandle dumH = 0;
    rval              = mb->tag_get_handle( CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag, MB_TAG_DENSE, &dumH );
    if( rval == MB_TAG_NOT_FOUND ) corrTag = 0;

    Tag gidTag;
    rval = mb->tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gidTag, MB_TAG_DENSE );MB_CHK_ERR( rval );

    std::vector< double > coords;
    coords.resize( 3 * MAXEDGES );  // at most 10 vertices per polygon
    // we should create a queue with new polygons that need processing for reflex angles
    //  (obtuse)
    std::queue< EntityHandle > newPolys;
    int brokenPolys     = 0;
    Range::iterator eit = inputRange.begin();
    while( eit != inputRange.end() || !newPolys.empty() )
    {
        EntityHandle eh;
        if( eit != inputRange.end() )
        {
            eh = *eit;
            ++eit;
        }
        else
        {
            eh = newPolys.front();
            newPolys.pop();
        }
        // get the nodes, then the coordinates
        const EntityHandle* verts;
        int num_nodes;
        rval = mb->get_connectivity( eh, verts, num_nodes );MB_CHK_ERR( rval );
        int nsides = num_nodes;
        // account for possible padded polygons
        while( verts[nsides - 2] == verts[nsides - 1] && nsides > 3 )
            nsides--;
        EntityHandle corrHandle = 0;
        if( corrTag )
        {
            rval = mb->tag_get_data( corrTag, &eh, 1, &corrHandle );MB_CHK_ERR( rval );
        }
        int gid = 0;
        rval    = mb->tag_get_data( gidTag, &eh, 1, &gid );MB_CHK_ERR( rval );
        coords.resize( 3 * nsides );
        if( nsides < 4 ) continue;  // if already triangles, don't bother
        // get coordinates
        rval = mb->get_coords( verts, nsides, &coords[0] );MB_CHK_ERR( rval );
        // compute each angle
        bool alreadyBroken = false;

        for( int i = 0; i < nsides; i++ )
        {
            double* A    = &coords[3 * i];
            double* B    = &coords[3 * ( ( i + 1 ) % nsides )];
            double* C    = &coords[3 * ( ( i + 2 ) % nsides )];
            double angle = IntxUtils::oriented_spherical_angle( A, B, C );
            if( angle - M_PI > 0. )  // even almost reflex is bad; break it!
            {
                if( alreadyBroken )
                {
                    mb->list_entities( &eh, 1 );
                    mb->list_entities( verts, nsides );
                    double* D = &coords[3 * ( ( i + 3 ) % nsides )];
                    std::cout << "ABC: " << angle << " \n";
                    std::cout << "BCD: " << IntxUtils::oriented_spherical_angle( B, C, D ) << " \n";
                    std::cout << "CDA: " << IntxUtils::oriented_spherical_angle( C, D, A ) << " \n";
                    std::cout << "DAB: " << IntxUtils::oriented_spherical_angle( D, A, B ) << " \n";
                    std::cout << " this cell has at least 2 angles > 180, it has serious issues\n";

                    return MB_FAILURE;
                }
                // the bad angle is at i+1;
                // create 1 triangle and one polygon; add the polygon to the input range, so
                // it will be processed too
                // also, add both to the set :) and remove the original polygon from the set
                // break the next triangle, even though not optimal
                // so create the triangle i+1, i+2, i+3; remove i+2 from original list
                // even though not optimal in general, it is good enough.
                EntityHandle conn3[3] = { verts[( i + 1 ) % nsides], verts[( i + 2 ) % nsides],
                                          verts[( i + 3 ) % nsides] };
                // create a polygon with num_nodes-1 vertices, and connectivity
                // verts[i+1], verts[i+3], (all except i+2)
                std::vector< EntityHandle > conn( nsides - 1 );
                for( int j = 1; j < nsides; j++ )
                {
                    conn[j - 1] = verts[( i + j + 2 ) % nsides];
                }
                EntityHandle newElement;
                rval = mb->create_element( MBTRI, conn3, 3, newElement );MB_CHK_ERR( rval );

                rval = mb->add_entities( lset, &newElement, 1 );MB_CHK_ERR( rval );
                if( corrTag )
                {
                    rval = mb->tag_set_data( corrTag, &newElement, 1, &corrHandle );MB_CHK_ERR( rval );
                }
                rval = mb->tag_set_data( gidTag, &newElement, 1, &gid );MB_CHK_ERR( rval );
                if( nsides == 4 )
                {
                    // create another triangle
                    rval = mb->create_element( MBTRI, &conn[0], 3, newElement );MB_CHK_ERR( rval );
                }
                else
                {
                    // create another polygon, and add it to the inputRange
                    rval = mb->create_element( MBPOLYGON, &conn[0], nsides - 1, newElement );MB_CHK_ERR( rval );
                    newPolys.push( newElement );  // because it has less number of edges, the
                    // reverse should work to find it.
                }
                rval = mb->add_entities( lset, &newElement, 1 );MB_CHK_ERR( rval );
                if( corrTag )
                {
                    rval = mb->tag_set_data( corrTag, &newElement, 1, &corrHandle );MB_CHK_ERR( rval );
                }
                rval = mb->tag_set_data( gidTag, &newElement, 1, &gid );MB_CHK_ERR( rval );
                rval = mb->remove_entities( lset, &eh, 1 );MB_CHK_ERR( rval );
                brokenPolys++;
                alreadyBroken = true;  // get out of the loop, element is broken
            }
        }
    }
    if( brokenPolys > 0 )
    {
        std::cout << "on local process " << my_rank << ", " << brokenPolys
                  << " concave polygons were decomposed in convex ones \n";
#ifdef VERBOSE
        std::stringstream fff;
        fff << "file_set" << mb->id_from_handle( lset ) << "rk_" << my_rank << ".h5m";
        rval = mb->write_file( fff.str().c_str(), 0, 0, &lset, 1 );MB_CHK_ERR( rval );
        std::cout << "wrote new file set: " << fff.str() << "\n";
#endif
    }
    return MB_SUCCESS;
}

// looking at quad connectivity, collapse to triangle if 2 nodes equal
// then delete the old quad
ErrorCode IntxUtils::fix_degenerate_quads( Interface* mb, EntityHandle set )
{
    Range quads;
    ErrorCode rval = mb->get_entities_by_type( set, MBQUAD, quads );MB_CHK_ERR( rval );
    Tag gid;
    gid = mb->globalId_tag();
    for( Range::iterator qit = quads.begin(); qit != quads.end(); ++qit )
    {
        EntityHandle quad         = *qit;
        const EntityHandle* conn4 = NULL;
        int num_nodes             = 0;
        rval                      = mb->get_connectivity( quad, conn4, num_nodes );MB_CHK_ERR( rval );
        for( int i = 0; i < num_nodes; i++ )
        {
            int next_node_index = ( i + 1 ) % num_nodes;
            if( conn4[i] == conn4[next_node_index] )
            {
                // form a triangle and delete the quad
                // first get the global id, to set it on triangle later
                int global_id = 0;
                rval          = mb->tag_get_data( gid, &quad, 1, &global_id );MB_CHK_ERR( rval );
                int i2                = ( i + 2 ) % num_nodes;
                int i3                = ( i + 3 ) % num_nodes;
                EntityHandle conn3[3] = { conn4[i], conn4[i2], conn4[i3] };
                EntityHandle tri;
                rval = mb->create_element( MBTRI, conn3, 3, tri );MB_CHK_ERR( rval );
                mb->add_entities( set, &tri, 1 );
                mb->remove_entities( set, &quad, 1 );
                mb->delete_entities( &quad, 1 );
                rval = mb->tag_set_data( gid, &tri, 1, &global_id );MB_CHK_ERR( rval );
            }
        }
    }
    return MB_SUCCESS;
}

ErrorCode IntxAreaUtils::positive_orientation( Interface* mb, EntityHandle set, double R )
{
    Range cells2d;
    ErrorCode rval = mb->get_entities_by_dimension( set, 2, cells2d );MB_CHK_ERR( rval );
    for( Range::iterator qit = cells2d.begin(); qit != cells2d.end(); ++qit )
    {
        EntityHandle cell        = *qit;
        const EntityHandle* conn = NULL;
        int num_nodes            = 0;
        rval                     = mb->get_connectivity( cell, conn, num_nodes );MB_CHK_ERR( rval );
        if( num_nodes < 3 ) return MB_FAILURE;

        double coords[9];
        rval = mb->get_coords( conn, 3, coords );MB_CHK_ERR( rval );

        double area;
        if( R > 0 )
            area = area_spherical_triangle_lHuiller( coords, coords + 3, coords + 6, R );
        else
            area = IntxUtils::area2D( coords, coords + 3, coords + 6 );
        if( area < 0 )
        {
            // compute all area, do not revert if total area is positive
            std::vector< double > coords2( 3 * num_nodes );
            // get coordinates
            rval = mb->get_coords( conn, num_nodes, &coords2[0] );MB_CHK_ERR( rval );
            double totArea = area_spherical_polygon_lHuiller( &coords2[0], num_nodes, R );
            if( totArea < 0 )
            {
                std::vector< EntityHandle > newconn( num_nodes );
                for( int i = 0; i < num_nodes; i++ )
                {
                    newconn[num_nodes - 1 - i] = conn[i];
                }
                rval = mb->set_connectivity( cell, &newconn[0], num_nodes );MB_CHK_ERR( rval );
            }
            else
            {
                std::cout << " nonconvex problem first area:" << area << " total area: " << totArea << std::endl;
            }
        }
    }
    return MB_SUCCESS;
}

// distance along a great circle on a sphere of radius 1
// page 4
double IntxUtils::distance_on_sphere( double la1, double te1, double la2, double te2 )
{
    return acos( sin( te1 ) * sin( te2 ) + cos( te1 ) * cos( te2 ) * cos( la1 - la2 ) );
}

/*
 * given 2 great circle arcs, AB and CD, compute the unique intersection point, if it exists
 *  in between
 */
ErrorCode IntxUtils::intersect_great_circle_arcs( double* A, double* B, double* C, double* D, double R, double* E )
{
    // first verify A, B, C, D are on the same sphere
    double R2              = R * R;
    const double Tolerance = 1.e-12 * R2;

    CartVect a( A ), b( B ), c( C ), d( D );

    if( fabs( a.length_squared() - R2 ) + fabs( b.length_squared() - R2 ) + fabs( c.length_squared() - R2 ) +
            fabs( d.length_squared() - R2 ) >
        10 * Tolerance )
        return MB_FAILURE;

    CartVect n1 = a * b;
    if( n1.length_squared() < Tolerance ) return MB_FAILURE;

    CartVect n2 = c * d;
    if( n2.length_squared() < Tolerance ) return MB_FAILURE;
    CartVect n3 = n1 * n2;
    n3.normalize();

    n3 = R * n3;
    // the intersection is either n3 or -n3
    CartVect n4 = a * n3, n5 = n3 * b;
    if( n1 % n4 >= -Tolerance && n1 % n5 >= -Tolerance )
    {
        // n3 is good for ab, see if it is good for cd
        n4 = c * n3;
        n5 = n3 * d;
        if( n2 % n4 >= -Tolerance && n2 % n5 >= -Tolerance )
        {
            E[0] = n3[0];
            E[1] = n3[1];
            E[2] = n3[2];
        }
        else
            return MB_FAILURE;
    }
    else
    {
        // try -n3
        n3 = -n3;
        n4 = a * n3, n5 = n3 * b;
        if( n1 % n4 >= -Tolerance && n1 % n5 >= -Tolerance )
        {
            // n3 is good for ab, see if it is good for cd
            n4 = c * n3;
            n5 = n3 * d;
            if( n2 % n4 >= -Tolerance && n2 % n5 >= -Tolerance )
            {
                E[0] = n3[0];
                E[1] = n3[1];
                E[2] = n3[2];
            }
            else
                return MB_FAILURE;
        }
        else
            return MB_FAILURE;
    }

    return MB_SUCCESS;
}

// verify that result is in between a and b on a great circle arc, and between c and d on a constant
// latitude arc
static bool verify( CartVect a, CartVect b, CartVect c, CartVect d, double x, double y, double z )
{
    // to check, the point has to be between a and b on a great arc, and between c and d on a const
    // lat circle
    CartVect s( x, y, z );
    CartVect n1 = a * b;
    CartVect n2 = a * s;
    CartVect n3 = s * b;
    if( n1 % n2 < 0 || n1 % n3 < 0 ) return false;

    // do the same for c, d, s, in plane z=0
    c[2] = d[2] = s[2] = 0.;  // bring everything in the same plane, z=0;

    n1 = c * d;
    n2 = c * s;
    n3 = s * d;
    if( n1 % n2 < 0 || n1 % n3 < 0 ) return false;

    return true;
}

ErrorCode IntxUtils::intersect_great_circle_arc_with_clat_arc( double* A,
                                                               double* B,
                                                               double* C,
                                                               double* D,
                                                               double R,
                                                               double* E,
                                                               int& np )
{
    const double distTol   = R * 1.e-6;
    const double Tolerance = R * R * 1.e-12;  // radius should be 1, usually
    np                     = 0;               // number of points in intersection
    CartVect a( A ), b( B ), c( C ), d( D );
    // check input first
    double R2 = R * R;
    if( fabs( a.length_squared() - R2 ) + fabs( b.length_squared() - R2 ) + fabs( c.length_squared() - R2 ) +
            fabs( d.length_squared() - R2 ) >
        10 * Tolerance )
        return MB_FAILURE;

    if( ( a - b ).length_squared() < Tolerance ) return MB_FAILURE;
    if( ( c - d ).length_squared() < Tolerance )  // edges are too short
        return MB_FAILURE;

    // CD is the const latitude arc
    if( fabs( C[2] - D[2] ) > distTol )  // cd is not on the same z (constant latitude)
        return MB_FAILURE;

    if( fabs( R - C[2] ) < distTol || fabs( R + C[2] ) < distTol ) return MB_FAILURE;  // too close to the poles

    // find the points on the circle P(teta) = (r*sin(teta), r*cos(teta), C[2]) that are on the
    // great circle arc AB normal to the AB circle:
    CartVect n1 = a * b;  // the normal to the great circle arc (circle)
    // solve the system of equations:
    /*
     *    n1%(x, y, z) = 0  // on the great circle
     *     z = C[2];
     *    x^2+y^2+z^2 = R^2
     */
    double z = C[2];
    if( fabs( n1[0] ) + fabs( n1[1] ) < 2 * Tolerance )
    {
        // it is the Equator; check if the const lat edge is Equator too
        if( fabs( C[2] ) > distTol )
        {
            return MB_FAILURE;  // no intx, too far from Eq
        }
        else
        {
            // all points are on the equator
            //
            CartVect cd = c * d;
            // by convention,  c<d, positive is from c to d
            // is a or b between c , d?
            CartVect ca = c * a;
            CartVect ad = a * d;
            CartVect cb = c * b;
            CartVect bd = b * d;
            bool agtc   = ( ca % cd >= -Tolerance );  // a>c?
            bool dgta   = ( ad % cd >= -Tolerance );  // d>a?
            bool bgtc   = ( cb % cd >= -Tolerance );  // b>c?
            bool dgtb   = ( bd % cd >= -Tolerance );  // d>b?
            if( agtc )
            {
                if( dgta )
                {
                    // a is for sure a point
                    E[0] = a[0];
                    E[1] = a[1];
                    E[2] = a[2];
                    np++;
                    if( bgtc )
                    {
                        if( dgtb )
                        {
                            // b is also in between c and d
                            E[3] = b[0];
                            E[4] = b[1];
                            E[5] = b[2];
                            np++;
                        }
                        else
                        {
                            // then order is c a d b, intx is ad
                            E[3] = d[0];
                            E[4] = d[1];
                            E[5] = d[2];
                            np++;
                        }
                    }
                    else
                    {
                        // b is less than c, so b c a d, intx is ac
                        E[3] = c[0];
                        E[4] = c[1];
                        E[5] = c[2];
                        np++;  // what if E[0] is E[3]?
                    }
                }
                else  // c < d < a
                {
                    if( dgtb )  // d is for sure in
                    {
                        E[0] = d[0];
                        E[1] = d[1];
                        E[2] = d[2];
                        np++;
                        if( bgtc )  // c<b<d<a
                        {
                            // another point is b
                            E[3] = b[0];
                            E[4] = b[1];
                            E[5] = b[2];
                            np++;
                        }
                        else  // b<c<d<a
                        {
                            // another point is c
                            E[3] = c[0];
                            E[4] = c[1];
                            E[5] = c[2];
                            np++;
                        }
                    }
                    else
                    {
                        // nothing, order is c, d < a, b
                    }
                }
            }
            else  // a < c < d
            {
                if( bgtc )
                {
                    // c is for sure in
                    E[0] = c[0];
                    E[1] = c[1];
                    E[2] = c[2];
                    np++;
                    if( dgtb )
                    {
                        // a < c < b < d; second point is b
                        E[3] = b[0];
                        E[4] = b[1];
                        E[5] = b[2];
                        np++;
                    }
                    else
                    {
                        // a < c < d < b; second point is d
                        E[3] = d[0];
                        E[4] = d[1];
                        E[5] = d[2];
                        np++;
                    }
                }
                else  // a, b < c < d
                {
                    // nothing
                }
            }
        }
        // for the 2 points selected, see if it is only one?
        // no problem, maybe it will be collapsed later anyway
        if( np > 0 ) return MB_SUCCESS;
        return MB_FAILURE;  // no intersection
    }
    {
        if( fabs( n1[0] ) <= fabs( n1[1] ) )
        {
            // resolve eq in x:  n0 * x + n1 * y +n2*z = 0; y = -n2/n1*z -n0/n1*x
            //  (u+v*x)^2+x^2=R2-z^2
            //  (v^2+1)*x^2 + 2*u*v *x + u^2+z^2-R^2 = 0
            //  delta = 4*u^2*v^2 - 4*(v^2-1)(u^2+z^2-R^2)
            // x1,2 =
            double u = -n1[2] / n1[1] * z, v = -n1[0] / n1[1];
            double a1 = v * v + 1, b1 = 2 * u * v, c1 = u * u + z * z - R2;
            double delta = b1 * b1 - 4 * a1 * c1;
            if( delta < -Tolerance ) return MB_FAILURE;  // no intersection
            if( delta > Tolerance )                      // 2 solutions possible
            {
                double x1 = ( -b1 + sqrt( delta ) ) / 2 / a1;
                double x2 = ( -b1 - sqrt( delta ) ) / 2 / a1;
                double y1 = u + v * x1;
                double y2 = u + v * x2;
                if( verify( a, b, c, d, x1, y1, z ) )
                {
                    E[0] = x1;
                    E[1] = y1;
                    E[2] = z;
                    np++;
                }
                if( verify( a, b, c, d, x2, y2, z ) )
                {
                    E[3 * np + 0] = x2;
                    E[3 * np + 1] = y2;
                    E[3 * np + 2] = z;
                    np++;
                }
            }
            else
            {
                // one solution
                double x1 = -b1 / 2 / a1;
                double y1 = u + v * x1;
                if( verify( a, b, c, d, x1, y1, z ) )
                {
                    E[0] = x1;
                    E[1] = y1;
                    E[2] = z;
                    np++;
                }
            }
        }
        else
        {
            // resolve eq in y, reverse
            //   n0 * x + n1 * y +n2*z = 0; x = -n2/n0*z -n1/n0*y = u+v*y
            //  (u+v*y)^2+y^2 -R2+z^2 =0
            //  (v^2+1)*y^2 + 2*u*v *y + u^2+z^2-R^2 = 0
            //
            // x1,2 =
            double u = -n1[2] / n1[0] * z, v = -n1[1] / n1[0];
            double a1 = v * v + 1, b1 = 2 * u * v, c1 = u * u + z * z - R2;
            double delta = b1 * b1 - 4 * a1 * c1;
            if( delta < -Tolerance ) return MB_FAILURE;  // no intersection
            if( delta > Tolerance )                      // 2 solutions possible
            {
                double y1 = ( -b1 + sqrt( delta ) ) / 2 / a1;
                double y2 = ( -b1 - sqrt( delta ) ) / 2 / a1;
                double x1 = u + v * y1;
                double x2 = u + v * y2;
                if( verify( a, b, c, d, x1, y1, z ) )
                {
                    E[0] = x1;
                    E[1] = y1;
                    E[2] = z;
                    np++;
                }
                if( verify( a, b, c, d, x2, y2, z ) )
                {
                    E[3 * np + 0] = x2;
                    E[3 * np + 1] = y2;
                    E[3 * np + 2] = z;
                    np++;
                }
            }
            else
            {
                // one solution
                double y1 = -b1 / 2 / a1;
                double x1 = u + v * y1;
                if( verify( a, b, c, d, x1, y1, z ) )
                {
                    E[0] = x1;
                    E[1] = y1;
                    E[2] = z;
                    np++;
                }
            }
        }
    }

    if( np <= 0 ) return MB_FAILURE;
    return MB_SUCCESS;
}

#if 0
ErrorCode set_edge_type_flag(Interface * mb, EntityHandle sf1)
{
  Range cells;
  ErrorCode rval = mb->get_entities_by_dimension(sf1, 2, cells);
  if (MB_SUCCESS!= rval)
  return rval;
  Range edges;
  rval = mb->get_adjacencies(cells, 1, true, edges, Interface::UNION);
  if (MB_SUCCESS!= rval)
  return rval;

  Tag edgeTypeTag;
  int default_int=0;
  rval = mb->tag_get_handle("edge_type", 1, MB_TYPE_INTEGER, edgeTypeTag,
      MB_TAG_DENSE | MB_TAG_CREAT, &default_int);
  if (MB_SUCCESS!= rval)
  return rval;
  // add edges to the set? not yet, maybe later
  // if edge horizontal, set value to 1
  int type_constant_lat=1;
  for (Range::iterator eit=edges.begin(); eit!=edges.end(); ++eit)
  {
    EntityHandle edge = *eit;
    const EntityHandle *conn=0;
    int num_n=0;
    rval = mb->get_connectivity(edge, conn, num_n );
    if (MB_SUCCESS!= rval)
    return rval;
    double coords[6];
    rval = mb->get_coords(conn, 2, coords);
    if (MB_SUCCESS!= rval)
    return rval;
    if (fabs( coords[2]-coords[5] )< 1.e-6 )
    {
      rval = mb->tag_set_data(edgeTypeTag, &edge, 1, &type_constant_lat);
      if (MB_SUCCESS!= rval)
      return rval;
    }
  }

  return MB_SUCCESS;
}
#endif

// decide in a different metric if the corners of CS quad are
// in the interior of an RLL quad
int IntxUtils::borderPointsOfCSinRLL( CartVect* redc,
                                      double* red2dc,
                                      int nsRed,
                                      CartVect* bluec,
                                      int nsBlue,
                                      int* blueEdgeType,
                                      double* P,
                                      int* side,
                                      double epsil )
{
    int extraPoints = 0;
    // first decide the blue z coordinates
    CartVect A( 0. ), B( 0. ), C( 0. ), D( 0. );
    for( int i = 0; i < nsBlue; i++ )
    {
        if( blueEdgeType[i] == 0 )
        {
            int iP1 = ( i + 1 ) % nsBlue;
            if( bluec[i][2] > bluec[iP1][2] )
            {
                A = bluec[i];
                B = bluec[iP1];
                C = bluec[( i + 2 ) % nsBlue];
                D = bluec[( i + 3 ) % nsBlue];  // it could be back to A, if triangle on top
                break;
            }
        }
    }
    if( nsBlue == 3 && B[2] < 0 )
    {
        // select D to be C
        D = C;
        C = B;  // B is the south pole then
    }
    // so we have A, B, C, D, with A at the top, b going down, then C, D, D > C, A > B
    // AB is const longitude, BC and DA constant latitude
    // check now each of the red points if they are inside this rectangle
    for( int i = 0; i < nsRed; i++ )
    {
        CartVect& X = redc[i];
        if( X[2] > A[2] || X[2] < B[2] ) continue;  // it is above or below the rectangle
        // now decide if it is between the planes OAB and OCD
        if( ( ( A * B ) % X >= -epsil ) && ( ( C * D ) % X >= -epsil ) )
        {
            side[i] = 1;  //
            // it means point X is in the rectangle that we want , on the sphere
            // pass the coords 2d
            P[extraPoints * 2]     = red2dc[2 * i];
            P[extraPoints * 2 + 1] = red2dc[2 * i + 1];
            extraPoints++;
        }
    }
    return extraPoints;
}

ErrorCode IntxUtils::deep_copy_set_with_quads( Interface* mb, EntityHandle source_set, EntityHandle dest_set )
{
    ReadUtilIface* read_iface;
    ErrorCode rval = mb->query_interface( read_iface );MB_CHK_ERR( rval );
    // create the handle tag for the corresponding element / vertex

    EntityHandle dum = 0;
    Tag corrTag      = 0;  // it will be created here
    rval             = mb->tag_get_handle( CORRTAGNAME, 1, MB_TYPE_HANDLE, corrTag, MB_TAG_DENSE | MB_TAG_CREAT, &dum );MB_CHK_ERR( rval );

    // give the same global id to new verts and cells created in the lagr(departure) mesh
    Tag gid = mb->globalId_tag();

    Range quads;
    rval = mb->get_entities_by_type( source_set, MBQUAD, quads );MB_CHK_ERR( rval );

    Range connecVerts;
    rval = mb->get_connectivity( quads, connecVerts );MB_CHK_ERR( rval );

    std::map< EntityHandle, EntityHandle > newNodes;

    std::vector< double* > coords;
    EntityHandle start_vert, start_elem, *connect;
    int num_verts = connecVerts.size();
    rval          = read_iface->get_node_coords( 3, num_verts, 0, start_vert, coords );
    if( MB_SUCCESS != rval ) return rval;
    // fill it up
    int i = 0;
    for( Range::iterator vit = connecVerts.begin(); vit != connecVerts.end(); ++vit, i++ )
    {
        EntityHandle oldV = *vit;
        CartVect posi;
        rval = mb->get_coords( &oldV, 1, &( posi[0] ) );MB_CHK_ERR( rval );

        int global_id;
        rval = mb->tag_get_data( gid, &oldV, 1, &global_id );MB_CHK_ERR( rval );
        EntityHandle new_vert = start_vert + i;
        // Cppcheck warning (false positive): variable coords is assigned a value that is never used
        coords[0][i] = posi[0];
        coords[1][i] = posi[1];
        coords[2][i] = posi[2];

        newNodes[oldV] = new_vert;
        // set also the correspondent tag :)
        rval = mb->tag_set_data( corrTag, &oldV, 1, &new_vert );MB_CHK_ERR( rval );

        // also the other side
        // need to check if we really need this; the new vertex will never need the old vertex
        // we have the global id which is the same
        rval = mb->tag_set_data( corrTag, &new_vert, 1, &oldV );MB_CHK_ERR( rval );
        // set the global id on the corresponding vertex the same as the initial vertex
        rval = mb->tag_set_data( gid, &new_vert, 1, &global_id );MB_CHK_ERR( rval );
    }
    // now create new quads in order (in a sequence)

    rval = read_iface->get_element_connect( quads.size(), 4, MBQUAD, 0, start_elem, connect );
    if( MB_SUCCESS != rval ) return rval;
    int ie = 0;
    for( Range::iterator it = quads.begin(); it != quads.end(); ++it, ie++ )
    {
        EntityHandle q = *it;
        int nnodes;
        const EntityHandle* conn;
        rval = mb->get_connectivity( q, conn, nnodes );MB_CHK_ERR( rval );
        int global_id;
        rval = mb->tag_get_data( gid, &q, 1, &global_id );MB_CHK_ERR( rval );

        for( int ii = 0; ii < nnodes; ii++ )
        {
            EntityHandle v1      = conn[ii];
            connect[4 * ie + ii] = newNodes[v1];
        }
        EntityHandle newElement = start_elem + ie;

        // set the corresponding tag; not sure we need this one, from old to new
        rval = mb->tag_set_data( corrTag, &q, 1, &newElement );MB_CHK_ERR( rval );
        rval = mb->tag_set_data( corrTag, &newElement, 1, &q );MB_CHK_ERR( rval );

        // set the global id
        rval = mb->tag_set_data( gid, &newElement, 1, &global_id );MB_CHK_ERR( rval );

        rval = mb->add_entities( dest_set, &newElement, 1 );MB_CHK_ERR( rval );
    }

    rval = read_iface->update_adjacencies( start_elem, quads.size(), 4, connect );MB_CHK_ERR( rval );

    return MB_SUCCESS;
}

ErrorCode IntxUtils::remove_duplicate_vertices( Interface* mb,
                                                EntityHandle file_set,
                                                double merge_tol,
                                                std::vector< Tag >& tagList )
{
    Range verts;
    ErrorCode rval = mb->get_entities_by_dimension( file_set, 0, verts );MB_CHK_ERR( rval );
    rval = mb->remove_entities( file_set, verts );MB_CHK_ERR( rval );

    MergeMesh mm( mb );

    // remove the vertices from the set, before merging

    rval = mm.merge_all( file_set, merge_tol );MB_CHK_ERR( rval );

    // now correct vertices that are repeated in polygons
    rval = remove_padded_vertices( mb, file_set, tagList );
    return MB_SUCCESS;
}

ErrorCode IntxUtils::remove_padded_vertices( Interface* mb, EntityHandle file_set, std::vector< Tag >& tagList )
{

    // now correct vertices that are repeated in polygons
    Range cells;
    ErrorCode rval = mb->get_entities_by_dimension( file_set, 2, cells );MB_CHK_ERR( rval );

    Range verts;
    rval = mb->get_connectivity( cells, verts );MB_CHK_ERR( rval );

    Range modifiedCells;  // will be deleted at the end; keep the gid
    Range newCells;

    for( Range::iterator cit = cells.begin(); cit != cells.end(); ++cit )
    {
        EntityHandle cell          = *cit;
        const EntityHandle* connec = NULL;
        int num_verts              = 0;
        rval                       = mb->get_connectivity( cell, connec, num_verts );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );

        std::vector< EntityHandle > newConnec;
        newConnec.push_back( connec[0] );  // at least one vertex
        int index    = 0;
        int new_size = 1;
        while( index < num_verts - 2 )
        {
            int next_index = ( index + 1 );
            if( connec[next_index] != newConnec[new_size - 1] )
            {
                newConnec.push_back( connec[next_index] );
                new_size++;
            }
            index++;
        }
        // add the last one only if different from previous and first node
        if( ( connec[num_verts - 1] != connec[num_verts - 2] ) && ( connec[num_verts - 1] != connec[0] ) )
        {
            newConnec.push_back( connec[num_verts - 1] );
            new_size++;
        }
        if( new_size < num_verts )
        {
            // cout << "new cell from " << cell << " has only " << new_size << " vertices \n";
            modifiedCells.insert( cell );
            // create a new cell with type triangle, quad or polygon
            EntityType type = MBTRI;
            if( new_size == 3 )
                type = MBTRI;
            else if( new_size == 4 )
                type = MBQUAD;
            else if( new_size > 4 )
                type = MBPOLYGON;

            // create new cell
            EntityHandle newCell;
            rval = mb->create_element( type, &newConnec[0], new_size, newCell );MB_CHK_SET_ERR( rval, "Failed to create new cell" );
            // set the old id to the new element
            newCells.insert( newCell );
            double value;  // use the same value to reset the tags, even if the tags are int (like Global ID)
            for( size_t i = 0; i < tagList.size(); i++ )
            {
                rval = mb->tag_get_data( tagList[i], &cell, 1, (void*)( &value ) );MB_CHK_SET_ERR( rval, "Failed to get tag value" );
                rval = mb->tag_set_data( tagList[i], &newCell, 1, (void*)( &value ) );MB_CHK_SET_ERR( rval, "Failed to set tag value on new cell" );
            }
        }
    }

    rval = mb->remove_entities( file_set, modifiedCells );MB_CHK_SET_ERR( rval, "Failed to remove old cells from file set" );
    rval = mb->delete_entities( modifiedCells );MB_CHK_SET_ERR( rval, "Failed to delete old cells" );
    rval = mb->add_entities( file_set, newCells );MB_CHK_SET_ERR( rval, "Failed to add new cells to file set" );
    rval = mb->add_entities( file_set, verts );MB_CHK_SET_ERR( rval, "Failed to add verts to the file set" );

    return MB_SUCCESS;
}

ErrorCode IntxUtils::max_diagonal( Interface* mb, Range cells, int max_edges, double& diagonal )
{
    diagonal = 0.0;
    std::vector< CartVect > coords( max_edges );  // maximum number of nodes in a cell? hard coded ?
    for( auto it = cells.begin(); it != cells.end(); ++it )
    {
        // get the connectivity, then the coordinates
        EntityHandle cell          = *it;
        const EntityHandle* connec = NULL;
        int num_verts              = 0;
        MB_CHK_SET_ERR( mb->get_connectivity( cell, connec, num_verts ), "Failed to get connectivity" );
        MB_CHK_SET_ERR( mb->get_coords( connec, num_verts, &( coords[0][0] ) ), "Failed to get coordinates" );
        // compute the max diagonal, in a double loop
        for( int i = 0; i < num_verts - 1; i++ )
        {
            for( int j = i + 1; j < num_verts; j++ )
            {
                double len_sq = ( coords[i] - coords[j] ).length_squared();
                if( len_sq > diagonal ) diagonal = len_sq;
            }
        }
    }

    // return the root of the diagonal
    diagonal = std::sqrt( diagonal );
    return MB_SUCCESS;
}

}  // namespace moab
