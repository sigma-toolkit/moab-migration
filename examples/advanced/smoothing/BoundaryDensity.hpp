/** @file BoundaryDensity.hpp
 *  \brief Distance-to-boundary mesh sizing field on the unit sphere.
 *
 * Reads a boundary geometry file - coastlines, country outlines, a basin
 * perimeter - and turns it into the two weights a centroidal Voronoi
 * tessellation needs in order to refine near that boundary:
 *
 *   accept_probability()  for drawing the initial generators
 *   lloyd_weight()        for the density-weighted Lloyd sweep
 *
 * The sizing law is the tanh blend used by MPAS-Tools and JIGSAW for coastal
 * refinement,
 *
 *     h(d) = hmin + (hmax - hmin) * 0.5 * (1 + tanh((d - d0) / Lt))
 *
 * with d the geodesic distance to the nearest boundary segment in km.  Note
 * that h(0) is NOT hmin: the tanh only approaches its limits, so the width
 * actually reached at the boundary is h(0) = hmin + (hmax-hmin)*0.5*(1 +
 * tanh(-d0/Lt)), which for the defaults here is about 43 km rather than 30.
 * min_cell_width_km() reports the value actually attained, and the two
 * weights below are normalized by it rather than by hmin.
 *
 * All input geometry is expected in unit-sphere Cartesian coordinates, which
 * is what a lon/lat shapefile converted for MOAB ends up as.
 */

#ifndef MOAB_BOUNDARY_DENSITY_HPP
#define MOAB_BOUNDARY_DENSITY_HPP

#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"

namespace moab
{
namespace meshopt
{

//! Mean Earth radius, km.  Only used to express the sizing law in physical
//! units; the geometry itself is always on the unit sphere.
const double EARTH_RADIUS_KM = 6371.0;

struct BoundaryDensityOptions
{
    double hmin;   //!< km, cell width approached at the boundary
    double hmax;   //!< km, far-field cell width
    double d0;     //!< km, midpoint of the transition
    double lt;     //!< km, transition length

    BoundaryDensityOptions() : hmin( 30.0 ), hmax( 240.0 ), d0( 200.0 ), lt( 150.0 ) {}
};

/** Distance-to-boundary sizing field, backed by a MOAB kd-tree. */
class BoundaryDensity
{
  public:
    BoundaryDensity() : mTree( &mCore ), mRoot( 0 ), mMaxSegChord( 0.0 ), mMinWidth( 0.0 ) {}

    /** Load boundary geometry from a file.
     *
     * Every edge of every element in the file becomes a boundary segment:
     * polygons contribute their closed rings, edges contribute themselves.
     * That covers both a coastline stored as MBPOLYGON rings (what a
     * converted shapefile gives) and one stored as MBEDGE polylines.
     */
    ErrorCode load( const std::string& filename, const BoundaryDensityOptions& opts );

    /** Synthetic equatorial ring, for self-testing without a data file.
     *
     * The geodesic distance from any point to the equator is R*|asin(z)|, so
     * this gives the distance query an analytic answer to be checked against.
     */
    ErrorCode set_equator_ring( int nseg, const BoundaryDensityOptions& opts );

    //! Geodesic distance in km from a unit-sphere point to the nearest segment.
    double distance_km( const double* xyz ) const;

    //! Target cell width in km at a unit-sphere point.
    double cell_width_km( const double* xyz ) const
    {
        return width_from_distance( distance_km( xyz ) );
    }

    //! The sizing law itself, exposed so it can be tested without geometry.
    double width_from_distance( double d_km ) const
    {
        return mOpts.hmin + ( mOpts.hmax - mOpts.hmin ) * 0.5 * ( 1.0 + std::tanh( ( d_km - mOpts.d0 ) / mOpts.lt ) );
    }

    /** Weight for the density-weighted Lloyd sweep.
     *
     * A CVT does not place generators proportionally to its weight.  Lloyd's
     * algorithm converges to the optimal quantizer, whose point density in d
     * dimensions goes as rho^(d/(d+2)) - so rho^(1/2) in the plane (Gersho's
     * conjecture; Du, Faber & Gunzburger, SIAM Rev. 41 (1999) 637, sec. 5).
     * Cell width h means a point density of 1/h^2, so the weight has to be
     *
     *     rho = 1/h^4,   giving   rho^(1/2) = 1/h^2.
     *
     * Using 1/h^2 here - the intuitive choice, and a common mistake - would
     * converge to a point density of 1/h instead, i.e. the square root of the
     * grading actually asked for.  check_against_analytic() pins this down.
     */
    double lloyd_weight( const double* xyz ) const
    {
        const double r = mMinWidth / cell_width_km( xyz );
        return r * r * r * r;
    }

    /** Acceptance probability for rejection-sampling the initial generators.
     *
     * Rejection sampling sets the point density directly rather than through
     * a quantizer limit, so this is 1/h^2 and not the 1/h^4 above.  Normalized
     * by the smallest attainable width so the result is in (0,1].
     */
    double accept_probability( const double* xyz ) const
    {
        const double r = mMinWidth / cell_width_km( xyz );
        return r * r;
    }

    //! Smallest cell width the field actually attains, km.  See the file comment.
    double min_cell_width_km() const
    {
        return mMinWidth;
    }
    double max_cell_width_km() const
    {
        return mOpts.hmax;
    }
    size_t num_segments() const
    {
        return mSegs.size() / 2;
    }

    /** Self-check against the analytic equator distance.
     *
     * Returns the largest error in km over a lat/lon sample.  Only meaningful
     * after set_equator_ring().
     */
    double check_against_analytic( int nsample ) const;

  private:
    //! Chord distance from p to the segment [a,b], both on the unit sphere.
    static double point_segment_chord( const double* p, const double* a, const double* b );

    //! Chord length to geodesic arc length on the unit sphere, in km.
    static double chord_to_km( double chord )
    {
        // Clamp: round-off can push a coincident point microscopically past 2.
        const double h = std::min( 1.0, 0.5 * chord );
        return EARTH_RADIUS_KM * 2.0 * std::asin( h );
    }

    ErrorCode finalize_from_arrays();
    void best_over_leaves( const double* xyz, const std::vector< EntityHandle >& leaves, double& best ) const;

    Core mCore;                    //!< private instance; the boundary is not part of the working mesh
    mutable AdaptiveKDTree mTree;  //!< distance_search is non-const in the MOAB API
    EntityHandle mRoot;

    std::vector< double > mVerts;              //!< 3*nv boundary vertex coordinates
    std::vector< int > mSegs;                  //!< 2*ns endpoint indices into mVerts
    std::vector< std::vector< int > > mVSegs;  //!< vertex index -> incident segment ids
    std::map< EntityHandle, int > mHandleIdx;  //!< tree entity -> index into mVerts

    double mMaxSegChord;  //!< longest segment; the search radius has to allow for it
    double mMinWidth;     //!< h at distance 0
    BoundaryDensityOptions mOpts;
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline double BoundaryDensity::point_segment_chord( const double* p, const double* a, const double* b )
{
    // Straight-line distance to the chord of the arc.  The arc bulges away
    // from its chord by O(theta^3), utterly negligible for the short segments
    // a digitized coastline is made of, and it never changes which segment is
    // nearest by more than that.
    double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
    double ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
    const double denom = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    double t           = 0.0;
    if( denom > 0.0 ) t = ( ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2] ) / denom;
    if( t < 0.0 ) t = 0.0;
    if( t > 1.0 ) t = 1.0;
    const double dx = ap[0] - t * ab[0];
    const double dy = ap[1] - t * ab[1];
    const double dz = ap[2] - t * ab[2];
    return std::sqrt( dx * dx + dy * dy + dz * dz );
}

inline ErrorCode BoundaryDensity::load( const std::string& filename, const BoundaryDensityOptions& opts )
{
    mOpts = opts;

    ErrorCode rval = mCore.load_file( filename.c_str() );MB_CHK_SET_ERR( rval, "Cannot read boundary file " << filename );

    Range elems;
    rval = mCore.get_entities_by_dimension( 0, 2, elems );MB_CHK_ERR( rval );
    Range edges;
    rval = mCore.get_entities_by_dimension( 0, 1, edges );MB_CHK_ERR( rval );
    if( elems.empty() && edges.empty() )
        MB_SET_ERR( MB_FAILURE, "Boundary file " << filename << " holds no edges or faces to measure distance to" );

    // Collect the unique boundary vertices first, so segments can be stored as
    // indices and each vertex can own a list of the segments touching it.
    Range bverts;
    rval = mCore.get_entities_by_dimension( 0, 0, bverts );MB_CHK_ERR( rval );
    mVerts.resize( 3 * bverts.size() );
    if( !bverts.empty() )
    {
        rval = mCore.get_coords( bverts, &mVerts[0] );MB_CHK_ERR( rval );
    }
    mHandleIdx.clear();
    {
        int i = 0;
        for( Range::iterator it = bverts.begin(); it != bverts.end(); ++it, ++i )
            mHandleIdx[*it] = i;
    }

    mSegs.clear();
    std::vector< EntityHandle > connStorage;
    Range all = elems;
    all.merge( edges );
    for( Range::iterator it = all.begin(); it != all.end(); ++it )
    {
        const EntityHandle* conn = NULL;
        int nnodes               = 0;
        rval                     = mCore.get_connectivity( *it, conn, nnodes, false, &connStorage );MB_CHK_ERR( rval );
        if( nnodes < 2 ) continue;
        // A face's ring closes; a bare edge does not.
        const int nseg = ( mCore.dimension_from_handle( *it ) == 2 ) ? nnodes : nnodes - 1;
        for( int k = 0; k < nseg; ++k )
        {
            std::map< EntityHandle, int >::const_iterator a = mHandleIdx.find( conn[k] );
            std::map< EntityHandle, int >::const_iterator b = mHandleIdx.find( conn[( k + 1 ) % nnodes] );
            if( a == mHandleIdx.end() || b == mHandleIdx.end() ) continue;
            if( a->second == b->second ) continue;  // duplicated point in the ring
            mSegs.push_back( a->second );
            mSegs.push_back( b->second );
        }
    }
    if( mSegs.empty() ) MB_SET_ERR( MB_FAILURE, "Boundary file " << filename << " yielded no usable segments" );

    return finalize_from_arrays();
}

inline ErrorCode BoundaryDensity::set_equator_ring( int nseg, const BoundaryDensityOptions& opts )
{
    mOpts = opts;
    if( nseg < 3 ) MB_SET_ERR( MB_FAILURE, "An equator ring needs at least 3 segments" );

    mVerts.resize( 3 * nseg );
    for( int i = 0; i < nseg; ++i )
    {
        const double t     = 2.0 * M_PI * (double)i / (double)nseg;
        mVerts[3 * i]      = std::cos( t );
        mVerts[3 * i + 1]  = std::sin( t );
        mVerts[3 * i + 2]  = 0.0;
    }
    mSegs.clear();
    for( int i = 0; i < nseg; ++i )
    {
        mSegs.push_back( i );
        mSegs.push_back( ( i + 1 ) % nseg );
    }

    // The tree indexes vertices, so they have to exist as entities.
    mCore.delete_mesh();
    mHandleIdx.clear();
    for( int i = 0; i < nseg; ++i )
    {
        EntityHandle vh;
        ErrorCode rval = mCore.create_vertex( &mVerts[3 * i], vh );MB_CHK_ERR( rval );
        mHandleIdx[vh] = i;
    }
    return finalize_from_arrays();
}

inline ErrorCode BoundaryDensity::finalize_from_arrays()
{
    const size_t nv = mVerts.size() / 3;
    const size_t ns = mSegs.size() / 2;

    mVSegs.assign( nv, std::vector< int >() );
    mMaxSegChord = 0.0;
    for( size_t s = 0; s < ns; ++s )
    {
        const int a = mSegs[2 * s], b = mSegs[2 * s + 1];
        mVSegs[a].push_back( (int)s );
        mVSegs[b].push_back( (int)s );
        const double dx = mVerts[3 * a] - mVerts[3 * b];
        const double dy = mVerts[3 * a + 1] - mVerts[3 * b + 1];
        const double dz = mVerts[3 * a + 2] - mVerts[3 * b + 2];
        mMaxSegChord    = std::max( mMaxSegChord, std::sqrt( dx * dx + dy * dy + dz * dz ) );
    }

    // AdaptiveKDTree::distance_search returns nothing at all when the query
    // point falls outside the tree's bounding box - it checks contains_point
    // and gives up before ever looking at the search radius.  Boundary geometry
    // is a curve, so its box is degenerate (exactly flat for a great circle,
    // a thin shell for coastlines) and almost every point on the sphere is
    // outside it, which made every query return the no-match fallback.
    //
    // Pad the tree with the eight corners of the cube enclosing the unit
    // sphere.  Every possible query point is then inside the box and the
    // traversal runs normally.  The padding vertices are deliberately absent
    // from mHandleIdx, so best_over_leaves() skips them and they can never be
    // mistaken for boundary geometry.
    for( int i = 0; i < 8; ++i )
    {
        const double c[3] = { ( i & 1 ) ? 1.0 : -1.0, ( i & 2 ) ? 1.0 : -1.0, ( i & 4 ) ? 1.0 : -1.0 };
        EntityHandle vh;
        ErrorCode prval = mCore.create_vertex( c, vh );MB_CHK_ERR( prval );
    }

    Range bverts;
    ErrorCode rval = mCore.get_entities_by_dimension( 0, 0, bverts );MB_CHK_ERR( rval );
    mRoot          = 0;
    rval           = mTree.build_tree( bverts, &mRoot );MB_CHK_SET_ERR( rval, "Cannot build the boundary kd-tree" );

    mMinWidth = width_from_distance( 0.0 );
    if( mMinWidth <= 0.0 ) MB_SET_ERR( MB_FAILURE, "Sizing law gives a non-positive cell width at the boundary" );
    if( mOpts.hmax < mMinWidth )
        MB_SET_ERR( MB_FAILURE, "--hmax (" << mOpts.hmax << " km) is below the width at the boundary (" << mMinWidth
                                           << " km); the field would refine away from the boundary instead" );
    return MB_SUCCESS;
}

inline void BoundaryDensity::best_over_leaves( const double* xyz,
                                               const std::vector< EntityHandle >& leaves,
                                               double& best ) const
{
    for( size_t l = 0; l < leaves.size(); ++l )
    {
        Range leafVerts;
        if( MB_SUCCESS != mCore.get_entities_by_handle( leaves[l], leafVerts ) ) continue;
        for( Range::iterator it = leafVerts.begin(); it != leafVerts.end(); ++it )
        {
            std::map< EntityHandle, int >::const_iterator f = mHandleIdx.find( *it );
            if( f == mHandleIdx.end() ) continue;
            const std::vector< int >& segs = mVSegs[f->second];
            for( size_t k = 0; k < segs.size(); ++k )
            {
                const int a = mSegs[2 * segs[k]], b = mSegs[2 * segs[k] + 1];
                best        = std::min( best, point_segment_chord( xyz, &mVerts[3 * a], &mVerts[3 * b] ) );
            }
        }
    }
}

inline double BoundaryDensity::distance_km( const double* xyz ) const
{
    // Expand until some boundary vertex is in range.  Seeding from the longest
    // segment guarantees the first probe is not absurdly small.
    double r = std::max( mMaxSegChord, 1.0e-3 );
    double best = std::numeric_limits< double >::max();
    std::vector< EntityHandle > leaves;

    for( int attempt = 0; attempt < 40 && best == std::numeric_limits< double >::max(); ++attempt )
    {
        leaves.clear();
        if( MB_SUCCESS != mTree.distance_search( xyz, r, leaves ) ) break;
        best_over_leaves( xyz, leaves, best );
        r *= 2.0;
        if( r > 4.0 ) break;  // larger than the sphere's diameter; nothing left to find
    }
    if( best == std::numeric_limits< double >::max() ) return chord_to_km( 2.0 );

    // The first hit is only a candidate: a nearer vertex may sit in a leaf just
    // outside the radius that found it, and a segment can come closer than
    // either endpoint by up to its own length.  One more search over
    // best+mMaxSegChord covers both and makes the result exact.
    leaves.clear();
    if( MB_SUCCESS == mTree.distance_search( xyz, best + mMaxSegChord, leaves ) )
        best_over_leaves( xyz, leaves, best );

    return chord_to_km( best );
}

inline double BoundaryDensity::check_against_analytic( int nsample ) const
{
    double worst = 0.0;
    for( int i = 0; i < nsample; ++i )
    {
        // Deterministic lat/lon lattice, avoiding the poles and the equator
        // itself (where the discretized ring is exactly the analytic answer).
        const double u   = ( i + 0.5 ) / (double)nsample;
        const double lat = -1.4 + 2.8 * u;
        const double lon = 2.0 * M_PI * std::fmod( u * 7.0, 1.0 );
        const double p[3] = { std::cos( lat ) * std::cos( lon ), std::cos( lat ) * std::sin( lon ), std::sin( lat ) };

        const double exact = EARTH_RADIUS_KM * std::fabs( lat );
        worst              = std::max( worst, std::fabs( distance_km( p ) - exact ) );
    }
    return worst;
}

}  // namespace meshopt
}  // namespace moab

#endif  // MOAB_BOUNDARY_DENSITY_HPP
