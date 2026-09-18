/** @example SphericalCVTOptimization.cpp
 * \brief Centroidal Voronoi tessellation on a sphere, then mesh optimization,
 *        with MOAB and Eigen3 only.
 *
 * A spherical centroidal Voronoi tessellation is the mesh family MPAS and the
 * E3SM atmosphere and ocean components are built on, so this doubles as an
 * earth-system workflow:
 *
 *   1. Scatter N generators on the unit sphere, either on a Fibonacci spiral or
 *      pseudo-randomly from a seed.  Given --refine-file, they are instead
 *      rejection-sampled from a sizing field built from real boundary geometry
 *      - coastlines, country outlines, a basin perimeter - so the mesh is fine
 *      near the coast and coarse in the open ocean, which is how production
 *      MPAS meshes are actually built.
 *   2. Triangulate them.  The Delaunay triangulation of points on a sphere is
 *      exactly the 3D convex hull of those points, so an incremental hull gives
 *      it directly; MOAB has no hull utility, so one is built here.
 *   3. Lloyd-iterate to a CVT: move every generator to the centroid of its
 *      spherical Voronoi cell and retriangulate, until the generators stop
 *      moving.  With a sizing field the centroid is density-weighted, which is
 *      what holds the refinement in place as the generators relax.
 *   4. Distribute the triangulation across ranks and run the same L-BFGS shape
 *      optimizer the plane and box example uses, with a spherical constraint
 *      that keeps vertices on the surface and confines the search to the
 *      tangent plane.
 *   5. Optionally write both the triangulation and its dual Voronoi polygonal
 *      mesh, which is the MPAS mesh form.
 *
 * A CVT is already a good mesh, so step 4 is a genuine test: the optimizer has
 * to improve on a strong starting point without destroying it.
 *
 * Run it:
 *     ./SphericalCVTOptimization -N 500 -l 50 -n 100
 *     mpiexec -np 4 ./SphericalCVTOptimization -N 2000 -o sphere.vtk
 *
 * Coastal refinement, writing the MPAS-form polygonal dual:
 *     ./SphericalCVTOptimization -N 20000 -l 60 \
 *         --refine-file world-continents.geo.h5m \
 *         --hmin 30 --hmax 240 --dist 200 --width 150 \
 *         -o tri.h5m --dual mpas.h5m
 *
 * Check the sizing field against its analytic form (needs no data file):
 *     ./SphericalCVTOptimization --check-density
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/verdict/VerdictWrapper.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab_mpi.h"
#endif

#include "MeshOptimizationMetrics.hpp"
#include "MeshOptimizerCore.hpp"
#include "BoundaryDensity.hpp"

using namespace moab;
using namespace moab::meshopt;

// ---------------------------------------------------------------------------
// Small vector helpers.  Three-component arrays rather than Eigen types, to
// keep the geometry code readable and free of temporaries.
// ---------------------------------------------------------------------------

static inline void v_cross( const double* a, const double* b, double* r )
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}
static inline double v_dot( const double* a, const double* b )
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline double v_norm( const double* a )
{
    return std::sqrt( v_dot( a, a ) );
}
static inline void v_normalize( double* a )
{
    const double n = v_norm( a );
    if( n > 0.0 )
    {
        a[0] /= n;
        a[1] /= n;
        a[2] /= n;
    }
}
static inline void v_sub( const double* a, const double* b, double* r )
{
    r[0] = a[0] - b[0];
    r[1] = a[1] - b[1];
    r[2] = a[2] - b[2];
}

// ---------------------------------------------------------------------------
// Spherical constraint
// ---------------------------------------------------------------------------

/** Keeps vertices on a sphere and the search in its tangent plane.
 *
 * Unlike the plane and box case, the constrained direction here rotates as a
 * vertex moves, so the projector has to be rebuilt from the current position
 * every time and the trial position has to be snapped back onto the surface.
 * That also means the L-BFGS history is only approximately valid - the
 * curvature pairs were measured in a tangent space that has since turned.  In
 * practice re-projecting each iteration costs a little of the superlinear
 * convergence near the solution and nothing else.
 */
class SphereConstraint : public Constraint
{
  public:
    SphereConstraint( double radius ) : mRadius( radius ) {}

    void projector( const double* x, double* P ) const override
    {
        double n[3] = { x[0], x[1], x[2] };
        v_normalize( n );
        for( int i = 0; i < 3; ++i )
            for( int j = 0; j < 3; ++j )
                P[3 * i + j] = ( i == j ? 1.0 : 0.0 ) - n[i] * n[j];
    }

    void project_point( double* x ) const override
    {
        const double r = v_norm( x );
        if( r <= 0.0 ) return;
        const double s = mRadius / r;
        x[0] *= s;
        x[1] *= s;
        x[2] *= s;
    }

    bool surface_normal( const double* c, double* n ) const override
    {
        // Outward radial direction at the element centroid.  Hull faces are
        // wound counter-clockwise seen from outside, so this makes the signed
        // area of an untangled triangle positive.
        n[0] = c[0];
        n[1] = c[1];
        n[2] = c[2];
        v_normalize( n );
        return true;
    }

  private:
    double mRadius;
};

// ---------------------------------------------------------------------------
// Generator placement
// ---------------------------------------------------------------------------

/** Fibonacci spiral: a near-uniform, fully deterministic starting point. */
static void fibonacci_sphere( int n, std::vector< double >& pts )
{
    pts.resize( 3 * n );
    const double golden = M_PI * ( 3.0 - std::sqrt( 5.0 ) );
    for( int i = 0; i < n; ++i )
    {
        const double z = 1.0 - 2.0 * ( i + 0.5 ) / n;
        const double r = std::sqrt( std::max( 0.0, 1.0 - z * z ) );
        const double t = golden * i;
        pts[3 * i]     = r * std::cos( t );
        pts[3 * i + 1] = r * std::sin( t );
        pts[3 * i + 2] = z;
    }
}

/** Seeded pseudo-random points, uniform on the sphere.
 *
 * A hash rather than a library RNG so the point set is identical on every rank
 * and every platform, which is what lets each rank build the same hull
 * independently instead of one rank building and broadcasting it.
 */
static void random_sphere( int n, unsigned seed, std::vector< double >& pts )
{
    pts.resize( 3 * n );
    for( int i = 0; i < n; ++i )
    {
        unsigned long h = (unsigned long)i * 6364136223846793005ULL + 1442695040888963407ULL;
        h ^= (unsigned long)seed * 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 33;
        h *= 0xFF51AFD7ED558CCDULL;
        h ^= h >> 33;
        const double u1 = (double)( h % 1000003UL ) / 1000003.0;
        h *= 0xC4CEB9FE1A85EC53ULL;
        h ^= h >> 33;
        const double u2 = (double)( h % 1000003UL ) / 1000003.0;

        const double z   = 2.0 * u1 - 1.0;
        const double r   = std::sqrt( std::max( 0.0, 1.0 - z * z ) );
        const double phi = 2.0 * M_PI * u2;
        pts[3 * i]       = r * std::cos( phi );
        pts[3 * i + 1]   = r * std::sin( phi );
        pts[3 * i + 2]   = z;
    }
}

/** Generators drawn to match a sizing field, by rejection sampling.
 *
 * Candidates come from the same hash stream random_sphere() uses, so the
 * result is still identical on every rank and every platform.  A candidate at
 * x survives with probability (hmin/h(x))^2: rejection sampling controls the
 * point density directly, and a cell of width h covers area ~h^2, so the
 * density we want is 1/h^2.  (The Lloyd weight needs a different power - see
 * BoundaryDensity::lloyd_weight.)
 *
 * Starting from the graded distribution rather than a uniform one matters:
 * density-weighted Lloyd moves generators only a fraction of a cell per sweep,
 * so reaching a strong grading from a uniform start would take far more sweeps
 * than the default budget.
 */
static void density_sphere( int n, unsigned seed, const BoundaryDensity& dens, std::vector< double >& pts )
{
    pts.clear();
    pts.reserve( 3 * n );

    // Unbounded in principle; in practice the acceptance rate is (hmin/hmax)^2
    // at worst, so this cap is many times what is needed and exists only so a
    // pathological sizing field cannot spin forever.
    const unsigned long maxDraws = 4000UL * (unsigned long)n + 100000UL;

    for( unsigned long i = 0; i < maxDraws && (int)( pts.size() / 3 ) < n; ++i )
    {
        unsigned long h = i * 6364136223846793005ULL + 1442695040888963407ULL;
        h ^= (unsigned long)seed * 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 33;
        h *= 0xFF51AFD7ED558CCDULL;
        h ^= h >> 33;
        const double u1 = (double)( h % 1000003UL ) / 1000003.0;
        h *= 0xC4CEB9FE1A85EC53ULL;
        h ^= h >> 33;
        const double u2 = (double)( h % 1000003UL ) / 1000003.0;
        h *= 0x9E3779B97F4A7C15ULL;
        h ^= h >> 31;
        const double u3 = (double)( h % 1000003UL ) / 1000003.0;

        const double z   = 2.0 * u1 - 1.0;
        const double r   = std::sqrt( std::max( 0.0, 1.0 - z * z ) );
        const double phi = 2.0 * M_PI * u2;
        const double c[3] = { r * std::cos( phi ), r * std::sin( phi ), z };

        if( u3 > dens.accept_probability( c ) ) continue;
        pts.push_back( c[0] );
        pts.push_back( c[1] );
        pts.push_back( c[2] );
    }
}

// ---------------------------------------------------------------------------
// Convex hull == spherical Delaunay triangulation
// ---------------------------------------------------------------------------

struct HullFace
{
    int v[3];
};

/** Incremental 3D convex hull of points on a sphere.
 *
 * Every point on a sphere is a hull vertex, so this needs no interior test and
 * no point ordering heuristics: insert each point, delete the faces it can
 * see, and stitch a fan from the horizon back to it.  Faces come out wound
 * counter-clockwise viewed from outside.
 *
 * O(n^2) in the worst case, which is fine at the scale an example runs at.
 * Exactly cospherical degeneracies are not handled; the two generators offered
 * here do not produce them.
 */
static bool convex_hull( const std::vector< double >& pts, std::vector< HullFace >& faces )
{
    const int n = (int)pts.size() / 3;
    faces.clear();
    if( n < 4 ) return false;

    // Seed tetrahedron: the first four points that enclose a non-zero volume.
    int seed[4] = { 0, 1, -1, -1 };
    for( int k = 2; k < n && seed[2] < 0; ++k )
    {
        double e1[3], e2[3], cr[3];
        v_sub( &pts[3 * k], &pts[0], e1 );
        v_sub( &pts[3], &pts[0], e2 );
        v_cross( e1, e2, cr );
        if( v_norm( cr ) > 1.0e-12 ) seed[2] = k;
    }
    if( seed[2] < 0 ) return false;
    {
        double e1[3], e2[3], e3[3], cr[3];
        v_sub( &pts[3 * seed[1]], &pts[0], e1 );
        v_sub( &pts[3 * seed[2]], &pts[0], e2 );
        v_cross( e1, e2, cr );
        for( int k = 2; k < n && seed[3] < 0; ++k )
        {
            if( k == seed[2] ) continue;
            v_sub( &pts[3 * k], &pts[0], e3 );
            if( std::fabs( v_dot( cr, e3 ) ) > 1.0e-12 ) seed[3] = k;
        }
    }
    if( seed[3] < 0 ) return false;

    // Orient each seed face outward, away from the tetrahedron's centroid.
    double centre[3] = { 0, 0, 0 };
    for( int i = 0; i < 4; ++i )
        for( int d = 0; d < 3; ++d )
            centre[d] += pts[3 * seed[i] + d] / 4.0;

    const int combos[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 } };
    for( int f = 0; f < 4; ++f )
    {
        HullFace hf;
        for( int j = 0; j < 3; ++j )
            hf.v[j] = seed[combos[f][j]];
        double e1[3], e2[3], nrm[3], toFace[3];
        v_sub( &pts[3 * hf.v[1]], &pts[3 * hf.v[0]], e1 );
        v_sub( &pts[3 * hf.v[2]], &pts[3 * hf.v[0]], e2 );
        v_cross( e1, e2, nrm );
        v_sub( &pts[3 * hf.v[0]], centre, toFace );
        if( v_dot( nrm, toFace ) < 0.0 ) std::swap( hf.v[1], hf.v[2] );
        faces.push_back( hf );
    }

    std::vector< bool > inHull( n, false );
    for( int i = 0; i < 4; ++i )
        inHull[seed[i]] = true;

    std::vector< char > visible;
    for( int p = 0; p < n; ++p )
    {
        if( inHull[p] ) continue;

        visible.assign( faces.size(), 0 );
        bool any = false;
        for( size_t f = 0; f < faces.size(); ++f )
        {
            double e1[3], e2[3], nrm[3], toP[3];
            v_sub( &pts[3 * faces[f].v[1]], &pts[3 * faces[f].v[0]], e1 );
            v_sub( &pts[3 * faces[f].v[2]], &pts[3 * faces[f].v[0]], e2 );
            v_cross( e1, e2, nrm );
            v_sub( &pts[3 * p], &pts[3 * faces[f].v[0]], toP );
            if( v_dot( nrm, toP ) > 1.0e-14 )
            {
                visible[f] = 1;
                any        = true;
            }
        }
        if( !any ) continue;  // numerically inside; skip rather than corrupt the hull

        // A directed edge of a visible face is on the horizon when its reverse
        // does not also belong to a visible face.
        std::set< std::pair< int, int > > visEdges;
        for( size_t f = 0; f < faces.size(); ++f )
        {
            if( !visible[f] ) continue;
            for( int j = 0; j < 3; ++j )
                visEdges.insert( std::make_pair( faces[f].v[j], faces[f].v[( j + 1 ) % 3] ) );
        }

        std::vector< std::pair< int, int > > horizon;
        for( std::set< std::pair< int, int > >::iterator e = visEdges.begin(); e != visEdges.end(); ++e )
            if( !visEdges.count( std::make_pair( e->second, e->first ) ) ) horizon.push_back( *e );

        std::vector< HullFace > kept;
        kept.reserve( faces.size() );
        for( size_t f = 0; f < faces.size(); ++f )
            if( !visible[f] ) kept.push_back( faces[f] );
        faces.swap( kept );

        // The new face inherits the visible face's edge direction, which leaves
        // it wound outward like everything else.
        for( size_t e = 0; e < horizon.size(); ++e )
        {
            HullFace hf;
            hf.v[0] = horizon[e].first;
            hf.v[1] = horizon[e].second;
            hf.v[2] = p;
            faces.push_back( hf );
        }
        inHull[p] = true;
    }

    return !faces.empty();
}

// ---------------------------------------------------------------------------
// Voronoi dual and Lloyd iteration
// ---------------------------------------------------------------------------

/** Circumcentre of each spherical triangle, on the unit sphere.
 *
 * For points on a sphere the circumcentre of the spherical triangle is just the
 * normalized face normal, so the Voronoi diagram falls out of the hull with no
 * extra geometry.
 */
static void face_circumcentres( const std::vector< double >& pts,
                                const std::vector< HullFace >& faces,
                                std::vector< double >& cc )
{
    cc.resize( 3 * faces.size() );
    for( size_t f = 0; f < faces.size(); ++f )
    {
        double e1[3], e2[3], nrm[3];
        v_sub( &pts[3 * faces[f].v[1]], &pts[3 * faces[f].v[0]], e1 );
        v_sub( &pts[3 * faces[f].v[2]], &pts[3 * faces[f].v[0]], e2 );
        v_cross( e1, e2, nrm );
        v_normalize( nrm );
        for( int d = 0; d < 3; ++d )
            cc[3 * f + d] = nrm[d];
    }
}

/** Faces incident on each generator. */
static void build_vertex_faces( int nverts,
                                const std::vector< HullFace >& faces,
                                std::vector< std::vector< int > >& vf )
{
    vf.assign( nverts, std::vector< int >() );
    for( size_t f = 0; f < faces.size(); ++f )
        for( int j = 0; j < 3; ++j )
            vf[faces[f].v[j]].push_back( (int)f );
}

/** Voronoi cell of generator i: its incident circumcentres in angular order.
 *
 * Ordering matters - the polygon is fan-triangulated afterwards, and an
 * unordered vertex list would give a self-intersecting polygon with a
 * meaningless area and centroid.
 */
static void ordered_cell( const std::vector< double >& pts,
                          const std::vector< double >& cc,
                          const std::vector< int >& incident,
                          int i,
                          std::vector< int >& order )
{
    double n[3] = { pts[3 * i], pts[3 * i + 1], pts[3 * i + 2] };
    v_normalize( n );

    // Any pair of axes spanning the tangent plane will do; take the smallest
    // component of n to build a vector that cannot be parallel to it.
    double t1[3] = { 0, 0, 0 };
    int smallest = 0;
    for( int d = 1; d < 3; ++d )
        if( std::fabs( n[d] ) < std::fabs( n[smallest] ) ) smallest = d;
    t1[smallest] = 1.0;
    double t2[3];
    v_cross( n, t1, t2 );
    v_normalize( t2 );
    v_cross( n, t2, t1 );

    std::vector< std::pair< double, int > > byAngle;
    byAngle.reserve( incident.size() );
    for( size_t k = 0; k < incident.size(); ++k )
    {
        const double* c = &cc[3 * incident[k]];
        byAngle.push_back( std::make_pair( std::atan2( v_dot( c, t2 ), v_dot( c, t1 ) ), incident[k] ) );
    }
    std::sort( byAngle.begin(), byAngle.end() );

    order.clear();
    for( size_t k = 0; k < byAngle.size(); ++k )
        order.push_back( byAngle[k].second );
}

/** One Lloyd sweep.  Returns the largest distance any generator moved.
 *
 * Each generator goes to the area-weighted centroid of its Voronoi polygon,
 * projected back onto the sphere.  Iterating this to a fixed point is the
 * definition of a centroidal Voronoi tessellation.
 *
 * With a sizing field the weight becomes area*rho instead of area, which is
 * what grades the mesh; dens = NULL reproduces the uniform arithmetic exactly,
 * down to the last bit, so an unrefined run is unaffected by this parameter.
 */
static double lloyd_sweep( std::vector< double >& pts,
                           const std::vector< HullFace >& faces,
                           const BoundaryDensity* dens )
{
    const int n = (int)pts.size() / 3;
    std::vector< double > cc;
    face_circumcentres( pts, faces, cc );
    std::vector< std::vector< int > > vf;
    build_vertex_faces( n, faces, vf );

    std::vector< double > next( pts );
    double maxmove = 0.0;
    std::vector< int > order;

    for( int i = 0; i < n; ++i )
    {
        if( vf[i].size() < 3 ) continue;
        ordered_cell( pts, cc, vf[i], i, order );

        // Fan-triangulate the cell and take the area-weighted centroid.
        double acc[3] = { 0, 0, 0 };
        double area   = 0.0;
        for( size_t k = 1; k + 1 < order.size(); ++k )
        {
            const double* p0 = &cc[3 * order[0]];
            const double* p1 = &cc[3 * order[k]];
            const double* p2 = &cc[3 * order[k + 1]];
            double e1[3], e2[3], cr[3];
            v_sub( p1, p0, e1 );
            v_sub( p2, p0, e2 );
            v_cross( e1, e2, cr );
            double a = 0.5 * v_norm( cr );
            if( dens )
            {
                // Sample the field at the fan-triangle centroid, pushed back out
                // to the sphere: the circumcentres are on the surface but their
                // average is not, and the sizing field is only defined there.
                double s[3] = { ( p0[0] + p1[0] + p2[0] ) / 3.0, ( p0[1] + p1[1] + p2[1] ) / 3.0,
                                ( p0[2] + p1[2] + p2[2] ) / 3.0 };
                v_normalize( s );
                a *= dens->lloyd_weight( s );
            }
            area += a;
            // Left exactly as it was: a*(x)/3 and a*(x/3) do not round alike, so
            // hoisting the centroid out of here would perturb an unrefined run.
            for( int d = 0; d < 3; ++d )
                acc[d] += a * ( p0[d] + p1[d] + p2[d] ) / 3.0;
        }
        if( area <= 0.0 ) continue;

        for( int d = 0; d < 3; ++d )
            acc[d] /= area;
        v_normalize( acc );

        double diff[3];
        v_sub( acc, &pts[3 * i], diff );
        maxmove = std::max( maxmove, v_norm( diff ) );
        for( int d = 0; d < 3; ++d )
            next[3 * i + d] = acc[d];
    }

    pts.swap( next );
    return maxmove;
}

/** Ratio of largest to smallest Voronoi cell area.
 *
 * A CVT on a sphere has near-uniform cells, so this heads towards 1 and is a
 * cheap scalar summary of how converged the tessellation is.
 */
static double cell_area_ratio( const std::vector< double >& pts, const std::vector< HullFace >& faces )
{
    const int n = (int)pts.size() / 3;
    std::vector< double > cc;
    face_circumcentres( pts, faces, cc );
    std::vector< std::vector< int > > vf;
    build_vertex_faces( n, faces, vf );

    double amin = 1.0e300, amax = 0.0;
    std::vector< int > order;
    for( int i = 0; i < n; ++i )
    {
        if( vf[i].size() < 3 ) continue;
        ordered_cell( pts, cc, vf[i], i, order );
        double area = 0.0;
        for( size_t k = 1; k + 1 < order.size(); ++k )
        {
            double e1[3], e2[3], cr[3];
            v_sub( &cc[3 * order[k]], &cc[3 * order[0]], e1 );
            v_sub( &cc[3 * order[k + 1]], &cc[3 * order[0]], e2 );
            v_cross( e1, e2, cr );
            area += 0.5 * v_norm( cr );
        }
        if( area > 0.0 )
        {
            amin = std::min( amin, area );
            amax = std::max( amax, area );
        }
    }
    return ( amin < 1.0e299 && amin > 0.0 ) ? amax / amin : 0.0;
}

// ---------------------------------------------------------------------------
// MOAB mesh construction
// ---------------------------------------------------------------------------

/** Create the triangulation, distributed across ranks.
 *
 * Every rank has already built the same hull from the same deterministic point
 * set, so instead of one rank building and scattering it, each rank simply
 * creates the slice of triangles it owns.  Vertices are tagged with their
 * global index and resolve_shared_ents matches them up across the partition
 * boundary, which is what makes the optimizer's reductions correct.
 */
static ErrorCode build_triangulation( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                      ParallelComm* pcomm,
#endif
                                      const std::vector< double >& pts,
                                      const std::vector< HullFace >& faces,
                                      int rank,
                                      int nprocs,
                                      Range& tris )
{
    ErrorCode rval;
    std::map< int, EntityHandle > vmap;
    tris.clear();

    for( size_t f = 0; f < faces.size(); ++f )
    {
        if( (int)( f % nprocs ) != rank ) continue;
        EntityHandle conn[3];
        for( int j = 0; j < 3; ++j )
        {
            const int gi                            = faces[f].v[j];
            std::map< int, EntityHandle >::iterator m = vmap.find( gi );
            if( m == vmap.end() )
            {
                EntityHandle vh;
                rval = mb->create_vertex( &pts[3 * gi], vh );MB_CHK_ERR( rval );
                vmap[gi] = vh;
                conn[j]  = vh;
            }
            else
                conn[j] = m->second;
        }
        EntityHandle tri;
        rval = mb->create_element( MBTRI, conn, 3, tri );MB_CHK_ERR( rval );
        tris.insert( tri );
    }

    Tag gidTag;
    rval = mb->tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gidTag, MB_TAG_DENSE | MB_TAG_CREAT );MB_CHK_ERR( rval );
    for( std::map< int, EntityHandle >::iterator m = vmap.begin(); m != vmap.end(); ++m )
    {
        const int gid = m->first + 1;  // GLOBAL_ID is 1-based
        rval          = mb->tag_set_data( gidTag, &m->second, 1, &gid );MB_CHK_ERR( rval );
    }

#ifdef MOAB_HAVE_MPI
    if( pcomm && nprocs > 1 )
    {
        // Pass the triangles explicitly.  The overload that takes this_set = 0
        // instead collects its entities from ParallelComm's registered partition
        // sets, and a mesh built in memory like this one has none - so it
        // resolved an empty range, reported success, and left every rank with
        // private copies of the partition-boundary vertices.
        rval = pcomm->resolve_shared_ents( 0, tris, 2, 0, NULL, &gidTag );MB_CHK_ERR( rval );

        // Verify rather than trust: a torn partition is silent.  The optimizer
        // still runs, but each rank relaxes its own copies of the shared
        // vertices, so the triangles decouple and converge to a perfect - and
        // meaningless - mesh that disagrees with the serial answer.  After a
        // correct resolve every vertex is owned by exactly one rank, so the
        // owned count must sum to the generator count.
        Range verts;
        rval = mb->get_entities_by_dimension( 0, 0, verts );MB_CHK_ERR( rval );
        Range ownedVerts;
        rval = pcomm->filter_pstatus( verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ownedVerts );MB_CHK_ERR( rval );
        long localOwned = (long)ownedVerts.size(), globalOwned = 0;
        MPI_Allreduce( &localOwned, &globalOwned, 1, MPI_LONG, MPI_SUM, pcomm->comm() );
        const long expected = (long)( pts.size() / 3 );
        if( globalOwned != expected )
        {
            MB_SET_ERR( MB_FAILURE, "resolve_shared_ents left the partition torn: "
                                        << globalOwned << " owned vertices across all ranks, expected " << expected
                                        << ".  The parallel result would not match serial." );
        }
    }
#endif
    return MB_SUCCESS;
}

/** Create the dual Voronoi polygonal mesh: the MPAS mesh form.
 *
 * Built only on rank 0 and only when an output file is requested, since it
 * exists to be looked at rather than optimized.
 */
static ErrorCode build_voronoi_mesh( Interface* mb,
                                     const std::vector< double >& pts,
                                     const std::vector< HullFace >& faces,
                                     EntityHandle& outSet )
{
    ErrorCode rval;
    const int n = (int)pts.size() / 3;

    std::vector< double > cc;
    face_circumcentres( pts, faces, cc );
    std::vector< std::vector< int > > vf;
    build_vertex_faces( n, faces, vf );

    rval = mb->create_meshset( MESHSET_SET, outSet );MB_CHK_ERR( rval );

    std::vector< EntityHandle > ccHandles( faces.size() );
    for( size_t f = 0; f < faces.size(); ++f )
    {
        rval = mb->create_vertex( &cc[3 * f], ccHandles[f] );MB_CHK_ERR( rval );
    }

    std::vector< int > order;
    std::vector< EntityHandle > conn;
    for( int i = 0; i < n; ++i )
    {
        if( vf[i].size() < 3 ) continue;
        ordered_cell( pts, cc, vf[i], i, order );
        conn.clear();
        for( size_t k = 0; k < order.size(); ++k )
            conn.push_back( ccHandles[order[k]] );
        EntityHandle poly;
        rval = mb->create_element( MBPOLYGON, &conn[0], (int)conn.size(), poly );MB_CHK_ERR( rval );
        rval = mb->add_entities( outSet, &poly, 1 );MB_CHK_ERR( rval );
    }
    return MB_SUCCESS;
}

/** Collect the optimized generator positions back into a flat array on rank 0.
 *
 * The generators in `pts` are the Lloyd output; the optimizer then moved them,
 * but it moved them in the distributed MOAB mesh, not in this array.  Writing
 * the dual straight from `pts` - which is what used to happen - therefore
 * discarded every bit of the optimization.  Read the coordinates back instead,
 * indexed by GLOBAL_ID.
 *
 * Only owned vertices are contributed, so no vertex is sent twice, and the
 * union across ranks is the whole generator set.
 */
static ErrorCode gather_optimized_points( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                          ParallelComm* pcomm,
#endif
                                          const Range& tris,
                                          int nprocs,
                                          int rank,
                                          std::vector< double >& pts )
{
    ErrorCode rval;
    Range verts;
    rval = mb->get_connectivity( tris, verts );MB_CHK_ERR( rval );

    Range owned = verts;
#ifdef MOAB_HAVE_MPI
    if( nprocs > 1 )
    {
        owned.clear();
        rval = pcomm->filter_pstatus( verts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned );MB_CHK_ERR( rval );
    }
#else
    (void)nprocs;
    (void)rank;
#endif

    Tag gidTag;
    rval = mb->tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gidTag );MB_CHK_ERR( rval );

    const int nloc = (int)owned.size();
    std::vector< int > gids( nloc );
    std::vector< double > coords( 3 * nloc );
    if( nloc )
    {
        rval = mb->tag_get_data( gidTag, owned, &gids[0] );MB_CHK_ERR( rval );
        rval = mb->get_coords( owned, &coords[0] );MB_CHK_ERR( rval );
    }

    const int npoints = (int)( pts.size() / 3 );

#ifdef MOAB_HAVE_MPI
    if( nprocs > 1 )
    {
        std::vector< int > counts( nprocs, 0 );
        MPI_Gather( &nloc, 1, MPI_INT, counts.empty() ? NULL : &counts[0], 1, MPI_INT, 0, pcomm->comm() );

        std::vector< int > gdisp( nprocs, 0 ), cCounts( nprocs, 0 ), cDisp( nprocs, 0 );
        int total = 0;
        if( !rank )
        {
            for( int p = 0; p < nprocs; ++p )
            {
                gdisp[p]   = total;
                cDisp[p]   = 3 * total;
                cCounts[p] = 3 * counts[p];
                total += counts[p];
            }
        }
        std::vector< int > allGids( rank ? 0 : total );
        std::vector< double > allCoords( rank ? 0 : 3 * total );

        MPI_Gatherv( nloc ? &gids[0] : NULL, nloc, MPI_INT, rank ? NULL : ( total ? &allGids[0] : NULL ),
                     rank ? NULL : &counts[0], rank ? NULL : &gdisp[0], MPI_INT, 0, pcomm->comm() );
        MPI_Gatherv( nloc ? &coords[0] : NULL, 3 * nloc, MPI_DOUBLE,
                     rank ? NULL : ( total ? &allCoords[0] : NULL ), rank ? NULL : &cCounts[0],
                     rank ? NULL : &cDisp[0], MPI_DOUBLE, 0, pcomm->comm() );

        if( !rank )
        {
            if( total != npoints )
                MB_SET_ERR( MB_FAILURE, "gathered " << total << " owned generators but expected " << npoints );
            for( int i = 0; i < total; ++i )
            {
                const int idx = allGids[i] - 1;  // GLOBAL_ID is 1-based
                if( idx < 0 || idx >= npoints )
                    MB_SET_ERR( MB_FAILURE, "GLOBAL_ID " << allGids[i] << " is outside 1.." << npoints );
                for( int d = 0; d < 3; ++d )
                    pts[3 * idx + d] = allCoords[3 * i + d];
            }
        }
        return MB_SUCCESS;
    }
#endif

    if( nloc != npoints )
        MB_SET_ERR( MB_FAILURE, "found " << nloc << " generators in the mesh but expected " << npoints );
    for( int i = 0; i < nloc; ++i )
    {
        const int idx = gids[i] - 1;
        if( idx < 0 || idx >= npoints )
            MB_SET_ERR( MB_FAILURE, "GLOBAL_ID " << gids[i] << " is outside 1.." << npoints );
        for( int d = 0; d < 3; ++d )
            pts[3 * idx + d] = coords[3 * i + d];
    }
    return MB_SUCCESS;
}

// ---------------------------------------------------------------------------
// Quality reporting
// ---------------------------------------------------------------------------

struct QualityStats
{
    double minq, maxq, mean;
    long count;
};

static ErrorCode report_quality( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                 ParallelComm* pcomm,
#endif
                                 const Range& elems,
                                 QualityStats& s )
{
    VerdictWrapper vw( mb );
    double localMin = 1.0e300, localMax = -1.0e300, localSum = 0.0;
    long localCount = 0;
    std::vector< EntityHandle > connStorage;

    for( Range::iterator it = elems.begin(); it != elems.end(); ++it )
    {
        const EntityType type = mb->type_from_handle( *it );
        const int nnodes      = metric_num_nodes( type );
        if( !nnodes ) continue;
        const EntityHandle* conn = NULL;
        int connlen              = 0;
        if( MB_SUCCESS != mb->get_connectivity( *it, conn, connlen, false, &connStorage ) ) continue;
        double c[24];
        if( MB_SUCCESS != mb->get_coords( conn, nnodes, c ) ) continue;

        double q = 0.0;
        if( MB_SUCCESS != vw.quality_measure( *it, MB_SHAPE, q, nnodes, type, c ) || !std::isfinite( q ) ) continue;
        localMin = std::min( localMin, q );
        localMax = std::max( localMax, q );
        localSum += q;
        ++localCount;
    }

    s.minq  = localMin;
    s.maxq  = localMax;
    s.mean  = localSum;
    s.count = localCount;
#ifdef MOAB_HAVE_MPI
    if( pcomm && pcomm->size() > 1 )
    {
        MPI_Allreduce( &localMin, &s.minq, 1, MPI_DOUBLE, MPI_MIN, pcomm->comm() );
        MPI_Allreduce( &localMax, &s.maxq, 1, MPI_DOUBLE, MPI_MAX, pcomm->comm() );
        MPI_Allreduce( &localSum, &s.mean, 1, MPI_DOUBLE, MPI_SUM, pcomm->comm() );
        MPI_Allreduce( &localCount, &s.count, 1, MPI_LONG, MPI_SUM, pcomm->comm() );
    }
#endif
    if( s.count > 0 ) s.mean /= (double)s.count;
    return MB_SUCCESS;
}

// ---------------------------------------------------------------------------
// Sizing-field self-test
// ---------------------------------------------------------------------------

/** Validate the sizing field without needing a geometry file.
 *
 * The boundary is the equator, whose distance function is known exactly
 * (R*|asin(z)|), so both halves can be checked against an analytic answer:
 *
 *   1. the kd-tree distance query, against R*|lat|;
 *   2. the sizing law, at distances where tanh is known;
 *   3. the exponent in lloyd_weight, which is the part most easily got wrong.
 *      A graded CVT ends up with cell areas spanning (hmax/hmin)^2.  Weighting
 *      by 1/h^2 instead of 1/h^4 would produce a span of only (hmax/hmin),
 *      and the two predictions are far enough apart to tell apart even from a
 *      short Lloyd run.
 *
 * Returns 0 on success.
 */
static int run_density_selftest( const BoundaryDensityOptions& dopts, int verbosity )
{
    BoundaryDensity dens;
    ErrorCode rval = dens.set_equator_ring( 2000, dopts );
    if( MB_SUCCESS != rval )
    {
        std::cerr << "check-density: could not build the equator ring" << std::endl;
        return 1;
    }

    int failures = 0;
    std::printf( "\nSizing-field self-test (analytic equator boundary)\n" );
    std::printf( "  segments    %lu\n", (unsigned long)dens.num_segments() );

    // 1. distance query -----------------------------------------------------
    // A 2000-segment ring is a polygon inscribed in the equator, so it sits
    // inside the true circle by R*(1-cos(pi/n)) = 7.9e-3 km.  Allow a hair more.
    const double ringErr = EARTH_RADIUS_KM * ( 1.0 - std::cos( M_PI / 2000.0 ) );
    const double tolDist = 10.0 * ringErr + 1.0e-6;
    const double worst   = dens.check_against_analytic( 997 );
    std::printf( "  distance    max error %.6e km (tolerance %.6e)  %s\n", worst, tolDist,
                 worst <= tolDist ? "ok" : "FAILED" );
    if( worst > tolDist ) ++failures;

    // 2. sizing law ---------------------------------------------------------
    // At d = d0 the tanh vanishes and the width is exactly the midpoint.
    const double mid      = dens.width_from_distance( dopts.d0 );
    const double expected = 0.5 * ( dopts.hmin + dopts.hmax );
    std::printf( "  law h(d0)   %.9f km (expected %.9f)             %s\n", mid, expected,
                 std::fabs( mid - expected ) < 1.0e-9 ? "ok" : "FAILED" );
    if( std::fabs( mid - expected ) >= 1.0e-9 ) ++failures;

    // Monotone increasing away from the boundary, and bounded by hmin/hmax.
    bool monotone = true, bounded = true;
    double prev = -1.0;
    for( int i = 0; i <= 400; ++i )
    {
        const double d = 40.0 * i;  // out to 16000 km, past the antipode
        const double h = dens.width_from_distance( d );
        if( h < prev ) monotone = false;
        if( h < dopts.hmin || h > dopts.hmax ) bounded = false;
        prev = h;
    }
    std::printf( "  law shape   monotone %s, within [hmin,hmax] %s\n", monotone ? "yes" : "NO", bounded ? "yes" : "NO" );
    if( !monotone || !bounded ) ++failures;

    std::printf( "  attained    h(0) = %.6f km, hmax = %.6f km (ratio %.4f)\n", dens.min_cell_width_km(),
                 dens.max_cell_width_km(), dens.max_cell_width_km() / dens.min_cell_width_km() );

    // 3. the exponent -------------------------------------------------------
    const int n = 3000;
    std::vector< double > pts;
    density_sphere( n, 12345, dens, pts );
    if( (int)( pts.size() / 3 ) < 4 )
    {
        std::cerr << "check-density: rejection sampling produced too few generators" << std::endl;
        return failures + 1;
    }

    std::vector< HullFace > faces;
    if( !convex_hull( pts, faces ) )
    {
        std::cerr << "check-density: hull failed on the graded point set" << std::endl;
        return failures + 1;
    }
    for( int s = 0; s < 40; ++s )
    {
        const double moved = lloyd_sweep( pts, faces, &dens );
        if( !convex_hull( pts, faces ) ) break;
        if( moved < 1.0e-8 ) break;
    }

    const double widthRatio = dens.max_cell_width_km() / dens.min_cell_width_km();
    const double predicted  = widthRatio * widthRatio;  // areas go as width^2
    const double wrongPower = widthRatio;               // what rho = 1/h^2 would give
    const double measured   = cell_area_ratio( pts, faces );

    // A short Lloyd run on a finite point set lands short of the asymptotic
    // prediction, so an absolute tolerance here would just be a tuned magic
    // number.  Ask the question that actually matters instead: of the two
    // candidate exponents, which one is the measurement nearer to?  Ratios are
    // compared in log space, where "twice as big" and "half as big" are the
    // same distance, so neither prediction is favoured by the metric.
    const double lm = std::log( measured ), l4 = std::log( predicted ), l2 = std::log( wrongPower );
    const double to4 = std::fabs( lm - l4 ), to2 = std::fabs( lm - l2 );
    const bool ok    = ( to4 < to2 );
    std::printf( "  grading     cell area max/min %.4f\n", measured );
    std::printf( "              rho=1/h^4 predicts %.4f (log distance %.4f)\n", predicted, to4 );
    std::printf( "              rho=1/h^2 would give %.4f (log distance %.4f)\n", wrongPower, to2 );
    std::printf( "              nearer the 1/h^4 prediction  %s\n", ok ? "ok" : "FAILED" );
    if( !ok ) ++failures;
    // Guard the guard: if the two predictions were close together the test
    // above would be meaningless, and that happens whenever hmax/hmin is small.
    if( std::fabs( l4 - l2 ) < 0.5 )
    {
        std::printf( "              INCONCLUSIVE: the predictions are too close to tell apart; "
                     "raise --hmax/--hmin\n" );
        ++failures;
    }

    // Euler characteristic: the grading must not have broken the triangulation.
    const int nv = (int)( pts.size() / 3 ), nf = (int)faces.size(), ne = 3 * nf / 2;
    const int euler = nv - ne + nf;
    std::printf( "  topology    V-E+F = %d (expected 2)                        %s\n", euler,
                 euler == 2 ? "ok" : "FAILED" );
    if( euler != 2 ) ++failures;

    if( verbosity > 1 )
        std::printf( "              %d generators, %d triangles\n", nv, nf );

    std::printf( "\n  %s\n\n", failures ? "SELF-TEST FAILED" : "all sizing-field checks passed" );
    return failures ? 1 : 0;
}

// ---------------------------------------------------------------------------

int main( int argc, char** argv )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    int npoints = 500, nlloyd = 40, maxiter = 100, verbosity = 1;
    double lloydTol = 1.0e-8, delta = 0.0, pnorm = 2.0, gtol = 1.0e-10;
    unsigned seed = 12345;
    bool useRandom = false;
    bool checkDensity = false;
    std::string metricName = "imr", outfile, dualfile, refineFile;
    BoundaryDensityOptions dopts;

    ProgOptions opts( "Spherical centroidal Voronoi tessellation followed by parallel mesh optimization, "
                      "using MOAB and Eigen3 only." );
    opts.addOpt< int >( std::string( "npoints,N" ), std::string( "Number of generators (default=500)" ), &npoints );
    opts.addOpt< int >( std::string( "lloyd,l" ), std::string( "Maximum Lloyd sweeps towards the CVT (default=40)" ),
                        &nlloyd );
    opts.addOpt< double >( std::string( "lloyd-tol" ),
                           std::string( "Stop Lloyd when the largest generator movement falls below this "
                                        "(default=1e-8)" ),
                           &lloydTol );
    opts.addOpt< void >( std::string( "random,r" ),
                         std::string( "Start from seeded random points instead of a Fibonacci spiral.  Random "
                                      "generators need far more Lloyd sweeps to reach a CVT" ),
                         &useRandom );
    opts.addOpt< int >( std::string( "seed,s" ), std::string( "Seed for --random (default=12345)" ), (int*)&seed );
    opts.addOpt< std::string >( std::string( "metric,m" ), std::string( "Shape metric: imr or condition "
                                                                        "(default=imr; identical for triangles)" ),
                                &metricName );
    opts.addOpt< int >( std::string( "niter,n" ), std::string( "Maximum L-BFGS iterations (default=100)" ), &maxiter );
    opts.addOpt< double >( std::string( "gtol,e" ), std::string( "Gradient norm convergence tolerance (default=1e-10)" ),
                           &gtol );
    opts.addOpt< double >( std::string( "pnorm,p" ), std::string( "Objective is the p-mean of the metric (default=2)" ),
                           &pnorm );
    opts.addOpt< double >( std::string( "delta" ), std::string( "Escobar regularization (default=0)" ), &delta );
    opts.addOpt< std::string >( std::string( "output,o" ), std::string( "Write the optimized triangulation here" ),
                                &outfile );
    opts.addOpt< std::string >( std::string( "dual" ),
                                std::string( "Write the dual Voronoi polygonal mesh here (the MPAS mesh form)" ),
                                &dualfile );
    opts.addOpt< std::string >( std::string( "refine-file" ),
                                std::string( "Boundary geometry (coastlines, basin outlines) on the unit sphere; "
                                             "generators are graded towards it" ),
                                &refineFile );
    opts.addOpt< double >( std::string( "hmin" ),
                           std::string( "Cell width approached at the boundary, km (default=30)" ), &dopts.hmin );
    opts.addOpt< double >( std::string( "hmax" ), std::string( "Far-field cell width, km (default=240)" ),
                           &dopts.hmax );
    opts.addOpt< double >( std::string( "dist" ),
                           std::string( "Distance from the boundary at which the width is halfway between hmin and "
                                        "hmax, km (default=200)" ),
                           &dopts.d0 );
    opts.addOpt< double >( std::string( "width" ), std::string( "Transition length of the tanh ramp, km "
                                                                "(default=150)" ),
                           &dopts.lt );
    opts.addOpt< void >( std::string( "check-density" ),
                         std::string( "Verify the sizing field against an analytic equator boundary and exit; "
                                      "needs no data file" ),
                         &checkDensity );
    opts.addOpt< int >( std::string( "verbose,v" ), std::string( "0 quiet, 1 normal, 2 per-iteration (default=1)" ),
                        &verbosity );
    opts.parseCommandLine( argc, argv );

    if( npoints < 4 )
    {
        std::cerr << "need at least 4 generators to triangulate a sphere" << std::endl;
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    Core moab;
    Interface* mb = &moab;
#ifdef MOAB_HAVE_MPI
    ParallelComm* pcomm = new ParallelComm( mb, MPI_COMM_WORLD );
    const int rank      = pcomm->rank();
    const int nprocs    = pcomm->size();
#else
    const int rank   = 0;
    const int nprocs = 1;
#endif

    if( checkDensity )
    {
        const int rc = rank ? 0 : run_density_selftest( dopts, verbosity );
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return rc;
    }

    // ---- 0. sizing field ---------------------------------------------------
    // Loaded redundantly on every rank, like the generators below: the Lloyd
    // phase is replicated, so every rank has to be able to evaluate the field.
    BoundaryDensity density;
    BoundaryDensity* densp = NULL;
    if( !refineFile.empty() )
    {
        ErrorCode drval = density.load( refineFile, dopts );
        if( MB_SUCCESS != drval )
        {
            if( !rank ) std::cerr << "could not load the refinement boundary " << refineFile << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return 1;
        }
        densp = &density;
    }

    if( !rank && verbosity > 0 )
    {
        std::printf( "\nSpherical CVT and mesh optimization (MOAB + Eigen3)\n" );
        std::printf( "  generators  %d (%s)\n", npoints,
                     densp ? "density-weighted rejection sampling"
                           : ( useRandom ? "seeded random" : "Fibonacci spiral" ) );
        std::printf( "  ranks       %d\n", nprocs );
        if( densp )
        {
            std::printf( "  refinement  %s (%lu boundary segments)\n", refineFile.c_str(),
                         (unsigned long)density.num_segments() );
            std::printf( "  cell width  %.2f km at the boundary -> %.2f km far field (d0 %.1f km, Lt %.1f km)\n",
                         density.min_cell_width_km(), density.max_cell_width_km(), dopts.d0, dopts.lt );
            if( density.min_cell_width_km() > 1.01 * dopts.hmin )
                std::printf( "              note: hmin=%.2f km is only approached, not reached; the tanh gives "
                             "%.2f km at d=0\n",
                             dopts.hmin, density.min_cell_width_km() );

            // The sizing law fixes the *shape* of the grading; -N fixes the
            // resolution.  Ask for a transition finer than a cell and the mesh
            // simply cannot represent it - the realized grading then falls well
            // short of hmax/hmin, and falls further the longer Lloyd runs,
            // because each sweep averages the density over cells too coarse to
            // resolve it.  Measured on real coastlines with hmax/h(0) = 5.50:
            // N=4000 (357 km cells) reached 2.86, N=16000 (179 km cells)
            // reached 3.76.  Say so rather than let it look like a target met.
            const double meanCell = std::sqrt( 4.0 * M_PI / (double)npoints ) * EARTH_RADIUS_KM;
            std::printf( "  resolution  %d generators give ~%.0f km mean cell width\n", npoints, meanCell );
            if( meanCell > dopts.lt )
            {
                const double nres = 4.0 * M_PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM / ( dopts.lt * dopts.lt );
                std::printf( "              WARNING: coarser than the %.0f km transition, so the grading will "
                             "fall short of\n"
                             "              the %.2fx the law asks for.  Resolving the transition needs "
                             "about %.0f generators.\n",
                             dopts.lt, dopts.hmax / density.min_cell_width_km(), nres );
            }
            if( meanCell > 2.0 * density.min_cell_width_km() )
            {
                const double nabs = 4.0 * M_PI * EARTH_RADIUS_KM * EARTH_RADIUS_KM /
                                    ( density.min_cell_width_km() * density.min_cell_width_km() );
                std::printf( "              note: the absolute widths in the law are out of reach at this -N; "
                             "%.0f km cells\n"
                             "              everywhere would need about %.0f generators.  The grading is still "
                             "applied, relative to\n"
                             "              the resolution -N does afford.\n",
                             density.min_cell_width_km(), nabs );
            }
        }
    }

    // ---- 1. generators -----------------------------------------------------
    std::vector< double > pts;
    if( densp )
    {
        density_sphere( npoints, seed, density, pts );
        const int got = (int)( pts.size() / 3 );
        if( got < 4 )
        {
            if( !rank )
                std::cerr << "rejection sampling produced only " << got << " generators; widen the sizing field or "
                          << "raise -N" << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return 1;
        }
        if( got < npoints && !rank && verbosity > 0 )
            std::printf( "  sampled     %d of %d requested (draw cap reached)\n", got, npoints );
        npoints = got;
    }
    else if( useRandom )
        random_sphere( npoints, seed, pts );
    else
        fibonacci_sphere( npoints, pts );

    // ---- 2 and 3. Delaunay + Lloyd to a CVT --------------------------------
    // Run redundantly on every rank.  The point set and the hull are both
    // deterministic, so all ranks reach bit-identical generators without a
    // single message; broadcasting would cost more than recomputing.
    std::vector< HullFace > faces;
    if( !convex_hull( pts, faces ) )
    {
        std::cerr << "convex hull failed; the generators may be degenerate" << std::endl;
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if( !rank && verbosity > 0 )
        std::printf( "  triangles   %lu (Euler check: expected %d)\n", (unsigned long)faces.size(), 2 * npoints - 4 );

    int sweeps    = 0;
    double moved  = 0.0;
    for( ; sweeps < nlloyd; ++sweeps )
    {
        moved = lloyd_sweep( pts, faces, densp );
        // Generators have moved, so the Delaunay triangulation has to be rebuilt
        // before the next sweep; the Voronoi diagram is its dual.
        if( !convex_hull( pts, faces ) ) break;
        if( verbosity > 1 && !rank ) std::printf( "  lloyd %4d  max movement %.6e\n", sweeps + 1, moved );
        if( moved < lloydTol ) break;
    }

    if( !rank && verbosity > 0 )
    {
        std::printf( "  lloyd       %d sweeps, final max movement %.6e%s\n", sweeps, moved,
                     ( sweeps >= nlloyd && moved >= lloydTol ) ? "  (sweep limit reached; raise -l for a truer CVT)"
                                                               : "" );
        std::printf( "  cell areas  max/min ratio %.6f\n", cell_area_ratio( pts, faces ) );
    }

    // ---- 4. distribute and optimize ----------------------------------------
    Range tris;
    ErrorCode rval = build_triangulation( mb,
#ifdef MOAB_HAVE_MPI
                                          pcomm,
#endif
                                          pts, faces, rank, nprocs, tris );MB_CHK_ERR( rval );

    Range ownedTris = tris;
#ifdef MOAB_HAVE_MPI
    if( nprocs > 1 )
    {
        ownedTris.clear();
        rval = pcomm->filter_pstatus( tris, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ownedTris );MB_CHK_ERR( rval );
    }
#endif

    SphereConstraint constraint( 1.0 );

    QualityStats before, after;
    rval = report_quality( mb,
#ifdef MOAB_HAVE_MPI
                           pcomm,
#endif
                           ownedTris, before );MB_CHK_ERR( rval );
    if( !rank && verbosity > 0 )
        std::printf( "\n  CVT mesh    MB_SHAPE  min %10.6f   mean %10.6f   max %10.6f   (%ld triangles)\n", before.minq,
                     before.mean, before.maxq, before.count );

    OptimizerOptions oo;
    oo.metric    = ( metricName == "condition" ? CONDITION_NUMBER : INVERSE_MEAN_RATIO );
    oo.delta     = delta;
    oo.pnorm     = pnorm;
    oo.maxiter   = maxiter;
    oo.gtol      = gtol;
    oo.verbosity = verbosity;

    MeshOptimizer optimizer( mb,
#ifdef MOAB_HAVE_MPI
                             pcomm,
#endif
                             oo, constraint );
    rval = optimizer.setup( ownedTris );MB_CHK_ERR( rval );

    if( !rank && verbosity > 0 )
        std::printf( "\nOptimizing %ld triangles on the sphere, metric %s\n", optimizer.num_global_elements(),
                     metric_name( oo.metric ) );

    OptimizerResult res;
    rval = optimizer.optimize( res );MB_CHK_ERR( rval );

    rval = report_quality( mb,
#ifdef MOAB_HAVE_MPI
                           pcomm,
#endif
                           ownedTris, after );MB_CHK_ERR( rval );

    if( !rank && verbosity > 0 )
    {
        std::printf( "\n  optimized   MB_SHAPE  min %10.6f   mean %10.6f   max %10.6f   (%ld triangles)\n", after.minq,
                     after.mean, after.maxq, after.count );
        std::printf( "\n  objective   %.12e -> %.12e\n", res.f_initial, res.f_final );
        std::printf( "  gradient    %.6e -> %.6e\n", res.gnorm_initial, res.gnorm_final );
        std::printf( "  iterations  %d (%d objective evaluations)%s\n", res.iterations, res.evaluations,
                     res.converged ? ", converged" : "" );
    }

    // Vertices must still be on the sphere; a drift here would mean the
    // projection is not being applied on every accepted step.
    {
        Range verts;
        rval = mb->get_connectivity( ownedTris, verts );MB_CHK_ERR( rval );
        std::vector< double > c( 3 * verts.size() );
        double localDev = 0.0;
        if( !verts.empty() )
        {
            rval = mb->get_coords( verts, &c[0] );MB_CHK_ERR( rval );
            for( size_t i = 0; i < verts.size(); ++i )
                localDev = std::max( localDev, std::fabs( v_norm( &c[3 * i] ) - 1.0 ) );
        }
        double dev = localDev;
#ifdef MOAB_HAVE_MPI
        if( nprocs > 1 ) MPI_Allreduce( &localDev, &dev, 1, MPI_DOUBLE, MPI_MAX, pcomm->comm() );
#endif
        if( !rank && verbosity > 0 ) std::printf( "  on-sphere   max radial deviation %.6e\n", dev );
    }

    // ---- 5. output ---------------------------------------------------------
    if( !outfile.empty() )
    {
        const char* wopts = ( nprocs > 1 ? "PARALLEL=WRITE_PART" : "" );
        rval              = mb->write_file( outfile.c_str(), 0, wopts );MB_CHK_ERR( rval );
        if( !rank && verbosity > 0 ) std::printf( "  wrote       %s\n", outfile.c_str() );
    }

    if( !dualfile.empty() )
    {
        // Refresh pts from the optimized mesh before taking the dual, gathering
        // to rank 0 when distributed.  `faces` is still the exact connectivity -
        // the optimizer moves vertices and never retriangulates - so the result
        // is the dual of the mesh that was just written, not of a re-derived
        // Delaunay triangulation of the same points.
        rval = gather_optimized_points( mb,
#ifdef MOAB_HAVE_MPI
                                        pcomm,
#endif
                                        tris, nprocs, rank, pts );MB_CHK_ERR( rval );

        if( !rank )
        {
            Core dualMoab;
            EntityHandle dualSet = 0;
            rval                 = build_voronoi_mesh( &dualMoab, pts, faces, dualSet );MB_CHK_ERR( rval );
            rval                 = dualMoab.write_file( dualfile.c_str(), 0, 0, &dualSet, 1 );MB_CHK_ERR( rval );
            if( verbosity > 0 ) std::printf( "  wrote       %s (dual Voronoi mesh)\n", dualfile.c_str() );
        }
    }

    if( !rank && verbosity > 0 ) std::printf( "\n" );

#ifdef MOAB_HAVE_MPI
    delete pcomm;
    MPI_Finalize();
#endif
    return 0;
}
