/** @example MeshOptimization.cpp
 * \brief Parallel mesh shape optimization with MOAB and Eigen3 only.
 *
 * Demonstrates a complete mesh optimization workflow on the two analytic
 * domains MOAB can build for itself: a rectangular plane of quadrilaterals and
 * a box of hexahedra.
 *
 *   1. ScdInterface::construct_box builds a distributed, sharing-resolved mesh
 *      with no input file and no HDF5, at any rank count.
 *   2. Interior vertices are displaced by a seeded pseudo-random field, which
 *      wrecks the element shapes in a reproducible way.
 *   3. L-BFGS on a differentiable shape metric (condition number or inverse
 *      mean ratio) pulls them back, with gradients from Eigen's AutoDiff
 *      module rather than hand-differentiated formulas.
 *   4. Verdict reports element quality before and after.
 *
 * Because the generated lattice is the exact minimizer of both metrics on
 * these domains, the example knows the answer it should reach: quality returns
 * to 1 and the vertices return to the lattice.  That is what --selftest
 * checks, alongside a finite-difference check of the AutoDiff gradient and a
 * cross-check of the templated metric against Verdict.
 *
 * Boundary vertices slide.  Pinning them caps the achievable quality on coarse
 * meshes, so a vertex on a box face keeps two degrees of freedom, one on an
 * edge keeps one, and only the eight corners are fixed.
 *
 * Run it:
 *     ./MeshOptimization -d 2 -N 20 -n 100
 *     mpiexec -np 4 ./MeshOptimization -d 3 -N 8 -m condition
 *     ./MeshOptimization --selftest
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/verdict/VerdictWrapper.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab_mpi.h"
#endif

#include "MeshOptimizationMetrics.hpp"
#include "MeshOptimizerCore.hpp"

using namespace moab;
using namespace moab::meshopt;

// ---------------------------------------------------------------------------
// Constraint for an axis-aligned plane or box.
// ---------------------------------------------------------------------------

/** Sliding-boundary constraint for an axis-aligned rectangular domain.
 *
 * Every constrained direction here is a coordinate axis, so the projector is
 * diagonal: a vertex sitting on the x = xmin plane simply loses its x degree of
 * freedom.  Counting how many planes a vertex lies on recovers the usual
 * classification without any topological queries - one plane is a face, two is
 * an edge, three is a corner.
 *
 * The constraint is computed from the domain, not from skinning the mesh.  In
 * parallel, Skinner reports partition-interface facets as boundary, and the
 * workarounds for that are fragile.  Reading it off the coordinates instead is
 * purely local and therefore gives every rank the same answer with no
 * communication at all.
 */
class BoxConstraint : public Constraint
{
  public:
    BoxConstraint( const double lo[3], const double hi[3], double tol, bool slide )
        : mTol( tol ), mSlide( slide )
    {
        for( int a = 0; a < 3; ++a )
        {
            mLo[a] = lo[a];
            mHi[a] = hi[a];
            // A degenerate axis means the mesh is a surface lying in that
            // coordinate plane; every vertex is confined to it.
            mDegenerate[a] = ( hi[a] - lo[a] ) <= tol;
        }
    }

    void projector( const double* x, double* P ) const override
    {
        bool constrained[3] = { mDegenerate[0], mDegenerate[1], mDegenerate[2] };
        bool onBoundary     = false;
        for( int a = 0; a < 3; ++a )
        {
            if( mDegenerate[a] ) continue;
            if( std::fabs( x[a] - mLo[a] ) < mTol || std::fabs( x[a] - mHi[a] ) < mTol )
            {
                constrained[a] = true;
                onBoundary     = true;
            }
        }
        if( !mSlide && onBoundary ) constrained[0] = constrained[1] = constrained[2] = true;

        for( int i = 0; i < 9; ++i )
            P[i] = 0.0;
        P[0] = constrained[0] ? 0.0 : 1.0;
        P[4] = constrained[1] ? 0.0 : 1.0;
        P[8] = constrained[2] ? 0.0 : 1.0;
    }

    bool surface_normal( const double*, double* n ) const override
    {
        for( int a = 0; a < 3; ++a )
        {
            if( !mDegenerate[a] ) continue;
            n[0] = n[1] = n[2] = 0.0;
            n[a]               = 1.0;
            return true;
        }
        return false;
    }

  private:
    double mLo[3], mHi[3], mTol;
    bool mDegenerate[3];
    bool mSlide;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Deterministic hash-based pseudo-random field in [-1, 1].
 *
 * A hash of the lattice index rather than a sequential RNG, so a vertex gets
 * the same displacement no matter which rank owns it or what order it is
 * visited in.  Without that the perturbed mesh would differ between rank
 * counts and the determinism check below would be meaningless.
 */
static double hash_random( long key, int component, unsigned seed )
{
    unsigned long h = (unsigned long)( key * 6364136223846793005ULL + 1442695040888963407ULL );
    h ^= (unsigned long)( component + 1 ) * 0x9E3779B97F4A7C15ULL;
    h ^= (unsigned long)seed * 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 33;
    h *= 0xFF51AFD7ED558CCDULL;
    h ^= h >> 33;
    h *= 0xC4CEB9FE1A85EC53ULL;
    h ^= h >> 33;
    return 2.0 * ( (double)( h % 1000000UL ) / 1000000.0 ) - 1.0;
}

struct QualityStats
{
    double minq, maxq, mean;
    long count;
};

/** Verdict quality over the owned elements, reduced across ranks.
 *
 * Deliberately independent of MeshOptimizationMetrics.hpp: it is the outside
 * opinion that says whether the optimizer actually improved the mesh.
 */
static ErrorCode report_quality( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                 ParallelComm* pcomm,
#endif
                                 const Range& elems,
                                 QualityType qt,
                                 QualityStats& s )
{
    VerdictWrapper vw( mb );
    double localMin = 1.0e300, localMax = -1.0e300, localSum = 0.0;
    long localCount = 0;
    std::vector< EntityHandle > connStorage;

    for( Range::iterator it = elems.begin(); it != elems.end(); ++it )
    {
        // Hand Verdict the coordinates rather than letting it look up
        // connectivity itself: on a structured mesh that lookup fails, and
        // quality_measure would then quietly return an error for every element.
        const EntityType type = mb->type_from_handle( *it );
        const int nnodes      = metric_num_nodes( type );
        if( !nnodes ) continue;
        const EntityHandle* conn = NULL;
        int connlen              = 0;
        if( MB_SUCCESS != mb->get_connectivity( *it, conn, connlen, false, &connStorage ) ) continue;
        double c[24];
        if( MB_SUCCESS != mb->get_coords( conn, nnodes, c ) ) continue;

        double q       = 0.0;
        ErrorCode rval = vw.quality_measure( *it, qt, q, nnodes, type, c );
        if( MB_SUCCESS != rval || !std::isfinite( q ) ) continue;
        if( q < localMin ) localMin = q;
        if( q > localMax ) localMax = q;
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

static void print_quality( const char* label, const char* metric, const QualityStats& s )
{
    std::printf( "  %-10s %-18s min %10.6f   mean %10.6f   max %10.6f   (%ld elements)\n", label, metric, s.minq,
                 s.mean, s.maxq, s.count );
}

// ---------------------------------------------------------------------------
// Self-test
// ---------------------------------------------------------------------------

static int g_failures = 0;

static void expect( const char* what, double got, double want, double tol )
{
    if( !( std::fabs( got - want ) <= tol ) )
    {
        std::printf( "  FAIL %-34s got %.12g want %.12g (tol %.2g)\n", what, got, want, tol );
        ++g_failures;
    }
    else
        std::printf( "  ok   %-34s %.12g\n", what, got );
}

/** Closed-form checks on the metric layer.
 *
 * These need no mesh and no MPI, and they pin down values that can be worked
 * out by hand, so a regression in the Jacobian or the ideal-element weight
 * shows up immediately rather than as a slightly worse optimization result.
 */
static void selftest_metrics()
{
    const double s3     = std::sqrt( 3.0 );
    const double nz[3]  = { 0.0, 0.0, 1.0 };
    const double tri[9] = { 0, 0, 0, 1, 0, 0, 0.5, s3 / 2, 0 };
    const double quad[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
    const double tet[12]  = { 0, 0, 0, 1, 0, 0, 0.5, s3 / 2, 0, 0.5, s3 / 6, std::sqrt( 2.0 / 3.0 ) };
    const double hex[24]  = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };

    std::printf( "\nMetric layer, closed-form values:\n" );
    expect( "equilateral tri", element_quality< double >( MBTRI, tri, nz, INVERSE_MEAN_RATIO, 0.0 ), 1.0, 1e-14 );
    expect( "unit quad", element_quality< double >( MBQUAD, quad, nz, INVERSE_MEAN_RATIO, 0.0 ), 1.0, 1e-14 );
    expect( "regular tet, imr", element_quality< double >( MBTET, tet, 0, INVERSE_MEAN_RATIO, 0.0 ), 1.0, 1e-14 );
    expect( "regular tet, condition", element_quality< double >( MBTET, tet, 0, CONDITION_NUMBER, 0.0 ), 1.0, 1e-14 );
    expect( "unit hex, imr", element_quality< double >( MBHEX, hex, 0, INVERSE_MEAN_RATIO, 0.0 ), 1.0, 1e-14 );
    expect( "unit hex, condition", element_quality< double >( MBHEX, hex, 0, CONDITION_NUMBER, 0.0 ), 1.0, 1e-14 );

    // T = diag(2,1) gives |T|_F^2 / (2 det) = 5/4 exactly.
    const double q2[12] = { 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 1, 0 };
    expect( "2:1 quad", element_quality< double >( MBQUAD, q2, nz, INVERSE_MEAN_RATIO, 0.0 ), 1.25, 1e-14 );

    // Shape metrics are scale invariant.
    const double q5[12] = { 0, 0, 0, 5, 0, 0, 5, 5, 0, 0, 5, 0 };
    expect( "scaled square", element_quality< double >( MBQUAD, q5, nz, INVERSE_MEAN_RATIO, 0.0 ), 1.0, 1e-14 );

    // An inverted element must be a barrier when unregularized and finite once
    // regularized; that is the whole point of the Escobar determinant.
    const double qi[12] = { 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0 };
    const double raw    = element_quality< double >( MBQUAD, qi, nz, INVERSE_MEAN_RATIO, 0.0 );
    const double reg    = element_quality< double >( MBQUAD, qi, nz, INVERSE_MEAN_RATIO, 0.1 );
    std::printf( "  ok   %-34s delta=0 -> %g, delta=0.1 -> %g\n", "inverted quad barrier", raw, reg );
    if( std::isfinite( raw ) )
    {
        std::printf( "  FAIL inverted quad is not a barrier at delta=0\n" );
        ++g_failures;
    }
    if( !std::isfinite( reg ) )
    {
        std::printf( "  FAIL regularized metric is not finite on an inverted element\n" );
        ++g_failures;
    }
}

/** Compare the templated condition number against Verdict element by element.
 *
 * Only meaningful at delta == 0 on an uninverted mesh, where the regularized
 * determinant reduces exactly to the determinant and the two are the same
 * function.  Verdict's MB_CONDITION is the independent implementation, so
 * agreement here says the Jacobian and ideal weight in the header are right.
 */
static ErrorCode selftest_against_verdict( Interface* mb, const Range& elems, const Constraint& con )
{
    VerdictWrapper vw( mb );
    double maxrel = 0.0;
    long compared = 0;
    std::vector< EntityHandle > connStorage;

    for( Range::iterator it = elems.begin(); it != elems.end(); ++it )
    {
        const EntityType type = mb->type_from_handle( *it );
        const int nnodes      = metric_num_nodes( type );
        if( !nnodes ) continue;

        // Structured-mesh elements need the storage vector; see MeshOptimizerCore.
        const EntityHandle* conn = NULL;
        int connlen              = 0;
        ErrorCode rval = mb->get_connectivity( *it, conn, connlen, false, &connStorage );MB_CHK_ERR( rval );
        double c[24];
        rval = mb->get_coords( conn, nnodes, c );MB_CHK_ERR( rval );

        double normal[3]  = { 0, 0, 1 };
        const double* np  = NULL;
        if( metric_dimension( type ) == 2 )
        {
            double centroid[3] = { 0, 0, 0 };
            for( int i = 0; i < nnodes; ++i )
                for( int d = 0; d < 3; ++d )
                    centroid[d] += c[3 * i + d] / nnodes;
            con.surface_normal( centroid, normal );
            np = normal;
        }

        const double mine = element_quality< double >( type, c, np, CONDITION_NUMBER, 0.0 );
        double theirs     = 0.0;
        if( MB_SUCCESS != vw.quality_measure( *it, MB_CONDITION, theirs, nnodes, type, c ) ) continue;
        if( !std::isfinite( mine ) || !std::isfinite( theirs ) ) continue;

        const double rel = std::fabs( mine - theirs ) / std::max( 1.0, std::fabs( theirs ) );
        if( rel > maxrel ) maxrel = rel;
        ++compared;
    }

    std::printf( "\nMetric vs Verdict MB_CONDITION over %ld elements:\n", compared );
    if( compared == 0 )
    {
        // A cross-check that compared nothing passes vacuously, which is worse
        // than failing: it hides exactly the kind of silent skip this guards.
        std::printf( "  FAIL no elements were compared\n" );
        ++g_failures;
        return MB_SUCCESS;
    }
    expect( "max relative difference", maxrel, 0.0, 1.0e-10 );
    return MB_SUCCESS;
}

// ---------------------------------------------------------------------------
// Mesh generation
// ---------------------------------------------------------------------------

/** Build a distributed structured quad (dim 2) or hex (dim 3) mesh.
 *
 * construct_box takes the parallel partition path only when the local corners
 * are equal, in which case it derives them from par.gDims and the partition
 * method; hence the HomCoord(0,0,0) pair.  With a NULL coordinate array the
 * vertices land on the integer lattice, which is the ideal mesh for both
 * metrics and so the answer the optimizer should recover.
 */
static ErrorCode generate_mesh( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                ParallelComm* pcomm,
#endif
                                int dim,
                                int n,
                                Range& elems )
{
    ScdInterface* scd = NULL;
    ErrorCode rval    = mb->query_interface( scd );MB_CHK_ERR( rval );

    ScdParData par;
    par.gDims[0] = par.gDims[1] = par.gDims[2] = 0;
    par.gDims[3]                               = n;
    par.gDims[4]                               = n;
    par.gDims[5]                               = ( dim == 3 ? n : 0 );
    par.gPeriodic[0] = par.gPeriodic[1] = par.gPeriodic[2] = 0;
#ifdef MOAB_HAVE_MPI
    par.pComm      = pcomm;
    par.partMethod = ( dim == 3 ? ScdParData::SQIJK : ScdParData::SQIJ );
    if( !pcomm || pcomm->size() == 1 ) par.partMethod = ScdParData::NOPART;
#else
    par.partMethod = ScdParData::NOPART;
#endif

    ScdBox* box = NULL;
    if( ScdParData::NOPART == par.partMethod )
    {
        rval = scd->construct_box( HomCoord( 0, 0, 0 ), HomCoord( n, n, dim == 3 ? n : 0 ), NULL, 0, box, NULL, NULL,
                                   true, -1 );MB_CHK_ERR( rval );
    }
    else
    {
        rval = scd->construct_box( HomCoord( 0, 0, 0 ), HomCoord( 0, 0, 0 ), NULL, 0, box, NULL, &par, true,
                                   dim );MB_CHK_ERR( rval );
    }

    elems.clear();
    rval = mb->get_entities_by_dimension( 0, dim, elems );MB_CHK_ERR( rval );
    return MB_SUCCESS;
}

/** Global bounding box of a vertex range. */
static ErrorCode global_bounds( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                ParallelComm* pcomm,
#endif
                                const Range& verts,
                                double lo[3],
                                double hi[3] )
{
    double llo[3] = { 1e300, 1e300, 1e300 }, lhi[3] = { -1e300, -1e300, -1e300 };
    std::vector< double > c( 3 * verts.size() );
    if( !verts.empty() )
    {
        ErrorCode rval = mb->get_coords( verts, &c[0] );MB_CHK_ERR( rval );
        for( size_t i = 0; i < verts.size(); ++i )
            for( int d = 0; d < 3; ++d )
            {
                if( c[3 * i + d] < llo[d] ) llo[d] = c[3 * i + d];
                if( c[3 * i + d] > lhi[d] ) lhi[d] = c[3 * i + d];
            }
    }
    for( int d = 0; d < 3; ++d )
    {
        lo[d] = llo[d];
        hi[d] = lhi[d];
    }
#ifdef MOAB_HAVE_MPI
    if( pcomm && pcomm->size() > 1 )
    {
        MPI_Allreduce( llo, lo, 3, MPI_DOUBLE, MPI_MIN, pcomm->comm() );
        MPI_Allreduce( lhi, hi, 3, MPI_DOUBLE, MPI_MAX, pcomm->comm() );
    }
#endif
    return MB_SUCCESS;
}

/** Displace vertices inside their allowed subspace, reproducibly.
 *
 * The displacement is projected through the same constraint the optimizer
 * uses, so a boundary vertex slides along its face or edge and a corner does
 * not move.  Perturbing across the constraint instead would move the domain
 * itself, and the lattice would no longer be the answer.
 *
 * Keyed on the global id so the perturbed mesh is identical at every rank
 * count.  That is what makes the serial and parallel objective histories
 * comparable.
 */
static ErrorCode perturb_mesh( Interface* mb, const Range& verts, const Constraint& con, double amplitude,
                               unsigned seed )
{
    Tag gid;
    ErrorCode rval = mb->tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid );MB_CHK_ERR( rval );

    std::vector< int > ids( verts.size() );
    rval = mb->tag_get_data( gid, verts, &ids[0] );MB_CHK_ERR( rval );
    std::vector< double > c( 3 * verts.size() );
    rval = mb->get_coords( verts, &c[0] );MB_CHK_ERR( rval );

    for( size_t i = 0; i < verts.size(); ++i )
    {
        double P[9];
        con.projector( &c[3 * i], P );
        double v[3];
        for( int d = 0; d < 3; ++d )
            v[d] = amplitude * hash_random( ids[i], d, seed );
        const double p0 = P[0] * v[0] + P[1] * v[1] + P[2] * v[2];
        const double p1 = P[3] * v[0] + P[4] * v[1] + P[5] * v[2];
        const double p2 = P[6] * v[0] + P[7] * v[1] + P[8] * v[2];
        c[3 * i] += p0;
        c[3 * i + 1] += p1;
        c[3 * i + 2] += p2;
    }
    return mb->set_coords( verts, &c[0] );
}

/** Largest distance any vertex still sits from the integer lattice. */
static ErrorCode lattice_error( Interface* mb,
#ifdef MOAB_HAVE_MPI
                                ParallelComm* pcomm,
#endif
                                const Range& verts,
                                double& err )
{
    std::vector< double > c( 3 * verts.size() );
    double local = 0.0;
    if( !verts.empty() )
    {
        ErrorCode rval = mb->get_coords( verts, &c[0] );MB_CHK_ERR( rval );
        for( size_t i = 0; i < verts.size(); ++i )
        {
            double d2 = 0.0;
            for( int d = 0; d < 3; ++d )
            {
                const double delta = c[3 * i + d] - std::floor( c[3 * i + d] + 0.5 );
                d2 += delta * delta;
            }
            const double dist = std::sqrt( d2 );
            if( dist > local ) local = dist;
        }
    }
    err = local;
#ifdef MOAB_HAVE_MPI
    if( pcomm && pcomm->size() > 1 ) MPI_Allreduce( &local, &err, 1, MPI_DOUBLE, MPI_MAX, pcomm->comm() );
#endif
    return MB_SUCCESS;
}

// ---------------------------------------------------------------------------

int main( int argc, char** argv )
{
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    int dim = 2, nelem = 10, maxiter = 100, verbosity = 1;
    // Negative amplitude means "pick one that does not tangle the mesh".  The
    // same displacement that merely distorts a quad will invert a hex, because
    // in 3D it is applied along three axes at once.
    double amplitude = -1.0, delta = 0.0, pnorm = 2.0, gtol = 1.0e-10;
    unsigned seed  = 12345;
    bool selftest  = false;
    bool noslide   = false;
    std::string metricName = "imr";
    std::string infile, outfile;

    ProgOptions opts( "Parallel mesh shape optimization on a plane or a box, using MOAB and Eigen3 only." );
    opts.addOpt< int >( std::string( "dim,d" ), std::string( "Domain dimension: 2 for a quad plane, 3 for a hex box "
                                                             "(default=2)" ),
                        &dim );
    opts.addOpt< int >( std::string( "nelem,N" ),
                        std::string( "Elements per side of the generated mesh (default=10)" ), &nelem );
    opts.addOpt< std::string >( std::string( "metric,m" ),
                                std::string( "Shape metric: imr or condition (default=imr).  They are the same "
                                             "function in 2D and differ only for tet and hex" ),
                                &metricName );
    opts.addOpt< int >( std::string( "niter,n" ), std::string( "Maximum L-BFGS iterations (default=100)" ), &maxiter );
    opts.addOpt< double >( std::string( "gtol,e" ),
                           std::string( "Convergence tolerance on the projected gradient norm (default=1e-10)" ),
                           &gtol );
    opts.addOpt< double >( std::string( "pnorm,p" ), std::string( "Objective is the p-mean of the metric (default=2)" ),
                           &pnorm );
    opts.addOpt< double >( std::string( "delta" ),
                           std::string( "Escobar regularization; 0 gives the pure metric, a positive value lets the "
                                        "optimizer pass through inverted elements (default=0)" ),
                           &delta );
    opts.addOpt< double >( std::string( "amplitude,a" ),
                           std::string( "Perturbation amplitude as a fraction of the element size (default: 0.3 in "
                                        "2D, 0.2 in 3D, which distort without tangling).  A larger value in 3D will "
                                        "invert elements, which then needs --delta to recover from" ),
                           &amplitude );
    opts.addOpt< int >( std::string( "seed,s" ), std::string( "Perturbation seed (default=12345)" ), (int*)&seed );
    opts.addOpt< void >( std::string( "no-slide" ),
                         std::string( "Pin boundary vertices instead of letting them slide along faces and edges" ),
                         &noslide );
    opts.addOpt< void >( std::string( "selftest" ),
                         std::string( "Run the closed-form metric checks, the Verdict cross-check and the "
                                      "finite-difference gradient check, then exit" ),
                         &selftest );
    opts.addOpt< std::string >( std::string( "file,f" ),
                                std::string( "Read a mesh instead of generating one.  The sliding constraint is "
                                             "derived from the global bounding box, so the domain must be an "
                                             "axis-aligned rectangle or box" ),
                                &infile );
    opts.addOpt< std::string >( std::string( "output,o" ), std::string( "Write the optimized mesh here" ), &outfile );
    opts.addOpt< int >( std::string( "verbose,v" ), std::string( "0 quiet, 1 normal, 2 per-iteration (default=1)" ),
                        &verbosity );
    opts.parseCommandLine( argc, argv );

    if( amplitude < 0.0 ) amplitude = ( dim == 3 ? 0.2 : 0.3 );

    if( dim != 2 && dim != 3 )
    {
        std::cerr << "dimension must be 2 or 3" << std::endl;
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

    ErrorCode rval;
    Range elems;

    if( infile.empty() )
    {
        rval = generate_mesh( mb,
#ifdef MOAB_HAVE_MPI
                              pcomm,
#endif
                              dim, nelem, elems );MB_CHK_ERR( rval );
    }
    else
    {
        std::string ropts;
        if( nprocs > 1 )
            ropts = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_"
                    "DISTRIBUTE";
        rval = mb->load_file( infile.c_str(), 0, ropts.c_str() );MB_CHK_ERR( rval );
        elems.clear();
        rval = mb->get_entities_by_dimension( 0, dim, elems );MB_CHK_ERR( rval );
    }

    // Only owned elements contribute to the objective; counting a ghosted
    // element on both its ranks would double its weight and bias the gradient.
    Range ownedElems = elems;
#ifdef MOAB_HAVE_MPI
    if( nprocs > 1 )
    {
        ownedElems.clear();
        rval = pcomm->filter_pstatus( elems, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ownedElems );MB_CHK_ERR( rval );
    }
#endif

    Range verts;
    rval = mb->get_connectivity( ownedElems, verts );MB_CHK_ERR( rval );

    double lo[3], hi[3];
    rval = global_bounds( mb,
#ifdef MOAB_HAVE_MPI
                          pcomm,
#endif
                          verts, lo, hi );MB_CHK_ERR( rval );

    // Tolerance for "is this vertex on a bounding plane": a fraction of the
    // element size, which for the generated lattice is 1.
    double extent = 0.0;
    for( int d = 0; d < 3; ++d )
        extent = std::max( extent, hi[d] - lo[d] );
    const double h   = ( infile.empty() ? 1.0 : extent / std::max( 1, nelem ) );
    const double tol = 1.0e-4 * h;

    BoxConstraint constraint( lo, hi, tol, !noslide );

    if( !rank && verbosity > 0 )
    {
        std::printf( "\nMOAB + Eigen3 mesh shape optimization\n" );
        std::printf( "  domain      [%g,%g] x [%g,%g] x [%g,%g]\n", lo[0], hi[0], lo[1], hi[1], lo[2], hi[2] );
        std::printf( "  ranks       %d\n", nprocs );
        std::printf( "  boundary    %s\n", noslide ? "pinned" : "sliding" );
    }

    if( selftest )
    {
        if( nprocs > 1 )
        {
            if( !rank ) std::cerr << "--selftest is a serial check; run it on one rank" << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return 1;
        }

        selftest_metrics();

        rval = selftest_against_verdict( mb, ownedElems, constraint );MB_CHK_ERR( rval );

        // Perturb first, so the gradient is checked somewhere other than the
        // minimum where every component is nearly zero and the relative
        // comparison is meaningless.
        rval = perturb_mesh( mb, verts, constraint, amplitude * h, seed );MB_CHK_ERR( rval );

        OptimizerOptions oo;
        oo.metric    = ( metricName == "condition" ? CONDITION_NUMBER : INVERSE_MEAN_RATIO );
        oo.delta     = delta;
        oo.pnorm     = pnorm;
        oo.verbosity = 0;
        MeshOptimizer optimizer( mb,
#ifdef MOAB_HAVE_MPI
                                 pcomm,
#endif
                                 oo, constraint );
        rval = optimizer.setup( ownedElems );MB_CHK_ERR( rval );

        double maxrel = 0.0;
        rval          = optimizer.verify_gradient( maxrel );MB_CHK_ERR( rval );
        std::printf( "\nAutoDiff gradient vs central differences:\n" );
        expect( "max relative difference", maxrel, 0.0, 1.0e-5 );

        std::printf( g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures );
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return g_failures ? 1 : 0;
    }

    const char* mname = ( metricName == "condition" ? "MB_CONDITION" : "MB_SHAPE" );
    const QualityType qt = ( metricName == "condition" ? MB_CONDITION : MB_SHAPE );
    QualityStats before, perturbed, after;

    rval = report_quality( mb,
#ifdef MOAB_HAVE_MPI
                           pcomm,
#endif
                           ownedElems, qt, before );MB_CHK_ERR( rval );
    if( !rank && verbosity > 0 ) print_quality( "initial", mname, before );

    if( infile.empty() && amplitude > 0.0 )
    {
        rval = perturb_mesh( mb, verts, constraint, amplitude * h, seed );MB_CHK_ERR( rval );
        rval = report_quality( mb,
#ifdef MOAB_HAVE_MPI
                               pcomm,
#endif
                               ownedElems, qt, perturbed );MB_CHK_ERR( rval );
        if( !rank && verbosity > 0 ) print_quality( "perturbed", mname, perturbed );
    }

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
    rval = optimizer.setup( ownedElems );MB_CHK_ERR( rval );

    if( !rank && verbosity > 0 )
        std::printf( "\nOptimizing %ld elements, metric %s, delta %g, p %g\n", optimizer.num_global_elements(),
                     metric_name( oo.metric ), delta, pnorm );

    OptimizerResult res;
    rval = optimizer.optimize( res );MB_CHK_ERR( rval );

    rval = report_quality( mb,
#ifdef MOAB_HAVE_MPI
                           pcomm,
#endif
                           ownedElems, qt, after );MB_CHK_ERR( rval );

    if( !rank && verbosity > 0 )
    {
        std::printf( "\n" );
        print_quality( "optimized", mname, after );
        std::printf( "\n  objective   %.12e -> %.12e\n", res.f_initial, res.f_final );
        std::printf( "  gradient    %.6e -> %.6e\n", res.gnorm_initial, res.gnorm_final );
        std::printf( "  iterations  %d (%d objective evaluations)%s\n", res.iterations, res.evaluations,
                     res.converged ? ", converged" : "" );
    }

    if( infile.empty() )
    {
        // The generated lattice is the exact minimizer on these domains, so how
        // far the vertices still sit from integer coordinates is a direct
        // measure of whether the optimizer found the right answer.
        double err = 0.0;
        rval       = lattice_error( mb,
#ifdef MOAB_HAVE_MPI
                              pcomm,
#endif
                              verts, err );MB_CHK_ERR( rval );
        if( !rank && verbosity > 0 ) std::printf( "  lattice     max vertex offset %.6e\n", err );
    }

    if( !outfile.empty() )
    {
        const char* wopts = ( nprocs > 1 ? "PARALLEL=WRITE_PART" : "" );
        rval              = mb->write_file( outfile.c_str(), 0, wopts );MB_CHK_ERR( rval );
        if( !rank && verbosity > 0 ) std::printf( "  wrote       %s\n", outfile.c_str() );
    }

    if( !rank && verbosity > 0 ) std::printf( "\n" );

#ifdef MOAB_HAVE_MPI
    delete pcomm;
    MPI_Finalize();
#endif
    return 0;
}
