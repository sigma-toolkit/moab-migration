/** @example MeshOptimizerCore.hpp
 * Distributed shape optimization of a MOAB mesh, using nothing but MOAB and
 * Eigen3.
 *
 * The objective is the p-mean of a differentiable element shape metric
 * (MeshOptimizationMetrics.hpp).  Gradients come from Eigen's unsupported
 * AutoDiff module, so there is no hand-differentiated metric anywhere, and the
 * minimization is L-BFGS with an Armijo backtracking line search.
 *
 * Everything the parallel case needs is here: each element is evaluated by
 * exactly one rank, vertex gradient contributions are summed across the
 * partition interface, and every inner product in the L-BFGS recursion is an
 * MPI_Allreduce.  The objective is therefore independent of the rank count,
 * which is the property the example's determinism check exercises.
 */

#ifndef MOAB_EXAMPLE_MESH_OPTIMIZER_CORE_HPP
#define MOAB_EXAMPLE_MESH_OPTIMIZER_CORE_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MeshOptimizationMetrics.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab_mpi.h"
#endif

namespace moab
{
namespace meshopt
{

/** Where a vertex is allowed to move.
 *
 * A subclass returns a symmetric 3x3 projector P for each vertex.  Search
 * directions and gradients are both premultiplied by P, so a vertex only ever
 * moves inside its permitted subspace:
 *
 *   rank(P) == 3  interior of a volume, free
 *   rank(P) == 2  on a plane or a box face, slides in that plane
 *   rank(P) == 1  on a box edge, slides along it
 *   rank(P) == 0  box corner, pinned
 *
 * project_point() additionally snaps a trial position back onto a curved
 * domain; it is a no-op for planar and box domains, where the projector alone
 * keeps vertices on the boundary exactly.
 */
class Constraint
{
  public:
    virtual ~Constraint() {}

    /** Symmetric projector at position x, row-major 3x3. */
    virtual void projector( const double* x, double* P ) const = 0;

    /** Snap a trial position back onto the domain.  Default: nothing to do. */
    virtual void project_point( double* /*x*/ ) const {}

    /** Reference normal for a 2D element centred at c.
     *
     * Surface metrics need this for the sign of the element area; without it an
     * inverted surface element looks valid.  Return false for volume meshes.
     */
    virtual bool surface_normal( const double* /*c*/, double* /*n*/ ) const
    {
        return false;
    }
};

/** Identity constraint: every vertex free in all three directions. */
class FreeConstraint : public Constraint
{
  public:
    void projector( const double*, double* P ) const override
    {
        for( int i = 0; i < 9; ++i )
            P[i] = 0.0;
        P[0] = P[4] = P[8] = 1.0;
    }
};

struct OptimizerOptions
{
    MetricType metric;
    double delta;      //!< Escobar regularization; 0 gives the pure metric
    double pnorm;      //!< objective is the p-mean of the element metric
    int maxiter;
    double gtol;       //!< converged when ||projected gradient|| falls below this
    int history;       //!< L-BFGS history pairs
    int verbosity;

    OptimizerOptions()
        : metric( INVERSE_MEAN_RATIO ), delta( 0.0 ), pnorm( 2.0 ), maxiter( 100 ), gtol( 1.0e-8 ), history( 7 ),
          verbosity( 1 )
    {
    }
};

struct OptimizerResult
{
    int iterations;
    int evaluations;
    double f_initial;
    double f_final;
    double gnorm_initial;
    double gnorm_final;
    bool converged;
};

class MeshOptimizer
{
  public:
    MeshOptimizer( Interface* mb,
#ifdef MOAB_HAVE_MPI
                   ParallelComm* pcomm,
#endif
                   const OptimizerOptions& opts,
                   const Constraint& constraint )
        : mMB( mb ),
#ifdef MOAB_HAVE_MPI
          mPcomm( pcomm ),
#endif
          mOpts( opts ), mConstraint( constraint ), mRank( 0 ), mSize( 1 ), mGlobalElems( 0 ), mCharLength( 1.0 ),
          mEvaluations( 0 )
    {
#ifdef MOAB_HAVE_MPI
        if( mPcomm )
        {
            mRank = mPcomm->rank();
            mSize = mPcomm->size();
        }
#endif
    }

    /** Cache connectivity, classify vertices and size the DOF vector.
     *
     * @param elements Elements this rank owns.  Passing elements owned by
     *                 another rank would double-count them in the objective.
     */
    ErrorCode setup( const Range& elements );

    /** Run L-BFGS to convergence or to the iteration limit. */
    ErrorCode optimize( OptimizerResult& result );

    /** Compare the AutoDiff gradient against central differences.
     *
     * Serial only: a finite difference here perturbs one global DOF and
     * re-evaluates the whole objective, which in parallel would need a
     * collective call per DOF inside a loop only one rank knows the length of.
     *
     * @param maxrel Largest relative discrepancy over the DOFs checked.
     * @param maxdof Cap on how many DOFs to check; <= 0 means all of them.
     */
    ErrorCode verify_gradient( double& maxrel, int maxdof = 60 );

    int num_free_vertices() const
    {
        return (int)mFreeVerts.size();
    }
    long num_global_elements() const
    {
        return mGlobalElems;
    }
    double characteristic_length() const
    {
        return mCharLength;
    }
    const std::vector< double >& objective_history() const
    {
        return mHistoryF;
    }

  private:
    // AutoDiff scalar.  The derivative vector is dynamically sized but capped
    // at 24 (a hexahedron's 8 nodes x 3), so Eigen keeps it on the stack and
    // the element loop does not allocate.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1, 0, 24, 1 > DerType;
    typedef Eigen::AutoDiffScalar< DerType > ADScalar;

    ErrorCode evaluate( double& f, bool need_gradient );
    ErrorCode push_positions();
    ErrorCode gather_gradient( std::vector< double >& g );
    ErrorCode set_free_positions( const std::vector< double >& x );
    ErrorCode get_free_positions( std::vector< double >& x );
    void apply_projector( const std::vector< double >& x, std::vector< double >& v ) const;
    double dot( const std::vector< double >& a, const std::vector< double >& b ) const;
    ErrorCode compute_characteristic_length();

    Interface* mMB;
#ifdef MOAB_HAVE_MPI
    ParallelComm* mPcomm;
#endif
    OptimizerOptions mOpts;
    const Constraint& mConstraint;
    int mRank, mSize;

    Range mElems;       //!< owned elements
    Range mAllVerts;    //!< every vertex touched by mElems, ghosts included
    Range mSharedVerts; //!< subset of mAllVerts on the partition interface
    std::vector< EntityHandle > mFreeVerts;  //!< owned, movable
    std::map< EntityHandle, int > mFreeIndex;

    Tag mGradTag, mPosTag;
    long mGlobalElems;
    double mCharLength;
    int mEvaluations;
    std::vector< double > mHistoryF;
};

// ---------------------------------------------------------------------------

inline ErrorCode MeshOptimizer::setup( const Range& elements )
{
    ErrorCode rval;
    mElems = elements;

    // Drop anything the metric layer does not know how to differentiate, so an
    // unsupported element type is reported rather than silently contributing
    // zero to the objective.
    Range supported;
    for( Range::iterator it = mElems.begin(); it != mElems.end(); ++it )
        if( metric_num_nodes( mMB->type_from_handle( *it ) ) > 0 ) supported.insert( *it );
    if( supported.size() != mElems.size() && !mRank )
        std::printf( "  warning: %lu element(s) of unsupported type are excluded from the objective\n",
                     (unsigned long)( mElems.size() - supported.size() ) );
    mElems = supported;

    rval = mMB->get_connectivity( mElems, mAllVerts );MB_CHK_ERR( rval );

    // Owned, movable vertices are the degrees of freedom.  A vertex is immovable
    // when its projector is rank 0 (a box corner), or when the mesh carries the
    // conventional "fixed" tag and has set it.
    Range ownedVerts = mAllVerts;
#ifdef MOAB_HAVE_MPI
    if( mPcomm && mSize > 1 )
    {
        ownedVerts.clear();
        rval = mPcomm->filter_pstatus( mAllVerts, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &ownedVerts );MB_CHK_ERR( rval );
        mSharedVerts.clear();
        rval = mPcomm->get_shared_entities( -1, mSharedVerts, 0, false, false );MB_CHK_ERR( rval );
        mSharedVerts = intersect( mSharedVerts, mAllVerts );
    }
#endif

    Tag fixedTag = 0;
    if( MB_SUCCESS != mMB->tag_get_handle( "fixed", 1, MB_TYPE_INTEGER, fixedTag ) ) fixedTag = 0;

    mFreeVerts.clear();
    mFreeIndex.clear();
    for( Range::iterator it = ownedVerts.begin(); it != ownedVerts.end(); ++it )
    {
        if( fixedTag )
        {
            int isfixed = 0;
            if( MB_SUCCESS == mMB->tag_get_data( fixedTag, &( *it ), 1, &isfixed ) && isfixed ) continue;
        }
        double x[3], P[9];
        rval = mMB->get_coords( &( *it ), 1, x );MB_CHK_ERR( rval );
        mConstraint.projector( x, P );
        double trace = P[0] + P[4] + P[8];
        if( trace < 1.0e-12 ) continue;  // pinned
        mFreeIndex[*it] = (int)mFreeVerts.size();
        mFreeVerts.push_back( *it );
    }

    const double zero3[3] = { 0, 0, 0 };
    rval = mMB->tag_get_handle( "MESHOPT_GRAD", 3, MB_TYPE_DOUBLE, mGradTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                zero3 );MB_CHK_ERR( rval );
    rval = mMB->tag_get_handle( "MESHOPT_POS", 3, MB_TYPE_DOUBLE, mPosTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                zero3 );MB_CHK_ERR( rval );

    long localElems = (long)mElems.size();
    mGlobalElems    = localElems;
#ifdef MOAB_HAVE_MPI
    if( mPcomm && mSize > 1 )
        MPI_Allreduce( &localElems, &mGlobalElems, 1, MPI_LONG, MPI_SUM, mPcomm->comm() );
#endif
    if( mGlobalElems == 0 ) MB_SET_ERR( MB_FAILURE, "No supported elements to optimize" );

    return compute_characteristic_length();
}

/** Mean edge length of the first element, reduced globally.
 *
 * Only used to scale the first trial step and the finite-difference
 * increment, so a rough value is enough - but it has to be the *same* rough
 * value everywhere, or ranks would take different first steps.
 */
inline ErrorCode MeshOptimizer::compute_characteristic_length()
{
    double localSum = 0.0;
    long localCount = 0;
    std::vector< EntityHandle > connStorage;
    for( Range::iterator it = mElems.begin(); it != mElems.end(); ++it )
    {
        const EntityHandle* conn = NULL;
        int nnodes               = 0;
        if( MB_SUCCESS != mMB->get_connectivity( *it, conn, nnodes, false, &connStorage ) || nnodes < 2 ) continue;
        double c[24 * 3];
        if( nnodes > 24 ) nnodes = 24;
        if( MB_SUCCESS != mMB->get_coords( conn, nnodes, c ) ) continue;
        for( int i = 1; i < nnodes; ++i )
        {
            const double dx = c[3 * i] - c[0], dy = c[3 * i + 1] - c[1], dz = c[3 * i + 2] - c[2];
            localSum += std::sqrt( dx * dx + dy * dy + dz * dz );
            ++localCount;
        }
    }
    double sum   = localSum;
    long count   = localCount;
#ifdef MOAB_HAVE_MPI
    if( mPcomm && mSize > 1 )
    {
        MPI_Allreduce( &localSum, &sum, 1, MPI_DOUBLE, MPI_SUM, mPcomm->comm() );
        MPI_Allreduce( &localCount, &count, 1, MPI_LONG, MPI_SUM, mPcomm->comm() );
    }
#endif
    mCharLength = ( count > 0 && sum > 0.0 ) ? sum / (double)count : 1.0;
    return MB_SUCCESS;
}

/** Publish owned vertex positions to the ranks that ghost them.
 *
 * MOAB coordinates are not a tag, so they travel through a mirror tag:
 * every rank writes its coordinates into MESHOPT_POS, exchange_tags overwrites
 * each sharer's copy with the owner's, and the result is read back into the
 * coordinates.  An owned vertex gets its own value back unchanged.
 */
inline ErrorCode MeshOptimizer::push_positions()
{
#ifdef MOAB_HAVE_MPI
    if( !mPcomm || mSize == 1 || mSharedVerts.empty() ) return MB_SUCCESS;

    ErrorCode rval;
    std::vector< double > buf( 3 * mAllVerts.size() );
    rval = mMB->get_coords( mAllVerts, &buf[0] );MB_CHK_ERR( rval );
    rval = mMB->tag_set_data( mPosTag, mAllVerts, &buf[0] );MB_CHK_ERR( rval );
    rval = mPcomm->exchange_tags( mPosTag, mSharedVerts );MB_CHK_ERR( rval );
    rval = mMB->tag_get_data( mPosTag, mAllVerts, &buf[0] );MB_CHK_ERR( rval );
    rval = mMB->set_coords( mAllVerts, &buf[0] );MB_CHK_ERR( rval );
#endif
    return MB_SUCCESS;
}

inline ErrorCode MeshOptimizer::evaluate( double& f, bool need_gradient )
{
    ErrorCode rval;
    ++mEvaluations;

    if( need_gradient )
    {
        std::vector< double > zeros( 3 * mAllVerts.size(), 0.0 );
        rval = mMB->tag_set_data( mGradTag, mAllVerts, &zeros[0] );MB_CHK_ERR( rval );
    }

    double localSum = 0.0;
    std::vector< EntityHandle > connStorage;

    for( Range::iterator it = mElems.begin(); it != mElems.end(); ++it )
    {
        const EntityType type = mMB->type_from_handle( *it );
        const int nnodes      = metric_num_nodes( type );
        const int nv          = 3 * nnodes;

        // The storage vector matters: a structured-mesh element has no stored
        // connectivity array, and without somewhere to materialize one this
        // call fails with MB_STRUCTURED_MESH.
        const EntityHandle* conn = NULL;
        int connlen              = 0;
        rval                     = mMB->get_connectivity( *it, conn, connlen, false, &connStorage );MB_CHK_ERR( rval );
        if( connlen < nnodes ) continue;

        double c[24];
        rval = mMB->get_coords( conn, nnodes, c );MB_CHK_ERR( rval );

        // 2D metrics need a reference normal to give the element area a sign.
        double normal[3] = { 0.0, 0.0, 1.0 };
        const double* np = NULL;
        if( metric_dimension( type ) == 2 )
        {
            double centroid[3] = { 0, 0, 0 };
            for( int i = 0; i < nnodes; ++i )
                for( int d = 0; d < 3; ++d )
                    centroid[d] += c[3 * i + d] / nnodes;
            if( !mConstraint.surface_normal( centroid, normal ) )
            {
                normal[0] = normal[1] = 0.0;
                normal[2]             = 1.0;
            }
            np = normal;
        }

        if( !need_gradient )
        {
            const double q = element_quality< double >( type, c, np, mOpts.metric, mOpts.delta );
            localSum += std::pow( q, mOpts.pnorm );
            continue;
        }

        ADScalar ac[24];
        for( int i = 0; i < nv; ++i )
            ac[i] = ADScalar( c[i], nv, i );

        const ADScalar q  = element_quality< ADScalar >( type, ac, np, mOpts.metric, mOpts.delta );
        const ADScalar qp = pow( q, mOpts.pnorm );
        localSum += qp.value();

        // Scatter the element gradient onto its vertices.  Ghost vertices
        // receive their share too; the reduction below moves it to the owner.
        double vg[3];
        for( int i = 0; i < nnodes; ++i )
        {
            rval = mMB->tag_get_data( mGradTag, conn + i, 1, vg );MB_CHK_ERR( rval );
            for( int d = 0; d < 3; ++d )
                vg[d] += qp.derivatives()( 3 * i + d );
            rval = mMB->tag_set_data( mGradTag, conn + i, 1, vg );MB_CHK_ERR( rval );
        }
    }

    double globalSum = localSum;
#ifdef MOAB_HAVE_MPI
    if( mPcomm && mSize > 1 ) MPI_Allreduce( &localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, mPcomm->comm() );
#endif
    f = globalSum / (double)mGlobalElems;

#ifdef MOAB_HAVE_MPI
    if( need_gradient && mPcomm && mSize > 1 && !mSharedVerts.empty() )
    {
        // Sum the partial contributions onto the owner, then hand the total
        // back to the sharers.  The exchange is not strictly needed for the
        // L-BFGS update, which only reads owned entries, but it keeps the
        // gradient tag consistent for anything that inspects it.
        rval = mPcomm->reduce_tags( mGradTag, MPI_SUM, mSharedVerts );MB_CHK_ERR( rval );
        rval = mPcomm->exchange_tags( mGradTag, mSharedVerts );MB_CHK_ERR( rval );
    }
#endif
    return MB_SUCCESS;
}

inline ErrorCode MeshOptimizer::gather_gradient( std::vector< double >& g )
{
    const int n = (int)mFreeVerts.size();
    g.assign( 3 * n, 0.0 );
    if( n == 0 ) return MB_SUCCESS;

    ErrorCode rval = mMB->tag_get_data( mGradTag, &mFreeVerts[0], n, &g[0] );MB_CHK_ERR( rval );
    // The objective is a mean, so the accumulated sum has to be scaled to match.
    const double scale = 1.0 / (double)mGlobalElems;
    for( int i = 0; i < 3 * n; ++i )
        g[i] *= scale;
    return MB_SUCCESS;
}

inline ErrorCode MeshOptimizer::get_free_positions( std::vector< double >& x )
{
    const int n = (int)mFreeVerts.size();
    x.assign( 3 * n, 0.0 );
    if( n == 0 ) return MB_SUCCESS;
    return mMB->get_coords( &mFreeVerts[0], n, &x[0] );
}

inline ErrorCode MeshOptimizer::set_free_positions( const std::vector< double >& x )
{
    const int n = (int)mFreeVerts.size();
    if( n > 0 )
    {
        ErrorCode rval = mMB->set_coords( &mFreeVerts[0], n, &x[0] );MB_CHK_ERR( rval );
    }
    return push_positions();
}

/** v <- P(x) v, vertex by vertex. */
inline void MeshOptimizer::apply_projector( const std::vector< double >& x, std::vector< double >& v ) const
{
    const int n = (int)mFreeVerts.size();
    for( int i = 0; i < n; ++i )
    {
        double P[9];
        mConstraint.projector( &x[3 * i], P );
        const double v0 = v[3 * i], v1 = v[3 * i + 1], v2 = v[3 * i + 2];
        v[3 * i]     = P[0] * v0 + P[1] * v1 + P[2] * v2;
        v[3 * i + 1] = P[3] * v0 + P[4] * v1 + P[5] * v2;
        v[3 * i + 2] = P[6] * v0 + P[7] * v1 + P[8] * v2;
    }
}

inline double MeshOptimizer::dot( const std::vector< double >& a, const std::vector< double >& b ) const
{
    double local = 0.0;
    for( size_t i = 0; i < a.size(); ++i )
        local += a[i] * b[i];
    double global = local;
#ifdef MOAB_HAVE_MPI
    if( mPcomm && mSize > 1 ) MPI_Allreduce( &local, &global, 1, MPI_DOUBLE, MPI_SUM, mPcomm->comm() );
#endif
    return global;
}

inline ErrorCode MeshOptimizer::optimize( OptimizerResult& result )
{
    ErrorCode rval;
    const int n = (int)mFreeVerts.size();

    std::vector< double > x, g, d, xtrial, gnew;
    rval = get_free_positions( x );MB_CHK_ERR( rval );

    double f = 0.0;
    rval     = evaluate( f, true );MB_CHK_ERR( rval );
    rval     = gather_gradient( g );MB_CHK_ERR( rval );
    apply_projector( x, g );

    double gnorm = std::sqrt( dot( g, g ) );

    result.f_initial     = f;
    result.gnorm_initial = gnorm;
    result.converged     = false;
    result.iterations    = 0;
    mHistoryF.clear();
    mHistoryF.push_back( f );

    if( mOpts.verbosity > 0 && !mRank )
        std::printf( "  iter %4d   f = %.12e   |g| = %.6e\n", 0, f, gnorm );

    if( !std::isfinite( f ) )
    {
        // An unregularized barrier metric is +inf on an inverted or degenerate
        // element and its gradient is NaN, so there is no descent direction to
        // find anywhere.  Say why rather than letting the line search fail 30
        // times and report zero iterations.
        if( !mRank )
            std::printf( "  the mesh contains inverted or degenerate elements, so the objective is infinite\n"
                         "  at the starting point.  Re-run with --delta > 0 (try %.3g) to regularize the\n"
                         "  determinant, which makes the metric finite through inversion and lets the\n"
                         "  optimizer untangle the mesh before improving its shape.\n",
                         0.1 * mCharLength );
        result.f_final     = f;
        result.gnorm_final = gnorm;
        result.evaluations = mEvaluations;
        return MB_SUCCESS;
    }

    // L-BFGS history.
    std::vector< std::vector< double > > S, Y;
    std::vector< double > rho;

    const double c1 = 1.0e-4;  // Armijo sufficient-decrease constant

    int iter = 0;
    for( ; iter < mOpts.maxiter; ++iter )
    {
        if( gnorm <= mOpts.gtol )
        {
            result.converged = true;
            break;
        }

        // Two-loop recursion for d = -H g.
        d.assign( 3 * n, 0.0 );
        for( int i = 0; i < 3 * n; ++i )
            d[i] = -g[i];

        const int m = (int)S.size();
        std::vector< double > alpha( m, 0.0 );
        for( int i = m - 1; i >= 0; --i )
        {
            alpha[i] = rho[i] * dot( S[i], d );
            for( int k = 0; k < 3 * n; ++k )
                d[k] -= alpha[i] * Y[i][k];
        }
        if( m > 0 )
        {
            const double yy = dot( Y[m - 1], Y[m - 1] );
            const double sy = dot( S[m - 1], Y[m - 1] );
            const double gamma = ( yy > 0.0 ) ? sy / yy : 1.0;
            for( int k = 0; k < 3 * n; ++k )
                d[k] *= gamma;
        }
        for( int i = 0; i < m; ++i )
        {
            const double beta = rho[i] * dot( Y[i], d );
            for( int k = 0; k < 3 * n; ++k )
                d[k] += ( alpha[i] - beta ) * S[i][k];
        }

        apply_projector( x, d );

        double gd = dot( g, d );
        if( gd >= 0.0 )
        {
            // The quasi-Newton direction is not a descent direction, which can
            // happen after a constraint projection rotates the subspace.  Drop
            // the history and fall back to steepest descent.
            S.clear();
            Y.clear();
            rho.clear();
            for( int k = 0; k < 3 * n; ++k )
                d[k] = -g[k];
            apply_projector( x, d );
            gd = dot( g, d );
            if( gd >= 0.0 ) break;  // projected gradient is zero; nothing to do
        }

        // First step is scaled to the mesh: a unit step along -g has no natural
        // length, and an overlong first step can inject tangling the
        // regularization then has to undo.
        double step = 1.0;
        if( S.empty() )
        {
            const double dnorm = std::sqrt( dot( d, d ) );
            if( dnorm > 0.0 ) step = std::min( 1.0, 0.1 * mCharLength / dnorm );
        }

        bool accepted = false;
        double fnew   = f;
        for( int ls = 0; ls < 30; ++ls )
        {
            xtrial.assign( x.begin(), x.end() );
            for( int k = 0; k < 3 * n; ++k )
                xtrial[k] += step * d[k];
            for( int i = 0; i < n; ++i )
                mConstraint.project_point( &xtrial[3 * i] );

            rval = set_free_positions( xtrial );MB_CHK_ERR( rval );
            rval = evaluate( fnew, false );MB_CHK_ERR( rval );

            if( std::isfinite( fnew ) && fnew <= f + c1 * step * gd )
            {
                accepted = true;
                break;
            }
            step *= 0.5;
        }

        if( !accepted )
        {
            // Restore the last good iterate before giving up.
            rval = set_free_positions( x );MB_CHK_ERR( rval );
            if( mOpts.verbosity > 0 && !mRank )
            {
                // Distinguish a stall at the precision limit from a real
                // failure.  On a curved domain the tangent space rotates as
                // vertices move, so the curvature pairs go stale and L-BFGS
                // runs out of progress while the gradient is already tiny -
                // that is expected, not an error.
                const double reduction = ( result.gnorm_initial > 0.0 ) ? gnorm / result.gnorm_initial : 1.0;
                if( reduction < 1.0e-3 )
                    std::printf( "  no further decrease available at this precision "
                                 "(gradient reduced %.1e-fold); stopping\n",
                                 1.0 / reduction );
                else
                    std::printf( "  line search failed with the gradient still at %.3e; stopping\n", gnorm );
            }
            break;
        }

        rval = evaluate( fnew, true );MB_CHK_ERR( rval );
        rval = gather_gradient( gnew );MB_CHK_ERR( rval );
        apply_projector( xtrial, gnew );

        // s is the realized displacement, not step*d: project_point may have
        // moved the trial position off the search ray, and feeding the
        // unrealized step into the history would corrupt the curvature pairs.
        std::vector< double > s( 3 * n ), y( 3 * n );
        for( int k = 0; k < 3 * n; ++k )
        {
            s[k] = xtrial[k] - x[k];
            y[k] = gnew[k] - g[k];
        }
        const double sy = dot( s, y );
        if( sy > 1.0e-14 * std::sqrt( dot( s, s ) * dot( y, y ) ) )
        {
            // Armijo alone does not guarantee the curvature condition, so pairs
            // that would make the inverse Hessian indefinite are skipped rather
            // than silently corrupting the recursion.
            S.push_back( s );
            Y.push_back( y );
            rho.push_back( 1.0 / sy );
            if( (int)S.size() > mOpts.history )
            {
                S.erase( S.begin() );
                Y.erase( Y.begin() );
                rho.erase( rho.begin() );
            }
        }

        x.swap( xtrial );
        g.swap( gnew );
        f     = fnew;
        gnorm = std::sqrt( dot( g, g ) );
        mHistoryF.push_back( f );

        if( mOpts.verbosity > 1 && !mRank )
            std::printf( "  iter %4d   f = %.12e   |g| = %.6e   step = %.3e\n", iter + 1, f, gnorm, step );
    }

    if( mOpts.verbosity > 0 && !mRank )
        std::printf( "  iter %4d   f = %.12e   |g| = %.6e%s\n", iter, f, gnorm,
                     result.converged ? "   (converged)" : "" );

    result.iterations  = iter;
    result.evaluations = mEvaluations;
    result.f_final     = f;
    result.gnorm_final = gnorm;
    return MB_SUCCESS;
}

inline ErrorCode MeshOptimizer::verify_gradient( double& maxrel, int maxdof )
{
    maxrel = 0.0;
    if( mSize > 1 ) MB_SET_ERR( MB_NOT_IMPLEMENTED, "--verify-gradient is a serial check" );

    ErrorCode rval;
    std::vector< double > x, g;
    rval = get_free_positions( x );MB_CHK_ERR( rval );

    double f0 = 0.0;
    rval      = evaluate( f0, true );MB_CHK_ERR( rval );
    rval      = gather_gradient( g );MB_CHK_ERR( rval );
    // Compare the raw gradient, not the projected one: the projector is exact
    // by construction and would mask a wrong derivative in a pinned direction.

    const int ndof = (int)x.size();
    const int nchk = ( maxdof > 0 && maxdof < ndof ) ? maxdof : ndof;
    const double h = 1.0e-6 * mCharLength;

    std::vector< double > xp( x );
    for( int i = 0; i < nchk; ++i )
    {
        double fp = 0.0, fm = 0.0;
        xp[i] = x[i] + h;
        rval  = set_free_positions( xp );MB_CHK_ERR( rval );
        rval  = evaluate( fp, false );MB_CHK_ERR( rval );
        xp[i] = x[i] - h;
        rval  = set_free_positions( xp );MB_CHK_ERR( rval );
        rval  = evaluate( fm, false );MB_CHK_ERR( rval );
        xp[i] = x[i];

        const double fd    = ( fp - fm ) / ( 2.0 * h );
        const double scale = std::max( 1.0e-8, std::fabs( fd ) );
        const double rel   = std::fabs( fd - g[i] ) / scale;
        if( rel > maxrel ) maxrel = rel;
    }

    rval = set_free_positions( x );MB_CHK_ERR( rval );
    return MB_SUCCESS;
}

}  // namespace meshopt
}  // namespace moab

#endif
