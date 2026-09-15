/** @example MeshOptimizationMetrics.hpp
 * Differentiable element shape metrics for the mesh optimization examples.
 *
 * Everything here is templated on the scalar type and written in plain
 * arithmetic, so the same code evaluates under \c double and under
 * \c Eigen::AutoDiffScalar.  That is the whole point of the header: MOAB's own
 * VerdictWrapper::quality_measure takes a hard \c double& for the result, so a
 * templated scalar cannot propagate through it and the metrics have to be
 * reimplemented here to get gradients.
 *
 * There is deliberately no Eigen dependency in this file.  Unqualified sqrt()
 * with "using std::sqrt" resolves to std::sqrt for double and, by ADL, to
 * Eigen's overload for AutoDiffScalar.
 *
 * Both metrics equal 1 for the ideal element and grow without bound as the
 * element degenerates, so they are minimized.
 */

#ifndef MOAB_EXAMPLE_MESH_OPTIMIZATION_METRICS_HPP
#define MOAB_EXAMPLE_MESH_OPTIMIZATION_METRICS_HPP

#include <cmath>
#include "moab/EntityType.hpp"

namespace moab
{
namespace meshopt
{

enum MetricType
{
    INVERSE_MEAN_RATIO = 0,
    CONDITION_NUMBER
};

inline const char* metric_name( MetricType m )
{
    return ( m == INVERSE_MEAN_RATIO ? "inverse-mean-ratio" : "condition-number" );
}

/** Nodes per supported element type; 0 means unsupported. */
inline int metric_num_nodes( EntityType type )
{
    switch( type )
    {
        case MBTRI:
            return 3;
        case MBQUAD:
            return 4;
        case MBTET:
            return 4;
        case MBHEX:
            return 8;
        default:
            return 0;
    }
}

/** Topological dimension of the element, 2 for tri/quad and 3 for tet/hex. */
inline int metric_dimension( EntityType type )
{
    return ( type == MBTRI || type == MBQUAD ) ? 2 : 3;
}

/** Number of corners the metric is sampled at.
 *
 * Tri and tet have a constant Jacobian, and with the equilateral/regular ideal
 * weight the metric comes out the same at every corner, so one sample is
 * exact.  Quad and hex have a Jacobian that varies over the element, so every
 * corner is sampled.
 */
inline int metric_num_samples( EntityType type )
{
    switch( type )
    {
        case MBTRI:
        case MBTET:
            return 1;
        case MBQUAD:
            return 4;
        case MBHEX:
            return 8;
        default:
            return 0;
    }
}

/** Escobar-style regularized determinant.
 *
 * Barrier metrics divide by the Jacobian determinant, so they are +inf on a
 * degenerate element and meaningless on an inverted one.  Replacing d by
 *
 *     (d + sqrt(d^2 + 4 delta^2)) / 2
 *
 * keeps the quotient finite and smoothly differentiable straight through
 * inversion, which lets one optimizer pass untangle a mesh and then improve
 * its shape.  See Escobar et al., "Simultaneous untangling and smoothing of
 * tetrahedral meshes", CMAME 192 (2003) 2775-2787.
 *
 * At delta == 0 this collapses to max(d, 0), so for a valid element it is
 * exactly d.  That is what makes the --verify-metric cross-check against
 * Verdict meaningful: with delta == 0 on an uninverted mesh the metric below
 * is the unregularized one.
 */
template < typename Scalar >
inline Scalar regularized_det( const Scalar& d, double delta )
{
    using std::sqrt;
    if( delta <= 0.0 ) return ( d + ( d < Scalar( 0.0 ) ? -d : d ) ) * Scalar( 0.5 );
    return ( d + sqrt( d * d + Scalar( 4.0 * delta * delta ) ) ) * Scalar( 0.5 );
}

namespace detail
{

/** Corner node triples for a hexahedron.
 *
 * Row c holds the three nodes the edge vectors at corner c point to, ordered so
 * that the triple (n0-c, n1-c, n2-c) is right-handed.  Corners 4-7 have their
 * first two entries swapped relative to 0-3 because the third edge there runs
 * along -zeta; without the swap those corners would report every element as
 * inverted.  This is Verdict's convention.
 */
const int HEX_CORNERS[8][3] = { { 1, 3, 4 }, { 2, 0, 5 }, { 3, 1, 6 }, { 0, 2, 7 },
                                { 7, 5, 0 }, { 4, 6, 1 }, { 5, 7, 2 }, { 6, 4, 3 } };

/** Corner node pairs for a quadrilateral, counter-clockwise at every corner. */
const int QUAD_CORNERS[4][2] = { { 1, 3 }, { 2, 0 }, { 3, 1 }, { 0, 2 } };

/** Inverse of the ideal-element weight W for a triangle.
 *
 * W maps the unit right triangle onto the equilateral one, so T = A * W^-1 is
 * the identity exactly when A is equilateral.  Stored column-major as the
 * multipliers applied to the two edge vectors.
 */
inline void tri_weight_inverse( double wi[2][2] )
{
    const double s3 = std::sqrt( 3.0 );
    wi[0][0]        = 1.0;
    wi[1][0]        = 0.0;
    wi[0][1]        = -1.0 / s3;
    wi[1][1]        = 2.0 / s3;
}

/** Inverse of the ideal-element weight W for a tetrahedron (regular tet). */
inline void tet_weight_inverse( double wi[3][3] )
{
    const double s3 = std::sqrt( 3.0 ), s6 = std::sqrt( 6.0 );
    wi[0][0] = 1.0;        wi[1][0] = 0.0;       wi[2][0] = 0.0;
    wi[0][1] = -1.0 / s3;  wi[1][1] = 2.0 / s3;  wi[2][1] = 0.0;
    wi[0][2] = -1.0 / s6;  wi[1][2] = -1.0 / s6; wi[2][2] = std::sqrt( 1.5 );
}

/** Shape metric of a 2D corner, from two edge vectors and a reference normal.
 *
 * The normal supplies the sign of the area, which a 3x2 Jacobian on its own
 * cannot: without it an inverted surface element is indistinguishable from a
 * valid one and the barrier never fires.  For a mesh in the z = const plane
 * the normal is (0,0,1) and this reduces to the ordinary 2x2 determinant; on
 * the sphere it is the local radial direction.
 *
 * Note that both metrics are the same function in 2D.  For a 2x2 matrix the
 * adjugate gives |T^-1|_F = |T|_F / det T, so the condition number
 * |T|_F |T^-1|_F / 2 equals |T|_F^2 / (2 det T), which is the inverse mean
 * ratio.  They differ only for tet and hex.  Both are computed by the same
 * expression here rather than pretending otherwise.
 */
template < typename Scalar >
inline Scalar corner_metric_2d( const Scalar t0[3], const Scalar t1[3], const double normal[3], double delta )
{
    const Scalar f = t0[0] * t0[0] + t0[1] * t0[1] + t0[2] * t0[2] + t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2];

    // Signed area scaling: normal . (t0 x t1).
    const Scalar cx = t0[1] * t1[2] - t0[2] * t1[1];
    const Scalar cy = t0[2] * t1[0] - t0[0] * t1[2];
    const Scalar cz = t0[0] * t1[1] - t0[1] * t1[0];
    const Scalar d  = cx * Scalar( normal[0] ) + cy * Scalar( normal[1] ) + cz * Scalar( normal[2] );

    return f / ( Scalar( 2.0 ) * regularized_det( d, delta ) );
}

/** Shape metric of a 3D corner, from three edge vectors. */
template < typename Scalar >
inline Scalar corner_metric_3d( const Scalar t0[3],
                                const Scalar t1[3],
                                const Scalar t2[3],
                                MetricType metric,
                                double delta )
{
    using std::pow;
    using std::sqrt;

    const Scalar f = t0[0] * t0[0] + t0[1] * t0[1] + t0[2] * t0[2] + t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2] +
                     t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2];

    const Scalar cx = t1[1] * t2[2] - t1[2] * t2[1];
    const Scalar cy = t1[2] * t2[0] - t1[0] * t2[2];
    const Scalar cz = t1[0] * t2[1] - t1[1] * t2[0];
    const Scalar d  = t0[0] * cx + t0[1] * cy + t0[2] * cz;

    const Scalar dr = regularized_det( d, delta );

    if( metric == INVERSE_MEAN_RATIO )
    {
        // |T|_F^2 / (3 det^(2/3)).  The exponent stays a plain double: Eigen's
        // AutoDiffScalar only overloads pow() for a raw scalar exponent.
        const Scalar c = pow( dr, 2.0 / 3.0 );
        return f / ( Scalar( 3.0 ) * c );
    }

    // Condition number |T|_F |T^-1|_F / 3.  T^-1 = adj(T)/det, and the adjugate
    // is polynomial in the entries, so this stays a single division and
    // AutoDiff never sees a reciprocal of a near-zero intermediate.
    const Scalar a[9] = { cx,
                          t2[1] * t0[2] - t2[2] * t0[1],
                          t0[1] * t1[2] - t0[2] * t1[1],
                          cy,
                          t2[2] * t0[0] - t2[0] * t0[2],
                          t0[2] * t1[0] - t0[0] * t1[2],
                          cz,
                          t2[0] * t0[1] - t2[1] * t0[0],
                          t0[0] * t1[1] - t0[1] * t1[0] };

    Scalar fa = a[0] * a[0];
    for( int i = 1; i < 9; ++i )
        fa += a[i] * a[i];

    return sqrt( f ) * sqrt( fa ) / ( Scalar( 3.0 ) * dr );
}

}  // namespace detail

/** Shape quality of one element, RMS-averaged over its corner samples.
 *
 * @param type    MBTRI, MBQUAD, MBTET or MBHEX.
 * @param coords  3 * metric_num_nodes(type) interleaved xyz values.
 * @param normal  Reference normal, required for MBTRI and MBQUAD, ignored
 *                otherwise.  May be NULL for the 3D types.
 * @param metric  Which metric; the two coincide for MBTRI and MBQUAD.
 * @param delta   Escobar regularization parameter, 0 for the pure metric.
 *
 * Corners are combined by RMS rather than by max because the optimizer needs a
 * differentiable objective and max is not.  RMS of values that are each >= 1 is
 * >= 1, with equality only when every corner is ideal, so the ideal element
 * still scores exactly 1.
 */
template < typename Scalar >
inline Scalar element_quality( EntityType type,
                               const Scalar* coords,
                               const double* normal,
                               MetricType metric,
                               double delta )
{
    using std::sqrt;

    const int nsamp = metric_num_samples( type );
    Scalar accum( 0.0 );

    if( type == MBTRI )
    {
        double wi[2][2];
        detail::tri_weight_inverse( wi );
        Scalar t0[3], t1[3];
        for( int d = 0; d < 3; ++d )
        {
            const Scalar e0 = coords[3 + d] - coords[d];
            const Scalar e1 = coords[6 + d] - coords[d];
            t0[d]           = e0 * Scalar( wi[0][0] ) + e1 * Scalar( wi[1][0] );
            t1[d]           = e0 * Scalar( wi[0][1] ) + e1 * Scalar( wi[1][1] );
        }
        const Scalar q = detail::corner_metric_2d( t0, t1, normal, delta );
        accum          = q * q;
    }
    else if( type == MBQUAD )
    {
        // The ideal quad is the unit square, so W is the identity and the edge
        // vectors are the Jacobian columns directly.
        for( int c = 0; c < 4; ++c )
        {
            Scalar t0[3], t1[3];
            const int n0 = detail::QUAD_CORNERS[c][0], n1 = detail::QUAD_CORNERS[c][1];
            for( int d = 0; d < 3; ++d )
            {
                t0[d] = coords[3 * n0 + d] - coords[3 * c + d];
                t1[d] = coords[3 * n1 + d] - coords[3 * c + d];
            }
            const Scalar q = detail::corner_metric_2d( t0, t1, normal, delta );
            accum += q * q;
        }
    }
    else if( type == MBTET )
    {
        double wi[3][3];
        detail::tet_weight_inverse( wi );
        Scalar t[3][3];
        for( int d = 0; d < 3; ++d )
        {
            const Scalar e[3] = { coords[3 + d] - coords[d], coords[6 + d] - coords[d], coords[9 + d] - coords[d] };
            for( int j = 0; j < 3; ++j )
                t[j][d] = e[0] * Scalar( wi[0][j] ) + e[1] * Scalar( wi[1][j] ) + e[2] * Scalar( wi[2][j] );
        }
        const Scalar q = detail::corner_metric_3d( t[0], t[1], t[2], metric, delta );
        accum          = q * q;
    }
    else if( type == MBHEX )
    {
        // The ideal hex is the unit cube, so again W is the identity.
        for( int c = 0; c < 8; ++c )
        {
            Scalar t[3][3];
            for( int j = 0; j < 3; ++j )
            {
                const int n = detail::HEX_CORNERS[c][j];
                for( int d = 0; d < 3; ++d )
                    t[j][d] = coords[3 * n + d] - coords[3 * c + d];
            }
            const Scalar q = detail::corner_metric_3d( t[0], t[1], t[2], metric, delta );
            accum += q * q;
        }
    }
    else
    {
        return Scalar( 0.0 );
    }

    return sqrt( accum / Scalar( (double)nsamp ) );
}

}  // namespace meshopt
}  // namespace moab

#endif
