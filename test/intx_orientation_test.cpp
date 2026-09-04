/*
 * intx_orientation_test.cpp
 *
 * Regression tests for IntxAreaUtils::positive_orientation() and for the
 * negative-area diagnostic.
 *
 * Background: TempestRemap-generated ".g" meshes store their cells wound
 * clockwise with respect to MOAB's convention, so positive_orientation() has to
 * reverse every single cell.  The orientation probe used to run through
 * area_spherical_triangle_lHuiller(), which carried an unconditional
 * "negative area" report.  Because the probe deliberately looks for negatively
 * oriented cells, that report fired once per fan sub-triangle on perfectly valid
 * input -- 7570 messages for the two meshes used by intx_rll_cs_sphere_test --
 * drowning out any report of a genuinely broken cell.
 *
 * These tests pin down the three properties that fix depends on:
 *   1. positive_orientation() actually repairs a fully-inverted mesh, and is
 *      idempotent: a second pass finds nothing left to do.
 *   2. The repaired mesh integrates to the exact surface area of the sphere.
 *   3. The negative-area diagnostic still fires for a genuinely inverted element,
 *      and stays silent for a correctly wound one.
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

#include "moab/Core.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "TestUtil.hpp"

using namespace moab;

// Count the cells whose stored winding is negative, without emitting anything.
// Van Oosterom & Strackee returns the signed area directly, so this is an exact
// orientation census rather than a heuristic.
static int count_inverted_cells( Interface* mb, EntityHandle set, double R )
{
    Range cells;
    if( MB_SUCCESS != mb->get_entities_by_dimension( set, 2, cells ) ) return -1;

    IntxAreaUtils vos( IntxAreaUtils::VanOosteromStrackee );
    int inverted = 0;
    for( Range::iterator it = cells.begin(); it != cells.end(); ++it )
    {
        const EntityHandle* conn = nullptr;
        int num_nodes            = 0;
        if( MB_SUCCESS != mb->get_connectivity( *it, conn, num_nodes ) ) return -1;

        std::vector< double > coords( 3 * num_nodes );
        if( MB_SUCCESS != mb->get_coords( conn, num_nodes, &coords[0] ) ) return -1;

        if( vos.area_spherical_polygon( &coords[0], num_nodes, R ) < 0 ) inverted++;
    }
    return inverted;
}

// Build a single spherical quad on the unit sphere, wound counter-clockwise as
// seen from outside (i.e. a valid, positively oriented cell).
static ErrorCode build_unit_quad( Interface* mb, EntityHandle& set, EntityHandle& quad )
{
    MB_CHK_ERR( mb->create_meshset( MESHSET_SET, set ) );

    double coords[12] = { 1.0, 0.0, 0.0, 0.9, 0.3, 0.0, 0.88, 0.28, 0.3, 0.95, 0.0, 0.3 };
    for( int i = 0; i < 4; i++ )
    {
        double* p    = coords + 3 * i;
        const double n = std::sqrt( p[0] * p[0] + p[1] * p[1] + p[2] * p[2] );
        p[0] /= n;
        p[1] /= n;
        p[2] /= n;
    }

    EntityHandle verts[4];
    for( int i = 0; i < 4; i++ )
        MB_CHK_ERR( mb->create_vertex( coords + 3 * i, verts[i] ) );

    MB_CHK_ERR( mb->create_element( MBQUAD, verts, 4, quad ) );
    MB_CHK_ERR( mb->add_entities( set, &quad, 1 ) );

    return MB_SUCCESS;
}

/*
 * positive_orientation() must repair a mesh whose cells are all wound backwards,
 * and must be a no-op when run a second time.
 *
 * Both input files are TempestRemap ".g" meshes, so *every* cell starts inverted:
 * 1200 cells in the RLL mesh, 1350 in the CS mesh.  Before the fix this pair
 * produced 7570 "negative area" lines while doing exactly the right thing.
 */
void test_positive_orientation_repairs_and_is_idempotent()
{
    const double R = 1.0;
    Core moab;
    Interface* mb = &moab;

    EntityHandle rll, cs;
    CHECK_ERR( mb->create_meshset( MESHSET_SET, rll ) );
    CHECK_ERR( mb->create_meshset( MESHSET_SET, cs ) );
    CHECK_ERR( mb->load_file( ( TestDir + "unittest/mbcslam/outRLLMesh.g" ).c_str(), &rll ) );
    CHECK_ERR( mb->load_file( ( TestDir + "unittest/mbcslam/outCSMesh.g" ).c_str(), &cs ) );

    CHECK_ERR( IntxUtils::fix_degenerate_quads( mb, rll ) );

    // Sanity check on the premise: these meshes really are stored fully inverted.
    Range rll_cells, cs_cells;
    CHECK_ERR( mb->get_entities_by_dimension( rll, 2, rll_cells ) );
    CHECK_ERR( mb->get_entities_by_dimension( cs, 2, cs_cells ) );
    CHECK( !rll_cells.empty() );
    CHECK( !cs_cells.empty() );
    CHECK_EQUAL( (int)rll_cells.size(), count_inverted_cells( mb, rll, R ) );
    CHECK_EQUAL( (int)cs_cells.size(), count_inverted_cells( mb, cs, R ) );

    IntxAreaUtils areaAdaptor;

    // First pass repairs everything.
    CHECK_ERR( areaAdaptor.positive_orientation( mb, rll, R ) );
    CHECK_ERR( areaAdaptor.positive_orientation( mb, cs, R ) );
    CHECK_EQUAL( 0, count_inverted_cells( mb, rll, R ) );
    CHECK_EQUAL( 0, count_inverted_cells( mb, cs, R ) );

    // Second pass must find nothing left to repair.
    CHECK_ERR( areaAdaptor.positive_orientation( mb, rll, R ) );
    CHECK_ERR( areaAdaptor.positive_orientation( mb, cs, R ) );
    CHECK_EQUAL( 0, count_inverted_cells( mb, rll, R ) );
    CHECK_EQUAL( 0, count_inverted_cells( mb, cs, R ) );

    // A correctly oriented closed surface integrates to 4*pi*R^2.  This is the
    // check that would catch an over-eager repair that flipped cells it shouldn't.
    const double exact = 4.0 * M_PI * R * R;
    CHECK_REAL_EQUAL( exact, areaAdaptor.area_on_sphere( mb, rll, R ), 1.0e-12 );
    CHECK_REAL_EQUAL( exact, areaAdaptor.area_on_sphere( mb, cs, R ), 1.0e-12 );
}

/*
 * The negative-area diagnostic must discriminate: silent for a valid cell, and
 * reporting for a genuinely inverted one.  A whole element carrying negative area
 * is always a defect, unlike an individual fan sub-triangle of a concave cell.
 */
void test_negative_area_reported_only_for_inverted_elements()
{
    Core moab;
    Interface* mb = &moab;

    EntityHandle set, quad;
    CHECK_ERR( build_unit_quad( mb, set, quad ) );

    IntxAreaUtils areaAdaptor;

    // Capture stdout so the diagnostic itself can be asserted on.
    std::ostringstream captured;
    std::streambuf* saved = std::cout.rdbuf( captured.rdbuf() );

    const double area_valid = areaAdaptor.area_spherical_element( mb, quad, 1.0 );
    const std::string out_valid = captured.str();

    // Invert the cell by reversing its connectivity.
    const EntityHandle* conn = nullptr;
    int num_nodes            = 0;
    ErrorCode rval           = mb->get_connectivity( quad, conn, num_nodes );
    if( MB_SUCCESS != rval ) std::cout.rdbuf( saved );
    CHECK_ERR( rval );

    std::vector< EntityHandle > reversed( num_nodes );
    for( int i = 0; i < num_nodes; i++ )
        reversed[num_nodes - 1 - i] = conn[i];
    rval = mb->set_connectivity( quad, &reversed[0], num_nodes );
    if( MB_SUCCESS != rval ) std::cout.rdbuf( saved );
    CHECK_ERR( rval );

    captured.str( "" );
    const double area_inverted   = areaAdaptor.area_spherical_element( mb, quad, 1.0 );
    const std::string out_inverted = captured.str();

    std::cout.rdbuf( saved );

    // The valid cell has positive area and produces no output at all.
    CHECK( area_valid > 0.0 );
    CHECK( out_valid.empty() );

    // Reversing connectivity negates the area exactly.
    CHECK_REAL_EQUAL( -area_valid, area_inverted, 1.0e-15 );

    // The inverted cell is reported, and the report identifies the element.
    CHECK( out_inverted.find( "negative area" ) != std::string::npos );
    std::ostringstream idtext;
    idtext << "element " << mb->id_from_handle( quad );
    CHECK( out_inverted.find( idtext.str() ) != std::string::npos );
}

int main( int /*argc*/, char** /*argv*/ )
{
    int failures = 0;

    failures += RUN_TEST( test_positive_orientation_repairs_and_is_idempotent );
    failures += RUN_TEST( test_negative_area_reported_only_for_inverted_elements );

    if( failures ) std::cout << failures << " tests failed\n";

    return failures;
}
