#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include <cmath>
#include <algorithm>

using namespace moab;

/* Input test file: rtttest_v100.rtt
 */

std::string example1 = TestDir + "unittest/io/rtttest_v100.rtt";
std::string example2 = TestDir + "unittest/io/rtttest_v101.rtt";
std::string example3 = TestDir + "unittest/io/rtttest_v101_2.rtt";
std::string example4 = TestDir + "unittest/io/rtttest_mcnp_mesh.rtt";
std::string example5 = TestDir + "unittest/io/rtttest_discontiguous.rtt";

void test_loadfile_1();
void test_meshset_tags_1();
void test_tets_1();
void test_tet_tags_1();
void test_triangles_1();
void test_triangles_tags_1();
void test_vertices_1();

void test_loadfile_2();
void test_meshset_tags_2();
void test_tets_2();
void test_tet_tags_2();
void test_triangles_2();
void test_triangles_tags_2();
void test_vertices_2();

void test_loadfile_3();
void test_tets_3();
void test_tet_tags_3();

void test_loadfile_4();
void test_tets_4();
void test_tet_tags_4();

void test_loadfile_5();
void test_tets_5();
void test_no_material_tag_5();
void test_surfaces_generated_5();
void test_volume_surface_links_5();

void read_file( Interface& moab, const char* input_file );

int main()
{
    int result = 0;
    // first batch
    result += RUN_TEST( test_loadfile_1 );
    result += RUN_TEST( test_meshset_tags_1 );
    result += RUN_TEST( test_tets_1 );
    result += RUN_TEST( test_tet_tags_1 );
    result += RUN_TEST( test_triangles_1 );
    result += RUN_TEST( test_vertices_1 );
    // second batch
    result += RUN_TEST( test_loadfile_2 );
    result += RUN_TEST( test_meshset_tags_2 );
    result += RUN_TEST( test_tets_2 );
    result += RUN_TEST( test_tet_tags_2 );
    result += RUN_TEST( test_triangles_2 );
    result += RUN_TEST( test_vertices_2 );
    // third batch
    result += RUN_TEST( test_loadfile_3 );
    result += RUN_TEST( test_tets_3 );
    result += RUN_TEST( test_tet_tags_3 );
    // fourth batch
    result += RUN_TEST( test_loadfile_4 );
    result += RUN_TEST( test_tets_4 );
    result += RUN_TEST( test_tet_tags_4 );
    // fifth batch - discontiguous mesh, no facets, no MATERIAL flag
    result += RUN_TEST( test_loadfile_5 );
    result += RUN_TEST( test_tets_5 );
    result += RUN_TEST( test_no_material_tag_5 );
    result += RUN_TEST( test_surfaces_generated_5 );
    result += RUN_TEST( test_volume_surface_links_5 );
    return result;
}

void read_file( Interface& moab, const char* input_file )
{
    ErrorCode rval;
    rval = moab.load_file( input_file );CHECK_ERR( rval );
}

void test_loadfile_1()
{
    Core moab;
    read_file( moab, example1.c_str() );
}

void test_meshset_tags_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBENTITYSET, entities );CHECK_ERR( rval );

    Tag id_tag = moab.globalId_tag();

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &id_tag, 0, 1, entities );CHECK_ERR( rval );

    // each tet should have the material tag
    int num_vols     = 10;
    int num_surfaces = 129;
    int num_sets     = entities.size();
    CHECK_EQUAL( num_sets, num_vols + num_surfaces + 1 + 1 );

    // from the entities, get the volume meshsets, get it from the dim tag
    Tag dim_tag;
    entities.clear();
    // get the tag handle
    rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, dim_tag, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged with dim 3
    int dim3                          = 3;
    const void* const tag_vals_dim3[] = { &dim3 };
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim3, 1, entities );CHECK_ERR( rval );

    // there should only be 10 meshsets of dimension 3
    num_sets = entities.size();
    CHECK_EQUAL( num_vols, num_sets );

    entities.clear();
    // get the entities that are tagged with dim 2
    int dim2                          = 2;
    const void* const tag_vals_dim2[] = { &dim2 };
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim2, 1, entities );CHECK_ERR( rval );

    // there should only be 129 meshsets of dimension 2
    num_sets = entities.size();
    CHECK_EQUAL( num_surfaces, num_sets );
}

void test_tets_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    // cells = 26710 - number of tets
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );
    int num_tets        = 26710;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );
}

void test_tet_tags_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );

    int num_tets        = 26710;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );

    // rtttest_v100.rtt has no MATERIAL cell flag, so MATERIAL_NUMBER should not be set
    entities.clear();
    Tag material_number;
    rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER, material_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTET, &material_number, 0, 1, entities );CHECK_ERR( rval );
    CHECK_EQUAL( (int)entities.size(), 0 );
}

void test_triangles_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTRI, entities );CHECK_ERR( rval );

    int num_tri         = 6383;  // num tris annotated in rtttest.rtt
    int num_tri_in_moab = entities.size();
    CHECK_EQUAL( num_tri_in_moab, num_tri );
}

void test_triangles_tags_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTRI, entities );CHECK_ERR( rval );

    int num_tri         = 6383;  // num tris annotated in rtttest.rtt
    int num_tri_in_moab = entities.size();
    CHECK_EQUAL( num_tri_in_moab, num_tri );

    // get the number of tris tagged with SURFACE_NUMBER
    entities.clear();
    Tag surface_number;
    // get the tag handle
    rval = moab.tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER, surface_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTRI, &surface_number, 0, 1, entities );
    // each tri should have the surface number tag
    int num_tri_tag = entities.size();
    CHECK_EQUAL( num_tri_tag, num_tri );CHECK_ERR( rval );

    // get the number of tris tagged with SIDEID_TAG
    entities.clear();
    Tag sideid_tag;
    // get the tag handle
    rval = moab.tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER, sideid_tag, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTRI, &sideid_tag, 0, 1, entities );
    // each tri should have the sideid tag
    num_tri_tag = entities.size();
    CHECK_EQUAL( num_tri_tag, num_tri );CHECK_ERR( rval );
}

void test_vertices_1()
{
    Core moab;
    // load the data into moab
    read_file( moab, example1.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBVERTEX, entities );CHECK_ERR( rval );

    int num_verts         = 5397;  // num verts annotated in rtttest.rtt
    int num_verts_in_moab = entities.size();
    CHECK_EQUAL( num_verts_in_moab, num_verts );
}

void test_loadfile_2()
{
    Core moab;
    read_file( moab, example2.c_str() );
}

void test_meshset_tags_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBENTITYSET, entities );CHECK_ERR( rval );

    Tag id_tag = moab.globalId_tag();

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &id_tag, 0, 1, entities );CHECK_ERR( rval );

    // each tet should have the material tag
    int num_vols     = 3;
    int num_surfaces = 24;
    int num_sets     = entities.size();
    CHECK_EQUAL( num_sets, num_vols + num_surfaces + 1 + 1 );

    // from the entities, get the volume meshsets, get it from the dim tag
    Tag dim_tag;
    entities.clear();
    // get the tag handle
    rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, dim_tag, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged with dim 3
    int dim3                          = 3;
    const void* const tag_vals_dim3[] = { &dim3 };
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim3, 1, entities );CHECK_ERR( rval );

    // there should only be 10 meshsets of dimension 3
    num_sets = entities.size();
    CHECK_EQUAL( num_vols, num_sets );

    entities.clear();
    // get the entities that are tagged with dim 2
    int dim2                          = 2;
    const void* const tag_vals_dim2[] = { &dim2 };
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim2, 1, entities );CHECK_ERR( rval );

    // there should only be 129 meshsets of dimension 2
    num_sets = entities.size();
    CHECK_EQUAL( num_surfaces, num_sets );
}

void test_tets_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    // cells = 26710 - number of tets
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );
    int num_tets        = 84;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );
}

void test_tet_tags_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );

    int num_tets        = 84;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );

    // rtttest_v101.rtt has no MATERIAL cell flag, so MATERIAL_NUMBER should not be set
    entities.clear();
    Tag material_number;
    rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER, material_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTET, &material_number, 0, 1, entities );CHECK_ERR( rval );
    CHECK_EQUAL( (int)entities.size(), 0 );
}

void test_triangles_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTRI, entities );CHECK_ERR( rval );

    int num_tri         = 60;  // num tris annotated in rtttest.rtt
    int num_tri_in_moab = entities.size();
    CHECK_EQUAL( num_tri_in_moab, num_tri );
}

void test_triangles_tags_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTRI, entities );CHECK_ERR( rval );

    int num_tri         = 60;  // num tris annotated in rtttest.rtt
    int num_tri_in_moab = entities.size();
    CHECK_EQUAL( num_tri_in_moab, num_tri );

    // get the number of tris tagged with SURFACE_NUMBER
    entities.clear();
    Tag surface_number;
    // get the tag handle
    rval = moab.tag_get_handle( "SURFACE_NUMBER", 1, MB_TYPE_INTEGER, surface_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTRI, &surface_number, 0, 1, entities );
    // each tri should have the surface number tag
    int num_tri_tag = entities.size();
    CHECK_EQUAL( num_tri_tag, num_tri );CHECK_ERR( rval );

    // get the number of tris tagged with SIDEID_TAG
    entities.clear();
    Tag sideid_tag;
    // get the tag handle
    rval = moab.tag_get_handle( "SIDEID_TAG", 1, MB_TYPE_INTEGER, sideid_tag, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTRI, &sideid_tag, 0, 1, entities );
    // each tri should have the sideid tag
    num_tri_tag = entities.size();
    CHECK_EQUAL( num_tri_tag, num_tri );CHECK_ERR( rval );
}

void test_vertices_2()
{
    Core moab;
    // load the data into moab
    read_file( moab, example2.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBVERTEX, entities );CHECK_ERR( rval );

    int num_verts         = 40;  // num verts annotated in rtttest.rtt
    int num_verts_in_moab = entities.size();
    CHECK_EQUAL( num_verts_in_moab, num_verts );
}

void test_loadfile_3()
{
    Core moab;
    read_file( moab, example3.c_str() );
}

void test_tets_3()
{
    Core moab;
    // load the data into moab
    read_file( moab, example3.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    // cells = 26710 - number of tets
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );
    int num_tets        = 84;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );
}

void test_tet_tags_3()
{
    Core moab;
    // load the data into moab
    read_file( moab, example3.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );

    int num_tets        = 84;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );

    // rtttest_v101_2.rtt has no MATERIAL cell flag, so MATERIAL_NUMBER should not be set
    entities.clear();
    Tag material_number;
    rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER, material_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTET, &material_number, 0, 1, entities );CHECK_ERR( rval );
    CHECK_EQUAL( (int)entities.size(), 0 );
}

void test_loadfile_4()
{
    Core moab;
    read_file( moab, example4.c_str() );
}

void test_tets_4()
{
    Core moab;
    // load the data into moab
    read_file( moab, example4.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    // cells = 26710 - number of tets
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );
    int num_tets        = 30;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );
}

void test_tet_tags_4()
{
    Core moab;
    // load the data into moab
    read_file( moab, example4.c_str() );
    // query the dataset to make sure that there are the correct number of cells
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );

    int num_tets        = 30;
    int num_tet_in_moab = entities.size();
    CHECK_EQUAL( num_tet_in_moab, num_tets );

    // get the number of tets tagged with
    entities.clear();
    Tag material_number;
    // get the tag handle
    rval = moab.tag_get_handle( "MATERIAL_NUMBER", 1, MB_TYPE_INTEGER, material_number, MB_TAG_SPARSE | MB_TAG_CREAT );CHECK_ERR( rval );

    // get the entities that are tagged
    rval = moab.get_entities_by_type_and_tag( 0, moab::MBTET, &material_number, 0, 1, entities );
    // each tet should have the material tag
    int num_tet_tag = entities.size();
    CHECK_EQUAL( num_tet_tag, num_tets );CHECK_ERR( rval );
}

// Fifth batch: discontiguous mesh with no sides/facets and no MATERIAL flag.
// Tests that the reader handles this case via tet skinning to produce surfaces.

void test_loadfile_5()
{
    Core moab;
    read_file( moab, example5.c_str() );
}

void test_tets_5()
{
    Core moab;
    read_file( moab, example5.c_str() );
    Range entities;
    ErrorCode rval = moab.get_entities_by_type( 0, moab::MBTET, entities );CHECK_ERR( rval );
    CHECK_EQUAL( (int)entities.size(), 2 );
}

void test_no_material_tag_5()
{
    // The file has no MATERIAL cell flag, so no MATERIAL_NUMBER tags should exist.
    Core moab;
    read_file( moab, example5.c_str() );

    Tag material_number;
    // this test loads a file with no MATERIAL_NUMBER tag, so the tag should not be found
    ErrorCode rval = moab.tag_get_handle( "MATERIAL_NUMBER", material_number );
    CHECK_EQUAL(rval, MB_TAG_NOT_FOUND);
}

void test_surfaces_generated_5()
{
    // With no facet data, surfaces should be generated by skinning each volume's tets.
    // Two volumes of 1 tet each -> 2 surface meshsets and 8 skin triangles (4 per tet).
    Core moab;
    read_file( moab, example5.c_str() );

    Tag dim_tag;
    CHECK_ERR(moab.tag_get_handle(GEOM_DIMENSION_TAG_NAME, dim_tag));

    int dim2                          = 2;
    const void* const tag_vals_dim2[] = { &dim2 };
    Range surf_sets;
    ErrorCode rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim2, 1, surf_sets );CHECK_ERR( rval );
    CHECK_EQUAL( (int)surf_sets.size(), 2 );

    // Each surface should contain 4 skin triangles (all faces of a single tet)
    Range all_tris;
    for( Range::iterator it = surf_sets.begin(); it != surf_sets.end(); ++it )
    {
        Range tris;
        rval = moab.get_entities_by_type( *it, moab::MBTRI, tris );CHECK_ERR( rval );
        CHECK_EQUAL( (int)tris.size(), 4 );
        all_tris.merge( tris );
    }
    CHECK_EQUAL( (int)all_tris.size(), 8 );
}

void test_volume_surface_links_5()
{
    // Each volume should have exactly one surface child (its skin surface).
    Core moab;
    read_file( moab, example5.c_str() );

    Tag dim_tag;
    CHECK_ERR(moab.tag_get_handle(GEOM_DIMENSION_TAG_NAME, dim_tag));

    int dim3                          = 3;
    const void* const tag_vals_dim3[] = { &dim3 };
    Range vol_sets;
    ErrorCode rval = moab.get_entities_by_type_and_tag( 0, moab::MBENTITYSET, &dim_tag, tag_vals_dim3, 1, vol_sets );CHECK_ERR( rval );
    CHECK_EQUAL( (int)vol_sets.size(), 2 );

    for( Range::iterator it = vol_sets.begin(); it != vol_sets.end(); ++it )
    {
        Range children;
        rval = moab.get_child_meshsets( *it, children );CHECK_ERR( rval );
        CHECK_EQUAL( (int)children.size(), 1 );
    }
}