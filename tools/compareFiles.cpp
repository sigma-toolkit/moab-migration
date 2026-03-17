/*
 * compareFiles.cpp
 * this tool will take two existing h5m files, for the same mesh;
 *  they will both have the same GLOBAL_IDs for the elements, but the entity handles can be
 *  very different (depending on how the mesh was partitioned, and saved in parallel)
 *
 *  will compare then the difference between tags, and store the result on one of the files (saved
 * again)
 *
 *
 * example of usage:
 * ./mbcmpfiles -i file1.h5m -j file2.h5m -n \c tag_name  -o out.file
 *
 * if no tag name is specified, it will try to compare all tags in the files
 *
 *
 * Basically, will output a new h5m file (out.file), which has an extra tag, corresponding to the
 * difference between the 2 values
 *
 */

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace moab;
using namespace std;

namespace
{

// Simple hashable key from coordinates with tolerance snapping
struct CoordKey
{
    std::array< long long, 3 > v;
    CoordKey() : v() {}
    CoordKey( long long x, long long y, long long z ) : v{{ x, y, z }} {}
    bool operator==( const CoordKey& other ) const { return v == other.v; }
};

struct CoordKeyHash
{
    std::size_t operator()( const CoordKey& k ) const
    {
        // A small hash combiner
        return std::hash< long long >{}( k.v[0] ) ^ ( std::hash< long long >{}( k.v[1] ) << 1 ) ^
               ( std::hash< long long >{}( k.v[2] ) << 2 );
    }
};

CoordKey make_key( const CartVect& c, double tol )
{
    double inv = 1.0 / tol;
    return CoordKey( llround( c[0] * inv ), llround( c[1] * inv ), llround( c[2] * inv ) );
}

CartVect entity_centroid( Interface* mb, EntityHandle eh )
{
    const EntityHandle* conn = nullptr;
    int nnodes               = 0;
    MB_CHK_ERR_RET_VAL( mb->get_connectivity( eh, conn, nnodes ), CartVect( 0.0 ) );
    std::vector< double > coords( 3 * nnodes );
    MB_CHK_ERR_RET_VAL( mb->get_coords( conn, nnodes, coords.data() ), CartVect( 0.0 ) );
    CartVect c( 0.0 );
    for( int i = 0; i < nnodes; ++i )
    {
        c[0] += coords[3 * i + 0];
        c[1] += coords[3 * i + 1];
        c[2] += coords[3 * i + 2];
    }
    if( nnodes > 0 ) c /= static_cast< double >( nnodes );
    return c;
}

// Build a coordinate-key map for either vertices (dim==0) or elements (by centroid)
std::unordered_map< CoordKey, EntityHandle, CoordKeyHash >
build_coordinate_map( Interface* mb, const Range& ents, bool use_centroid, double tol )
{
    std::unordered_map< CoordKey, EntityHandle, CoordKeyHash > key_to_ent;
    key_to_ent.reserve( ents.size() );
    for( Range::const_iterator it = ents.begin(); it != ents.end(); ++it )
    {
        CartVect c;
        if( use_centroid )
            c = entity_centroid( mb, *it );
        else
        {
            ErrorCode rval = mb->get_coords( &( *it ), 1, c.array() );
            if( MB_SUCCESS != rval )
            {
                std::cerr << "Warning: failed to get coords for entity " << mb->id_from_handle( *it ) << "\n";
                continue;
            }
        }

        CoordKey key = make_key( c, tol );
        if( key_to_ent.find( key ) != key_to_ent.end() )
        {
            std::cerr << "Warning: duplicate coordinate key encountered; keeping first occurrence\n";
            continue;
        }
        key_to_ent[key] = *it;
    }
    return key_to_ent;
}

template < typename T >
void update_minmax( const T val, T& minv, T& maxv )
{
    if( val < minv ) minv = val;
    if( val > maxv ) maxv = val;
}

}  // namespace

int main( int argc, char* argv[] )
{
    ProgOptions opts;

    std::string inputfile1, inputfile2, outfile;

    std::string tag_name;
    int dim = 2;
    bool coord_match = false;  // when true, match by coordinates/centroids instead of GLOBAL_ID
    double tol        = 1e-12; // tolerance for coordinate snapping

    opts.addOpt< std::string >( "input1,i", "input mesh filename 1", &inputfile1 );
    opts.addOpt< std::string >( "input2,j", "input mesh filename 2", &inputfile2 );
    opts.addOpt< std::string >( "tagname,n", "tag to compare (compare all tags if not specified)", &tag_name );
    opts.addOpt< int >( "dimension,d", "topological dimension of entities to compare", &dim );
    opts.addOpt< std::string >( "outfile,o", "output file with differences", &outfile );
    opts.addOpt< void >( "coord,c", "match entities by coordinates (verts) or centroids (elements)", &coord_match );
    opts.addOpt< double >( "tolerance,t", "tolerance for coordinate matching", &tol );

    opts.parseCommandLine( argc, argv );

    ErrorCode rval;
    Core* mb = new Core();

    MB_CHK_SET_ERR( mb->load_file( inputfile1.c_str() ), "can't load input file 1" );

    Core* mb2 = new Core();
    MB_CHK_SET_ERR( mb2->load_file( inputfile2.c_str() ), "can't load input file 2" );

    std::cout << " opened " << inputfile1 << " and " << inputfile2 << " with initial h5m data.\n";
    // open the netcdf file, and see if it has that variable we are looking for

    Range nodes;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( 0, 0, nodes ), "can't get nodes" );

    Range edges;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( 0, 1, edges ), "can't get edges" );

    Range cells;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( 0, 2, cells ), "can't get cells" );

    Range solids;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( 0, 3, solids ), "can't get cells" );

    std::cout << inputfile1 << " has " << nodes.size() << " vertices " << edges.size() << " edges " << cells.size()
              << " cells " << solids.size() << " solids \n";

    // construct maps between global id and handles

    std::map< int, EntityHandle > cGidHandle;
    std::vector< int > gids;
    Tag gid = mb->globalId_tag();

    Range ents = cells;
    if( dim == 0 ) ents = nodes;
    if( dim == 1 ) ents = edges;
    if( dim == 3 ) ents = solids;
    gids.resize( ents.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( gid, ents, &gids[0] ), "can't get global id on entities" );

    int i = 0;
    for( Range::iterator vit = ents.begin(); vit != ents.end(); ++vit )
    {
        cGidHandle[gids[i++]] = *vit;
    }

    Range nodes2;
    MB_CHK_SET_ERR( mb2->get_entities_by_dimension( 0, 0, nodes2 ), "can't get nodes2" );

    Range edges2;
    MB_CHK_SET_ERR( mb2->get_entities_by_dimension( 0, 1, edges2 ), "can't get edges2" );

    Range cells2;
    MB_CHK_SET_ERR( mb2->get_entities_by_dimension( 0, 2, cells2 ), "can't get cells2" );

    Range solids2;
    MB_CHK_SET_ERR( mb2->get_entities_by_dimension( 0, 3, solids2 ), "can't get solids2" );

    std::cout << inputfile2 << " has " << nodes2.size() << " vertices " << edges2.size() << " edges " << cells2.size()
              << " cells " << solids2.size() << " solids \n";

    Range ents2 = cells2;
    if( dim == 0 ) ents2 = nodes2;
    if( dim == 1 ) ents2 = edges2;
    if( dim == 3 ) ents2 = solids2;
    // construct maps between global id and handles
    std::map< int, EntityHandle > cGidHandle2;

    Tag gid2 = mb2->globalId_tag();
    std::vector< int > gids2( ents2.size() );
    MB_CHK_SET_ERR( mb2->tag_get_data( gid2, ents2, &gids2[0] ), "can't get global id on second entities" );

    i = 0;
    for( Range::iterator vit = ents2.begin(); vit != ents2.end(); ++vit )
    {
        cGidHandle2[gids2[i++]] = *vit;
    }

    // Build aligned match lists based on either GLOBAL_ID or coordinate keys
    std::vector< EntityHandle > match1, match2;
    match1.reserve( ents.size() );
    match2.reserve( ents2.size() );

    if( coord_match )
    {
        bool use_centroid = ( dim != 0 );
        auto map2         = build_coordinate_map( mb2, ents2, use_centroid, tol );
        for( Range::iterator it = ents.begin(); it != ents.end(); ++it )
        {
            CartVect c;
            if( use_centroid )
                c = entity_centroid( mb, *it );
            else
                MB_CHK_SET_ERR( mb->get_coords( &( *it ), 1, c.array() ), "can't get coords for matching" );
            CoordKey key = make_key( c, tol );
            auto f       = map2.find( key );
            if( f == map2.end() )
            {
                std::cerr << "No coordinate match found for entity " << mb->id_from_handle( *it ) << "\n";
                continue;
            }
            match1.push_back( *it );
            match2.push_back( f->second );
        }
        if( match1.size() != match2.size() )
        {
            std::cerr << "Coordinate matching produced unequal pair counts: " << match1.size() << " vs "
                      << match2.size() << "\n";
        }
    }
    else
    {
        if( ents.size() != ents2.size() )
        {
            std::cout << "cannot compare tags, because number of entities is different:" << ents.size() << " "
                      << ents2.size() << "\n";
            exit( 1 );
        }
        // Align via GLOBAL_ID (original behavior)
        std::map< int, int > gidMap2;
        for( size_t j = 0; j < ents2.size(); j++ )
            gidMap2[gids2[j]] = static_cast< int >( j );
        for( size_t j = 0; j < ents.size(); j++ )
        {
            int gid_val = gids[j];
            if( gidMap2.find( gid_val ) == gidMap2.end() )
            {
                std::cerr << "No GLOBAL_ID match for " << gid_val << "\n";
                continue;
            }
            match1.push_back( ents[j] );
            match2.push_back( ents2[gidMap2[gid_val]] );
        }
    }

    if( match1.empty() )
    {
        std::cerr << "No matching entities to compare.\n";
        return 1;
    }
    if( match1.size() != match2.size() )
    {
        std::cerr << "Mismatched match counts between files: " << match1.size() << " vs " << match2.size() << "\n";
        return 1;
    }

    if( tag_name.length() > 0 )  // old tool
    {
        Tag tag;
        MB_CHK_SET_ERR( mb->tag_get_handle( tag_name.c_str(), tag ), "can't get tag on file 1" );

        int len_tag = 0;
        MB_CHK_SET_ERR( mb->tag_get_length( tag, len_tag ), "can't get tag length on tag" );
        std::cout << "length tag : " << len_tag << "\n";

        moab::DataType dtype;
        MB_CHK_SET_ERR( mb->tag_get_data_type( tag, dtype ), "can't get tag data type" );
        if( dtype != MB_TYPE_INTEGER && dtype != MB_TYPE_DOUBLE )
        {
            std::cout << "tag data type is not integer or double, do not compare \n";
            exit( 1 );
        }

        bool doubleType = ( dtype == MB_TYPE_DOUBLE );
        std::vector< double > vals;
        std::vector< int > ivals;
        if( doubleType )
        {
            vals.resize( len_tag * match1.size() );
            MB_CHK_SET_ERR( mb->tag_get_data( tag, &match1[0], match1.size(), &vals[0] ),
                            "can't get tag data on double tag" );
        }
        else
        {
            ivals.resize( len_tag * match1.size() );
            MB_CHK_SET_ERR( mb->tag_get_data( tag, &match1[0], match1.size(), &ivals[0] ),
                            "can't get tag data on integer tag" );
        }

        Tag tag2;
        MB_CHK_SET_ERR( mb2->tag_get_handle( tag_name.c_str(), tag2 ), "can't get tag on file 2" );
        std::vector< double > vals2;
        std::vector< int > ivals2;
        if( doubleType )
        {
            vals2.resize( len_tag * match2.size() );
            MB_CHK_SET_ERR( mb2->tag_get_data( tag2, &match2[0], match2.size(), &vals2[0] ),
                            "can't get tag data on file 2" );
        }
        else
        {
            ivals2.resize( len_tag * match2.size() );
            MB_CHK_SET_ERR( mb2->tag_get_data( tag2, &match2[0], match2.size(), &ivals2[0] ),
                            "can't get tag data on file 2" );
        }
        std::string new_tag_name = tag_name + "_2";
        Tag newTag, newTagDiff;
        std::string tag_name_diff = tag_name + "_diff";
        if( doubleType )
        {
            std::vector< double > def_vald( len_tag, 0.0 );
            MB_CHK_SET_ERR( mb->tag_get_default_value( tag, &def_vald[0] ), "can't get default double tag value" );
            MB_CHK_SET_ERR( mb->tag_get_handle( new_tag_name.c_str(), len_tag, dtype, newTag,
                                                MB_TAG_CREAT | MB_TAG_DENSE, &def_vald[0] ),
                            "can't define new double tag" );
            MB_CHK_SET_ERR( mb->tag_get_handle( tag_name_diff.c_str(), len_tag, dtype, newTagDiff,
                                                MB_TAG_CREAT | MB_TAG_DENSE | MB_TAG_DFTOK, &def_vald[0] ),
                            "can't define new double tag diff" );
        }
        else
        {
            std::vector< int > def_vali( len_tag, 0 );
            MB_CHK_SET_ERR( mb->tag_get_default_value( tag, &def_vali[0] ), "can't get default integer tag value" );
            // the difference should be the same size tag
            MB_CHK_SET_ERR( mb->tag_get_handle( new_tag_name.c_str(), len_tag, dtype, newTag,
                                                MB_TAG_CREAT | MB_TAG_DENSE, &def_vali[0] ),
                            "can't define new integer tag" );
            MB_CHK_SET_ERR( mb->tag_get_handle( tag_name_diff.c_str(), len_tag, dtype, newTagDiff,
                                                MB_TAG_CREAT | MB_TAG_DENSE | MB_TAG_DFTOK, &def_vali[0] ),
                            "can't define new integer tag diff" );
        }

        double l2norm = 0;
        for( size_t idx = 0; idx < match1.size(); ++idx )
        {
            EntityHandle h1        = match1[idx];
            const double* val2d    = doubleType ? &vals2[idx * len_tag] : nullptr;
            const int* val2i       = doubleType ? nullptr : &ivals2[idx * len_tag];

            if( doubleType )
                MB_CHK_SET_ERR( mb->tag_set_data( newTag, &h1, 1, val2d ), "can't set new tag" );
            else
                MB_CHK_SET_ERR( mb->tag_set_data( newTag, &h1, 1, val2i ), "can't set new tag" );

            int indx = static_cast< int >( idx );  // aligned order
            if( doubleType )
            {
                double* diff = &vals[indx * len_tag];
                for( int k = 0; k < len_tag; k++ )
                {
                    diff[k] -= val2d[k];
                    l2norm += diff[k] * diff[k];
                }
                MB_CHK_SET_ERR( mb->tag_set_data( newTagDiff, &h1, 1, diff ), "can't set diff double tag" );
            }
            else
            {
                int* diffi = &ivals[indx * len_tag];
                for( int k = 0; k < len_tag; k++ )
                {
                    diffi[k] -= val2i[k];
                    l2norm += static_cast< double >( diffi[k] ) * diffi[k];
                }
                MB_CHK_SET_ERR( mb->tag_set_data( newTagDiff, &h1, 1, diffi ), "can't set diff int tag" );
            }
        }
        l2norm = sqrt( l2norm );

        if( !outfile.empty() )
        {
            MB_CHK_SET_ERR( mb->write_file( outfile.c_str() ), "can't write file" );
            std::cout << " wrote file " << outfile << "\n";
        }
        std::cout << " l2norm of the diff: " << l2norm << "\n";
    }
    else  // look at all tags that can be compared on ents
    {
        // compare all tags
        std::vector< Tag > list1;
        MB_CHK_SET_ERR( mb->tag_get_tags( list1 ), "can't get tags 1" );
        int k  = 0;  // number of different fields
        int k1 = 0;  // number of exactly the same fields

        std::vector< std::string > same_fields;
        std::vector< std::string > skipped_fields;
        std::vector< Tag > diffTags;
        for( size_t i = 0; i < list1.size(); i++ )
        {
            Tag tag = list1[i];
            if( tag == gid ) continue;  // do not compare global id tag
            std::string name;
            MB_CHK_SET_ERR( mb->tag_get_name( tag, name ), "can't get tag name" );
            DataType type;
            MB_CHK_SET_ERR( mb->tag_get_data_type( tag, type ), "can't get tag data type" );
            if( MB_TYPE_DOUBLE != type && MB_TYPE_INTEGER != type ) continue;
            bool doubleType = ( type == MB_TYPE_DOUBLE );
            TagType tag_type;
            MB_CHK_SET_ERR( mb->tag_get_type( tag, tag_type ), "can't get tag type" );
            if( MB_TAG_DENSE != tag_type ) continue;
            int length = 0;
            MB_CHK_SET_ERR( mb->tag_get_length( tag, length ), "can't get tag length" );
            if( 1 != length ) continue;
            Tag tag2;
            //std::cout <<" tag : " << name << "\n";
            MB_CHK_SET_ERR( mb2->tag_get_handle( name.c_str(), tag2 ), "can't get tag on second model" );
            std::vector< double > vals1;
            std::vector< int > ivals1;
            if( doubleType )
            {
                vals1.resize( ents.size() );
                rval = mb->tag_get_data( tag, ents, &vals1[0] );
            }
            else
            {
                ivals1.resize( ents.size() );
                rval = mb->tag_get_data( tag, ents, &ivals1[0] );
            }

            if( MB_SUCCESS != rval )
            {
                std::cout << " can't get values for tag " << name << " on model 1; skip it in comparison \n";
                skipped_fields.push_back( name );
                continue;
            }
            std::vector< double > vals2;
            std::vector< int > ivals2;
            if( doubleType )
            {
                vals2.resize( ents2.size() );
                rval = mb2->tag_get_data( tag2, ents2, &vals2[0] );
            }
            else
            {
                ivals2.resize( ents2.size() );
                rval = mb2->tag_get_data( tag2, ents2, &ivals2[0] );
            }

            if( MB_SUCCESS != rval )
            {
                std::cout << " can't get values for tag " << name << " on model 2; skip it in comparison \n";
                skipped_fields.push_back( name );
            }

            double minv1 = std::numeric_limits< double >::max();
            double maxv1 = std::numeric_limits< double >::lowest();
            double minv2 = std::numeric_limits< double >::max();
            double maxv2 = std::numeric_limits< double >::lowest();
            // compute the difference
            double sum = 0;
            for( size_t j = 0; j < match1.size(); j++ )
            {
                double value1 = doubleType ? vals1[j] : static_cast< double >( ivals1[j] );
                double value2 = doubleType ? vals2[j] : static_cast< double >( ivals2[j] );

                sum += fabs( value1 - value2 );
                update_minmax( value1, minv1, maxv1 );
                update_minmax( value2, minv2, maxv2 );
            }

            if( sum > 0. )
            {
                std::cout << " tag: " << name << " \t difference : " << sum << " \t min/max (" << minv1 << "/" << maxv1
                          << ") \t (" << minv2 << "/" << maxv2 << ") \n";
                k++;

                for( size_t j = 0; j < match1.size(); j++ )
                {
                    if( doubleType )
                        vals1[j] -= vals2[j];
                    else
                        ivals1[j] -= ivals2[j];
                }

                std::string diffTagName = name + "_diff";
                Tag newTag;
                MB_CHK_ERR( mb->tag_get_handle( diffTagName.c_str(), 1, type, newTag, MB_TAG_CREAT | MB_TAG_DENSE ) );
                if( doubleType )
                {
                    MB_CHK_ERR( mb->tag_set_data( newTag, &match1[0], match1.size(), &vals1[0] ) );
                }
                else
                {
                    MB_CHK_ERR( mb->tag_set_data( newTag, &match1[0], match1.size(), &ivals1[0] ) );
                }
                diffTags.push_back( newTag );
            }
            else
            {
                same_fields.push_back( name );
                k1++;
            }
        }
        if( k > 0 && !outfile.empty() )
        {
            MB_CHK_ERR( mb->write_file( outfile.c_str(), 0, 0, 0, 0, &diffTags[0], diffTags.size() ) );
            std::cout << " wrote difference file: " << outfile << "\n";
        }
        std::cout << " different fields:" << k << " \n exactly the same fields:" << k1 << "\n";
        std::cout << " number of skipped fields: " << skipped_fields.size() << "\n";
        for( size_t i = 0; i < same_fields.size(); i++ )
        {
            std::cout << " " << same_fields[i];
        }
        std::cout << "\n";
    }
    delete mb;
    delete mb2;
    return 0;
}
