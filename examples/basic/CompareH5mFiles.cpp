/*
 * CompareH5mFiles.cpp
 *
 *  Created on: Sep 26, 2023
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReadUtilIface.hpp"
#include "math.h"
#include <map>
#include <iostream>
#include <cassert>

using namespace moab;
using namespace std;

int main( int argc, char** argv )
{

    std::string file1( "NE4pg2_imp/run/OcnCplAftMm25.h5m" );
    std::string file2( "NE4pg2_imp2/run/OcnCplAftMm25.h5m" );

    ProgOptions opts;

    int dim = 2;
    opts.addOpt< std::string >( "file1,f", "first file", &file1 );
    opts.addOpt< std::string >( "file2,g", "second file", &file2 );

    opts.addOpt<void>( "global,G", "use global id for correspondence" );


    opts.parseCommandLine( argc, argv );

    bool use_global_id = opts.numOptSet("global") > 0;

    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( file1.c_str() );MB_CHK_SET_ERR( rval, "can't load file1" );

    Core* mb2 = new Core();
    rval      = mb2->load_file( file2.c_str() );MB_CHK_SET_ERR( rval, "can't load file2" );

    std::vector< Tag > list1;
    rval = mb->tag_get_tags( list1 );MB_CHK_SET_ERR( rval, "can't get tags 1" );

    Tag gidTag = mb->globalId_tag();
    Tag gidTag2 = mb2->globalId_tag();

    Range cells1;
    rval = mb->get_entities_by_dimension( 0, dim, cells1 );MB_CHK_SET_ERR( rval, "can't get cells 1" );
    if( cells1.size() == 0 )
    {
        dim  = 0;
        rval = mb->get_entities_by_dimension( 0, dim, cells1 );MB_CHK_SET_ERR( rval, "can't get cells 1" );
    }


    Range cells2;
    rval = mb2->get_entities_by_dimension( 0, dim, cells2 );MB_CHK_SET_ERR( rval, "can't get cells 2" );

    if( cells1.size() != cells2.size() ) MB_CHK_SET_ERR( MB_FAILURE, "different size models " );
    vector< double > vals1, vals2;
    vals1.resize( cells1.size() );
    vals2.resize( cells2.size() );

    std::vector<int> gids, gids2;
    std::map<int, int> gidMap2;
    if (use_global_id)
    {
        gids.resize(cells1.size());
        rval = mb->tag_get_data( gidTag, cells1, &gids[0] );MB_CHK_SET_ERR( rval, "can't get ids for cells1" );
        gids2.resize(cells2.size());
        rval = mb2->tag_get_data( gidTag2, cells2, &gids2[0] );MB_CHK_SET_ERR( rval, "can't get ids for cells2" );
        for (int i=0; i<cells2.size(); i++)
        {
            gidMap2[ gids2[i] ] = i;
        }
    }

    int k  = 0;  // number of different fields
    int k1 = 0;  // number of exactly the same fields
    std::cout << " compare files: " << file1 << " and " << file2 << " dimension entity: " << dim << "\n";
    std::vector< std::string > same_fields;
    std::vector<Tag> diffTags;
    for( size_t i = 0; i < list1.size(); i++ )
    {
        Tag tag = list1[i];
        std::string name;
        rval = mb->tag_get_name( tag, name );MB_CHK_SET_ERR( rval, "can't get tag name" );
        DataType type;
        rval = mb->tag_get_data_type( tag, type );MB_CHK_SET_ERR( rval, "can't get tag data type" );
        if( MB_TYPE_DOUBLE != type ) continue;
        TagType tag_type;
        rval = mb->tag_get_type( tag, tag_type );MB_CHK_SET_ERR( rval, "can't get tag type" );
        if( MB_TAG_DENSE != tag_type ) continue;
        int length = 0;
        rval       = mb->tag_get_length( tag, length );MB_CHK_SET_ERR( rval, "can't get tag length" );
        if( 1 != length ) continue;
        Tag tag2;
        rval = mb2->tag_get_handle( name.c_str(), tag2 );MB_CHK_SET_ERR( rval, "can't get tag on second model" );
        rval = mb->tag_get_data( tag, cells1, &vals1[0] );MB_CHK_SET_ERR( rval, "can't get values on tag on model 1" );
        rval = mb2->tag_get_data( tag2, cells2, &vals2[0] );MB_CHK_SET_ERR( rval, "can't get values on tag on model 2" );
        double minv1, maxv1, minv2, maxv2;
        if( vals1.size() > 0 )
        {
            minv1 = maxv1 = vals1[0];
        }
        if( vals2.size() > 0 )
        {
            minv2 = maxv2 = vals2[0];
        }
        // compute the difference
        double sum = 0, value1, value2;
        for( int j = 0; j < vals1.size(); j++ )
        {
            value1 = vals1[j];
            if (use_global_id)
            {

                int index2 = gidMap2[ gids[j] ];
                value2 = vals2[index2];

            }
            else
            {
                value2 = vals2[j];


            }
            sum += fabs( value1 - value2 );
            if( value1 < minv1 ) minv1 = value1;
            if( value1 > maxv1 ) maxv1 = value1;
            if( value2 < minv2 ) minv2 = value2;
            if( value2 > maxv2 ) maxv2 = value2;
        }

        if( sum > 0. )
        {
            std::cout << " tag: " << name << " \t difference : " << sum << " \t min/max (" << minv1 << "/" << maxv1
                      << ") \t (" << minv2 << "/" << maxv2 << ") \n";
            k++;
            if (use_global_id)
            {
                for (int j=0; j< vals1.size(); j++)
                {
                    int index2 = gidMap2[ gids[j] ];
                    vals1[j] -= vals2[index2];
                }
            }
            else
            {
                for (int j=0; j< vals1.size(); j++)
                {
                    vals1[j] -= vals2[j];
                }
            }
            std::string diffTagName = name+"_diff";
            Tag newTag;
            rval = mb->tag_get_handle( diffTagName.c_str(), 1, MB_TYPE_DOUBLE, newTag, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_ERR( rval );
            rval = mb->tag_set_data(newTag, cells1, &vals1[0]);MB_CHK_ERR( rval );
            diffTags.push_back(newTag);

        }
        else
        {
            same_fields.push_back( name );
            k1++;
        }
    }
    if (k>0)
    {
        rval = mb->write_file("diff_tags.h5m", 0, 0, 0, 0, &diffTags[0], diffTags.size() ); MB_CHK_ERR( rval );
    }
    std::cout << " different fields:" << k << " \n exactly the same fields:" << k1 << "\n";
    for( size_t i = 0; i < same_fields.size(); i++ )
    {
        std::cout << " " << same_fields[i];
    }
    std::cout << "\n";
    // Cr
    delete mb;
    delete mb2;
}
