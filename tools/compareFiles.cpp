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
 * ./mbcmpfiles -i file1.h5m -j file2.h5m -n <tag_name>  -o out.file
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

#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile1, inputfile2, outfile;

    std::string tag_name;
    int dim = 2;

    opts.addOpt< std::string >( "input1,i", "input mesh filename 1", &inputfile1 );
    opts.addOpt< std::string >( "input2,j", "input mesh filename 2", &inputfile2 );
    opts.addOpt< std::string >( "tagname,n", "tag to compare (compare all tags if not specified)", &tag_name );
    opts.addOpt< int >( "dimension,d", "topological dimension of entities to compare", &dim );
    opts.addOpt< std::string >( "outfile,o", "output file with differences", &outfile );

    opts.parseCommandLine( argc, argv );

    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( inputfile1.c_str() );MB_CHK_SET_ERR( rval, "can't load input file 1" );

    Core* mb2 = new Core();
    rval      = mb2->load_file( inputfile2.c_str() );MB_CHK_SET_ERR( rval, "can't load input file 2" );

    std::cout << " opened " << inputfile1 << " and " << inputfile2 << " with initial h5m data.\n";
    // open the netcdf file, and see if it has that variable we are looking for

    Range nodes;
    rval = mb->get_entities_by_dimension( 0, 0, nodes );MB_CHK_SET_ERR( rval, "can't get nodes" );

    Range edges;
    rval = mb->get_entities_by_dimension( 0, 1, edges );MB_CHK_SET_ERR( rval, "can't get edges" );

    Range cells;
    rval = mb->get_entities_by_dimension( 0, 2, cells );MB_CHK_SET_ERR( rval, "can't get cells" );

    Range solids;
    rval = mb->get_entities_by_dimension( 0, 3, solids );MB_CHK_SET_ERR( rval, "can't get cells" );

    std::cout << inputfile1 << " has " << nodes.size() << " vertices " << edges.size() << " edges " << cells.size()
              << " cells " << solids.size() << " solids \n";

    // construct maps between global id and handles

    std::map< int, EntityHandle > cGidHandle;
    std::vector< int > gids;
    Tag gid;
    rval = mb->tag_get_handle( "GLOBAL_ID", gid );MB_CHK_SET_ERR( rval, "can't get global id tag" );

    Range ents = cells;
    if( dim == 0 ) ents = nodes;
    if( dim == 1 ) ents = edges;
    if( dim == 3 ) ents = solids;
    gids.resize( ents.size() );
    rval = mb->tag_get_data( gid, ents, &gids[0] );MB_CHK_SET_ERR( rval, "can't get global id on entities" );

    int i = 0;
    for( Range::iterator vit = ents.begin(); vit != ents.end(); ++vit )
    {
        cGidHandle[gids[i++]] = *vit;
    }

    Range nodes2;
    rval = mb2->get_entities_by_dimension( 0, 0, nodes2 );MB_CHK_SET_ERR( rval, "can't get nodes2" );

    Range edges2;
    rval = mb2->get_entities_by_dimension( 0, 1, edges2 );MB_CHK_SET_ERR( rval, "can't get edges2" );

    Range cells2;
    rval = mb2->get_entities_by_dimension( 0, 2, cells2 );MB_CHK_SET_ERR( rval, "can't get cells2" );

    Range solids2;
    rval = mb2->get_entities_by_dimension( 0, 3, solids2 );MB_CHK_SET_ERR( rval, "can't get solids2" );

    std::cout << inputfile2 << " has " << nodes2.size() << " vertices " << edges2.size() << " edges " << cells2.size()
              << " cells " << solids2.size() << " solids \n";

    Range ents2 = cells2;
    if( dim == 0 ) ents2 = nodes2;
    if( dim == 1 ) ents2 = edges2;
    if( dim == 3 ) ents2 = solids2;
    // construct maps between global id and handles
    std::map< int, EntityHandle > cGidHandle2;

    Tag gid2;
    rval = mb2->tag_get_handle( "GLOBAL_ID", gid2 );MB_CHK_SET_ERR( rval, "can't get global id tag2" );

    std::vector< int > gids2( ents2.size() );
    rval = mb2->tag_get_data( gid2, ents2, &gids2[0] );MB_CHK_SET_ERR( rval, "can't get global id on second entities" );
    i = 0;
    for( Range::iterator vit = ents2.begin(); vit != ents2.end(); ++vit )
    {
        cGidHandle2[gids2[i++]] = *vit;
    }

    if( ents.size() != ents2.size() )
    {
        std::cout << "cannot compare tags, because number of entities is different:" << ents.size() << " "
                  << ents2.size() << "\n";
        exit( 1 );
    }
    if( tag_name.length() > 0 )  // old tool
    {
        Tag tag;
        rval = mb->tag_get_handle( tag_name.c_str(), tag );MB_CHK_SET_ERR( rval, "can't get tag on file 1" );

        int len_tag = 0;
        rval        = mb->tag_get_length( tag, len_tag );MB_CHK_SET_ERR( rval, "can't get tag length on tag" );
        std::cout << "length tag : " << len_tag << "\n";
        moab::DataType dtype;
        rval = mb->tag_get_data_type( tag, dtype );MB_CHK_SET_ERR( rval, "can't get tag data type" );
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
            vals.resize( len_tag * ents.size() );
            rval = mb->tag_get_data( tag, ents, &vals[0] );MB_CHK_SET_ERR( rval, "can't get tag data on double tag" );
        }
        else
        {
            ivals.resize( len_tag * ents.size() );
            rval = mb->tag_get_data( tag, ents, &ivals[0] );MB_CHK_SET_ERR( rval, "can't get tag data on integer tag" );
        }

        Tag tag2;
        rval = mb2->tag_get_handle( tag_name.c_str(), tag2 );MB_CHK_SET_ERR( rval, "can't get tag on file 2" );
        std::vector< double > vals2;
        std::vector< int > ivals2;
        if( doubleType )
        {
            vals2.resize( len_tag * ents2.size() );
            rval = mb2->tag_get_data( tag2, ents2, &vals2[0] );MB_CHK_SET_ERR( rval, "can't get tag data on file 2" );
        }
        else
        {
            ivals2.resize( len_tag * ents2.size() );
            rval = mb2->tag_get_data( tag2, ents2, &ivals2[0] );MB_CHK_SET_ERR( rval, "can't get tag data on file 2" );
        }
        std::string new_tag_name = tag_name + "_2";
        Tag newTag, newTagDiff;
        std::string tag_name_diff = tag_name + "_diff";
        if( doubleType )
        {
            double def_val2 = 0.;
            double def_vald[len_tag];
            rval = mb->tag_get_default_value( tag, def_vald );MB_CHK_SET_ERR( rval, "can't get default" );

            rval = mb->tag_get_handle( new_tag_name.c_str(), 1, dtype, newTag, MB_TAG_CREAT | MB_TAG_DENSE, def_vald );MB_CHK_SET_ERR( rval, "can't define new tag" );
            rval = mb->tag_get_handle( tag_name_diff.c_str(), 1, dtype, newTagDiff,
                                       MB_TAG_CREAT | MB_TAG_DENSE | MB_TAG_DFTOK, &def_val2 );MB_CHK_SET_ERR( rval, "can't define new tag diff" );
        }
        else
        {
            int def_vali[len_tag];
            int dev_vali2             = 0;
            rval = mb->tag_get_default_value( tag, def_vali );MB_CHK_SET_ERR( rval, "can't get default" );

            rval = mb->tag_get_handle( new_tag_name.c_str(), 1, dtype, newTag, MB_TAG_CREAT | MB_TAG_DENSE, def_vali );MB_CHK_SET_ERR( rval, "can't define new tag" );
            rval = mb->tag_get_handle( tag_name_diff.c_str(), 1, dtype, newTagDiff,
                                       MB_TAG_CREAT | MB_TAG_DENSE | MB_TAG_DFTOK, &dev_vali2 );MB_CHK_SET_ERR( rval, "can't define new tag diff" );
        }

        i             = 0;
        double l2norm = 0;
        for( Range::iterator c2it = ents2.begin(); c2it != ents2.end(); ++c2it )
        {
            double val2;
            int ival2;
            if( doubleType )
                val2 = vals2[i];
            else
                ival2 = ivals2[i];

            int id2 = gids2[i];
            i++;
            EntityHandle c1 = cGidHandle[id2];
            if( doubleType )
            {
                rval = mb->tag_set_data( newTag, &c1, 1, &val2 );MB_CHK_SET_ERR( rval, "can't set new tag" );
            }
            else
            {
                rval = mb->tag_set_data( newTag, &c1, 1, &ival2 );MB_CHK_SET_ERR( rval, "can't set new tag" );
            }
            int indx    = ents.index( c1 );
            if( doubleType )
            {
                double diff = vals[indx] - val2;
                rval = mb->tag_set_data( newTagDiff, &c1, 1, &diff );MB_CHK_SET_ERR( rval, "can't set new tag" );
                l2norm += diff * diff;
            }
            else
            {
                int diffi   = ivals[indx] - ival2;
                rval  = mb->tag_set_data( newTagDiff, &c1, 1, &diffi );MB_CHK_SET_ERR( rval, "can't set new tag" );
                l2norm += diffi * diffi;
            }
        }
        l2norm = sqrt( l2norm );

        if (!outfile.empty())
        {
            rval = mb->write_file( outfile.c_str() );MB_CHK_SET_ERR( rval, "can't write file" );
            std::cout << " wrote file " << outfile << "\n";
        }
        std::cout << " l2norm of the diff: " << l2norm << "\n";

    }
    else  // look at all tags that can be compared on ents
    {
        // compare all tags
        std::vector< Tag > list1;
        rval = mb->tag_get_tags( list1 );MB_CHK_SET_ERR( rval, "can't get tags 1" );

        std::map< int, int > gidMap2;
        for( int i = 0; i < ents2.size(); i++ )
        {
            gidMap2[gids2[i]] = i;
        }
        int k  = 0;  // number of different fields
        int k1 = 0;  // number of exactly the same fields

        std::vector< std::string > same_fields;
        std::vector< std::string > skipped_fields;
        std::vector< Tag > diffTags;
        for( size_t i = 0; i < list1.size(); i++ )
        {
            Tag tag = list1[i];
            if (tag == gid) continue; // do not compare global id tag
            std::string name;
            rval = mb->tag_get_name( tag, name );MB_CHK_SET_ERR( rval, "can't get tag name" );
            DataType type;
            rval = mb->tag_get_data_type( tag, type );MB_CHK_SET_ERR( rval, "can't get tag data type" );
            if( MB_TYPE_DOUBLE != type && MB_TYPE_INTEGER != type ) continue;
            bool doubleType = ( type == MB_TYPE_DOUBLE );
            TagType tag_type;
            rval = mb->tag_get_type( tag, tag_type );MB_CHK_SET_ERR( rval, "can't get tag type" );
            if( MB_TAG_DENSE != tag_type ) continue;
            int length = 0;
            rval       = mb->tag_get_length( tag, length );MB_CHK_SET_ERR( rval, "can't get tag length" );
            if( 1 != length ) continue;
            Tag tag2;
            //std::cout <<" tag : " << name << "\n";
            rval = mb2->tag_get_handle( name.c_str(), tag2 );MB_CHK_SET_ERR( rval, "can't get tag on second model" );
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

            double minv1, maxv1, minv2, maxv2;
            if( vals1.size() > 0 )
            {
                minv1 = maxv1 = vals1[0];
            }
            if( vals2.size() > 0 )
            {
                minv2 = maxv2 = vals2[0];
            }
            if( ivals1.size() > 0 )
            {
                minv1 = maxv1 = ivals1[0];
            }
            if( ivals2.size() > 0 )
            {
                minv2 = maxv2 = ivals2[0];
            }
            // compute the difference
            double sum = 0;
            for( int j = 0; j < ents.size(); j++ )
            {
                double value1, value2;
                if( doubleType )
                    value1 = vals1[j];
                else
                    value1 = ivals1[j];

                int index2 = gidMap2[gids[j]];
                if( doubleType )
                    value2 = vals2[index2];
                else
                    value2 = ivals2[j];

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

                for( int j = 0; j < ents.size(); j++ )
                {
                    int index2 = gidMap2[gids[j]];
                    if( doubleType )
                        vals1[j] -= vals2[index2];
                    else
                        ivals1[j] -= ivals2[index2];
                }

                std::string diffTagName = name + "_diff";
                Tag newTag;
                rval = mb->tag_get_handle( diffTagName.c_str(), 1, type, newTag, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_ERR( rval );
                if( doubleType )
                {
                    rval = mb->tag_set_data( newTag, ents, &vals1[0] );MB_CHK_ERR( rval );
                }
                else
                {
                    rval = mb->tag_set_data( newTag, ents, &ivals1[0] );MB_CHK_ERR( rval );
                }
                diffTags.push_back( newTag );
            }
            else
            {
                same_fields.push_back( name );
                k1++;
            }
        }
        if( k > 0 && !outfile.empty())
        {
            rval = mb->write_file( outfile.c_str(), 0, 0, 0, 0, &diffTags[0], diffTags.size() );MB_CHK_ERR( rval );
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
