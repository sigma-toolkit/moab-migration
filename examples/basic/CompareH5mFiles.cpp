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


    std::string file1("NE4pg2_imp/run/OcnCplAftMm25.h5m");
    std::string file2("NE4pg2_imp2/run/OcnCplAftMm25.h5m");

    ProgOptions opts;

    int dim = 2;
    opts.addOpt< std::string >( "file1,f", "first file", &file1 );
    opts.addOpt< std::string >( "file2,g", "second file", &file2 );
    opts.addOpt< int  >( "dimension,d", "dimension for entities", &dim );

    opts.parseCommandLine( argc, argv );

    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( file1.c_str() );MB_CHK_SET_ERR( rval, "can't load file1" );

    Core* mb2 = new Core();
    rval      = mb2->load_file( file2.c_str() );MB_CHK_SET_ERR( rval, "can't load file2" );

    std::vector<Tag> list1;
    rval = mb->tag_get_tags(list1);MB_CHK_SET_ERR( rval, "can't get tags 1" );

    Range cells1;
    rval = mb->get_entities_by_dimension(0, dim, cells1);MB_CHK_SET_ERR( rval, "can't get cells 1" );

    Range cells2;
    rval = mb2->get_entities_by_dimension(0, dim, cells2);MB_CHK_SET_ERR( rval, "can't get cells 2" );

    if (cells1.size() != cells2.size()) MB_CHK_SET_ERR( MB_FAILURE, "different size models " );
    vector<double> vals1, vals2;
    vals1.resize(cells1.size());
    vals2.resize(cells2.size());
    int k=0; // number of different fields
    int k1=0; // number of exactly the same fields
    std::cout << " compare files: " << file1 << " and " << file2 << " dimension entity: " << dim << "\n";
    for (size_t i=0; i< list1.size(); i++)
    {
        Tag tag=list1[i];
        std::string name;
        rval = mb->tag_get_name(tag, name);MB_CHK_SET_ERR( rval, "can't get tag name" );
        DataType type;
        rval = mb->tag_get_data_type(tag, type);MB_CHK_SET_ERR( rval, "can't get tag data type" );
        if (MB_TYPE_DOUBLE != type) continue;
        TagType tag_type;
        rval = mb->tag_get_type(tag, tag_type);MB_CHK_SET_ERR( rval, "can't get tag type" );
        if (MB_TAG_DENSE != tag_type) continue;
        int length = 0;
        rval = mb->tag_get_length(tag, length);MB_CHK_SET_ERR( rval, "can't get tag length" );
        if (1 != length) continue;
        Tag tag2;
        rval = mb2->tag_get_handle(name.c_str(), tag2); MB_CHK_SET_ERR( rval, "can't get tag on second model" );
        rval = mb->tag_get_data(tag, cells1, &vals1[0]);MB_CHK_SET_ERR( rval, "can't get values on tag on model 1" );
        rval = mb2->tag_get_data(tag2, cells2, &vals2[0]);MB_CHK_SET_ERR( rval, "can't get values on tag on model 2" );
        // compute the difference
        double sum = 0;
        for (int j=0; j<vals1.size(); j++)
        {
            sum += fabs(vals1[j] - vals2[j]);
        }
        if (sum > 0.)
        {
            std::cout<<" tag: " << name << " \t difference : "<< sum <<"\n";
            k++;
        }
	else
        {
            k1++;
        }
    }
    std::cout<<" different fields:" << k << " \n exactly the same fields:" << k1 << "\n";
    // Cr
    delete mb;
    delete mb2;
}
