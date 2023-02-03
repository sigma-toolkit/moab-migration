/*
 *
 *
/** @example ExtractFewTags.cpp  extract only few tags from a file, to reduce its size
 *
 * example of usage:
 * ./ExtractFewTags -i OcnCplAftMm02.h5m -o outfile.h5m -t Faxa_swnet:Foxx_swnet
 *
 * Basically, it will read the file, and write out only the tags in the list (and global id tag) *
 *
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <cmath>
#include <sstream>

using namespace moab;

// Utility function
static void split_tag_names( std::string input_names,
                             std::string& separator,
                             std::vector< std::string >& list_tag_names )
{
    size_t pos = 0;
    std::string token;
    while( ( pos = input_names.find( separator ) ) != std::string::npos )
    {
        token = input_names.substr( 0, pos );
        if( !token.empty() ) list_tag_names.push_back( token );
        // std::cout << token << std::endl;
        input_names.erase( 0, pos + separator.length() );
    }
    if( !input_names.empty() )
    {
        // if leftover something, or if not ended with delimiter
        list_tag_names.push_back( input_names );
    }
    return;
}

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile, outfile( "out.h5m" ), sourcefile, tagname;

    opts.addOpt< std::string >( "input,i", "input mesh filename", &inputfile );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );
    opts.addOpt< std::string >( "tag,t", "tag name", &tagname );

    opts.parseCommandLine( argc, argv );

    std::cout << "input mesh file: " << inputfile << "\n";
    std::cout << "output file: " << outfile << "\n";
    std::cout << "tagname: " << tagname << "\n";

    if( inputfile.empty() )
    {
        opts.printHelp();
        return 0;
    }

    ErrorCode rval;
    Core* mb = new Core();

    rval = mb->load_file( inputfile.c_str() );MB_CHK_SET_ERR( rval, "can't load input file" );

    std::vector< std::string > tagNames;
    std::string separator( ":" );
    split_tag_names( tagname, separator, tagNames );
    if (tagNames.size() < 1)
    {
        std::cout << " no tags identified in " << tagname << "\n";
        return 0;
    }
    std::vector< Tag > tagHandles;
    for (size_t i=0; i<tagNames.size(); i++ )
    {
        Tag tagHandle;
        rval = mb->tag_get_handle( tagNames[i].c_str(), tagHandle );
        if( MB_SUCCESS != rval || NULL == tagHandle )
        {
            std::cout << " can't get tag handle for tag named:" << tagNames[i].c_str() << " at index " << i << "\n";MB_CHK_SET_ERR( rval, "can't get tag handle" );
        }
        tagHandles.push_back( tagHandle );
    }
    Tag gid = mb->globalId_tag();
    tagHandles.push_back( gid );
    rval = mb->write_file( outfile.c_str(), 0, 0, 0, 0, &tagHandles[0],
                                                  (int)tagHandles.size() );MB_CHK_ERR( rval );

    delete mb;
    return 0;
}
