/**
 * @file AddFieldtoPG2.cpp
 * @brief Example demonstrating addition of field data from phys grid to PG2 mesh
 *
 * This example shows how to:
 * - Load PG2 mesh and phys grid solution files
 * - Extract variable data from phys grid
 * - Match entities between phys grid and PG2 mesh using global IDs
 * - Copy variable data to PG2 mesh cells
 * - Write enhanced PG2 mesh files for visualization
 *
 * This tool is useful for transferring field data from phys grid
 * solutions to PG2 mesh representations for climate model analysis.
 *
 * @author MOAB Development Team
 * @date 2024
 *

 * Description: Add field data from phys grid to PG2 mesh for visualization and analysis
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return 0 on success, 1 on failure
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"
#include <iostream>

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile, outfile( "out.h5m" ), physgridfile, variable_name;

    opts.addOpt< std::string >( "input,i", "input mesh filename", &inputfile );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );
    opts.addOpt< std::string >( "phys,p", "phys grid solution filename", &physgridfile );
    opts.addOpt< std::string >( "var,v", "variable to extract and add to output file", &variable_name );

    opts.parseCommandLine( argc, argv );

    if( inputfile.empty() )
    {
        opts.printHelp();
        return 0;
    }
    ErrorCode rval;
    Core* mb = new Core();

    EntityHandle fset1, fset2;
    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, fset1 ), "can't create mesh set" );
    MB_CHK_SET_ERR( mb->load_file( inputfile.c_str(), &fset1 ), "can't load input file" );

    cout << " opened " << inputfile << " with initial h5m data.\n";

    MB_CHK_SET_ERR( mb->create_meshset( MESHSET_SET, fset2 ), "can't create mesh set" );
    MB_CHK_SET_ERR( mb->load_file( physgridfile.c_str(), &fset2 ), "can't load phys grid file" );

    Tag tagv;
    MB_CHK_SET_ERR( mb->tag_get_handle( variable_name.c_str(), tagv ), "can't get tag handle" );

    Tag gitag = mb->globalId_tag();

    Range verts;  // from phys grid
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( fset2, 0, verts ), "can't get vertices" );
    std::vector< int > gids;
    gids.resize( verts.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( gitag, verts, &gids[0] ), "can't get gi tag values" );
    std::vector< double > valsTag;
    valsTag.resize( verts.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( tagv, verts, &valsTag[0] ), "can't get tag vals" );
    Range cells;

    MB_CHK_SET_ERR( mb->get_entities_by_dimension( fset1, 2, cells ), "can't get cells" );

    std::map< int, double > valsByID;
    for( int i = 0; i < (int)gids.size(); i++ )
        valsByID[gids[i]] = valsTag[i];

    // set now cells values
    std::vector< int > cellsIds;
    cellsIds.resize( cells.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( gitag, cells, &cellsIds[0] ), "can't get cells ids" );
    for( int i = 0; i < (int)cells.size(); i++ )
    {
        valsTag[i] = valsByID[cellsIds[i]];
    }
    MB_CHK_SET_ERR( mb->tag_set_data( tagv, cells, &valsTag[0] ), "can't set  cells tags" );

    MB_CHK_SET_ERR( mb->write_file( outfile.c_str(), 0, 0, &fset1, 1 ), "can't write file" );

    return 0;
}
