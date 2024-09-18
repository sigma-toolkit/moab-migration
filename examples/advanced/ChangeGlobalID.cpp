/*
 * change global id , reset the old one to 1->size cells
 *
 */
#include <iostream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/ProgOptions.hpp"

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{
    string filein  = "mpas_2d_source.h5m";
    string fileout = "mpas_modified_id.h5m";

    ProgOptions opts;
    opts.addOpt< std::string >( "model,m", "input file ", &filein );

    opts.addOpt< std::string >( "output,o", "output filename", &fileout );

    opts.parseCommandLine( argc, argv );

    Core moab;
    Interface* mb = &moab;
    EntityHandle sf;
    ErrorCode rval = mb->create_meshset( MESHSET_SET, sf );MB_CHK_ERR( rval );

    rval = mb->load_file( filein.c_str(), &sf );MB_CHK_ERR( rval );
    // get all 2d cells, and reset global id
    Range cells;
    rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_ERR( rval );
    Tag global_id_tag = mb->globalId_tag();
    Tag oldIdTag;
    int dum_id = -1;
    rval = mb->tag_get_handle( "OLD_GLOBAL_ID", 1, MB_TYPE_INTEGER, oldIdTag,  MB_TAG_CREAT | MB_TAG_DENSE, &dum_id );MB_CHK_SET_ERR( rval, "Can't get parallel partition tag" );
    std::vector<int> globalIds(cells.size());
    rval = mb->tag_get_data(global_id_tag, cells, &globalIds[0]);MB_CHK_SET_ERR( rval, "Can't get global id vals" );

    rval = mb->tag_set_data(oldIdTag, cells, &globalIds[0]);MB_CHK_SET_ERR( rval, "Can't set old global id vals" );
    // set the new vals in order
	for (size_t i=0; i<globalIds.size(); i++)
		globalIds[i] = i+1;
	rval = mb->tag_set_data(global_id_tag, cells, &globalIds[0]);MB_CHK_SET_ERR( rval, "Can't set new global id vals" );

    std::cout <<"writing output file " << fileout.c_str() << "\n";
    rval = mb->write_file( fileout.c_str(), 0, 0, &sf, 1 );MB_CHK_ERR( rval );

    return 0;
}
