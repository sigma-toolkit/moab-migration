/** @example ExtractLand.cpp
 * this tool will take an existing h5m  pg2 mesh (fv)  file and a land point cloud mesh
 *
 * will extract the corresponding land mesh
 *
 * example of usage:
 * ./ExtractLand -p wholeATM_PG2.h5m -l wholeLnd.h5m -o LandMesh.h5m
 *
 * the *PG2" style atm mesh is available in E3SM only for pg2 runs, something like
 *  --res ne30pg2_r05_oECv3_ICG --compset A_WCYCL1850S_CMIP6
 *  or
 *  --res ne4pg2_ne4pg2 --compset FC5AV1C-L
 */
#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

#include <cmath>
#include <sstream>

using namespace moab;

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string pg2file, lndfile, outfile;

    opts.addOpt< std::string >( "land,l", "phys grid filename", &lndfile );
    opts.addOpt< std::string >( "pg2file,p", "pg2 mesh file", &pg2file );
    opts.addOpt< std::string >( "output,o", "output mesh filename", &outfile );

    opts.parseCommandLine( argc, argv );

    std::cout << " land file " << lndfile << "\n";
    std::cout << "pg2 mesh file: " << pg2file << "\n";
    std::cout << "output file: " << outfile << "\n";

    if( lndfile.empty() || pg2file.empty() || outfile.empty() )
    {
        opts.printHelp();
        return 0;
    }
    ErrorCode rval;
    Core* mb = new Core();

    MB_CHK_SET_ERR( mb->load_file( lndfile.c_str() ), "can't load land pc file" );

    Core* mb2 = new Core();
    MB_CHK_SET_ERR( mb2->load_file( pg2file.c_str() ), "can't load pg2 mesh file" );

    Tag globalIDTag1 = mb->globalId_tag();

    Tag globalIDTag2 = mb2->globalId_tag();

    Range verts1;
    MB_CHK_SET_ERR( mb->get_entities_by_dimension( 0, 0, verts1 ), "can't get vertices " );

    Range cells;
    MB_CHK_SET_ERR( mb2->get_entities_by_dimension( 0, 2, cells ), "can't get 2d cells " );

    std::vector< int > globalIdsCells;
    globalIdsCells.resize( cells.size() );
    MB_CHK_SET_ERR( mb2->tag_get_data( globalIDTag2, cells, &globalIdsCells[0] ), "can't get global ids cells " );

    std::vector< int > globalIdsVerts;
    globalIdsVerts.resize( verts1.size() );
    MB_CHK_SET_ERR( mb->tag_get_data( globalIDTag1, verts1, &globalIdsVerts[0] ), "can't get global ids cells " );

    // now, every cell will be put into one set, by looking at the global id of cell

    std::map< int, EntityHandle > gidToCell;
    int i = 0;
    for( Range::iterator it = cells.begin(); it != cells.end(); ++it, i++ )
    {
        gidToCell[globalIdsCells[i]] = *it;
    }
    // create a new file set with land cells
    EntityHandle fileSet;
    MB_CHK_SET_ERR( mb2->create_meshset( MESHSET_SET, fileSet ), "Error creating file set" );
    // empty all sets

    Range landCells;
    // look now at gid values for vertices
    for( i = 0; i < (int)verts1.size(); i++ )
    {
        int gid = globalIdsVerts[i];
        landCells.insert( gidToCell[gid] );
    }

    MB_CHK_SET_ERR( mb2->add_entities( fileSet, landCells ), "can't add land cells" );

    MB_CHK_SET_ERR( mb2->write_file( outfile.c_str(), 0, 0, &fileSet, 1 ), "can't write file" );

    // write the original with mask 0/1, default -1
    Tag mask;
    double def_val = -1.;
    MB_CHK_SET_ERR( mb2->tag_get_handle( "mask", 1, MB_TYPE_DOUBLE, mask, MB_TAG_CREAT | MB_TAG_DENSE, &def_val ),
                    "can't create mask tag" );
    for( Range::iterator it = cells.begin(); it != cells.end(); ++it, i++ )
    {
        EntityHandle cell = *it;
        // set to 0
        double val = 0.;
        MB_CHK_SET_ERR( mb2->tag_set_data( mask, &cell, 1, &val ), "can't set mask tag" );
    }

    for( Range::iterator it = landCells.begin(); it != landCells.end(); ++it, i++ )
    {
        EntityHandle cell = *it;
        // set to 0
        double val = 1.;
        MB_CHK_SET_ERR( mb2->tag_set_data( mask, &cell, 1, &val ), "can't set mask tag" );
    }
    MB_CHK_SET_ERR( mb2->delete_entities( &fileSet, 1 ), "can't delete set" );
    MB_CHK_SET_ERR( mb2->write_file( "AtmWithLandMask.h5m" ), "can't write file" );
    delete mb;
    delete mb2;

    return 0;
}
