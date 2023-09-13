/*
 * NCHelperESMF.cpp
 *
 *  Created on: Sep. 12, 2023
 */

#include "NCHelperESMF.h"
#include "moab/ReadUtilIface.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"

#ifdef MOAB_HAVE_ZOLTAN
#include "moab/ZoltanPartitioner.hpp"
#endif

#include <cmath>

namespace moab
{

const int DEFAULT_MAX_EDGES_PER_CELL = 10;


NCHelperESMF::NCHelperESMF( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
: UcdNCHelper( readNC, fileId, opts, fileSet ), maxEdgesPerCell( DEFAULT_MAX_EDGES_PER_CELL ), numCellGroups( 0 )
{
}

NCHelperESMF::~NCHelperESMF() {
}

bool NCHelperMPAS::can_read_file( ReadNC* readNC )
{
    std::vector< std::string >& dimNames = readNC->dimNames;
    if( ( std::find( dimNames.begin(), dimNames.end(), std::string( "nodeCount" ) ) != dimNames.end() ) &&
            ( std::find( dimNames.begin(), dimNames.end(), std::string( "elementCount" ) ) != dimNames.end() ) &&
            ( std::find( dimNames.begin(), dimNames.end(), std::string( "connectionCount" ) ) != dimNames.end() ) &&
            ( std::find( dimNames.begin(), dimNames.end(), std::string( "coordDim" ) ) != dimNames.end() ) )
    {
        return true;
    }

    return false;
}
