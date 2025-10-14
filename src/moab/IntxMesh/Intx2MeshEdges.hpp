/*
 * Intx2MeshEdges.hpp
 *
 *  Created on: Aug. 4, 2025
 *      Author: iulian
 */

#ifndef INTX2MESHEDGES_HPP_
#define INTX2MESHEDGES_HPP_

#include "Intx2MeshOnSphere.hpp"
#include <iomanip>

#ifdef MOAB_HAVE_NETCDF
#include <netcdf.h>
#define ERRCODE 2
#define ERRN( e )                                  \
    {                                              \
        printf( "Error: %s\n", nc_strerror( e ) ); \
        exit( ERRCODE );                           \
    }
#endif

#ifdef MOAB_HAVE_PNETCDF
#include <pnetcdf.h>
#define ERRCODE 2
#define ERR( e )                                      \
    {                                                 \
        printf( "Error: %s\n", ncmpi_strerror( e ) ); \
        exit( ERRCODE );                              \
    }
#endif

namespace moab
{

class Intx2MeshEdges : public moab::Intx2MeshOnSphere
{
  public:
    Intx2MeshEdges( Interface* mbimpl, IntxAreaUtils::AreaMethod amethod = IntxAreaUtils::lHuiller );
    virtual ~Intx2MeshEdges();

    ErrorCode EdgeSplits( double areaTolerance );

#ifdef MOAB_HAVE_PNETCDF
#ifdef MOAB_HAVE_MPI
    ErrorCode write_edge_map_parallel( const char* fileName );
#endif
#endif

#ifdef MOAB_HAVE_NETCDF
    ErrorCode write_edge_map( const char* fileName );
#endif

  private:
    // local helper
    ErrorCode orderSubEdges( std::vector< EntityHandle >& subEdges,
                             std::vector< EntityHandle >& VerticesSubEdges,
                             const EntityHandle* connEdge,
                             std::vector< EntityHandle >& chainVertices,
                             std::vector< int >& polygonIds,
                             Tag otherParentTag );

    // internal map / arrays
    std::map< EntityHandle, std::vector< EntityHandle > >
        edgeVertices;  // for each recovered edge, the chain of vertices that form subedges
    std::map< EntityHandle, std::vector< int > >
        edgePolygons;  // for each recovered edge, the list of intersected polygons;
    Range recoveredCells;
};

#endif /* INTX2MESHEDGES_HPP_ */
}
