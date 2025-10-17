/*
 * Intx2MeshEdges.cpp
 *
 *  Created on: Aug. 4, 2025
 *      Author: iulian
 */

#include "moab/IntxMesh/Intx2MeshEdges.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

namespace moab
{

Intx2MeshEdges::Intx2MeshEdges( Interface* mbimpl, IntxAreaUtils::AreaMethod amethod )
    : Intx2MeshOnSphere( mbimpl, amethod )
{
}

Intx2MeshEdges::~Intx2MeshEdges()
{
    // TODO Auto-generated destructor stub
}

// local method
ErrorCode Intx2MeshEdges::orderSubEdges( std::vector< EntityHandle >& subEdges,
                                         std::vector< EntityHandle >& VerticesSubEdges,
                                         const EntityHandle* connEdge,
                                         std::vector< EntityHandle >& chainVertices,
                                         std::vector< int >& polygonIds,
                                         Tag otherParentTag )
{
    int numEdges = (int)subEdges.size();
    if( numEdges == 1 ) return MB_SUCCESS;  // nothing to do

    EntityHandle currentVertex = connEdge[0];  // start vertex
    chainVertices.push_back( currentVertex );
    EntityHandle endVertex = connEdge[1];
    std::vector< EntityHandle > chain;
    std::vector< int > markedEdge( numEdges, 0 );  // 0 means not found yet; -1 or +1 for orientation
    // start finding the current vertex, until we close the chain; double loop, we could be smarter :)
    for( int i = 0; i < numEdges; i++ )
    {
        for( int j = 0; j < numEdges; j++ )
        {
            if( 0 != markedEdge[j] ) continue;  // do not use it anymore
            if( VerticesSubEdges[j * 2] == currentVertex )
            {
                currentVertex = VerticesSubEdges[j * 2 + 1];
                chainVertices.push_back( currentVertex );
                chain.push_back( subEdges[j] );
                markedEdge[j] = 1;  // positive
                break;              // break the j loop
            }
            if( VerticesSubEdges[j * 2 + 1] == currentVertex )
            {
                currentVertex = VerticesSubEdges[j * 2];
                chainVertices.push_back( currentVertex );
                chain.push_back( subEdges[j] );
                markedEdge[j] = -1;  // reversed
                break;               // break the j loop
            }
        }
    }
    if( (int)chain.size() == numEdges && currentVertex == endVertex )
    {
        subEdges = chain;  // reordered list, no orientation saved; maybe we should ?
        // from chain, form the list of original polygons that contain each subedge
        // get the parent tag of 2 intx polys connected to each edge in chain
        for( int j = 0; j < (int)chain.size(); j++ )
        {
            EntityHandle sEdge = chain[j];
            std::vector< EntityHandle > intxPolys;
            ErrorCode rval = mb->get_adjacencies( &sEdge, 1, 2, false, intxPolys, Interface::UNION );MB_CHK_ERR( rval );
            EntityHandle onePolygon = intxPolys[0];
            int global_id           = 0;
            rval                    = mb->tag_get_data( otherParentTag, &onePolygon, 1, &global_id );MB_CHK_ERR( rval );
            polygonIds.push_back( global_id );
        }
        return MB_SUCCESS;
    }
    else
        return MB_FAILURE;  // we did not find a chain, do not change anything
}

ErrorCode Intx2MeshEdges::EdgeSplits( double areaTolerance )
{

    Tag parentTag, otherParentTag;
    ErrorCode rval;

    rval = mb->tag_get_handle( "TargetParent", parentTag );MB_CHK_SET_ERR( rval, "can't get parent target tag in edge map" );
    rval = mb->tag_get_handle( "SourceParent", otherParentTag );MB_CHK_SET_ERR( rval, "can't get parent source tag in edge map" );

    Tag fractionTag;
    rval = mb->tag_get_handle( "EdgeRecoveryFraction", fractionTag );MB_CHK_SET_ERR( rval, "can't get tag for recovery fraction" );
    Tag subTag;
    rval = mb->tag_get_handle( "NumSubEnts", subTag );MB_CHK_SET_ERR( rval, "can't get tag NumSubEdges" );
    Tag areaDiffTag;
    rval = mb->tag_get_handle( "AreaDiff", areaDiffTag );MB_CHK_SET_ERR( rval, "can't get tag AreaDiff" );
    Tag areaTag;
    rval = mb->tag_get_handle( "Area", areaTag );MB_CHK_SET_ERR( rval, "can't get tag Area" );

    Range parentCells;
    rval = mb->get_entities_by_dimension( mbs2, 2, parentCells );MB_CHK_SET_ERR( rval, "can't get parent cells" );
    Range initialEdges;
    // we will assume the edges are already included in the target mesh
    // we do not want to create them
    rval = mb->get_adjacencies( parentCells, 1, false, initialEdges, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get parent edges" );
    std::cout << " number of original input edges:" << initialEdges.size() << "\n";

    // get all polygons in the intx set
    Range cells;
    rval = mb->get_entities_by_dimension( outSet, 2, cells );MB_CHK_SET_ERR( rval, "can't get intersection cells" );
    // create all edges adjacent to the cells in the intx set
    // some might be original edges from edge or target meshes
    Range intxEdges;
    rval = mb->get_adjacencies( cells, 1, true, intxEdges, Interface::UNION );MB_CHK_SET_ERR( rval, "can't create adjacent intx edges " );
    std::cout << " number of intx edges:" << intxEdges.size() << "\n";

    // first, identify input polygons that are recovered fully
    // get global ids of the initial cells
    std::vector< int > parentGids( parentCells.size() );
    Tag gidTag = mb->globalId_tag();
    rval       = mb->tag_get_data( gidTag, parentCells, &parentGids[0] );MB_CHK_SET_ERR( rval, "can't get parent global ids" );
    std::map< int, double > initAreas;  // get areas of those initial cells
    int i = 0;
    std::vector< double > coords( 3 * 20 );  // enough vertices, at most 30, good for intx polygons too
    int num_nodes = 0;                       // num nodes in cells
    const EntityHandle* verts;
    IntxAreaUtils areaAdaptor( IntxAreaUtils::lHuiller );  // GaussQuadrature , lHuiller
    std::map< int, double > recoveredAreas;                // get areas of from intersection cells
    std::map< int, EntityHandle > mapFromGIDToParent;
    std::map< int, std::vector< EntityHandle > > mapFromParentGIDToIntxCells;
    for( Range::iterator it = parentCells.begin(); it != parentCells.end(); ++it, i++ )
    {
        EntityHandle parentCell = *it;
        rval                    = mb->get_connectivity( parentCell, verts, num_nodes );MB_CHK_SET_ERR( rval, "can't get connectivity of parent cell" );
        // get coordinates
        rval = mb->get_coords( verts, num_nodes, &coords[0] );MB_CHK_SET_ERR( rval, "can't get coords of parent cell" );

        double area = areaAdaptor.area_spherical_polygon( &coords[0], num_nodes, 1. );
        rval        = mb->tag_set_data( areaTag, &parentCell, 1, &area );MB_CHK_SET_ERR( rval, "can't set area Tag on parent cell" );
        int parentID                 = parentGids[i];
        initAreas[parentID]          = area;
        recoveredAreas[parentID]     = 0.;
        mapFromGIDToParent[parentID] = parentCell;
        mapFromParentGIDToIntxCells[parentID];  // just initialize it with empty vector
    }

    for( Range::iterator it = cells.begin(); it != cells.end(); ++it )
    {
        EntityHandle cell = *it;
        rval              = mb->get_connectivity( cell, verts, num_nodes );MB_CHK_SET_ERR( rval, "can't get connectivity of intx cell" );
        // get coordinates
        rval = mb->get_coords( verts, num_nodes, &coords[0] );MB_CHK_SET_ERR( rval, "can't get coods of intx cell" );
        double intx_area = areaAdaptor.area_spherical_polygon( &coords[0], num_nodes, 1. );
        int parentID     = 0;
        rval             = mb->tag_get_data( parentTag, &cell, 1, &parentID );MB_CHK_SET_ERR( rval, "can't get parent Tag" );
        recoveredAreas[parentID] += intx_area;
        mapFromParentGIDToIntxCells[parentID].push_back( cell );
    }
    int recovered = 0, notRecovered = 0;
    //Range recoveredCells;
    for( size_t j = 0; j < parentGids.size(); j++ )
    {
        int parentID            = parentGids[j];
        EntityHandle parentCell = parentCells[j];
        double areaDiff         = initAreas[parentID] - recoveredAreas[parentID];
        rval                    = mb->tag_set_data( areaDiffTag, &parentCell, 1, &areaDiff );MB_CHK_SET_ERR( rval, "can't set diff area Tag" );
        double numIntxcells = static_cast< double >( mapFromParentGIDToIntxCells[parentID].size() );
        rval                = mb->tag_set_data( subTag, &parentCell, 1, &numIntxcells );MB_CHK_SET_ERR( rval, "can't set num sub ents Tag" );
        if( fabs( areaDiff ) < areaTolerance )
        {
            recovered++;
            recoveredCells.insert( parentCells[j] );  // should we use a std::vector, that will be ordered already ?
        }
        else
            notRecovered++;
    }
    std::cout << "recovered initial cells: " << recovered << " vs:" << notRecovered << " not recovered \n";
    // initial edges that should be decomposable from intx edges
    Range recoverableEdges;
    rval = mb->get_adjacencies( recoveredCells, 1, false, recoverableEdges, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get recoverable edges" );

    // now recover each edge, looking at intx polygons forming one of the adjacent cells, that is in initial recoverable set
    int recoveredEdges = 0;
    int unrecovered    = 0;
    int identity_edges = 0;  // edges that are formed by one intx edge, itself, the original
    std::map< EntityHandle, std::vector< EntityHandle > > mapEdges;
    for( Range::iterator eit = recoverableEdges.begin(); eit != recoverableEdges.end(); ++eit )
    {
        EntityHandle initialEdge = *eit;
        // if this edge is among intxEdges, we are done
        int index = intxEdges.index( initialEdge );
        if( index >= 0 )
        {
            mapEdges[initialEdge].push_back( initialEdge );
            identity_edges++;
            recoveredEdges++;
            double numSubEdges = 1.;
            rval               = mb->tag_set_data( subTag, &initialEdge, 1, &numSubEdges );MB_CHK_SET_ERR( rval, "can't set num sub ents Tag" );
            // get vertices and intx poly attached to it
            int nve = 0;
            const EntityHandle* connCell;
            rval = mb->get_connectivity( initialEdge, connCell, nve );MB_CHK_SET_ERR( rval, "can't get connectivity of parent cell" );
            edgeVertices[initialEdge].push_back( connCell[0] );
            edgeVertices[initialEdge].push_back( connCell[1] );
            // find cells in intx set adjacent to it, and get the other tag parent
            Range adjPolys;
            rval = mb->get_adjacencies( &initialEdge, 1, 2, false, adjPolys, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get adj polys" );
            adjPolys = subtract( adjPolys, parentCells );
            if( adjPolys.empty() ) MB_CHK_SET_ERR( MB_FAILURE, "no adjacent intx cells" );
            EntityHandle intxPoly = adjPolys[0];
            // get its parent tag
            int gid;
            rval = mb->tag_get_data( otherParentTag, &intxPoly, 1, &gid );MB_CHK_SET_ERR( rval, "can't get global id from other tag" );
            edgePolygons[initialEdge].push_back( gid );
            continue;  // no need to sweat it anymore
        }
        // get adjacent polygons; if there is an adj polygon in intx cells, we are done; if not, get one in the initial recoveredCells range
        Range initialAdjCells;
        rval = mb->get_adjacencies( &initialEdge, 1, 2, false, initialAdjCells, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get adjacent initial cells" );
        if( initialAdjCells.size() == 0 ) MB_CHK_SET_ERR( MB_FAILURE, "no adjacent initial cells" );
        // intersect with recoveredCells
        Range adjRecoveredInitialCell = intersect( initialAdjCells, recoveredCells );
        if( adjRecoveredInitialCell.size() == 0 )
            MB_CHK_SET_ERR( MB_FAILURE, "no adjacent initial cells that are recovered" );
        // get all edges adjacent to intersection polygons that form one of the recovered initial cell
        EntityHandle recoveredCell = adjRecoveredInitialCell[0];  // first one
        int nv                     = 0;
        const EntityHandle* connCell;
        rval = mb->get_connectivity( recoveredCell, connCell, nv );MB_CHK_SET_ERR( rval, "can't get connectivity of parent cell" );
        int cellGlobalId = 0;
        rval             = mb->tag_get_data( gidTag, &recoveredCell, 1, &cellGlobalId );MB_CHK_SET_ERR( rval, "can't get id of parent cell" );
        // list of intx polygons in it:
        std::vector< EntityHandle >& listIntxCells = mapFromParentGIDToIntxCells[cellGlobalId];
        std::vector< EntityHandle > candidateEdges;
        rval =
            mb->get_adjacencies( &listIntxCells[0], listIntxCells.size(), 1, false, candidateEdges, Interface::UNION );MB_CHK_SET_ERR( rval, "can't get adj edges" );
        // add all edges that have both vertices on the original edge
        CartVect verticesOriginal[2];
        const EntityHandle* connEdge;
        int numVerts = 0;
        rval         = mb->get_connectivity( initialEdge, connEdge, numVerts );MB_CHK_SET_ERR( rval, "can't get connectivity of initial edge" );

        rval = mb->get_coords( connEdge, numVerts, verticesOriginal[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of initial edge" );
        // get gnomonic plane for the start of the edge
        int gnomonicPlane = 0;
        IntxUtils::decide_gnomonic_plane( verticesOriginal[0], gnomonicPlane );
        double coords2D[6];  // coords in gnomonic plane, for interior point decision
        rval = IntxUtils::gnomonic_projection( verticesOriginal[0], 1.0, gnomonicPlane, coords2D[0], coords2D[1] );MB_CHK_SET_ERR( rval, "can't get gnomonic coords" );
        rval = IntxUtils::gnomonic_projection( verticesOriginal[1], 1.0, gnomonicPlane, coords2D[2], coords2D[3] );MB_CHK_SET_ERR( rval, "can't get gnomonic coords" );

        double edgeLength =
            angle_robust( verticesOriginal[0], verticesOriginal[1] );  // distance on sphere of radius 1 is the angle
        // loop now over candidateEdges
        double recoveredLength = 0.;
        std::vector< EntityHandle > VerticesSubEdges;  // push all vertices here, so we can order the subedges later
        for( size_t i = 0; i < candidateEdges.size(); i++ )
        {
            EntityHandle subedge = candidateEdges[i];
            // get its vertex coordinates:
            const EntityHandle* connEdge2;
            rval = mb->get_connectivity( subedge, connEdge2, numVerts );MB_CHK_SET_ERR( rval, "can't get connectivity of candidate edge" );
            CartVect verticesSubEdge[2];
            rval = mb->get_coords( connEdge2, numVerts, verticesSubEdge[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of initial edge" );
            // decide if they are on the initial edge (form an angle)
            bool onEdge = true;
            for( int j = 0; j < 2 && onEdge; j++ )
            {
                rval =
                    IntxUtils::gnomonic_projection( verticesSubEdge[j], 1.0, gnomonicPlane, coords2D[4], coords2D[5] );MB_CHK_SET_ERR( rval, "can't get gnomonic coords of subedge" );
                // area of triangle in gnomonic plane should be 0
                double areaAbs = fabs( IntxUtils::area2D( &coords2D[0], &coords2D[2], &coords2D[4] ) );

                if( areaAbs > 1.e-14 ) onEdge = false;
            }
            if( onEdge )
            {
                double subEdgeLen = angle_robust( verticesSubEdge[0], verticesSubEdge[1] );
                recoveredLength += subEdgeLen;
                mapEdges[initialEdge].push_back( subedge );
                VerticesSubEdges.push_back( connEdge2[0] );
                VerticesSubEdges.push_back( connEdge2[1] );
            }
        }
        double fraction = recoveredLength / edgeLength;
        rval            = mb->tag_set_data( fractionTag, &initialEdge, 1, &fraction );MB_CHK_SET_ERR( rval, "can't set fraction on initial edge" );
        double numSubEdge = double( mapEdges[initialEdge].size() );
        rval              = mb->tag_set_data( subTag, &initialEdge, 1, &numSubEdge );MB_CHK_SET_ERR( rval, "can't set number of subedges on initial edge" );
        // now , reorder the subedges to create a chain along the original edge
        // if we cannot form a chain, it means we have a problem
        // order subedges on the original edge, and find out their orientation
        // set also the parent tag on the edge, either source or target parent
        // in some cases, the parents can be both
        rval = orderSubEdges( mapEdges[initialEdge], VerticesSubEdges, connEdge, edgeVertices[initialEdge],
                              edgePolygons[initialEdge], otherParentTag );
        if( fabs( edgeLength - recoveredLength ) < 1.e-10 || rval == MB_SUCCESS )
            recoveredEdges++;
        else
        {
            mb->list_entity( initialEdge );
            std::vector< EntityHandle >& listSubEdges = mapEdges[initialEdge];
            double newCheckLength                     = 0;
            CartVect vertices[2];
            const EntityHandle* conn2;
            int numVerts2 = 0;
            rval          = mb->get_connectivity( initialEdge, conn2, numVerts2 );MB_CHK_SET_ERR( rval, "can't get connectivity of initial edge " );
            rval = mb->get_coords( conn2, numVerts2, vertices[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of initial edge" );
            double length2 = angle_robust( vertices[0], vertices[1] );
            std::cout << std::setprecision( 14 );
            std::cout << " initial edge:" << mb->id_from_handle( initialEdge ) << " v: " << conn2[0] << ", " << conn2[1]
                      << " len: " << length2 << "\n";
            for( size_t j = 0; j < listSubEdges.size(); j++ )
            {
                EntityHandle subEdge = listSubEdges[j];
                rval                 = mb->get_connectivity( subEdge, conn2, numVerts2 );MB_CHK_SET_ERR( rval, "can't get connectivity of subedge " );
                rval = mb->get_coords( conn2, numVerts2, vertices[0].array() );MB_CHK_SET_ERR( rval, "can't get coordinates of vertices of subedge" );
                length2 = angle_robust( vertices[0], vertices[1] );
                newCheckLength += length2;
                std::cout << "     sub edge:" << mb->id_from_handle( subEdge ) << " v: " << conn2[0] << ", " << conn2[1]
                          << " len: " << length2 << "\n";
            }
            std::cout << " edge length:" << edgeLength << " newCheckLength:" << newCheckLength
                      << " diff:" << edgeLength - newCheckLength << " fraction:" << fraction
                      << " subedges:" << numSubEdge << "\n";
            unrecovered++;
            std::cout << std::setprecision( 7 );
        }
    }

    std::cout << " recoveredEdges:" << recoveredEdges << " identity edges:" << identity_edges
              << " unrecovered edges:" << unrecovered << "\n";
    return MB_SUCCESS;
}

#ifdef MOAB_HAVE_PNETCDF
#ifdef MOAB_HAVE_MPI
ErrorCode Intx2MeshEdges::write_edge_map_parallel( const char* fileName )
{
    // open for writing the pnetcdf nc file
    int ncid;  // file id
    int ncell_dimid, max_edge_dimid, max_sub_edge_dimid, max_sub_edgeP1_dimid;
    int retval;  // return val for nc
    //Below is an example code fragment that sets the file header hint to 1MB and pass it to PnetCDF when creating a file.

    if( ( retval = ncmpi_create( parcomm->comm(), fileName, NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid ) ) )
        ERR( retval );

    int num_cells_local;
    Range polys;
    ErrorCode rval = mb->get_entities_by_dimension( mbs2, 2, polys );MB_CHK_SET_ERR( rval, "Failed to get polygons" );
    num_cells_local = (int)polys.size();

    int num_cells = 0;
    // do  mpi reduce all, to find out the total num_cells;
    MPI_Allreduce( &num_cells_local, &num_cells, 1, MPI_INTEGER, MPI_SUM, parcomm->comm() );

    if( ( retval = ncmpi_def_dim( ncid, "num_cells", num_cells, &ncell_dimid ) ) ) ERR( retval );

    Tag gid = mb->globalId_tag();
    std::vector< int > global_ids_polys( num_cells_local );
    rval = mb->tag_get_data( gid, polys, global_ids_polys.data() );MB_CHK_SET_ERR( rval, "Failed to get ids of cells" );
    Range localGidCells;
    std::copy( global_ids_polys.rbegin(), global_ids_polys.rend(), range_inserter( localGidCells ) );

    // find max_edges and max subedges
    int max_edge = -1, max_sub_edge = -1, max_subedge1 = -1;

    for( auto it = polys.begin(); it != polys.end(); ++it )
    {
        EntityHandle polygon     = *it;
        const EntityHandle* conn = nullptr;
        int nv;
        rval = mb->get_connectivity( polygon, conn, nv );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );
        if( max_edge < nv ) max_edge = nv;
    }
    int max_edgeg = 0;
    // do  mpi reduce all, to find out the max number of edges;
    MPI_Allreduce( &max_edge, &max_edgeg, 1, MPI_INTEGER, MPI_MAX, parcomm->comm() );
    if( ( retval = ncmpi_def_dim( ncid, "max_edges", max_edgeg, &max_edge_dimid ) ) ) ERR( retval );
    max_edge = max_edgeg;
    for( auto mapit = edgePolygons.begin(); mapit != edgePolygons.end(); ++mapit )
    {
        int nsb = (int)mapit->second.size();
        if( max_sub_edge < nsb ) max_sub_edge = nsb;
    }
    int max_sub_edgeg = 0;
    // do  mpi reduce all, to find out the max number of subedges;
    MPI_Allreduce( &max_sub_edge, &max_sub_edgeg, 1, MPI_INTEGER, MPI_MAX, parcomm->comm() );
    if( ( retval = ncmpi_def_dim( ncid, "max_sub_edges", max_sub_edgeg, &max_sub_edge_dimid ) ) ) ERR( retval );
    max_sub_edge = max_sub_edgeg;
    max_subedge1 = max_sub_edge + 1;
    if( ( retval = ncmpi_def_dim( ncid, "max_sub_edges1", max_subedge1, &max_sub_edgeP1_dimid ) ) ) ERR( retval );

    int dimids_nbs[2];

    dimids_nbs[0] = ncell_dimid;
    dimids_nbs[1] = max_edge_dimid;
    int varid_nsub;
    if( ( retval = ncmpi_def_var( ncid, "nb_sub_edge", NC_INT, 2, dimids_nbs, &varid_nsub ) ) ) ERR( retval );

    int dimids_cell_assoc[3];
    dimids_cell_assoc[0] = ncell_dimid;
    dimids_cell_assoc[1] = max_edge_dimid;
    dimids_cell_assoc[2] = max_sub_edge_dimid;
    int varid_cell_assoc;
    if( ( retval = ncmpi_def_var( ncid, "cells_assoc", NC_INT, 3, dimids_cell_assoc, &varid_cell_assoc ) ) )
        ERR( retval );

    dimids_cell_assoc[2] = max_sub_edgeP1_dimid;
    int varid_lat, varid_lon;
    if( ( retval = ncmpi_def_var( ncid, "lat_sub_edge", NC_DOUBLE, 3, dimids_cell_assoc, &varid_lat ) ) ) ERR( retval );
    if( ( retval = ncmpi_def_var( ncid, "lon_sub_edge", NC_DOUBLE, 3, dimids_cell_assoc, &varid_lon ) ) ) ERR( retval );

    ncmpi_enddef( ncid );

    // need to find out the owned local cells on this task
    // write only those owned local cells, at the correct location, given by their global id
    // we are repeating the shared edges anyway;
    /// so basically, all local cells are owned; no need to worry, just write them at the correct location in the file, based on
    // their global id
    std::vector< int > nb_sub_edge_per_edge( num_cells_local * max_edge, -9999 );
    std::vector< int > cells_assoc_per_edge( num_cells_local * max_edge * max_sub_edge, -9999 );

    std::vector< double > latvals( num_cells_local * max_edge * max_subedge1, -9999 );
    std::vector< double > lonvals( num_cells_local * max_edge * max_subedge1, -9999 );

    for( auto it = recoveredCells.begin(); it != recoveredCells.end(); ++it )
    {
        EntityHandle polygon = *it;
        int gidPoly          = 0;
        rval                 = mb->tag_get_data( gid, &polygon, 1, &gidPoly );MB_CHK_SET_ERR( rval, "Failed to get id of poly" );
        int indexGidPoly = localGidCells.index( gidPoly );

        const EntityHandle* conn = nullptr;
        int nv;
        rval = mb->get_connectivity( polygon, conn, nv );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );
        for( int i = 0; i < nv; i++ )
        {
            EntityHandle v[2];
            v[0] = conn[i];
            v[1] = conn[( i + 1 ) % nv];
            Range edges;
            rval = mb->get_adjacencies( v, 2, 1, false, edges, Interface::INTERSECT );MB_CHK_SET_ERR( rval, "Failed to get edge" );
            EntityHandle edge = edges[0];
            // reverse or not?
            std::vector< int > poly_assoc                     = edgePolygons[edge];
            nb_sub_edge_per_edge[indexGidPoly * max_edge + i] = (int)poly_assoc.size();
            std::vector< EntityHandle > vertexEdges           = edgeVertices[edge];
            std::vector< CartVect > coords( vertexEdges.size() );
            rval = mb->get_coords( &vertexEdges[0], vertexEdges.size(), &( coords[0][0] ) );MB_CHK_SET_ERR( rval, "can't get coordinates" );
            // convert to lat/lon
            std::vector< double > latv( vertexEdges.size() ), lonv( vertexEdges.size() );
            for( int j = 0; j < (int)vertexEdges.size(); j++ )
            {
                IntxUtils::SphereCoords sph1 = IntxUtils::cart_to_spherical( coords[j] );
                lonv[j]                      = sph1.lon;
                latv[j]                      = sph1.lat;
            }
            // reversed edge or not?
            const EntityHandle* edgeconn = nullptr;
            int nve;
            rval = mb->get_connectivity( edge, edgeconn, nve );MB_CHK_SET_ERR( rval, "Failed to get edge connectivity" );
            bool reverse = false;
            if( v[0] == edgeconn[1] ) reverse = true;
            if( reverse )
            {
                // fill the arrays in reverse order
                int sizep = (int)poly_assoc.size();
                for( int j = 0; j < sizep; j++ )
                {
                    cells_assoc_per_edge[indexGidPoly * max_edge * max_sub_edge + i * max_sub_edge + j] =
                        poly_assoc[sizep - 1 - j];
                }
                sizep = (int)vertexEdges.size();
                for( int j = 0; j < sizep; j++ )
                {
                    latvals[indexGidPoly * max_edge * max_subedge1 + i * max_subedge1 + j] = latv[sizep - 1 - j];
                    lonvals[indexGidPoly * max_edge * max_subedge1 + i * max_subedge1 + j] = lonv[sizep - 1 - j];
                }
            }
            else
            {
                // max_sub_edges
                for( int j = 0; j < (int)poly_assoc.size(); j++ )
                {
                    cells_assoc_per_edge[indexGidPoly * max_edge * max_sub_edge + i * max_sub_edge + j] = poly_assoc[j];
                }
                for( int j = 0; j < (int)vertexEdges.size(); j++ )
                {
                    latvals[indexGidPoly * max_edge * max_subedge1 + i * max_subedge1 + j] = latv[j];
                    lonvals[indexGidPoly * max_edge * max_subedge1 + i * max_subedge1 + j] = lonv[j];
                }
            }
        }
    }

    size_t idxReq    = 0;
    size_t nb_writes = localGidCells.psize();
    std::vector< int > requests( nb_writes * 4 );
    std::vector< int > statuss( nb_writes * 4 );

    size_t indexInArray1 = 0;
    size_t indexInArray2 = 0;
    size_t indexInArray3 = 0;

    for( Range::pair_iterator pair_iter = localGidCells.pair_begin(); pair_iter != localGidCells.pair_end();
         ++pair_iter )  // these will be the size of nb_writes
    {
        EntityHandle starth = pair_iter->first;
        EntityHandle endh   = pair_iter->second;
        MPI_Offset write_start[3], write_count[3];  // max dimension
        write_start[0] = static_cast< MPI_Offset >( starth - 1 );
        write_start[1] = static_cast< MPI_Offset >( 0 );
        write_start[2] = static_cast< MPI_Offset >( 0 );
        write_count[0] = static_cast< MPI_Offset >( endh - starth + 1 );
        write_count[1] = static_cast< MPI_Offset >( max_edge );
        write_count[2] = static_cast< MPI_Offset >( max_sub_edge );

        if( ( retval = ncmpi_iput_vara_int( ncid, varid_nsub, write_start, write_count,  // first 2 are used
                                            &nb_sub_edge_per_edge[indexInArray1], &requests[idxReq++] ) ) )
            ERR( retval );
        indexInArray1 += ( endh - starth + 1 ) * max_edge;

        if( ( retval = ncmpi_iput_vara_int( ncid, varid_cell_assoc, write_start, write_count,  // all 3 are used
                                            &cells_assoc_per_edge[indexInArray2], &requests[idxReq++] ) ) )
            ERR( retval );
        indexInArray2 += ( endh - starth + 1 ) * max_edge * max_sub_edge;

        write_count[2] = max_subedge1;

        if( ( retval = ncmpi_iput_vara_double( ncid, varid_lat, write_start, write_count,  // all 3 are used
                                               &latvals[indexInArray3], &requests[idxReq++] ) ) )
            ERR( retval );

        if( ( retval = ncmpi_iput_vara_double( ncid, varid_lon, write_start, write_count,  // all 3 are used
                                               &lonvals[indexInArray3], &requests[idxReq++] ) ) )
            ERR( retval );

        indexInArray3 += ( endh - starth + 1 ) * max_edge * max_subedge1;
    }
    // Wait outside the loop
    if( ( retval = ncmpi_wait_all( ncid, requests.size(), &requests[0], &statuss[0] ) ) ) ERR( retval );
    if( ( retval = ncmpi_close( ncid ) ) ) ERR( retval );
    return MB_SUCCESS;
}
#endif
#endif

#ifdef MOAB_HAVE_NETCDF
ErrorCode Intx2MeshEdges::write_edge_map( const char* filename )
/*Interface * mb, EntityHandle sf1,
        std::map<EntityHandle, std::vector<EntityHandle>>  & edgeVertices,
        std::map<EntityHandle, std::vector<int>> & edgePolygons,
        moab::Range & recoveredCells)*/
{
    // open for writing the nc file
    int ncid;  // file id
    int ncell_dimid, max_edge_dimid, max_sub_edge_dimid, max_sub_edgeP1_dimid;
    int retval;  // return val for nc
    if( ( retval = nc_create( filename, NC_CLASSIC_MODEL | NC_CLOBBER, &ncid ) ) ) ERRN( retval );
    int num_cells;
    Range polys;
    ErrorCode rval = mb->get_entities_by_dimension( mbs2, 2, polys );MB_CHK_SET_ERR( rval, "Failed to get polygons" );
    num_cells = (int)polys.size();

    if( ( retval = nc_def_dim( ncid, "num_cells", num_cells, &ncell_dimid ) ) ) ERRN( retval );

    Tag gid = mb->globalId_tag();
    // find max_edges and max subedges
    int max_edge = -1, max_sub_edge = -1, max_subedge1 = -1;

    for( auto it = polys.begin(); it != polys.end(); ++it )
    {
        EntityHandle polygon     = *it;
        const EntityHandle* conn = nullptr;
        int nv;
        rval = mb->get_connectivity( polygon, conn, nv );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );
        if( max_edge < nv ) max_edge = nv;
    }
    if( ( retval = nc_def_dim( ncid, "max_edges", max_edge, &max_edge_dimid ) ) ) ERRN( retval );

    for( auto mapit = edgePolygons.begin(); mapit != edgePolygons.end(); ++mapit )
    {
        int nsb = (int)mapit->second.size();
        if( max_sub_edge < nsb ) max_sub_edge = nsb;
    }
    if( ( retval = nc_def_dim( ncid, "max_sub_edges", max_sub_edge, &max_sub_edge_dimid ) ) ) ERRN( retval );
    max_subedge1 = max_sub_edge + 1;
    if( ( retval = nc_def_dim( ncid, "max_sub_edges1", max_subedge1, &max_sub_edgeP1_dimid ) ) ) ERRN( retval );

    int dimids_nbs[2];

    dimids_nbs[0] = ncell_dimid;
    dimids_nbs[1] = max_edge_dimid;
    int varid_nsub;
    if( ( retval = nc_def_var( ncid, "nb_sub_edge", NC_INT, 2, dimids_nbs, &varid_nsub ) ) ) ERRN( retval );

    int dimids_cell_assoc[3];
    dimids_cell_assoc[0] = ncell_dimid;
    dimids_cell_assoc[1] = max_edge_dimid;
    dimids_cell_assoc[2] = max_sub_edge_dimid;
    int varid_cell_assoc;
    if( ( retval = nc_def_var( ncid, "cells_assoc", NC_INT, 3, dimids_cell_assoc, &varid_cell_assoc ) ) )
        ERRN( retval );

    dimids_cell_assoc[2] = max_sub_edgeP1_dimid;
    int varid_lat, varid_lon;
    if( ( retval = nc_def_var( ncid, "lat_sub_edge", NC_DOUBLE, 3, dimids_cell_assoc, &varid_lat ) ) ) ERRN( retval );
    if( ( retval = nc_def_var( ncid, "lon_sub_edge", NC_DOUBLE, 3, dimids_cell_assoc, &varid_lon ) ) ) ERRN( retval );

    nc_enddef( ncid );

    std::vector< int > nb_sub_edge_per_edge( num_cells * max_edge, -9999 );
    std::vector< int > cells_assoc_per_edge( num_cells * max_edge * max_sub_edge, -9999 );

    std::vector< double > latvals( num_cells * max_edge * max_subedge1, -9999 );
    std::vector< double > lonvals( num_cells * max_edge * max_subedge1, -9999 );

    for( auto it = recoveredCells.begin(); it != recoveredCells.end(); ++it )
    {
        EntityHandle polygon = *it;
        int gidPoly          = 0;
        rval                 = mb->tag_get_data( gid, &polygon, 1, &gidPoly );MB_CHK_SET_ERR( rval, "Failed to get id of poly" );
        const EntityHandle* conn = nullptr;
        int nv;
        rval = mb->get_connectivity( polygon, conn, nv );MB_CHK_SET_ERR( rval, "Failed to get connectivity" );
        for( int i = 0; i < nv; i++ )
        {
            EntityHandle v[2];
            v[0] = conn[i];
            v[1] = conn[( i + 1 ) % nv];
            Range edges;
            rval = mb->get_adjacencies( v, 2, 1, false, edges, Interface::INTERSECT );MB_CHK_SET_ERR( rval, "Failed to get edge" );
            EntityHandle edge = edges[0];
            // reverse or not?
            std::vector< int > poly_assoc                        = edgePolygons[edge];
            nb_sub_edge_per_edge[( gidPoly - 1 ) * max_edge + i] = (int)poly_assoc.size();
            std::vector< EntityHandle > vertexEdges              = edgeVertices[edge];
            std::vector< CartVect > coords( vertexEdges.size() );
            rval = mb->get_coords( &vertexEdges[0], vertexEdges.size(), &( coords[0][0] ) );MB_CHK_SET_ERR( rval, "can't get coordinates" );
            // convert to lat/lon
            std::vector< double > latv( vertexEdges.size() ), lonv( vertexEdges.size() );
            for( int j = 0; j < (int)vertexEdges.size(); j++ )
            {
                IntxUtils::SphereCoords sph1 = IntxUtils::cart_to_spherical( coords[j] );
                lonv[j]                      = sph1.lon;
                latv[j]                      = sph1.lat;
            }
            // reversed edge or not?
            const EntityHandle* edgeconn = nullptr;
            int nve;
            rval = mb->get_connectivity( edge, edgeconn, nve );MB_CHK_SET_ERR( rval, "Failed to get edge connectivity" );
            bool reverse = false;
            if( v[0] == edgeconn[1] ) reverse = true;
            if( reverse )
            {
                // fill the arrays in reverse order
                int sizep = (int)poly_assoc.size();
                for( int j = 0; j < sizep; j++ )
                {
                    cells_assoc_per_edge[( gidPoly - 1 ) * max_edge * max_sub_edge + i * max_sub_edge + j] =
                        poly_assoc[sizep - 1 - j];
                }
                sizep = (int)vertexEdges.size();
                for( int j = 0; j < sizep; j++ )
                {
                    latvals[( gidPoly - 1 ) * max_edge * max_subedge1 + i * max_subedge1 + j] = latv[sizep - 1 - j];
                    lonvals[( gidPoly - 1 ) * max_edge * max_subedge1 + i * max_subedge1 + j] = lonv[sizep - 1 - j];
                }
            }
            else
            {
                // max_sub_edges
                for( int j = 0; j < (int)poly_assoc.size(); j++ )
                {
                    cells_assoc_per_edge[( gidPoly - 1 ) * max_edge * max_sub_edge + i * max_sub_edge + j] =
                        poly_assoc[j];
                }
                for( int j = 0; j < (int)vertexEdges.size(); j++ )
                {
                    latvals[( gidPoly - 1 ) * max_edge * max_subedge1 + i * max_subedge1 + j] = latv[j];
                    lonvals[( gidPoly - 1 ) * max_edge * max_subedge1 + i * max_subedge1 + j] = lonv[j];
                }
            }
        }
    }

    if( ( retval = nc_put_var_int( ncid, varid_nsub, &nb_sub_edge_per_edge[0] ) ) ) ERRN( retval );

    if( ( retval = nc_put_var_int( ncid, varid_cell_assoc, &cells_assoc_per_edge[0] ) ) ) ERRN( retval );

    if( ( retval = nc_put_var_double( ncid, varid_lat, &latvals[0] ) ) ) ERRN( retval );

    if( ( retval = nc_put_var_double( ncid, varid_lon, &lonvals[0] ) ) ) ERRN( retval );

    if( ( retval = nc_close( ncid ) ) ) ERRN( retval );
    return MB_SUCCESS;
}
#endif

}  // namespace moab
