/*
 * AddEdgeFluxVector.cpp
 * will take an edge tag (baro-flux) and convert it to a 3d vector, associated to edges
 *
 */
#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"

using namespace moab;
using namespace std;

int main( int argc, char* argv[] )
{
    string filein  = "source_1_angle.h5m";
    string fileout = "source_with_vec.h5m";

    ProgOptions opts;
    opts.addOpt< std::string >( "model,m", "input file ", &filein );

    opts.addOpt< std::string >( "output,o", "output filename", &fileout );

    opts.parseCommandLine( argc, argv );

    Core moab;
    Interface* mb = &moab;

    ErrorCode rval = mb->load_file( filein.c_str() );MB_CHK_ERR( rval );

    // get the angle tag, and the baro tag; compute a 3d vector normal on the edge
    // (or look at the angle) compute also the latitude, long at mid edge, draw the normal and multiply with the edge length?
    //
    Tag angle, bflux;
    rval = mb->tag_get_handle( "barotropicThicknessFlux0", bflux );MB_CHK_ERR( rval );
    rval = mb->tag_get_handle( "angleEdge", angle );MB_CHK_ERR( rval );
    MB_CHK_ERR( rval );

    Range edges;
    rval = mb->get_entities_by_dimension( 0, 1, edges );MB_CHK_ERR( rval );
    std::vector< double > angles( edges.size() ), fluxes( edges.size() );
    std::cout << " file:" << filein << " nb edges:" << edges.size() << "\n";

    rval = mb->tag_get_data( angle, edges, &angles[0] );MB_CHK_ERR( rval );
    rval = mb->tag_get_data( bflux, edges, &fluxes[0] );MB_CHK_ERR( rval );

    Tag vectorTag;
    rval = mb->tag_get_handle( "FluxVector", 3, MB_TYPE_DOUBLE, vectorTag, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_SET_ERR( rval, "Couldn't get tag handle" );

    Tag vectorTagU;
    rval = mb->tag_get_handle( "FluxVectorU", 3, MB_TYPE_DOUBLE, vectorTagU, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_SET_ERR( rval, "Couldn't get tag handle" );
    Tag vectorTagV;
    rval = mb->tag_get_handle( "FluxVectorV", 3, MB_TYPE_DOUBLE, vectorTagV, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_SET_ERR( rval, "Couldn't get tag handle" );

    Tag edgeLenTag;
    rval = mb->tag_get_handle( "EdgeLength", 1, MB_TYPE_DOUBLE, edgeLenTag, MB_TAG_CREAT | MB_TAG_DENSE );MB_CHK_SET_ERR( rval, "Couldn't get tag handle" );

    // for each edge. get v1, v2, mid edge. compute lat, lon,
    int i = 0;
    for( Range::iterator eit = edges.begin(); eit != edges.end(); eit++, i++ )
    {
        EntityHandle eh = *eit;
        const EntityHandle* conn;
        int nv = 2;
        rval   = mb->get_connectivity( eh, conn, nv );MB_CHK_SET_ERR( rval, "Couldn't get edge conn" );
        CartVect verts[2];
        if( nv != 2 )
        {
            MB_CHK_SET_ERR( MB_FAILURE, "wrong number of vertices on edge" );
        }
        rval = mb->get_coords( conn, 2, &verts[0][0] );MB_CHK_SET_ERR( rval, "Couldn't get vertex coordinates" );
        CartVect mid   = 0.5 * ( verts[0] + verts[1] );
        double edgeLen = ( verts[1] - verts[0] ).length();
        rval           = mb->tag_set_data( edgeLenTag, &eh, 1, &edgeLen );MB_CHK_SET_ERR( rval, "Couldn't set edge length" );
        // utility to convert to lat/lon
        IntxUtils::SphereCoords sph = IntxUtils::cart_to_spherical( mid );
        double lat = sph.lat, lon = sph.lon;
        // vectors along par latitude, and meridian:
        CartVect u, v;  // these are tangent at the mid to the sphere u = d p / d lon
        // v = dp / d lat
        // CartVect res;
        // res[0] = sc.R * cos( sc.lat ) * cos( sc.lon );  // x coordinate
        // res[1] = sc.R * cos( sc.lat ) * sin( sc.lon );  // y
        // res[2] = sc.R * sin( sc.lat );                  // z
        u[0]          = -sin( lon );  // x coordinate
        u[1]          = cos( lon );
        u[2]          = 0.;  // no z component
        v[0]          = -sin( lat ) * cos( lon );
        v[1]          = -sin( lat ) * sin( lon );
        v[2]          = cos( lat );
        CartVect flux = cos( angles[i] ) * u + sin( angles[i] ) * v;  // unit vector
        flux          = fluxes[i] * flux;                             //

        rval = mb->tag_set_data( vectorTag, &eh, 1, &( flux[0] ) );MB_CHK_SET_ERR( rval, "Couldn't set vector tag" );
        double projU, projV;
        projU       = flux % u;  // this is dot product between vectors
        CartVect pU = projU * u;
        projV       = flux % v;
        CartVect pV = projV * v;
        rval        = mb->tag_set_data( vectorTagU, &eh, 1, &( pU[0] ) );MB_CHK_SET_ERR( rval, "Couldn't set vector tag" );
        rval = mb->tag_set_data( vectorTagV, &eh, 1, &( pV[0] ) );MB_CHK_SET_ERR( rval, "Couldn't set vector tag" );
    }
    // after computing each edge flux, add all vectors in a cell, to see the total flux along edges
    Range cells;
    rval = mb->get_entities_by_dimension( 0, 2, cells );MB_CHK_ERR( rval );

    // get all edges adjacent to a cell, then multiply by edge length and compute the total flux at
    // center of cell
    for( Range::iterator cit = cells.begin(); cit != cells.end(); ++cit )
    {
        EntityHandle cell = *cit;
        // get edges adjacent to the cell
        std::vector< EntityHandle > adjEdges;
        rval = mb->get_adjacencies( &cell, 1, 1, false, adjEdges, Interface::UNION );MB_CHK_ERR( rval );
        CartVect totalFlux = CartVect( 0. );
        for( size_t i = 0; i < adjEdges.size(); i++ )
        {
            EntityHandle edge = adjEdges[i];
            CartVect flux;
            rval = mb->tag_get_data( vectorTag, &edge, 1, &flux[0] );MB_CHK_ERR( rval );
            double edgeLen;
            rval = mb->tag_get_data( edgeLenTag, &edge, 1, &edgeLen );MB_CHK_ERR( rval );
            totalFlux = totalFlux + edgeLen * flux;
        }
        rval = mb->tag_set_data( vectorTag, &cell, 1, &totalFlux[0] );MB_CHK_ERR( rval );
    }
    std::cout << " writing " << fileout << "\n";
    rval = mb->write_file( fileout.c_str() );MB_CHK_ERR( rval );

    return 0;
}
