/*
 * =====================================================================================
 *
 *       Filename:  Remapper.hpp
 *
 *    Description:  Interface to the a general remapping capability on arbitrary topology
 *                  that performs both mesh intersection between a source and target grid,
 *                  with arbitrary decompositions. The intersections can then be used to
 *                  either evaluate interpolation weights or to perform high-order
 *                  conservative remapping of solutions defined on the source grid.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#ifndef MB_REMAPPER_HPP
#define MB_REMAPPER_HPP

#include <string>

#include "moab/Interface.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

namespace moab
{

class Remapper
{
  public:
#ifdef MOAB_HAVE_MPI
    Remapper( moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL ) : m_interface( mbInt ), m_pcomm( pcomm )
#else
    Remapper( moab::Interface* mbInt ) : m_interface( mbInt )
#endif
    {
    }

    virtual ~Remapper()
    {
#ifdef MOAB_HAVE_MPI
        m_pcomm = NULL;
#endif
        m_interface = NULL;
    }

    enum IntersectionContext
    {
        DEFAULT              = -1,  // default context
        SourceMesh           = 0,   // source mesh
        TargetMesh           = 1,   // target mesh
        OverlapMesh          = 2,   // overlap/intersection mesh
        CoveringMesh         = 3,   // source mesh covering target mesh
        SourceMeshWithGhosts = 4,   // mesh with extra ghost layers to compute coverage in high order case or bilin
        TargetMeshWithGhosts = 5    // mesh with extra ghost layers to impose target data limiting
    };

    moab::Interface* get_interface()
    {
        return m_interface;
    }
//#define MOAB_DBG
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* get_parallel_communicator()
    {
        return m_pcomm;
    }

    /// <summary>
    ///     ghost layers
    /// </summary>
    moab::ErrorCode GhostLayers( moab::ParallelComm * pcomm, moab::EntityHandle& meshset, const int ngh_layers, moab::EntityHandle& set_with_ghosts )
    {
        // meshset contains the mesh set distributed already
        //
        moab::ErrorCode rval = m_interface->create_meshset( MESHSET_SET, set_with_ghosts );MB_CHK_ERR( rval );
        // copy original content of mesh set here; we will use it later for local area, for example
        // it will not have any ghosts in it
        moab::Range orgEnts;
        rval = m_interface->get_entities_by_handle( meshset, orgEnts );MB_CHK_ERR( rval );
        rval = m_interface->add_entities( set_with_ghosts, orgEnts );MB_CHK_ERR( rval );
        bool global_id_filter = true;// only time this is true, so far
        rval = pcomm->exchange_ghost_cells( 2, 0, 1, 0, true, true, &set_with_ghosts, global_id_filter );MB_CHK_ERR( rval );
        for( int i = 2; i <= ngh_layers; i++ )
        {
            rval = pcomm->correct_thin_ghost_layers();MB_CHK_ERR( rval );
            rval = pcomm->exchange_ghost_cells( 2, 0, i, 0, true, true, &set_with_ghosts, global_id_filter );MB_CHK_ERR( rval );
        }

        // need to set global id tags
        // need also to propagate global id to ghost cells; it is not done by default :(
        moab::Tag gtag = m_interface->globalId_tag();
        moab::Range entities;
        rval = m_interface->get_entities_by_dimension( set_with_ghosts, 2, entities );MB_CHK_ERR( rval );

        moab::Tag doftag;
        rval = m_interface->tag_get_handle( "GLOBAL_DOFS", doftag );
        if ( rval == MB_SUCCESS )
        {
            moab::Range quads = entities.subset_by_type(moab::MBQUAD);
            rval = pcomm->exchange_tags( doftag, quads );MB_CHK_ERR( rval );
        }

        // get all vertices too, need to exchange global ids for vertices too
        moab::Range vertices;
        rval = m_interface->get_connectivity( entities, vertices );MB_CHK_ERR( rval );
        entities.merge( vertices );
        rval = pcomm->exchange_tags( gtag, entities );MB_CHK_ERR( rval );
#ifdef MOAB_DBG
        std::stringstream filename1;
        filename1 << "set_with_ghosts" << m_pcomm->rank() << ".h5m";
        rval = m_interface->write_file( filename1.str().c_str(), 0, 0, &set_with_ghosts, 1 );MB_CHK_ERR( rval );
        // dump global ids of 2d entities, in order
        moab::Range cells= entities.subset_by_dimension(2);
        std::vector<int> gids(cells.size());
        rval = m_interface->tag_get_data(gtag, cells, &gids[0]);MB_CHK_ERR( rval );
        std::ofstream id_file;
        std::stringstream filename2;
        filename2 << "fileIds_" << m_pcomm->rank() << ".txt";
        id_file.open (filename2.str());
        for (size_t k=0; k<gids.size(); k++)
        {
            id_file << " " << gids[k];
            if(k%10==9) id_file << "\n";
        }
        id_file.close();

#endif
        return rval;
    }

#endif

#undef MOAB_DBG

    ErrorCode LoadNativeMesh( std::string filename,
                              moab::EntityHandle& meshset,
                              std::vector< int >& metadata,
                              const char* readopts = 0 )
    {
#ifdef MOAB_HAVE_MPI
        std::string opts = "";
        if( readopts )
        {
            if( opts.size() )
                opts = opts + ";" + std::string( readopts );
            else
                opts = std::string( readopts );
        }

        if( !m_pcomm->rank() ) std::cout << "Reading file (" << filename << ") with options = [" << opts << "]\n";
#else
        const std::string opts = std::string( ( readopts ? readopts : "" ) );
        std::cout << "Reading file (" << filename << ") with options = [" << opts << "]\n";
#endif
        moab::ErrorCode rval = m_interface->load_file( filename.c_str(), &meshset, opts.c_str() );MB_CHK_ERR( rval );

        Tag rectilinearTag;
        rval = m_interface->tag_get_handle( "ClimateMetadata", rectilinearTag );

        if( rval != MB_FAILURE && rval != MB_TAG_NOT_FOUND && rval != MB_ALREADY_ALLOCATED &&
            rectilinearTag != nullptr )
        {
            int dimSizes[3];
            rval = m_interface->tag_get_data( rectilinearTag, &meshset, 1,
                                              dimSizes );  // MB_CHK_SET_ERR( rval, "Error geting tag data" );
            metadata.clear();
            metadata.push_back( dimSizes[0] );
            metadata.push_back( dimSizes[1] );
            metadata.push_back( dimSizes[2] );
        }

        return MB_SUCCESS;
    }

  protected:
    // member data
    Interface* m_interface;

#ifdef MOAB_HAVE_MPI
    ParallelComm* m_pcomm;
#endif
};

}  // namespace moab

#endif /* MB_REMAPPER_HPP */
