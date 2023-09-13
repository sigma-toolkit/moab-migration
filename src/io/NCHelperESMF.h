/*
 * NCHelperESMF.h
 *
 *  Created on: Sep. 12, 2023
 */

#ifndef SRC_IO_NCHELPERESMF_H_
#define SRC_IO_NCHELPERESMF_H_

#include "NCHelper.hpp"

namespace moab
{
#ifdef MOAB_HAVE_MPI
class ParallelComm;
#endif

class NCHelperESMF: public moab::UcdNCHelper {
public:
    NCHelperESMF(( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet );
    static bool can_read_file( ReadNC* readNC );

private:
    //! Implementation of NCHelper::init_mesh_vals()
    virtual ErrorCode init_mesh_vals();
    //! Implementation of NCHelper::check_existing_mesh()
    virtual ErrorCode check_existing_mesh();
    //! Implementation of NCHelper::create_mesh()
    virtual ErrorCode create_mesh( Range& faces );
    //! Implementation of NCHelper::get_mesh_type_name()
    virtual std::string get_mesh_type_name()
    {
        return "ESMF";
    }


#ifdef MOAB_HAVE_MPI
    //! Redistribute local cells after trivial partition (e.g. Zoltan partition, if applicable)
    ErrorCode redistribute_local_cells( int start_cell_index, ParallelComm* pco );
#endif

    //! Create local vertices
    ErrorCode create_local_vertices( const std::vector< int >& vertices_on_local_cells, EntityHandle& start_vertex );

    //! Create local edges (optional)
    ErrorCode create_local_edges( EntityHandle start_vertex, const std::vector< int >& num_edges_on_local_cells );

    //! Create local cells without padding (cells are divided into groups based on the number of
    //! edges)
    ErrorCode create_local_cells( const std::vector< int >& vertices_on_local_cells,
                                  const std::vector< int >& num_edges_on_local_cells,
                                  EntityHandle start_vertex,
                                  Range& faces );

    //! Create local cells with padding (padded cells will have the same number of edges)
    ErrorCode create_padded_local_cells( const std::vector< int >& vertices_on_local_cells,
                                         EntityHandle start_vertex,
                                         Range& faces );

  private:
    int maxEdgesPerCell;
    int numCellGroups;
    std::map< EntityHandle, int > cellHandleToGlobalID;
    Range facesOwned;
};

}  // namespace moab

#endif /* SRC_IO_NCHELPERESMF_H_ */
