/*
 * NCHelperTOPO.h
 *
 *  Created on: Aug. 30, 2025
 */

#ifndef SRC_IO_NCHELPERTOPO_H_
#define SRC_IO_NCHELPERTOPO_H_

#include "NCHelper.hpp"

namespace moab
{
#ifdef MOAB_HAVE_MPI
class ParallelComm;
#endif

class NCHelperTOPO : public UcdNCHelper
{
public:
    NCHelperTOPO( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet );
    static bool can_read_file( ReadNC* readNC );


private:

    //! Implementation of NCHelper::init_mesh_vals()
    virtual ErrorCode init_mesh_vals();

    //! Implementation of NCHelper::create_mesh()
    virtual ErrorCode create_mesh( Range& faces );
    //! Implementation of NCHelper::get_mesh_type_name()
    virtual std::string get_mesh_type_name()
    {
        return "TOPO";
    }


    //! Read lat/lon coordinate variables and convert to 3D coordinates
    ErrorCode read_coordinate_variables( int startLatIdx, int endLatIdx, int startLonIdx, int endLonIdx,
                                        std::vector< double* >& arrays );

    //! Read data variables (htopo, landfract) and store on vertices
    ErrorCode read_data_variables( int startLatIdx, int endLatIdx, int startLonIdx, int endLonIdx,
                                  const Range& local_verts );

    virtual ErrorCode check_existing_mesh() {return MB_SUCCESS;}

    //! Implementation of UcdNCHelper::read_ucd_variables_to_nonset_allocate()
    virtual ErrorCode read_ucd_variables_to_nonset_allocate( std::vector< ReadNC::VarData >& vdatas,
                                                             std::vector< int >& tstep_nums );

#ifdef MOAB_HAVE_PNETCDF
    //! Implementation of UcdNCHelper::read_ucd_variables_to_nonset_async()
    virtual ErrorCode read_ucd_variables_to_nonset_async( std::vector< ReadNC::VarData >& vdatas,
                                                          std::vector< int >& tstep_nums );
#else
    //! Implementation of UcdNCHelper::read_ucd_variables_to_nonset()
    virtual ErrorCode read_ucd_variables_to_nonset( std::vector< ReadNC::VarData >& vdatas,
                                                    std::vector< int >& tstep_nums );
#endif

  private:
    std::string latDimName;
    std::string lonDimName;
    int coordDim;
    bool degrees;
    int nLatVals;
    int nLonVals;
};

}  // namespace moab

#endif /* SRC_IO_NCHELPERTOPO_H_ */