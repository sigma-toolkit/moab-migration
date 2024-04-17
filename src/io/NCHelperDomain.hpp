//-------------------------------------------------------------------------
// Filename      : NCHelperFV.hpp
//
// Purpose       : Climate NC file helper for Domain grid

//-----------------------
#ifndef NCHELPERDOMAIN_HPP
#define NCHELPERDOMAIN_HPP

#include "NCHelper.hpp"

namespace moab
{
#ifdef MOAB_HAVE_MPI
class ParallelComm;
#endif
//! Child helper class for Domain grid
class NCHelperDomain : public ScdNCHelper
{
  public:
    NCHelperDomain( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
        : ScdNCHelper( readNC, fileId, opts, fileSet )
    {
    }
    static bool can_read_file( ReadNC* readNC, int fileId );

    ErrorCode create_mesh( Range& faces );

  private:
    virtual ErrorCode init_mesh_vals();
    virtual std::string get_mesh_type_name()
    {
        return "DOMAIN";
    }
#ifdef MOAB_HAVE_MPI
    ErrorCode redistribute_cells( moab::ParallelComm * myPcomm,
                                  std::vector<double> & xc,
                                  std::vector<double> & yc,
                                  std::vector<double> & xv,
                                  std::vector<double> & yv,
                                  std::vector<double> & frac,
                                  std::vector<int> & mask,
                                  std::vector<double> & area,
                                  std::vector<int> & gids,
                                  int nv,                     // number of vertices per cell
                                  bool nv_last);
#endif
    int nv;     // number of vertices per cell
    int nvDim;  // index of nv dim
};

#endif /* NCHELPERDOMAIN_HPP */
}  // namespace moab
