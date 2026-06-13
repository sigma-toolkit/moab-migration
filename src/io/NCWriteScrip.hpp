//-------------------------------------------------------------------------
// Filename      : NCWriteScrip.hpp
//
// Purpose       : Write a generic polygonal mesh out in SCRIP grid format.
//
//   Counterpart of NCHelperScrip — same wire format on disk, inverse
//   direction. Drives off any 2-D unstructured mesh in MOAB; not tied
//   to a particular source format. Selected via the WriteNC option
//   "WRITE_FORMAT=SCRIP", which lets us convert from arbitrary inputs
//   (MPAS, HOMME, h5m, exodus, ...) into SCRIP without first having
//   to go through the climate read-modify-write template.
//
//   On-disk schema (matches NCHelperScrip's reader):
//     dimensions:
//       grid_size    = number of cells
//       grid_corners = max vertices per cell (smaller cells are padded
//                      with the first corner, per SCRIP convention)
//       grid_rank    = 1 (unstructured)
//     variables:
//       int    grid_dims(grid_rank)
//       double grid_center_lat(grid_size)   ; units = "degrees"
//       double grid_center_lon(grid_size)   ; units = "degrees"
//       double grid_corner_lat(grid_size, grid_corners) ; units = "degrees"
//       double grid_corner_lon(grid_size, grid_corners) ; units = "degrees"
//       int    grid_imask(grid_size)
//       (optional) double grid_area(grid_size) ; units = "radians^2"
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#ifndef NCWRITESCRIP_HPP_
#define NCWRITESCRIP_HPP_

#include "NCWriteHelper.hpp"

namespace moab
{

class NCWriteScrip : public NCWriteHelper
{
  public:
    NCWriteScrip( WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
        : NCWriteHelper( writeNC, fileId, opts, fileSet ),
          mLocalCells( 0 ), mGlobalCells( 0 ), mMaxCornersGlobal( 0 ), mHasAreas( false ),
          mDimGridSize( -1 ), mDimGridCorners( -1 ), mDimGridRank( -1 ),
          mVarGridDims( -1 ), mVarCenterLon( -1 ), mVarCenterLat( -1 ), mVarCornerLon( -1 ),
          mVarCornerLat( -1 ), mVarImask( -1 ), mVarArea( -1 )
    {
    }

    virtual ~NCWriteScrip();

    // ----- NCWriteHelper interface --------------------------------------

    //! Collect cells + compute SCRIP-ready center/corner lat/lon arrays.
    virtual ErrorCode collect_mesh_info();

    //! Override: synthesize the SCRIP schema directly — no need for
    //! the climate var_names / dim plumbing that the base init_file
    //! walks through.
    virtual ErrorCode init_file( std::vector< std::string >& var_names,
                                 std::vector< std::string >& desired_names,
                                 bool _append );

    //! Override: write SCRIP grid arrays (cells gathered to rank 0
    //! and written serially via the dispatch layer).
    virtual ErrorCode write_values( std::vector< std::string >& var_names, std::vector< int >& tstep_nums );

  protected:
    //! Required by the abstract base; no climate "nonset" variables
    //! exist in a SCRIP grid file, so this is a no-op.
    virtual ErrorCode write_nonset_variables( std::vector< WriteNC::VarData >& vdatas,
                                              std::vector< int >& tstep_nums );

  private:
    //! Per-rank cell count and the global total (set in collect_mesh_info).
    long mLocalCells;
    long mGlobalCells;

    //! Global max corners per cell. Used to size the corner arrays.
    int mMaxCornersGlobal;

    //! True if a cell-area tag was present on the input mesh and we'll
    //! emit grid_area; false otherwise.
    bool mHasAreas;

    //! Per-cell local data (size mLocalCells; corners are
    //! mLocalCells * mMaxCornersGlobal, row-major).
    std::vector< int >    mLocalGids;
    std::vector< double > mCenterLon;
    std::vector< double > mCenterLat;
    std::vector< double > mCornerLon;
    std::vector< double > mCornerLat;
    std::vector< int >    mImask;
    std::vector< double > mAreas;  // empty if !mHasAreas

    //! NetCDF dimension and variable handles populated by init_file().
    int mDimGridSize;
    int mDimGridCorners;
    int mDimGridRank;
    int mVarGridDims;
    int mVarCenterLon;
    int mVarCenterLat;
    int mVarCornerLon;
    int mVarCornerLat;
    int mVarImask;
    int mVarArea;
};

}  // namespace moab

#endif  // NCWRITESCRIP_HPP_
