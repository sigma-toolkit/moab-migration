//-------------------------------------------------------------------------
// Filename      : NCWriteDomain.hpp
//
// Purpose       : Write a generic polygonal mesh out as a CESM domain
//                 NetCDF file. Pairs with NCHelperDomain on the read side.
//
//   On-disk schema (matches NCHelperDomain's reader):
//     dimensions:
//       n  = ni * nj (number of cells)
//       ni = i-extent (set to the total cell count for unstructured input)
//       nj = j-extent (set to 1 for unstructured input)
//       nv = max vertices per cell
//     variables:
//       double xc(nj, ni)         ; long_name = "longitude of grid cell center"
//                                 ; units = "degrees_east"
//                                 ; bounds = "xv"
//       double yc(nj, ni)         ; long_name = "latitude of grid cell center"
//                                 ; units = "degrees_north"
//                                 ; bounds = "yv"
//       double xv(nj, ni, nv)     ; long_name = "longitude of grid cell verticies"
//                                 ; units = "degrees_east"
//       double yv(nj, ni, nv)     ; long_name = "latitude of grid cell verticies"
//                                 ; units = "degrees_north"
//       int    mask(nj, ni)       ; long_name = "domain mask"
//       double area(nj, ni)       ; units = "radian2"
//       double frac(nj, ni)       ; fraction of grid cell that is active
//
//   Layout note: this first cut treats any input mesh as a flat
//   ni = total-cells, nj = 1 layout. The sample CESM domain files
//   for unstructured inputs (e.g. ne4np4_oQU240 with 866 cells, ni=866,
//   nj=1) follow exactly this convention, so reloading the output via
//   NCHelperDomain produces the original cell count. Genuine
//   rectangular (ni > 1, nj > 1) output is a follow-up that needs
//   structured-mesh topology detection from the input.
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#ifndef NCWRITEDOMAIN_HPP_
#define NCWRITEDOMAIN_HPP_

#include "NCWriteHelper.hpp"

namespace moab
{

class NCWriteDomain : public NCWriteHelper
{
  public:
    NCWriteDomain( WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
        : NCWriteHelper( writeNC, fileId, opts, fileSet ),
          mLocalCells( 0 ), mGlobalCells( 0 ), mMaxCornersGlobal( 0 ),
          mNi( 0 ), mNj( 1 ), mHasAreas( false ), mHasMask( false ), mHasFrac( false ),
          mDimN( -1 ), mDimNi( -1 ), mDimNj( -1 ), mDimNv( -1 ),
          mVarXc( -1 ), mVarYc( -1 ), mVarXv( -1 ), mVarYv( -1 ),
          mVarMask( -1 ), mVarArea( -1 ), mVarFrac( -1 )
    {
    }

    virtual ~NCWriteDomain();

    virtual ErrorCode collect_mesh_info();
    virtual ErrorCode init_file( std::vector< std::string >& var_names,
                                 std::vector< std::string >& desired_names,
                                 bool _append );
    virtual ErrorCode write_values( std::vector< std::string >& var_names, std::vector< int >& tstep_nums );

  protected:
    virtual ErrorCode write_nonset_variables( std::vector< WriteNC::VarData >& vdatas,
                                              std::vector< int >& tstep_nums );

  private:
    long mLocalCells;
    long mGlobalCells;
    int  mMaxCornersGlobal;

    //! On-disk grid extents. For now we always emit ni = mGlobalCells,
    //! nj = 1 — i.e. a flat layout that the CESM domain reader handles
    //! identically to the rectangular case but with j collapsed.
    long mNi;
    long mNj;

    bool mHasAreas;
    bool mHasMask;
    bool mHasFrac;

    std::vector< int >    mLocalGids;
    std::vector< double > mXc;
    std::vector< double > mYc;
    std::vector< double > mXv;  // size mLocalCells * mMaxCornersGlobal, row-major
    std::vector< double > mYv;
    std::vector< int >    mMask;
    std::vector< double > mAreas;
    std::vector< double > mFrac;

    int mDimN;
    int mDimNi;
    int mDimNj;
    int mDimNv;
    int mVarXc;
    int mVarYc;
    int mVarXv;
    int mVarYv;
    int mVarMask;
    int mVarArea;
    int mVarFrac;
};

}  // namespace moab

#endif  // NCWRITEDOMAIN_HPP_
