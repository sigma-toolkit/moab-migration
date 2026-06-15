//-------------------------------------------------------------------------
// Filename      : NCWriteESMF.hpp
//
// Purpose       : Write a generic polygonal mesh out as an ESMF unstructured
//                 grid file. Pairs with NCHelperESMF on the read side.
//
//   On-disk schema (matches NCHelperESMF's reader):
//     dimensions:
//       nodeCount        = number of unique vertices
//       elementCount     = number of cells
//       maxNodePElement  = global max vertices per cell
//       coordDim         = 2 (lat/lon) — sphere meshes
//     variables:
//       double nodeCoords(nodeCount, coordDim)    ; units = "degrees"
//                                                 ; ordering: [lon, lat]
//       int    elementConn(elementCount, maxNodePElement)  ; 1-based vertex indices
//       int    numElementConn(elementCount)       ; actual #nodes per element
//       double centerCoords(elementCount, coordDim) ; units = "degrees", [lon,lat]
//       (optional) double elementArea(elementCount) ; units = "radians^2"
//       (optional) int    elementMask(elementCount)
//
//   The 1-based connectivity and [lon, lat] ordering match the
//   reader's expectations (NCHelperESMF.cpp:334 / :441).
//
//   Creator     : Vijay Mahadevan, 2026-06-13
//-------------------------------------------------------------------------

#ifndef NCWRITEESMF_HPP_
#define NCWRITEESMF_HPP_

#include "NCWriteHelper.hpp"

namespace moab
{

class NCWriteESMF : public NCWriteHelper
{
  public:
    NCWriteESMF( WriteNC* writeNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
        : NCWriteHelper( writeNC, fileId, opts, fileSet ),
          mLocalCells( 0 ), mGlobalCells( 0 ), mGlobalNodes( 0 ), mMaxCornersGlobal( 0 ), mCoordDim( 2 ),
          mHasAreas( false ), mHasMask( false ),
          mDimNodeCount( -1 ), mDimElementCount( -1 ), mDimMaxNodePElement( -1 ), mDimCoordDim( -1 ),
          mVarNodeCoords( -1 ), mVarElementConn( -1 ), mVarNumElementConn( -1 ),
          mVarCenterCoords( -1 ), mVarElementArea( -1 ), mVarElementMask( -1 )
    {
    }

    virtual ~NCWriteESMF() override;

    ErrorCode collect_mesh_info() override;
    ErrorCode init_file( std::vector< std::string >& var_names,
                         std::vector< std::string >& desired_names,
                         bool _append ) override;
    ErrorCode write_values( std::vector< std::string >& var_names, std::vector< int >& tstep_nums ) override;

  protected:
    ErrorCode write_nonset_variables( std::vector< WriteNC::VarData >& vdatas,
                                      std::vector< int >& tstep_nums ) override;

  private:
    long mLocalCells;
    long mGlobalCells;
    long mGlobalNodes;
    int  mMaxCornersGlobal;
    int  mCoordDim;     // always 2 for our writer (lat/lon sphere meshes)
    bool mHasAreas;
    bool mHasMask;

    // Per-rank node and cell arrays (all unique, owned-only).
    std::vector< int >    mLocalCellGids;
    std::vector< double > mNodeLon;       // size mLocalNodes
    std::vector< double > mNodeLat;       // size mLocalNodes
    std::vector< int >    mNodeGids;      // global vertex id for dedup at rank 0

    std::vector< int >    mElementConn;   // size mLocalCells * mMaxCornersGlobal; global vertex GIDs (1-based)
    std::vector< int >    mElementNumNodes; // actual node count per cell
    std::vector< double > mCenterLon;
    std::vector< double > mCenterLat;
    std::vector< double > mAreas;         // empty if !mHasAreas
    std::vector< int >    mMask;          // empty if !mHasMask

    int mDimNodeCount;
    int mDimElementCount;
    int mDimMaxNodePElement;
    int mDimCoordDim;
    int mVarNodeCoords;
    int mVarElementConn;
    int mVarNumElementConn;
    int mVarCenterCoords;
    int mVarElementArea;
    int mVarElementMask;
};

}  // namespace moab

#endif  // NCWRITEESMF_HPP_
