/*
 * =====================================================================================
 *
 *       Filename:  TempestOnlineMap.hpp
 *
 *    Description:  Interface to the TempestRemap library to compute the consistent,
 *                  and accurate high-order conservative remapping weights for overlap
 *                  grids on the sphere in climate simulations.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#include "Announce.h"
#include "DataArray3D.h"
#include "FiniteVolumeTools.h"
#include "FiniteElementTools.h"
#include "TriangularQuadrature.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"
#include "LinearRemapFV.h"

#include "LinearRemapSE0.h"
#include "LinearRemapFV.h"

#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/Remapping/IntegerReprosum.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "DebugOutput.hpp"
#include "moab/TupleList.hpp"
#include "moab/MeshTopoUtil.hpp"

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <algorithm>

#ifdef MOAB_HAVE_NETCDFPAR
#include "netcdfcpp_par.hpp"
#else
#include "netcdfcpp.h"
#endif

// #define USE_NATIVE_TEMPESTREMAP_ROUTINES

///////////////////////////////////////////////////////////////////////////////

// #define VERBOSE
// #define VVERBOSE
// #define CHECK_INCREASING_DOF

///////////////////////////////////////////////////////////////////////////////

#define MPI_CHK_ERR( err )                                          \
    if( err )                                                       \
    {                                                               \
        std::cout << "MPI Failure. ErrorCode (" << ( err ) << ") "; \
        std::cout << "\nMPI Aborting... \n";                        \
        return moab::MB_FAILURE;                                    \
    }

moab::TempestOnlineMap::TempestOnlineMap( moab::TempestRemapper* remapper ) : OfflineMap(), m_remapper( remapper )
{
    // Get the references for the MOAB core objects
    m_interface = m_remapper->get_interface();
#ifdef MOAB_HAVE_MPI
    m_pcomm = m_remapper->get_parallel_communicator();
#endif

    // now let us re-update the reference to the input source mesh
    m_meshInput = m_remapper->GetMesh( moab::Remapper::SourceMesh );
    // now let us re-update the reference to the covering mesh
    m_meshInputCov = m_remapper->GetCoveringMesh();
    // now let us re-update the reference to the output target mesh
    m_meshOutput = m_remapper->GetMesh( moab::Remapper::TargetMesh );
    // now let us re-update the reference to the output target mesh
    m_meshOverlap = m_remapper->GetMesh( moab::Remapper::OverlapMesh );

    is_parallel = remapper->is_parallel;
    is_root     = remapper->is_root;
    rank        = remapper->rank;
    size        = remapper->size;

    // set default order
    m_input_order = m_output_order = 1;

    // unknown until a map file is read (ReadParallelMap sets it to n_a)
    m_nTotDofs_SrcGlobal = -1;

    // Initialize dimension information from file
    this->setup_sizes_dimensions();
}

void moab::TempestOnlineMap::setup_sizes_dimensions()
{
    if( m_meshInputCov )
    {
        std::vector< std::string > dimNames;
        std::vector< int > dimSizes;
        dimNames.push_back( "num_elem" );
        dimSizes.push_back( m_meshInputCov->faces.size() );

        this->InitializeSourceDimensions( dimNames, dimSizes );
    }

    if( m_meshOutput )
    {
        std::vector< std::string > dimNames;
        std::vector< int > dimSizes;
        dimNames.push_back( "num_elem" );
        dimSizes.push_back( m_meshOutput->faces.size() );

        this->InitializeTargetDimensions( dimNames, dimSizes );
    }
}

///////////////////////////////////////////////////////////////////////////////

moab::TempestOnlineMap::~TempestOnlineMap()
{
    m_interface = nullptr;
#ifdef MOAB_HAVE_MPI
    m_pcomm = nullptr;
#endif
    m_meshInput   = nullptr;
    m_meshOutput  = nullptr;
    m_meshOverlap = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::SetDOFmapTags( const std::string srcDofTagName,
                                                       const std::string tgtDofTagName )
{
    moab::ErrorCode rval;

    int tagSize = 0;
    tagSize     = ( m_eInputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Src * m_nDofsPEl_Src );
    rval =
        m_interface->tag_get_handle( srcDofTagName.c_str(), tagSize, MB_TYPE_INTEGER, this->m_dofTagSrc, MB_TAG_ANY );

    if( rval == moab::MB_TAG_NOT_FOUND && m_eInputType != DiscretizationType_FV )
    {
        MB_CHK_SET_ERR( MB_FAILURE, "DoF tag is not set correctly for source mesh." );
    }
    else
        MB_CHK_ERR( rval );

    tagSize = ( m_eOutputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Dest * m_nDofsPEl_Dest );
    rval =
        m_interface->tag_get_handle( tgtDofTagName.c_str(), tagSize, MB_TYPE_INTEGER, this->m_dofTagDest, MB_TAG_ANY );
    if( rval == moab::MB_TAG_NOT_FOUND && m_eOutputType != DiscretizationType_FV )
    {
        MB_CHK_SET_ERR( MB_FAILURE, "DoF tag is not set correctly for target mesh." );
    }
    else
        MB_CHK_ERR( rval );

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::SetDOFmapAssociation( DiscretizationType srcType,
                                                              int srcOrder,
                                                              bool isSrcContinuous,
                                                              DataArray3D< int >* srcdataGLLNodes,
                                                              DataArray3D< int >* srcdataGLLNodesSrc,
                                                              DiscretizationType destType,
                                                              int destOrder,
                                                              bool isTgtContinuous,
                                                              DataArray3D< int >* tgtdataGLLNodes )
{
    std::vector< bool > dgll_cgll_row_ldofmap, dgll_cgll_col_ldofmap, dgll_cgll_covcol_ldofmap;
    std::vector< int > src_soln_gdofs, locsrc_soln_gdofs, tgt_soln_gdofs;

    // We are assuming that these are element based tags that are sized: np * np
    m_srcDiscType  = srcType;
    m_destDiscType = destType;
    m_input_order  = srcOrder;
    m_output_order = destOrder;

    bool vprint = is_root && false;

    // Compute and store the total number of source and target DoFs corresponding
    // to number of rows and columns in the mapping.
    // Now compute the mapping and store it for the covering mesh
    int srcTagSize = ( m_eInputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Src * m_nDofsPEl_Src );
    if( m_remapper->point_cloud_source )
    {
        assert( m_nDofsPEl_Src == 1 );
        col_gdofmap.resize( m_remapper->m_covering_source_vertices.size(), UINT_MAX );
        col_dtoc_dofmap.resize( m_remapper->m_covering_source_vertices.size(), -1 );
        src_soln_gdofs.resize( m_remapper->m_covering_source_vertices.size(), -1 );
        MB_CHK_ERR(
            m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_covering_source_vertices, &src_soln_gdofs[0] ) );
        srcTagSize = 1;
    }
    else
    {
        col_gdofmap.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX );
        col_dtoc_dofmap.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, -1 );
        src_soln_gdofs.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, -1 );
        MB_CHK_ERR(
            m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_covering_source_entities, &src_soln_gdofs[0] ) );
    }

    m_nTotDofs_SrcCov = 0;
    if( srcdataGLLNodes == nullptr )
    {
        /* we only have a mapping for elements as DoFs */
        for( unsigned i = 0; i < col_gdofmap.size(); ++i )
        {
            auto gdof = src_soln_gdofs[i];
            assert( gdof > 0 );
            col_gdofmap[i]     = gdof - 1;
            col_dtoc_dofmap[i] = i;
            if( vprint ) std::cout << "Col: " << i << ", " << col_gdofmap[i] << "\n";
            m_nTotDofs_SrcCov++;
        }
    }
    else
    {
        if( isSrcContinuous )
            dgll_cgll_covcol_ldofmap.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, false );
        // Put these remap coefficients into the SparseMatrix map
        for( unsigned j = 0; j < m_remapper->m_covering_source_entities.size(); j++ )
        {
            for( int p = 0; p < m_nDofsPEl_Src; p++ )
            {
                for( int q = 0; q < m_nDofsPEl_Src; q++ )
                {
                    const int localDOF  = ( *srcdataGLLNodes )[p][q][j] - 1;
                    const int offsetDOF = j * srcTagSize + p * m_nDofsPEl_Src + q;
                    if( isSrcContinuous && !dgll_cgll_covcol_ldofmap[localDOF] )
                    {
                        m_nTotDofs_SrcCov++;
                        dgll_cgll_covcol_ldofmap[localDOF] = true;
                    }
                    if( !isSrcContinuous ) m_nTotDofs_SrcCov++;
                    assert( src_soln_gdofs[offsetDOF] > 0 );
                    // For CGLL: weight matrix uses localDOF (continuous shared node index)
                    //   as column index → col_gdofmap must be indexed by localDOF
                    // For DGLL: weight matrix uses offsetDOF (= elem*nP*nP + p*nP + q)
                    //   as column index → col_gdofmap must be indexed by offsetDOF
                    if( isSrcContinuous )
                    {
                        col_gdofmap[localDOF]      = src_soln_gdofs[offsetDOF] - 1;
                        col_dtoc_dofmap[offsetDOF] = localDOF;
                    }
                    else
                    {
                        col_gdofmap[offsetDOF]     = src_soln_gdofs[offsetDOF] - 1;
                        col_dtoc_dofmap[offsetDOF] = offsetDOF;
                    }
                }
            }
        }
    }

    if( m_remapper->point_cloud_source )
    {
        assert( m_nDofsPEl_Src == 1 );
        srccol_gdofmap.resize( m_remapper->m_source_vertices.size(), UINT_MAX );
        srccol_dtoc_dofmap.resize( m_remapper->m_covering_source_vertices.size(), -1 );
        locsrc_soln_gdofs.resize( m_remapper->m_source_vertices.size(), -1 );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_source_vertices, &locsrc_soln_gdofs[0] ) );
    }
    else
    {
        srccol_gdofmap.resize( m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX );
        srccol_dtoc_dofmap.resize( m_remapper->m_source_entities.size() * srcTagSize, -1 );
        locsrc_soln_gdofs.resize( m_remapper->m_source_entities.size() * srcTagSize, -1 );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_source_entities, &locsrc_soln_gdofs[0] ) );
    }

    // Now compute the mapping and store it for the original source mesh
    m_nTotDofs_Src = 0;
    if( srcdataGLLNodesSrc == nullptr )
    {
        /* we only have a mapping for elements as DoFs */
        for( unsigned i = 0; i < srccol_gdofmap.size(); ++i )
        {
            auto gdof = locsrc_soln_gdofs[i];
            assert( gdof > 0 );
            srccol_gdofmap[i]     = gdof - 1;
            srccol_dtoc_dofmap[i] = i;
            m_nTotDofs_Src++;
        }
    }
    else
    {
        if( isSrcContinuous ) dgll_cgll_col_ldofmap.resize( m_remapper->m_source_entities.size() * srcTagSize, false );
        // Put these remap coefficients into the SparseMatrix map
        for( unsigned j = 0; j < m_remapper->m_source_entities.size(); j++ )
        {
            for( int p = 0; p < m_nDofsPEl_Src; p++ )
            {
                for( int q = 0; q < m_nDofsPEl_Src; q++ )
                {
                    const int localDOF  = ( *srcdataGLLNodesSrc )[p][q][j] - 1;
                    const int offsetDOF = j * srcTagSize + p * m_nDofsPEl_Src + q;
                    if( isSrcContinuous && !dgll_cgll_col_ldofmap[localDOF] )
                    {
                        m_nTotDofs_Src++;
                        dgll_cgll_col_ldofmap[localDOF] = true;
                    }
                    if( !isSrcContinuous ) m_nTotDofs_Src++;
                    assert( locsrc_soln_gdofs[offsetDOF] > 0 );
                    if( isSrcContinuous )
                    {
                        srccol_gdofmap[localDOF]      = locsrc_soln_gdofs[offsetDOF] - 1;
                        srccol_dtoc_dofmap[offsetDOF] = localDOF;
                    }
                    else
                    {
                        srccol_gdofmap[offsetDOF]     = locsrc_soln_gdofs[offsetDOF] - 1;
                        srccol_dtoc_dofmap[offsetDOF] = offsetDOF;
                    }
                }
            }
        }
    }

    int tgtTagSize = ( m_eOutputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Dest * m_nDofsPEl_Dest );
    if( m_remapper->point_cloud_target )
    {
        assert( m_nDofsPEl_Dest == 1 );
        row_gdofmap.resize( m_remapper->m_target_vertices.size(), UINT_MAX );
        row_dtoc_dofmap.resize( m_remapper->m_target_vertices.size(), -1 );
        tgt_soln_gdofs.resize( m_remapper->m_target_vertices.size(), -1 );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagDest, m_remapper->m_target_vertices, &tgt_soln_gdofs[0] ) );
        tgtTagSize = 1;
    }
    else
    {
        row_gdofmap.resize( m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX );
        row_dtoc_dofmap.resize( m_remapper->m_target_entities.size() * tgtTagSize, -1 );
        tgt_soln_gdofs.resize( m_remapper->m_target_entities.size() * tgtTagSize, -1 );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagDest, m_remapper->m_target_entities, &tgt_soln_gdofs[0] ) );
    }

    // Now compute the mapping and store it for the target mesh
    // To access the GID for each row: row_gdofmap [ row_ldofmap [ 0 : local_ndofs ] ] = GDOF
    m_nTotDofs_Dest = 0;
    if( tgtdataGLLNodes == nullptr )
    {
        /* we only have a mapping for elements as DoFs */
        for( unsigned i = 0; i < row_gdofmap.size(); ++i )
        {
            auto gdof = tgt_soln_gdofs[i];
            assert( gdof > 0 );
            row_gdofmap[i]     = gdof - 1;
            row_dtoc_dofmap[i] = i;
            if( vprint ) std::cout << "Row: " << i << ", " << row_gdofmap[i] << "\n";
            m_nTotDofs_Dest++;
        }
    }
    else
    {
        if( isTgtContinuous ) dgll_cgll_row_ldofmap.resize( m_remapper->m_target_entities.size() * tgtTagSize, false );
        // Put these remap coefficients into the SparseMatrix map
        for( unsigned j = 0; j < m_remapper->m_target_entities.size(); j++ )
        {
            for( int p = 0; p < m_nDofsPEl_Dest; p++ )
            {
                for( int q = 0; q < m_nDofsPEl_Dest; q++ )
                {
                    const int localDOF  = ( *tgtdataGLLNodes )[p][q][j] - 1;
                    const int offsetDOF = j * tgtTagSize + p * m_nDofsPEl_Dest + q;
                    if( isTgtContinuous && !dgll_cgll_row_ldofmap[localDOF] )
                    {
                        m_nTotDofs_Dest++;
                        dgll_cgll_row_ldofmap[localDOF] = true;
                    }
                    if( !isTgtContinuous ) m_nTotDofs_Dest++;
                    assert( tgt_soln_gdofs[offsetDOF] > 0 );
                    if( isTgtContinuous )
                    {
                        row_gdofmap[localDOF]      = tgt_soln_gdofs[offsetDOF] - 1;
                        row_dtoc_dofmap[offsetDOF] = localDOF;
                    }
                    else
                    {
                        row_gdofmap[offsetDOF]     = tgt_soln_gdofs[offsetDOF] - 1;
                        row_dtoc_dofmap[offsetDOF] = offsetDOF;
                    }
                    if( vprint )
                        std::cout << "Row: " << offsetDOF << ", " << localDOF << ", " << row_gdofmap[offsetDOF] << ", "
                                  << m_nTotDofs_Dest << "\n";
                }
            }
        }
    }

    // Let us also allocate the local representation of the sparse matrix
#if defined( MOAB_HAVE_EIGEN3 ) && defined( VERBOSE )
    if( is_root )
    {
        std::cout << "[" << rank << "] DoFs: row = " << m_nTotDofs_Dest << " (gdofmap.size=" << row_gdofmap.size()
                  << "), col_src = " << m_nTotDofs_Src << ", col_cov = " << m_nTotDofs_SrcCov
                  << " (gdofmap.size=" << col_gdofmap.size() << ")\n";
    }
#endif

    // check monotonicity of row_gdofmap and col_gdofmap
#ifdef CHECK_INCREASING_DOF
    for( size_t i = 0; i < row_gdofmap.size() - 1; i++ )
    {
        if( row_gdofmap[i] > row_gdofmap[i + 1] )
            std::cout << " on rank " << rank << " in row_gdofmap[" << i << "]=" << row_gdofmap[i] << " > row_gdofmap["
                      << i + 1 << "]=" << row_gdofmap[i + 1] << " \n";
    }
    for( size_t i = 0; i < col_gdofmap.size() - 1; i++ )
    {
        if( col_gdofmap[i] > col_gdofmap[i + 1] )
            std::cout << " on rank " << rank << " in col_gdofmap[" << i << "]=" << col_gdofmap[i] << " > col_gdofmap["
                      << i + 1 << "]=" << col_gdofmap[i + 1] << " \n";
    }
#endif

    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::set_col_dc_dofs( std::vector< int >& values_entities )
{
    // col_gdofmap has global dofs , that should be in the list of values, such that
    // row_dtoc_dofmap[offsetDOF] = localDOF;
    // we need to find col_dtoc_dofmap such that: col_gdofmap[ col_dtoc_dofmap[i] ] == values_entities [i];
    // we know that col_gdofmap[0..(nbcols-1)] = global_col_dofs -> in values_entities
    // form first inverse
    //
    // resize and initialize to -1 to signal that this value should not be used, if not set below
    col_dtoc_dofmap.resize( values_entities.size(), -1 );
    for( size_t j = 0; j < values_entities.size(); j++ )
    {
        // values are 1 based, but rowMap, colMap are not
        const auto it = colMap.find( values_entities[j] - 1 );
        if( it != colMap.end() ) col_dtoc_dofmap[j] = it->second;
    }
    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::set_row_dc_dofs( std::vector< int >& values_entities )
{
    //  we need to find row_dtoc_dofmap such that: row_gdofmap[ row_dtoc_dofmap[i] ] == values_entities [i];
    // resize and initialize to -1 to signal that this value should not be used, if not set below
    row_dtoc_dofmap.resize( values_entities.size(), -1 );
    for( size_t j = 0; j < values_entities.size(); j++ )
    {
        // values are 1 based, but rowMap, colMap are not
        const auto it = rowMap.find( values_entities[j] - 1 );
        if( it != rowMap.end() ) row_dtoc_dofmap[j] = it->second;
    }
    return moab::MB_SUCCESS;
}

// Compute which weight-matrix columns the migrated coverage covers.
// delivered[mc] == true iff some covering cell maps to matrix column mc.
static void compute_delivered_columns( int ncols, const std::vector< int >& col_dtoc_dofmap,
                                       std::vector< bool >& delivered )
{
    delivered.assign( ncols, false );
    for( size_t k = 0; k < col_dtoc_dofmap.size(); k++ )
    {
        const int mc = col_dtoc_dofmap[k];
        if( mc >= 0 && mc < ncols ) delivered[mc] = true;
    }
}

int moab::TempestOnlineMap::CountAbsentColumns( int& first_absent_gid ) const
{
    first_absent_gid = -1;
    const int ncols = m_nTotDofs_SrcCov;
    std::vector< bool > delivered;
    compute_delivered_columns( ncols, col_dtoc_dofmap, delivered );
    int cnt = 0;
    for( int mc = 0; mc < ncols; mc++ )
    {
        if( !delivered[mc] )
        {
            cnt++;
            if( first_absent_gid < 0 && mc < (int)col_gdofmap.size() )
                first_absent_gid = (int)col_gdofmap[mc] + 1;  // col_gdofmap is 0-based
        }
    }
    return cnt;
}

int moab::TempestOnlineMap::DropAbsentColumns()
{
    const int ncols = m_nTotDofs_SrcCov;
    std::vector< bool > delivered;
    compute_delivered_columns( ncols, col_dtoc_dofmap, delivered );

    // Zero every stored coefficient whose column was not supplied by coverage.
    // The projection already treats these columns as zero-source (ApplyWeights
    // leaves m_colVector at 0 for them), so this changes no projected value; it
    // only lets the dual-map CAAS bounds loop skip them via its |w|<1e-50 test.
    int dropped = 0;
    for( int r = 0; r < m_weightMatrix.outerSize(); r++ )
    {
        for( WeightMatrix::InnerIterator it( m_weightMatrix, r ); it; ++it )
        {
            const int mc = (int)it.col();
            if( mc < 0 || mc >= ncols || !delivered[mc] )
            {
                if( it.value() != 0.0 ) dropped++;
                it.valueRef() = 0.0;
            }
        }
    }
    // Remove the explicit zeros so iterators no longer visit them.
    m_weightMatrix.prune( []( const Eigen::Index&, const Eigen::Index&, const double& v ) { return v != 0.0; } );
    return dropped;
}
///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::GenerateRemappingWeights( std::string strInputType,
                                                                  std::string strOutputType,
                                                                  const GenerateOfflineMapAlgorithmOptions& mapOptions,
                                                                  const std::string& srcDofTagName,
                                                                  const std::string& tgtDofTagName )
{
    NcError error( NcError::silent_nonfatal );

    moab::DebugOutput dbgprint( std::cout, rank, 0 );
    dbgprint.set_prefix( "[TempestOnlineMap]: " );
    moab::ErrorCode rval;

    const bool m_bPointCloudSource = ( m_remapper->point_cloud_source );
    const bool m_bPointCloudTarget = ( m_remapper->point_cloud_target );
    const bool m_bPointCloud       = m_bPointCloudSource || m_bPointCloudTarget;

    // Build a matrix of source and target discretization so that we know how
    // to assign the global DoFs in parallel for the mapping weights.
    // For example,
    //   for FV->FV: the rows represented target DoFs and cols represent source DoFs
    try
    {
        // Check command line parameters (data type arguments)
        STLStringHelper::ToLower( strInputType );
        STLStringHelper::ToLower( strOutputType );

        DiscretizationType eInputType;
        DiscretizationType eOutputType;

        if( strInputType == "fv" )
        {
            eInputType = DiscretizationType_FV;
        }
        else if( strInputType == "cgll" )
        {
            eInputType = DiscretizationType_CGLL;
        }
        else if( strInputType == "dgll" )
        {
            eInputType = DiscretizationType_DGLL;
        }
        else if( strInputType == "pcloud" )
        {
            eInputType = DiscretizationType_PCLOUD;
        }
        else
        {
            _EXCEPTION1( "Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]", strInputType.c_str() );
        }

        if( strOutputType == "fv" )
        {
            eOutputType = DiscretizationType_FV;
        }
        else if( strOutputType == "cgll" )
        {
            eOutputType = DiscretizationType_CGLL;
        }
        else if( strOutputType == "dgll" )
        {
            eOutputType = DiscretizationType_DGLL;
        }
        else if( strOutputType == "pcloud" )
        {
            eOutputType = DiscretizationType_PCLOUD;
        }
        else
        {
            _EXCEPTION1( "Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]", strOutputType.c_str() );
        }

        // set all required input params
        m_bConserved  = !mapOptions.fNoConservation;
        m_eInputType  = eInputType;
        m_eOutputType = eOutputType;

        // Method flags
        std::string strMapAlgorithm( "" );
        int nMonotoneType = ( mapOptions.fMonotone ) ? ( 1 ) : ( 0 );

        // Make an index of method arguments
        std::set< std::string > setMethodStrings;
        {
            int iLast = 0;
            for( size_t i = 0; i <= mapOptions.strMethod.length(); i++ )
            {
                if( ( i == mapOptions.strMethod.length() ) || ( mapOptions.strMethod[i] == ';' ) )
                {
                    std::string strMethodString = mapOptions.strMethod.substr( iLast, i - iLast );
                    STLStringHelper::RemoveWhitespaceInPlace( strMethodString );
                    if( strMethodString.length() > 0 )
                    {
                        setMethodStrings.insert( strMethodString );
                    }
                    iLast = i + 1;
                }
            }
        }

        for( const auto& it : setMethodStrings )
        {
            // Piecewise constant monotonicity
            if( it == "mono2" )
            {
                if( ( m_eInputType == DiscretizationType_FV ) && ( m_eOutputType == DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"mono2\" is only used when remapping to/from CGLL or DGLL grids" );
                }
                nMonotoneType = 2;

                // Piecewise linear monotonicity
            }
            else if( it == "mono3" )
            {
                if( ( m_eInputType == DiscretizationType_FV ) && ( m_eOutputType == DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"mono3\" is only used when remapping to/from CGLL or DGLL grids" );
                }
                nMonotoneType = 3;

                // Volumetric remapping from FV to GLL
            }
            else if( it == "volumetric" )
            {
                if( ( m_eInputType != DiscretizationType_FV ) || ( m_eOutputType == DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"volumetric\" may only be used for FV->CGLL or FV->DGLL remapping" );
                }
                strMapAlgorithm = "volumetric";

                // Inverse distance mapping
            }
            else if( it == "invdist" )
            {
                if( ( m_eInputType != DiscretizationType_FV ) || ( m_eOutputType != DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"invdist\" may only be used for FV->FV remapping" );
                }
                strMapAlgorithm = "invdist";

                // Delaunay triangulation mapping
            }
            else if( it == "delaunay" )
            {
                if( ( m_eInputType != DiscretizationType_FV ) || ( m_eOutputType != DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"delaunay\" may only be used for FV->FV remapping" );
                }
                strMapAlgorithm = "delaunay";

                // Bilinear
            }
            else if( it == "bilin" )
            {
                if( ( m_eInputType != DiscretizationType_FV ) || ( m_eOutputType != DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"bilin\" may only be used for FV->FV remapping" );
                }
                strMapAlgorithm = "fvbilin";

                // Integrated bilinear (same as mono3 when source grid is CGLL/DGLL)
            }
            else if( it == "intbilin" )
            {
                if( m_eOutputType != DiscretizationType_FV )
                {
                    _EXCEPTIONT( "--method \"intbilin\" may only be used when mapping to FV." );
                }
                if( m_eInputType == DiscretizationType_FV )
                {
                    strMapAlgorithm = "fvintbilin";
                }
                else
                {
                    strMapAlgorithm = "mono3";
                }

                // Integrated bilinear with generalized Barycentric coordinates
            }
            else if( it == "intbilingb" )
            {
                if( ( m_eInputType != DiscretizationType_FV ) || ( m_eOutputType != DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"intbilingb\" may only be used for FV->FV remapping" );
                }
                strMapAlgorithm = "fvintbilingb";
            }
            else
            {
                _EXCEPTION1( "Invalid --method argument \"%s\"", it.c_str() );
            }
        }

        m_nDofsPEl_Src =
            ( m_eInputType == DiscretizationType_FV || m_eInputType == DiscretizationType_PCLOUD ? 1
                                                                                                 : mapOptions.nPin );
        m_nDofsPEl_Dest =
            ( m_eOutputType == DiscretizationType_FV || m_eOutputType == DiscretizationType_PCLOUD ? 1
                                                                                                   : mapOptions.nPout );

        // Set the source and target mesh objects
        MB_CHK_ERR( SetDOFmapTags( srcDofTagName, tgtDofTagName ) );

        ///   the tag should be created already in the e3sm workflow; if not, create it here
        Tag areaTag;
        rval = m_interface->tag_get_handle( "aream", 1, MB_TYPE_DOUBLE, areaTag,
                                            MB_TAG_DENSE | MB_TAG_EXCL | MB_TAG_CREAT );
        if( MB_ALREADY_ALLOCATED == rval )
        {
            if( is_root ) dbgprint.printf( 0, "aream tag already defined \n" );
        }

        double local_areas[3] = { 0.0, 0.0, 0.0 }, global_areas[3] = { 0.0, 0.0, 0.0 };
        if( !m_bPointCloudSource )
        {
            // Calculate Input Mesh Face areas
            if( is_root ) dbgprint.printf( 0, "Calculating input mesh Face areas\n" );
            local_areas[0] = m_meshInput->CalculateFaceAreas( mapOptions.fSourceConcave );
            // Set source element areas as tag on the source mesh
            MB_CHK_ERR( m_interface->tag_set_data( areaTag, m_remapper->m_source_entities, m_meshInput->vecFaceArea ) );

            // Update coverage source mesh areas as well.
            m_meshInputCov->CalculateFaceAreas( mapOptions.fSourceConcave );
        }

        if( !m_bPointCloudTarget )
        {
            // Calculate Output Mesh Face areas
            if( is_root ) dbgprint.printf( 0, "Calculating output mesh Face areas\n" );
            local_areas[1] = m_meshOutput->CalculateFaceAreas( mapOptions.fTargetConcave );
            // Set target element areas as tag on the target mesh
            MB_CHK_ERR(
                m_interface->tag_set_data( areaTag, m_remapper->m_target_entities, m_meshOutput->vecFaceArea ) );
        }

        if( !m_bPointCloud )
        {
            // Calculate Face areas
            if (m_meshOverlap)
            {
                // Verify that overlap mesh is in the correct order (sanity check)
                assert( m_meshOverlap->vecSourceFaceIx.size() == m_meshOverlap->vecTargetFaceIx.size() );

                if( is_root ) dbgprint.printf( 0, "Calculating overlap mesh Face areas\n" );
                local_areas[2] =
                    m_meshOverlap->CalculateFaceAreas( mapOptions.fSourceConcave || mapOptions.fTargetConcave );
            }

            // store it as global output for now - used later in reduction
            std::copy( local_areas, local_areas + 3, global_areas );
#ifdef MOAB_HAVE_MPI
            // reduce the local source, target and overlap mesh areas to global areas
            if( m_pcomm && is_parallel )
                MPI_Reduce( local_areas, global_areas, 3, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
#endif
            if( is_root )
            {
                dbgprint.printf( 0, "Input Mesh Geometric Area: %1.15e\n", global_areas[0] );
                dbgprint.printf( 0, "Output Mesh Geometric Area: %1.15e\n", global_areas[1] );
                if (m_meshOverlap) dbgprint.printf( 0, "Overlap Mesh Recovered Area: %1.15e\n", global_areas[2] );
            }

            // Correct areas to match the areas calculated in the overlap mesh
            constexpr bool fCorrectAreas = true;
            if( fCorrectAreas && m_meshOverlap )  // In MOAB-TempestRemap, we will always keep this to be true
            {
                if( is_root ) dbgprint.printf( 0, "Correcting source/target areas to overlap mesh areas\n" );
                DataArray1D< double > dSourceArea( m_meshInputCov->faces.size() );
                DataArray1D< double > dTargetArea( m_meshOutput->faces.size() );

                assert( m_meshOverlap->vecSourceFaceIx.size() == m_meshOverlap->faces.size() );
                assert( m_meshOverlap->vecTargetFaceIx.size() == m_meshOverlap->faces.size() );
                assert( m_meshOverlap->vecFaceArea.GetRows() == m_meshOverlap->faces.size() );

                assert( m_meshInputCov->vecFaceArea.GetRows() == m_meshInputCov->faces.size() );
                assert( m_meshOutput->vecFaceArea.GetRows() == m_meshOutput->faces.size() );

                for( size_t i = 0; i < m_meshOverlap->faces.size(); i++ )
                {
                    if( m_meshOverlap->vecSourceFaceIx[i] < 0 || m_meshOverlap->vecTargetFaceIx[i] < 0 )
                        continue;  // skip this cell since it is ghosted

                    // let us recompute the source/target areas based on overlap mesh areas
                    assert( static_cast< size_t >( m_meshOverlap->vecSourceFaceIx[i] ) < m_meshInputCov->faces.size() );
                    dSourceArea[m_meshOverlap->vecSourceFaceIx[i]] += m_meshOverlap->vecFaceArea[i];
                    assert( static_cast< size_t >( m_meshOverlap->vecTargetFaceIx[i] ) < m_meshOutput->faces.size() );
                    dTargetArea[m_meshOverlap->vecTargetFaceIx[i]] += m_meshOverlap->vecFaceArea[i];
                }

                for( size_t i = 0; i < m_meshInputCov->faces.size(); i++ )
                {
                    if( fabs( dSourceArea[i] - m_meshInputCov->vecFaceArea[i] ) < 1.0e-10 )
                    {
                        m_meshInputCov->vecFaceArea[i] = dSourceArea[i];
                    }
                }
                for( size_t i = 0; i < m_meshOutput->faces.size(); i++ )
                {
                    if( fabs( dTargetArea[i] - m_meshOutput->vecFaceArea[i] ) < 1.0e-10 )
                    {
                        m_meshOutput->vecFaceArea[i] = dTargetArea[i];
                    }
                }
            }

            // Set source mesh areas in map
            if( !m_bPointCloudSource && eInputType == DiscretizationType_FV )
            {
                this->SetSourceAreas( m_meshInputCov->vecFaceArea );
                if( m_meshInputCov->vecMask.size() )
                {
                    this->SetSourceMask( m_meshInputCov->vecMask );
                }
            }

            // Set target mesh areas in map
            if( !m_bPointCloudTarget && eOutputType == DiscretizationType_FV )
            {
                this->SetTargetAreas( m_meshOutput->vecFaceArea );
                if( m_meshOutput->vecMask.size() )
                {
                    this->SetTargetMask( m_meshOutput->vecMask );
                }
            }

            /*
                // Recalculate input mesh area from overlap mesh
                if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
                    dbgprint.printf(0, "Overlap mesh only covers a sub-area of the sphere\n");
                    dbgprint.printf(0, "Recalculating source mesh areas\n");
                    dTotalAreaInput = m_meshInput->CalculateFaceAreasFromOverlap(m_meshOverlap);
                    dbgprint.printf(0, "New Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput);
                }
            */
        }

        // Finite volume input / Finite volume output
        if( ( eInputType == DiscretizationType_FV ) && ( eOutputType == DiscretizationType_FV ) )
        {
            // Generate reverse node array and edge map
            if( m_meshInputCov->revnodearray.size() == 0 ) m_meshInputCov->ConstructReverseNodeArray();
            if( m_meshInputCov->edgemap.size() == 0 ) m_meshInputCov->ConstructEdgeMap( false );

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFV( *m_meshInputCov );
            this->InitializeTargetCoordinatesFromMeshFV( *m_meshOutput );

            this->m_pdataGLLNodesIn  = nullptr;
            this->m_pdataGLLNodesOut = nullptr;

            // Finite volume input / Finite element output
            MB_CHK_ERR( this->SetDOFmapAssociation( eInputType, mapOptions.nPin, false, nullptr, nullptr, eOutputType,
                                                    mapOptions.nPout, false, nullptr ) );

            // Construct remap for FV-FV
            if( is_root ) dbgprint.printf( 0, "Calculating remap weights\n" );

            // Construct OfflineMap
            if( strMapAlgorithm == "invdist" )
            {
                if( m_meshInputCov->faces.size() )
                {
                    if( is_root ) dbgprint.printf( 0, "Calculating map (invdist)\n" );
                    LinearRemapFVtoFVInvDist( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
                }
            }
            else if( strMapAlgorithm == "delaunay" ) // does not need intersection mesh
            {
                if( m_meshInputCov->faces.size() )
                {
                    if( is_root ) dbgprint.printf( 0, "Calculating map (delaunay)\n" );
                    if (m_meshOverlap) LinearRemapTriangulation( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
                    else
                    {
                        Mesh dummy;
                        LinearRemapTriangulation( *m_meshInputCov, *m_meshOutput, dummy, *this );
                    }
                }
            }
            else if( strMapAlgorithm == "fvintbilin" )
            {
                if( m_meshInputCov->faces.size() )
                {
                    if( is_root ) dbgprint.printf( 0, "Calculating map (intbilin)\n" );
                    LinearRemapIntegratedBilinear( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
                }
            }
            else if( strMapAlgorithm == "fvintbilingb" )
            {
                if( m_meshInputCov->faces.size() )
                {
                    if( is_root ) dbgprint.printf( 0, "Calculating map (intbilingb)\n" );
                    LinearRemapIntegratedGeneralizedBarycentric( *m_meshInputCov, *m_meshOutput, *m_meshOverlap,
                                                                 *this );
                }
            }
            else if( strMapAlgorithm == "fvbilin" ) // does not need intersection mesh
            {
#ifdef VERBOSE
                if( is_root )
                {
                    m_meshInputCov->Write( "SourceMeshMBTR.g" );
                    m_meshOutput->Write( "TargetMeshMBTR.g" );
                }
                else
                {
                    m_meshInputCov->Write( "SourceMeshMBTR" + std::to_string( rank ) + ".g" );
                    m_meshOutput->Write( "TargetMeshMBTR" + std::to_string( rank ) + ".g" );
                }
#endif

                if( m_meshInputCov->faces.size() )
                {
                    if( is_root ) dbgprint.printf( 0, "Calculating map (bilin)\n" );
                    if (m_meshOverlap) LinearRemapBilinear( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
                    else
                    {
                        Mesh dummy;
                        LinearRemapBilinear( *m_meshInputCov, *m_meshOutput, dummy, *this );
                    }
                }
            }
            else
            {
                if( is_root ) dbgprint.printf( 0, "Calculating conservative FV-FV map\n" );
                if( m_meshInputCov->faces.size() )
                {
#ifdef USE_NATIVE_TEMPESTREMAP_ROUTINES
                    LinearRemapFVtoFV( *m_meshInputCov, *m_meshOutput, *m_meshOverlap,
                                       ( mapOptions.fMonotone ) ? ( 1 ) : ( mapOptions.nPin ), *this );
#else
                    LinearRemapFVtoFV_Tempest_MOAB( ( mapOptions.fMonotone ? 1 : mapOptions.nPin ) );
#endif
                }
            }
        }
        else if( eInputType == DiscretizationType_FV )
        {
            DataArray3D< double > dataGLLJacobian;

            if( is_root ) dbgprint.printf( 0, "Generating output mesh meta data\n" );
            double dNumericalArea_loc = GenerateMetaData( *m_meshOutput, mapOptions.nPout, mapOptions.fNoBubble,
                                                          dataGLLNodesDest, dataGLLJacobian );

            double dNumericalArea = dNumericalArea_loc;
#ifdef MOAB_HAVE_MPI
            if( m_pcomm )
                MPI_Reduce( &dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
#endif
            if( is_root ) dbgprint.printf( 0, "Output Mesh Numerical Area: %1.15e\n", dNumericalArea );

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFV( *m_meshInputCov );
            this->InitializeTargetCoordinatesFromMeshFE( *m_meshOutput, mapOptions.nPout, dataGLLNodesDest );

            this->m_pdataGLLNodesIn  = nullptr;
            this->m_pdataGLLNodesOut = &dataGLLNodesDest;

            // Generate the continuous Jacobian
            bool fContinuous = ( eOutputType == DiscretizationType_CGLL );

            if( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian( dataGLLNodesDest, dataGLLJacobian, this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian( dataGLLJacobian, this->GetTargetAreas() );
            }

            // Generate reverse node array and edge map
            if( m_meshInputCov->revnodearray.size() == 0 ) m_meshInputCov->ConstructReverseNodeArray();
            if( m_meshInputCov->edgemap.size() == 0 ) m_meshInputCov->ConstructEdgeMap( false );

            // Finite volume input / Finite element output
            MB_CHK_ERR( this->SetDOFmapAssociation( eInputType, mapOptions.nPin, false, nullptr, nullptr, eOutputType,
                                                    mapOptions.nPout, ( eOutputType == DiscretizationType_CGLL ),
                                                    &dataGLLNodesDest ) );

            // Generate remap weights
            if( strMapAlgorithm == "volumetric" )
            {
                if( is_root ) dbgprint.printf( 0, "Calculating remapping weights for FV->GLL (volumetric)\n" );
                LinearRemapFVtoGLL_Volumetric( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, dataGLLNodesDest,
                                               dataGLLJacobian, this->GetTargetAreas(), mapOptions.nPin, *this,
                                               nMonotoneType, fContinuous, mapOptions.fNoConservation );
            }
            else
            {
                if( is_root ) dbgprint.printf( 0, "Calculating remapping weights for FV->GLL\n" );
                LinearRemapFVtoGLL( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, dataGLLNodesDest, dataGLLJacobian,
                                    this->GetTargetAreas(), mapOptions.nPin, *this, nMonotoneType, fContinuous,
                                    mapOptions.fNoConservation );
            }
        }
        else if( ( eInputType == DiscretizationType_PCLOUD ) || ( eOutputType == DiscretizationType_PCLOUD ) )
        {
            DataArray3D< double > dataGLLJacobian;
            if( !m_bPointCloudSource )
            {
                // Generate reverse node array and edge map
                if( m_meshInputCov->revnodearray.size() == 0 ) m_meshInputCov->ConstructReverseNodeArray();
                if( m_meshInputCov->edgemap.size() == 0 ) m_meshInputCov->ConstructEdgeMap( false );

                // Initialize coordinates for map
                if( eInputType == DiscretizationType_FV )
                {
                    this->InitializeSourceCoordinatesFromMeshFV( *m_meshInputCov );
                }
                else
                {
                    if( is_root ) dbgprint.printf( 0, "Generating input mesh meta data\n" );
                    DataArray3D< double > dataGLLJacobianSrc;
                    GenerateMetaData( *m_meshInputCov, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrcCov,
                                      dataGLLJacobian );
                    GenerateMetaData( *m_meshInput, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrc,
                                      dataGLLJacobianSrc );
                }
            }
            // else { /* Source is a point cloud dataset */ }

            if( !m_bPointCloudTarget )
            {
                // Generate reverse node array and edge map
                if( m_meshOutput->revnodearray.size() == 0 ) m_meshOutput->ConstructReverseNodeArray();
                if( m_meshOutput->edgemap.size() == 0 ) m_meshOutput->ConstructEdgeMap( false );

                // Initialize coordinates for map
                if( eOutputType == DiscretizationType_FV )
                {
                    this->InitializeSourceCoordinatesFromMeshFV( *m_meshOutput );
                }
                else
                {
                    if( is_root ) dbgprint.printf( 0, "Generating output mesh meta data\n" );
                    GenerateMetaData( *m_meshOutput, mapOptions.nPout, mapOptions.fNoBubble, dataGLLNodesDest,
                                      dataGLLJacobian );
                }
            }
            // else { /* Target is a point cloud dataset */ }

            // Finite volume input / Finite element output
            MB_CHK_ERR( this->SetDOFmapAssociation(
                eInputType, mapOptions.nPin, ( eInputType == DiscretizationType_CGLL ),
                ( m_bPointCloudSource || eInputType == DiscretizationType_FV ? nullptr : &dataGLLNodesSrcCov ),
                ( m_bPointCloudSource || eInputType == DiscretizationType_FV ? nullptr : &dataGLLNodesSrc ),
                eOutputType, mapOptions.nPout, ( eOutputType == DiscretizationType_CGLL ),
                ( m_bPointCloudTarget ? nullptr : &dataGLLNodesDest ) ) );

            // Construct remap
            if( is_root ) dbgprint.printf( 0, "Calculating remap weights with Nearest-Neighbor method\n" );
            MB_CHK_ERR( LinearRemapNN_MOAB( true /*use_GID_matching*/, false /*strict_check*/ ) );
        }
        else if( ( eInputType != DiscretizationType_FV ) && ( eOutputType == DiscretizationType_FV ) )
        {
            DataArray3D< double > dataGLLJacobianSrc, dataGLLJacobian;

            if( is_root ) dbgprint.printf( 0, "Generating input mesh meta data\n" );
            // generate metadata for the input meshes (both source and covering source)
            GenerateMetaData( *m_meshInput, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrc,
                              dataGLLJacobianSrc );
            GenerateMetaData( *m_meshInputCov, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrcCov,
                              dataGLLJacobian );

            if( dataGLLNodesSrcCov.GetSubColumns() != m_meshInputCov->faces.size() )
            {
                _EXCEPTIONT( "Number of element does not match between metadata and "
                             "input mesh" );
            }

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE( *m_meshInputCov, mapOptions.nPin, dataGLLNodesSrcCov );
            this->InitializeTargetCoordinatesFromMeshFV( *m_meshOutput );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian( dataGLLNodesSrcCov, dataGLLJacobian, this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian( dataGLLJacobian, this->GetSourceAreas() );
            }

            // Finite element input / Finite volume output
            MB_CHK_ERR( this->SetDOFmapAssociation( eInputType, mapOptions.nPin,
                                                    ( eInputType == DiscretizationType_CGLL ), &dataGLLNodesSrcCov,
                                                    &dataGLLNodesSrc, eOutputType, mapOptions.nPout, false, nullptr ) );

            // Generate remap
            if( is_root ) dbgprint.printf( 0, "Calculating remap weights\n" );

            if( strMapAlgorithm == "volumetric" )
            {
                _EXCEPTIONT( "Unimplemented: Volumetric currently unavailable for"
                             "GLL input mesh" );
            }

            this->m_pdataGLLNodesIn  = &dataGLLNodesSrcCov;
            this->m_pdataGLLNodesOut = nullptr;

#ifdef USE_NATIVE_TEMPESTREMAP_ROUTINES
            LinearRemapSE4( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, dataGLLNodesSrcCov, dataGLLJacobian,
                            nMonotoneType, fContinuousIn, mapOptions.fNoConservation, mapOptions.fSparseConstraints,
                            *this );
#else
            LinearRemapSE4_Tempest_MOAB( dataGLLNodesSrcCov, dataGLLJacobian, nMonotoneType, fContinuousIn,
                                         mapOptions.fNoConservation, mapOptions.fSparseConstraints );
#endif
        }
        else if( ( eInputType != DiscretizationType_FV ) && ( eOutputType != DiscretizationType_FV ) )
        {
            DataArray3D< double > dataGLLJacobianIn, dataGLLJacobianSrc;
            DataArray3D< double > dataGLLJacobianOut;

            // Input metadata
            if( is_root ) dbgprint.printf( 0, "Generating input mesh meta data\n" );
            // generate metadata for the input meshes (both source and covering source)
            GenerateMetaData( *m_meshInput, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrc,
                              dataGLLJacobianSrc );
            // now coverage
            GenerateMetaData( *m_meshInputCov, mapOptions.nPin, mapOptions.fNoBubble, dataGLLNodesSrcCov,
                              dataGLLJacobianIn );
            // Output metadata
            if( is_root ) dbgprint.printf( 0, "Generating output mesh meta data\n" );
            GenerateMetaData( *m_meshOutput, mapOptions.nPout, mapOptions.fNoBubble, dataGLLNodesDest,
                              dataGLLJacobianOut );

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE( *m_meshInputCov, mapOptions.nPin, dataGLLNodesSrcCov );
            this->InitializeTargetCoordinatesFromMeshFE( *m_meshOutput, mapOptions.nPout, dataGLLNodesDest );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian( dataGLLNodesSrcCov, dataGLLJacobianIn, this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian( dataGLLJacobianIn, this->GetSourceAreas() );
            }

            // Generate the continuous Jacobian for output mesh
            bool fContinuousOut = ( eOutputType == DiscretizationType_CGLL );

            if( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian( dataGLLNodesDest, dataGLLJacobianOut, this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian( dataGLLJacobianOut, this->GetTargetAreas() );
            }

            // Input Finite Element to Output Finite Element
            MB_CHK_ERR( this->SetDOFmapAssociation( eInputType, mapOptions.nPin,
                                                    ( eInputType == DiscretizationType_CGLL ), &dataGLLNodesSrcCov,
                                                    &dataGLLNodesSrc, eOutputType, mapOptions.nPout,
                                                    ( eOutputType == DiscretizationType_CGLL ), &dataGLLNodesDest ) );

            this->m_pdataGLLNodesIn  = &dataGLLNodesSrcCov;
            this->m_pdataGLLNodesOut = &dataGLLNodesDest;

            // Generate remap
            if( is_root ) dbgprint.printf( 0, "Calculating remap weights\n" );

#ifdef USE_NATIVE_TEMPESTREMAP_ROUTINES
            LinearRemapGLLtoGLL_Integrated( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, dataGLLNodesSrcCov,
                                            dataGLLJacobianIn, dataGLLNodesDest, dataGLLJacobianOut,
                                            this->GetTargetAreas(), mapOptions.nPin, mapOptions.nPout, nMonotoneType,
                                            fContinuousIn, fContinuousOut, mapOptions.fSparseConstraints, *this );
#else
            LinearRemapGLLtoGLL2_MOAB( dataGLLNodesSrcCov, dataGLLJacobianIn, dataGLLNodesDest, dataGLLJacobianOut,
                                       this->GetTargetAreas(), mapOptions.nPin, mapOptions.nPout, nMonotoneType,
                                       fContinuousIn, fContinuousOut, mapOptions.fNoConservation );
#endif
        }
        else
        {
            _EXCEPTIONT( "Not implemented" );
        }

#ifdef MOAB_HAVE_EIGEN3
        copy_tempest_sparsemat_to_eigen3();
#endif

#ifdef MOAB_HAVE_MPI
        if (m_meshOverlap)
        {
            // Remove ghosted entities from overlap set
            moab::Range ghostedEnts;
            MB_CHK_ERR( m_remapper->GetOverlapAugmentedEntities( ghostedEnts ) );
            moab::EntityHandle m_meshOverlapSet = m_remapper->GetMeshSet( moab::Remapper::OverlapMesh );
            MB_CHK_SET_ERR( m_interface->remove_entities( m_meshOverlapSet, ghostedEnts ),
                            "Deleting ghosted entities failed" );
        }
#endif
        // Verify consistency, conservation and monotonicity, globally
        if( !mapOptions.fNoCheck )
        {
            if( is_root ) dbgprint.printf( 0, "Verifying map" );
            this->IsConsistent( 1.0e-8 );
            if( !mapOptions.fNoConservation ) this->IsConservative( 1.0e-8 );

            if( nMonotoneType != 0 )
            {
                this->IsMonotone( 1.0e-12 );
            }
        }
    }
    catch( Exception& e )
    {
        dbgprint.printf( 0, "%s", e.ToString().c_str() );
        return ( moab::MB_FAILURE );
    }
    catch( ... )
    {
        return ( moab::MB_FAILURE );
    }
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsConsistent( double dTolerance )
{
#ifndef MOAB_HAVE_MPI

    return OfflineMap::IsConsistent( dTolerance );

#else

    // Get map entries
    DataArray1D< int > dataRows;
    DataArray1D< int > dataCols;
    DataArray1D< double > dataEntries;

    // Calculate row sums
    DataArray1D< double > dRowSums;
    m_mapRemap.GetEntries( dataRows, dataCols, dataEntries );
    dRowSums.Allocate( m_mapRemap.GetRows() );

    for( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        dRowSums[dataRows[i]] += dataEntries[i];
    }

    // Verify all row sums are equal to 1
    int fConsistent = 0;
    for( unsigned i = 0; i < dRowSums.GetRows(); i++ )
    {
        if( fabs( dRowSums[i] - 1.0 ) > dTolerance )
        {
            fConsistent++;
            int rowGID = row_gdofmap[i];
            Announce( "TempestOnlineMap is not consistent in row %i (%1.15e)", rowGID, dRowSums[i] );
        }
    }

    int ierr;
    int fConsistentGlobal = 0;
    ierr                  = MPI_Allreduce( &fConsistent, &fConsistentGlobal, 1, MPI_INT, MPI_SUM, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;

    return fConsistentGlobal;
#endif
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsConservative( double dTolerance )
{
#ifndef MOAB_HAVE_MPI

    return OfflineMap::IsConservative( dTolerance );

#else
    // return OfflineMap::IsConservative(dTolerance);

    int ierr;
    // Get map entries
    DataArray1D< int > dataRows;
    DataArray1D< int > dataCols;
    DataArray1D< double > dataEntries;
    const DataArray1D< double >& dTargetAreas = this->GetTargetAreas();
    const DataArray1D< double >& dSourceAreas = this->GetSourceAreas();

    // Calculate column sums
    std::vector< int > dColumnsUnique;
    std::vector< double > dColumnSums;

    int nColumns = m_mapRemap.GetColumns();
    m_mapRemap.GetEntries( dataRows, dataCols, dataEntries );
    dColumnSums.resize( m_nTotDofs_SrcCov, 0.0 );
    dColumnsUnique.resize( m_nTotDofs_SrcCov, -1 );

    for( unsigned i = 0; i < dataEntries.GetRows(); i++ )
    {
        dColumnSums[dataCols[i]] += dataEntries[i] * dTargetAreas[dataRows[i]] / dSourceAreas[dataCols[i]];

        assert( dataCols[i] < m_nTotDofs_SrcCov );

        // GID for column DoFs: col_gdofmap[ col_ldofmap [ dataCols[i] ] ]
        int colGID = this->GetColGlobalDoF( dataCols[i] );  // col_gdofmap[ col_ldofmap [ dataCols[i] ] ];
        // int colGID = col_gdofmap[ col_ldofmap [ dataCols[i] ] ];
        dColumnsUnique[dataCols[i]] = colGID;

        // std::cout << "Column dataCols[i]=" << dataCols[i] << " with GID = " << colGID <<
        // std::endl;
    }

    int rootProc = 0;
    std::vector< int > nElementsInProc;
    const int nDATA = 3;
    nElementsInProc.resize( size * nDATA );
    int senddata[nDATA] = { nColumns, m_nTotDofs_SrcCov, m_nTotDofs_Src };
    ierr = MPI_Gather( senddata, nDATA, MPI_INT, nElementsInProc.data(), nDATA, MPI_INT, rootProc, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;

    int nTotVals = 0, nTotColumns = 0; // nTotColumnsUnq = 0;
    std::vector< int > dColumnIndices;
    std::vector< double > dColumnSumsTotal;
    std::vector< int > displs, rcount;
    if( rank == rootProc )
    {
        displs.resize( size + 1, 0 );
        rcount.resize( size, 0 );
        int gsum = 0;
        for( int ir = 0; ir < size; ++ir )
        {
            nTotVals += nElementsInProc[ir * nDATA];
            nTotColumns += nElementsInProc[ir * nDATA + 1];
            // nTotColumnsUnq += nElementsInProc[ir * nDATA + 2];

            displs[ir] = gsum;
            rcount[ir] = nElementsInProc[ir * nDATA + 1];
            gsum += rcount[ir];

            // printf( "%d: nTotColumns: %d, Displs: %d, rcount: %d, gsum = %d\n", ir, nTotColumns, displs[ir], rcount[ir], gsum );
        }

        printf( "Total nnz: %d, global source elements = %d\n", nTotVals, gsum );

        dColumnIndices.resize( nTotColumns, -1 );
        dColumnSumsTotal.resize( nTotColumns, 0.0 );
        // dColumnSourceAreas.resize ( nTotColumns, 0.0 );
    }

    // Gather all ColumnSums to root process and accumulate
    // We expect that the sums of all columns equate to 1.0 within user specified tolerance
    // Need to do a gatherv here since different processes have different number of elements
    // MPI_Reduce(&dColumnSums[0], &dColumnSumsTotal[0], m_mapRemap.GetColumns(), MPI_DOUBLE,
    // MPI_SUM, 0, m_pcomm->comm());
    // Use .data() rather than &vec[0] -- on non-root ranks dColumnIndices /
    // dColumnSumsTotal are empty (only resized on root, see ~10 lines above),
    // and &vec[0] indexing into an empty vector is undefined behavior. The
    // .data() form returns nullptr for an empty vector, which MPI_Gatherv
    // ignores since recvcount on non-root paths is effectively zero.
    ierr = MPI_Gatherv( dColumnsUnique.data(), m_nTotDofs_SrcCov, MPI_INT, dColumnIndices.data(), rcount.data(),
                        displs.data(), MPI_INT, rootProc, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;
    ierr = MPI_Gatherv( dColumnSums.data(), m_nTotDofs_SrcCov, MPI_DOUBLE, dColumnSumsTotal.data(), rcount.data(),
                        displs.data(), MPI_DOUBLE, rootProc, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;
    // ierr = MPI_Gatherv ( &dSourceAreas[0], m_nTotDofs_SrcCov, MPI_DOUBLE, &dColumnSourceAreas[0],
    // rcount.data(), displs.data(), MPI_DOUBLE, rootProc, m_pcomm->comm() ); if ( ierr !=
    // MPI_SUCCESS ) return -1;

    // Clean out unwanted arrays now
    dColumnSums.clear();
    dColumnsUnique.clear();

    // Verify all column sums equal the input Jacobian
    int fConservative = 0;
    if( rank == rootProc )
    {
        displs[size] = ( nTotColumns );
        // std::vector<double> dColumnSumsOnRoot(nTotColumnsUnq, 0.0);
        std::map< int, double > dColumnSumsOnRoot;
        // std::map<int, double> dColumnSourceAreasOnRoot;
        for( int ir = 0; ir < size; ir++ )
        {
            for( int ips = displs[ir]; ips < displs[ir + 1]; ips++ )
            {
                if( dColumnIndices[ips] < 0 ) continue;
                // printf("%d, %d: dColumnIndices[ips]: %d\n", ir, ips, dColumnIndices[ips]);
                // assert( dColumnIndices[ips] < nTotColumnsUnq );
                dColumnSumsOnRoot[dColumnIndices[ips]] += dColumnSumsTotal[ips];  // / dColumnSourceAreas[ips];
                // dColumnSourceAreasOnRoot[ dColumnIndices[ips] ] = dColumnSourceAreas[ips];
                // dColumnSourceAreas[ dColumnIndices[ips] ]
            }
        }

        for( std::map< int, double >::iterator it = dColumnSumsOnRoot.begin(); it != dColumnSumsOnRoot.end(); ++it )
        {
            // if ( fabs ( it->second - dColumnSourceAreasOnRoot[it->first] ) > dTolerance )
            if( fabs( it->second - 1.0 ) > dTolerance )
            {
                fConservative++;
                Announce( "TempestOnlineMap is not conservative in column "
                          // "%i (%1.15e)", it->first, it->second );
                          "%i (%1.15e)",
                          it->first, it->second /* / dColumnSourceAreasOnRoot[it->first] */ );
            }
        }
    }

    // TODO: Just do a broadcast from root instead of a reduction
    ierr = MPI_Bcast( &fConservative, 1, MPI_INT, rootProc, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;

    return fConservative;
#endif
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsMonotone( double dTolerance )
{
#ifndef MOAB_HAVE_MPI

    return OfflineMap::IsMonotone( dTolerance );

#else

    // Get map entries
    DataArray1D< int > dataRows;
    DataArray1D< int > dataCols;
    DataArray1D< double > dataEntries;

    m_mapRemap.GetEntries( dataRows, dataCols, dataEntries );

    // Verify all entries are in the range [0,1]
    int fMonotone = 0;
    for( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        if( ( dataEntries[i] < -dTolerance ) || ( dataEntries[i] > 1.0 + dTolerance ) )
        {
            fMonotone++;

            Announce( "TempestOnlineMap is not monotone in entry (%i): %1.15e", i, dataEntries[i] );
        }
    }

    int ierr;
    int fMonotoneGlobal = 0;
    ierr                = MPI_Allreduce( &fMonotone, &fMonotoneGlobal, 1, MPI_INT, MPI_SUM, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;

    return fMonotoneGlobal;
#endif
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::ComputeAdjacencyRelations( std::vector< std::unordered_set< int > >& vecAdjFaces,
                                                        int nrings,
                                                        const Range& entities,
                                                        bool useMOABAdjacencies,
                                                        Mesh* trMesh )
{
    assert( nrings > 0 );
    assert( useMOABAdjacencies || trMesh != nullptr );

    const size_t nrows = vecAdjFaces.size();
    moab::MeshTopoUtil mtu( m_interface );
    for( size_t index = 0; index < nrows; index++ )
    {
        vecAdjFaces[index].insert( index );  // add self target face first
        {
            // Compute the adjacent faces to the target face
            if( useMOABAdjacencies )
            {
                moab::Range ents;
                // ents.insert( entities.index( entities[index] ) );
                ents.insert( entities[index] );
                moab::Range adjEnts;
                moab::ErrorCode rval = mtu.get_bridge_adjacencies( ents, 0, 2, adjEnts, nrings );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
                for( moab::Range::iterator it = adjEnts.begin(); it != adjEnts.end(); ++it )
                {
                    // int adjIndex = m_interface->id_from_handle(*it)-1;
                    int adjIndex = entities.index( *it );
                    // printf("rank: %d, Element %lu, entity: %lu, adjIndex %d\n", rank, index, *it, adjIndex);
                    if( adjIndex >= 0 ) vecAdjFaces[index].insert( adjIndex );
                }
            }
            else
            {
                ///  Vector storing adjacent Faces.
                typedef std::pair< int, int > FaceDistancePair;
                typedef std::vector< FaceDistancePair > AdjacentFaceVector;
                AdjacentFaceVector adjFaces;
                Face& face = trMesh->faces[index];
                GetAdjacentFaceVectorByEdge( *trMesh, index, nrings * face.edges.size(), adjFaces );

                // Add the adjacent faces to the target face list
                for( auto adjFace : adjFaces )
                    if( adjFace.first >= 0 )
                        vecAdjFaces[index].insert( adjFace.first );  // map target face to source face
            }
        }
    }
}

moab::ErrorCode moab::TempestOnlineMap::ApplyWeights( moab::Tag srcSolutionTag,
                                                      moab::Tag tgtSolutionTag,
                                                      bool transpose,
                                                      CAASType caasType,
                                                      double default_projection )
{
    std::vector< double > solSTagVals;
    std::vector< double > solTTagVals;

    moab::Range sents, tents;
    if( m_remapper->point_cloud_source || m_remapper->point_cloud_target )
    {
        if( m_remapper->point_cloud_source )
        {
            moab::Range& covSrcEnts = m_remapper->GetMeshVertices( moab::Remapper::CoveringMesh );
            solSTagVals.resize( covSrcEnts.size(), default_projection );
            sents = covSrcEnts;
        }
        else
        {
            moab::Range& covSrcEnts = m_remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
            solSTagVals.resize( covSrcEnts.size() * this->GetSourceNDofsPerElement() * this->GetSourceNDofsPerElement(),
                                default_projection );
            sents = covSrcEnts;
        }
        if( m_remapper->point_cloud_target )
        {
            moab::Range& tgtEnts = m_remapper->GetMeshVertices( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size(), default_projection );
            tents = tgtEnts;
        }
        else
        {
            moab::Range& tgtEnts = m_remapper->GetMeshEntities( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size() * this->GetDestinationNDofsPerElement() *
                                    this->GetDestinationNDofsPerElement(),
                                default_projection );
            tents = tgtEnts;
        }
    }
    else
    {
        moab::Range& covSrcEnts = m_remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
        moab::Range& tgtEnts    = m_remapper->GetMeshEntities( moab::Remapper::TargetMesh );
        solSTagVals.resize( covSrcEnts.size() * this->GetSourceNDofsPerElement() * this->GetSourceNDofsPerElement(),
                            default_projection );
        solTTagVals.resize( tgtEnts.size() * this->GetDestinationNDofsPerElement() *
                                this->GetDestinationNDofsPerElement(),
                            default_projection );

        sents = covSrcEnts;
        tents = tgtEnts;
    }

    // The tag data is np*np*n_el_src
    MB_CHK_SET_ERR( m_interface->tag_get_data( srcSolutionTag, sents, &solSTagVals[0] ),
                    "Getting local tag data failed" );

    // Compute the application of weights on the suorce solution data and store it in the
    // destination solution vector data Optionally, can also perform the transpose application of
    // the weight matrix. Set the 3rd argument to true if this is needed
    MB_CHK_SET_ERR( this->ApplyWeights( solSTagVals, solTTagVals, transpose ),
                    "Applying remap operator onto source vector data failed" );

    // The tag data is np*np*n_el_dest
    MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                    "Setting target tag data failed" );

    if( caasType != CAAS_NONE )
    {
        std::string tgtSolutionTagName;
        MB_CHK_SET_ERR( m_interface->tag_get_name( tgtSolutionTag, tgtSolutionTagName ), "Getting tag name failed" );

        // Perform CAAS iterations iteratively until convergence
        constexpr int nmax_caas_iterations = 10;
        double mismatch                    = 1.0;
        int caasIteration                  = 0;
        double initialMismatch             = 0.0;
        while( ( fabs( mismatch / initialMismatch ) > 1e-15 && fabs( mismatch ) > 1e-15 ) &&
               caasIteration++ < nmax_caas_iterations )  // iterate until convergence or a maximum of 5 iterations
        {
            double dMassDiffPostGlobal;
            std::pair< double, double > mDefect =
                this->ApplyBoundsLimiting( solSTagVals, solTTagVals, caasType, caasIteration, mismatch );
#ifdef MOAB_HAVE_MPI
            double dMassDiffPost = mDefect.second;
            MPI_Allreduce( &dMassDiffPost, &dMassDiffPostGlobal, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#else
            dMassDiffPostGlobal = mDefect.second;
#endif
            if( caasIteration == 1 ) initialMismatch = mDefect.first;
            if( m_remapper->verbose && is_root )
            {
                printf( "Field {%s} -> CAAS iteration: %d, mass defect: %3.4e, post-CAAS: %3.4e\n",
                        tgtSolutionTagName.c_str(), caasIteration, mDefect.first, dMassDiffPostGlobal );
            }
            mismatch = dMassDiffPostGlobal;

            // The tag data is np*np*n_el_dest
            MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                            "Setting local tag data failed" );
        }
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::ApplyWeightsWithDualMap( moab::Tag srcSolutionTag,
                                                                  moab::Tag tgtSolutionTag,
                                                                  TempestOnlineMap* loWeightMap,
                                                                  CAASType caasType )
{
    // Setup entity ranges (same pattern as ApplyWeights(Tag, Tag))
    std::vector< double > solSTagVals, solTTagVals;
    moab::Range sents, tents;

    if( m_remapper->point_cloud_source || m_remapper->point_cloud_target )
    {
        if( m_remapper->point_cloud_source )
        {
            moab::Range& covSrcEnts = m_remapper->GetMeshVertices( moab::Remapper::CoveringMesh );
            solSTagVals.resize( covSrcEnts.size(), 0.0 );
            sents = covSrcEnts;
        }
        else
        {
            moab::Range& covSrcEnts = m_remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
            solSTagVals.resize( covSrcEnts.size() * this->GetSourceNDofsPerElement() *
                                    this->GetSourceNDofsPerElement(),
                                0.0 );
            sents = covSrcEnts;
        }
        if( m_remapper->point_cloud_target )
        {
            moab::Range& tgtEnts = m_remapper->GetMeshVertices( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size(), 0.0 );
            tents = tgtEnts;
        }
        else
        {
            moab::Range& tgtEnts = m_remapper->GetMeshEntities( moab::Remapper::TargetMesh );
            solTTagVals.resize( tgtEnts.size() * this->GetDestinationNDofsPerElement() *
                                    this->GetDestinationNDofsPerElement(),
                                0.0 );
            tents = tgtEnts;
        }
    }
    else
    {
        moab::Range& covSrcEnts = m_remapper->GetMeshEntities( moab::Remapper::CoveringMesh );
        moab::Range& tgtEnts    = m_remapper->GetMeshEntities( moab::Remapper::TargetMesh );
        solSTagVals.resize( covSrcEnts.size() * this->GetSourceNDofsPerElement() * this->GetSourceNDofsPerElement(),
                            0.0 );
        solTTagVals.resize(
            tgtEnts.size() * this->GetDestinationNDofsPerElement() * this->GetDestinationNDofsPerElement(), 0.0 );
        sents = covSrcEnts;
        tents = tgtEnts;
    }

    // Read source tag data from coverage mesh
    MB_CHK_SET_ERR( m_interface->tag_get_data( srcSolutionTag, sents, &solSTagVals[0] ),
                    "Getting source tag data failed" );

    // Apply high-order projection only (no CAAS — bounds come from the low-order map)
    MB_CHK_SET_ERR( this->ApplyWeights( solSTagVals, solTTagVals, false ),
                    "High-order projection failed" );

    // Write initial projection to target tag
    MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                    "Setting target tag data failed" );

    if( caasType == CAAS_NONE || loWeightMap == nullptr ) return moab::MB_SUCCESS;

    // =====================================================================
    // Dual-map CAAS (Clip-And-Assured-Sum) — bit-for-bit port of MCT's
    // seq_nlmap_avNormArr (driver-mct/main/seq_nlmap_mod.F90).
    //
    // PURPOSE
    //   Conservative, bounds-preserving remap of a source field x onto a
    //   target mesh using TWO weight matrices: a high-order
    //   non-monotone map (A, = `this`) and a low-order monotone &
    //   conservative map (Am, = `loWeightMap`). The high-order map gives
    //   accuracy; the low-order map gives the conservation reference and
    //   the bounds-preservation safety net. This routine wires them
    //   together using the Clip-And-Assured-Sum scheme of
    //     Bradley, Bosler & Guba, "Conservation with bounded variation
    //     and limiters in semi-Lagrangian transport schemes",
    //     SIAM J. Sci. Comput. 41(5), 2019, doi:10.1137/18M1165414.
    //
    // NOTATION (matching the reference)
    //   x          source field values on the coverage mesh (solSTagVals)
    //   A          high-order map  (this->m_weightMatrix)
    //   Am         low-order map   (loWeightMap)
    //   y_hi       = A  * x        high-order projection (in solTTagVals)
    //   y_lo       = Am * x        low-order projection (mass reference)
    //   [lo, hi]   per-row source-value bounds taken over A's stencil
    //   gmins/gmaxs unscaled global min/max of the per-row [lo, hi] —
    //              used as a final safety clip
    //   norm8wt    fractional-coverage weight (one scalar per source cell)
    //              propagated by the E3SM driver as a side-channel tag;
    //              when present, all bounds & redistribution arithmetic is
    //              rescaled to match MCT's lnorm=.true. branch exactly
    //
    // ALGORITHM (one pass, FP-order-preserved vs MCT)
    //   1) y_hi  = A * x                            (done above, in solTTagVals)
    //   2) y_lo  = Am * x                                            [Step 2]
    //   2b) Pull source norm8wt side-channel tag if present          [Step 2b]
    //   3) Per-row bounds [lo, hi] over A's stencil columns          [Step 3]
    //      Divide source value by srcNorm8wt before tracking
    //      min/max so bounds are over RECOVERED x, not (frac*x).
    //   4) y_lo == 0 mask: where the low-order projection is zero,
    //      force y_hi = lo = hi = 0 to drop the cell.               [Step 4]
    //   4b) Snapshot UNSCALED global extrema gmins/gmaxs from the
    //       masked, but not-yet-norm-scaled, per-row [lo, hi].      [Step 4b]
    //   4c) mappedNorm8wt = Am * srcNorm8wt (or Am * 1 if absent).  [Step 4c]
    //   4d) Scale per-row bounds: lo *= mappedNorm8wt, hi *= ...    [Step 4d]
    //   5) Per-cell CAAS quantities (clipping defect, room to lower/raise) [Step 5]
    //   6) Reproducible global reductions of the per-cell quantities [Step 7]
    //   7) dM_total = dM_clip + (M_low - M_hi_unclipped)             [Step 8]
    //   8) Redistribute the deficit across cells with room.          [Step 9]
    //   9) Final hard clip to gmins/gmaxs (skipping yLow==0 cells).
    //
    // REPRODUCIBILITY MODEL
    //   "BfB with MCT" means: for the same inputs, this routine produces
    //   the same target values MCT produces, BIT FOR BIT, regardless of
    //   MPI rank count or mesh decomposition. This requires three things
    //   that the code below enforces explicitly:
    //
    //   (i)  Same area values. MCT uses 'aream' = area_b from the netcdf
    //        map file. We read the same MOAB 'aream' tag (loaded by
    //        iMOAB_LoadMapFile). Recomputing spherical-polygon areas via
    //        lHuiller from mesh geometry is FP-different and is only used
    //        as a final fallback for online-computed maps with no aream.
    //
    //   (ii) Same FP operation order in the per-cell arithmetic. The
    //        redistribute step computes `(hi - yc)/cap_g * dM_total`, NOT
    //        the algebraically-equivalent `(hi - yc) * (dM_total/cap_g)`.
    //        See Step 9 below for why. The dM_total computation also uses
    //        MCT's exact two-step form (subtract, then add), preserving
    //        the catastrophic-cancellation residual MCT carries.
    //
    //   (iii) Order-independent global reductions. We use Worley's
    //         IntegerReprosum (the MOAB port of shr_reprosum_int), which
    //         is MCT's default reprosum path. It is decomposition- and
    //         order-independent by construction (integer-vector MPI sum).
    //         A Kahan + sort-by-gid summation lambda is also defined
    //         below as a reference alternative but is not the active
    //         reducer — using two different algorithms would defeat BfB.
    //
    // GUARDRAILS
    //   Bounds extraction (Step 3) hard-aborts the run if any owned
    //   high-order row references a coverage column not present on this
    //   rank. Silently dropping such columns would produce
    //   decomposition-dependent bounds and break BfB. The error message
    //   tells the caller exactly which row/column/weight failed and
    //   recommends widening the ghost-layer count (nghlay_cov in the
    //   E3SM coupler driver). See lines below the bounds loop.
    //
    // EARLY RETURN
    //   If caasType == CAAS_NONE or loWeightMap is null, we keep the raw
    //   high-order projection that was already written to tgtSolutionTag
    //   above. The dual-map machinery only runs when the caller explicitly
    //   activates it with a non-null low-order map and a non-CAAS_NONE
    //   filter type.

    const size_t nTargetDofs = solTTagVals.size();
    const size_t nSourceDofs = solSTagVals.size();

    // Map from target tag index to matrix row index. Both A and Am must share
    // the same row layout (same target mesh, same partitioning); this is true
    // because both maps are loaded onto the same intersection application.
    if( row_dtoc_dofmap.size() < nTargetDofs )
    {
        MB_CHK_SET_ERR( moab::MB_FAILURE, "row_dtoc_dofmap smaller than target tag size" );
    }

    // ----- Step 2: low-order projection y_lo = Am * x ---------------------
    std::vector< double > yLow( nTargetDofs, 0.0 );
    MB_CHK_SET_ERR( loWeightMap->ApplyWeights( solSTagVals, yLow, false ),
                    "Low-order projection failed" );

    // ----- Step 2b: pull source norm8wt side-channel (if present) ---------
    //
    // BACKGROUND
    //   When the E3SM driver coupler asks for a normalized projection
    //   (lnorm=.true.), it does NOT send raw source values x to the
    //   remapper. Instead, in seq_map_avNormArr it pre-multiplies each
    //   source data field by a fractional-coverage weight `frac`, and
    //   sends the products (frac * x) over to the intersection app
    //   together with a parallel single-component tag named "norm8wt"
    //   that carries (frac) on each source coverage cell.
    //
    //   The driver later UN-DOES this pre-norm on the target side by
    //   dividing each mapped data field by the mapped norm8wt — so the
    //   final value on the target is (Am*(frac*x)) / (Am*frac). That
    //   per-cell weighted average is the conservative answer when source
    //   cells are only partially covered (e.g. land/ocean coastlines).
    //
    // WHY THE CAAS KERNEL NEEDS TO SEE norm8wt
    //   To match MCT's seq_nlmap_avNormArr bit-for-bit, two things have
    //   to happen INSIDE the CAAS kernel — neither can be done by the
    //   driver as a post-pass:
    //
    //     (a) Per-row bounds [lo, hi] must be the min/max of RECOVERED x
    //         over the high-order stencil, not the min/max of (frac*x).
    //         MCT does
    //             tmp = solSTagVals[srcIdx]
    //             tmp = tmp / xPrimeAV(natt+1, col)   ! divide by frac
    //         in sMat_avMult_and_calc_bounds before extending bounds, and
    //         skips columns where frac == 0 (the field can't say anything
    //         meaningful at a cell with no source coverage). Without this
    //         divide, bounds would be 0-suppressed in coastal regions and
    //         the CAAS clip would lose accuracy.
    //
    //     (b) The mapped-norm8wt scale factor used in Step 4d must be the
    //         LOW-ORDER projection of the ACTUAL source `frac`, not the
    //         low-order projection of constant-1. MCT computes this in
    //         the same mct_sMat_avMult call that produces avp_o data — the
    //         natt+1 column gets sum_l w_lo[j,l] * frac(l), and that's
    //         what the bounds get scaled by.
    //
    // FALLBACK
    //   If no "norm8wt" tag exists on the intersection-side mesh (callers
    //   that never pre-normed), we set hasNorm8wt=false. In that branch
    //   srcNorm8wt is treated as constant-1 for both (a) and (b), which is
    //   mathematically correct: with no pre-norm, frac would have been
    //   1.0 everywhere and the divide / scale are no-ops.
    //
    // SHAPE CONSTRAINT
    //   norm8wt is single-component (one double per source coverage cell).
    //   For the FV-FV configuration on the active CAAS path,
    //   sents.size() == nSourceDofs and the tag values map directly to
    //   solSTagVals indices. If a future caller wires an SE source layout
    //   where nSourceDofs > sents.size() (multi-DOF per cell), the
    //   per-cell norm8wt cannot be unambiguously expanded to per-DOF
    //   values here — we deliberately fall back to the constant-1 path
    //   rather than guess an expansion that would silently break BfB.
    std::vector< double > srcNorm8wt;
    bool hasNorm8wt = false;
    {
        moab::Tag normTag = nullptr;
        moab::ErrorCode rvalN = m_interface->tag_get_handle( "norm8wt", normTag );
        if( MB_SUCCESS == rvalN && normTag != nullptr )
        {
            // Single-component tag (one double per source coverage entity).
            // For FV-FV (the only configuration on the active CAAS path)
            // sents.size() == nSourceDofs. If the source layout is multi-DOF
            // (e.g. SE) the per-cell norm8wt cannot be unambiguously expanded
            // to per-DOF values here; bail to the constant-1 fallback rather
            // than guess.
            srcNorm8wt.resize( sents.size(), 0.0 );
            moab::ErrorCode rvalD = m_interface->tag_get_data( normTag, sents, &srcNorm8wt[0] );
            if( MB_SUCCESS == rvalD && srcNorm8wt.size() == nSourceDofs )
            {
                hasNorm8wt = true;
            }
            else
            {
                srcNorm8wt.clear();
            }
        }
    }

    // ----- Step 3: per-row bounds from HIGH-ORDER stencil -----------------
    // bounds(A, x): for each target row r, [lo, hi] = [min, max] of x over
    // A(r,:)'s nonzero columns. The proof in the reference paper requires
    // bounds to come from the larger (high-order) stencil so the constraint
    // set is provably nonempty.
    std::vector< double > lcl_lo( nTargetDofs, 1e308 );
    std::vector< double > lcl_hi( nTargetDofs, -1e308 );

    WeightMatrix& hiW = this->m_weightMatrix;
    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        int r = row_dtoc_dofmap[i];
        if( r < 0 || r >= hiW.outerSize() ) continue;
        for( WeightMatrix::InnerIterator it( hiW, r ); it; ++it )
        {
            // it.col() is a matrix column index; map it to source vector index
            // by inverting col_dtoc_dofmap. For FV-FV with cell-based DOFs
            // and one-to-one mapping, this is the identity for owned columns.
            int mc = (int)it.col();
            // Search col_dtoc_dofmap[k]==mc; for typical FV cases the mapping
            // is dense and contiguous, so a linear scan over solSTagVals is
            // avoided by precomputing an inverse (below). For correctness we
            // fall back to scanning if the inverse is not available.
            // Build inverse once outside the loop (see below).
            (void)mc;
        }
    }

    // Precompute matrix-col -> source-vector-index inverse (cached per call;
    // O(nSourceDofs) construction). col_dtoc_dofmap has size nSourceDofs and
    // maps source-vector-index -> matrix-col.
    int maxMatCol = -1;
    for( size_t k = 0; k < nSourceDofs && k < col_dtoc_dofmap.size(); k++ )
        if( col_dtoc_dofmap[k] > maxMatCol ) maxMatCol = col_dtoc_dofmap[k];
    std::vector< int > col_inv( maxMatCol + 1, -1 );
    for( size_t k = 0; k < nSourceDofs && k < col_dtoc_dofmap.size(); k++ )
        if( col_dtoc_dofmap[k] >= 0 ) col_inv[col_dtoc_dofmap[k]] = (int)k;

    // Track whether any owned row's high-order stencil column failed to
    // resolve into the local coverage source vector. If that happens the
    // [lcl_lo, lcl_hi] bounds are computed over an INCOMPLETE stencil and
    // the CAAS clip + redistribute will produce decomposition-dependent
    // values — exactly the symptom seen as 1-2 ULP cross-rank-count drift
    // on file-loaded maps. Fail loudly with the offending coordinates so
    // the coverage layout (nghlay_cov in the calling code) can be widened.
    int    bndsLocalErr   = 0;
    int    bndsFirstRowG  = -1;   // global target row id where the first failure happened
    int    bndsFirstMc    = -1;   // matrix col index that failed to resolve
    double bndsFirstWgt   = 0.0;  // the dropped (nonzero) weight value
    int    bndsFirstKind  = 0;    // 1 = mc out of maxMatCol; 2 = col_inv -> -1; 3 = srcIdx OOB

    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        int r = row_dtoc_dofmap[i];
        if( r < 0 || r >= hiW.outerSize() ) continue;
        for( WeightMatrix::InnerIterator it( hiW, r ); it; ++it )
        {
            // Skip explicit-zero entries. Eigen's InnerIterator visits any
            // stored coefficient regardless of value; TempestRemap offline
            // maps routinely emit explicit zeros. MCT's
            // sMat_avMult_and_calc_bounds explicitly does
            //     if (wgt == 0) cycle
            // before extending bounds (seq_nlmap_mod.F90:855). We can't use
            // an exact-zero compare here (FP-fragile), but 1e-50 is below
            // any physically meaningful map weight while still robust to
            // sign and denormal noise — entries this small can't shift the
            // per-row [lo, hi] enough to cross a clip threshold either.
            if( fabs( it.value() ) < 1e-50 ) continue;
            const int mc = (int)it.col();

            // Hard checks: a nonzero high-order weight at column mc means
            // this owned row genuinely depends on source-coverage column mc.
            // If we cannot resolve mc to a local source-vector index, the
            // 3-ring coverage on this rank is too narrow for the high-order
            // stencil. Either the caller asked for too few ghost layers,
            // or the map file references columns not present in any rank's
            // coverage (a catastrophic mismatch). Either way, silently
            // skipping corrupts the bounds and breaks BFB.
            if( mc < 0 || mc > maxMatCol )
            {
                if( !bndsLocalErr )
                {
                    bndsLocalErr  = 1;
                    bndsFirstRowG = (r >= 0 && r < (int)row_gdofmap.size()) ? (int)row_gdofmap[r] : -1;
                    bndsFirstMc   = mc;
                    bndsFirstWgt  = it.value();
                    bndsFirstKind = 1;
                }
                continue;
            }
            const int srcIdx = col_inv[mc];
            if( srcIdx < 0 )
            {
                if( !bndsLocalErr )
                {
                    bndsLocalErr  = 1;
                    bndsFirstRowG = (r >= 0 && r < (int)row_gdofmap.size()) ? (int)row_gdofmap[r] : -1;
                    bndsFirstMc   = mc;
                    bndsFirstWgt  = it.value();
                    bndsFirstKind = 2;
                }
                continue;
            }
            if( srcIdx >= (int)nSourceDofs )
            {
                if( !bndsLocalErr )
                {
                    bndsLocalErr  = 1;
                    bndsFirstRowG = (r >= 0 && r < (int)row_gdofmap.size()) ? (int)row_gdofmap[r] : -1;
                    bndsFirstMc   = mc;
                    bndsFirstWgt  = it.value();
                    bndsFirstKind = 3;
                }
                continue;
            }
            // solSTagVals[srcIdx] holds (frac * x) when the driver pre-normed
            // (hasNorm8wt true); divide by frac to recover x for bounds, and
            // skip the source cell when frac == 0 (matches MCT
            // seq_nlmap_mod.F90:857 "if xPrimeAV(natt+1,col) == 0 cycle").
            // When hasNorm8wt is false, solSTagVals already holds raw x.
            double v = solSTagVals[srcIdx];
            if( hasNorm8wt )
            {
                const double n = srcNorm8wt[srcIdx];
                if( fabs(n) < 1E-20 ) continue;
                v /= n;
            }
            if( v < lcl_lo[i] ) lcl_lo[i] = v;
            if( v > lcl_hi[i] ) lcl_hi[i] = v;
        }
        // If row had no nonzero columns in the high-order stencil, set bounds
        // to 0 (matching MCT's sMat_avMult_and_calc_bounds: "lop = 0; hip = 0").
        // Together with the y_lo == 0 masking step below, this forces the cell
        // to 0 — a "rare, local reduction in order to one" per the reference.
        if( lcl_lo[i] > lcl_hi[i] )
        {
            lcl_lo[i] = 0.0;
            lcl_hi[i] = 0.0;
        }
    }

    // Globalize: any rank with bndsLocalErr triggers a collective failure.
#ifdef MOAB_HAVE_MPI
    {
        MPI_Comm comm = m_pcomm ? m_pcomm->comm() : MPI_COMM_SELF;
        int bndsGlobalErr = 0;
        MPI_Allreduce( &bndsLocalErr, &bndsGlobalErr, 1, MPI_INT, MPI_MAX, comm );
        if( bndsGlobalErr )
        {
            int myRank = 0;
            MPI_Comm_rank( comm, &myRank );
            if( bndsLocalErr )
            {
                static const char* kindStr[4] = { "?", "mc>maxMatCol", "col_inv[mc]==-1", "srcIdx>=nSourceDofs" };
                fprintf( stderr,
                         "FATAL: ApplyWeightsWithDualMap bounds extraction dropped a nonzero "
                         "high-order stencil column on rank %d.\n"
                         "       global_target_row=%d  matrix_col=%d  weight=%.17e  reason=%s\n"
                         "       This means the source coverage on this rank does NOT contain a "
                         "column the owned high-order row references — the 3-ring (or whatever) "
                         "ghost layer setting is too narrow, or the map file was generated against "
                         "a different mesh. Bounds computed over an incomplete stencil break BFB; "
                         "aborting rather than silently producing wrong CAAS output.\n",
                         myRank, bndsFirstRowG, bndsFirstMc, bndsFirstWgt, kindStr[bndsFirstKind] );
                fflush( stderr );
            }
            MPI_Abort( comm, 1 );
        }
    }
#else
    if( bndsLocalErr )
    {
        static const char* kindStr[4] = { "?", "mc>maxMatCol", "col_inv[mc]==-1", "srcIdx>=nSourceDofs" };
        fprintf( stderr,
                 "FATAL: ApplyWeightsWithDualMap bounds extraction dropped a nonzero "
                 "high-order stencil column.\n"
                 "       global_target_row=%d  matrix_col=%d  weight=%.17e  reason=%s\n",
                 bndsFirstRowG, bndsFirstMc, bndsFirstWgt, kindStr[bndsFirstKind] );
        fflush( stderr );
        return moab::MB_FAILURE;
    }
#endif

    // ----- Step 4: mask -- where y_lo == 0, zero out y_hi and bounds ------
    // (Per reference: "An exact 0 in the low-order field will mask the
    //  high-order field unnecessarily, but that's OK: it's a rare, local
    //  reduction in order to one, not a wrong value.")
    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        if( yLow[i] == 0.0 )
        {
            solTTagVals[i] = 0.0;
            lcl_lo[i]      = 0.0;
            lcl_hi[i]      = 0.0;
        }
    }

    // ----- Step 4b: compute UNSCALED global extrema for the final safety
    // clip (Item 4). MCT's seq_nlmap_avNormArr clips the redistributed result
    // against gmins/gmaxs = global min/max of the per-row (post-mask) bounds,
    // not against the per-row bounds themselves. Snapshot here, BEFORE the
    // bounds get scaled by the mapped norm in Step 4d below.
    double g_lo = 1e308, g_hi = -1e308;
    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        int r = row_dtoc_dofmap[i];
        if( r < 0 || r >= (int)row_gdofmap.size() ) continue;  // not owned
        if( lcl_lo[i] < g_lo ) g_lo = lcl_lo[i];
        if( lcl_hi[i] > g_hi ) g_hi = lcl_hi[i];
    }
#ifdef MOAB_HAVE_MPI
    {
        MPI_Comm comm = m_pcomm ? m_pcomm->comm() : MPI_COMM_SELF;
        double tmp_min = g_lo, tmp_max = g_hi;
        MPI_Allreduce( &tmp_min, &g_lo, 1, MPI_DOUBLE, MPI_MIN, comm );
        MPI_Allreduce( &tmp_max, &g_hi, 1, MPI_DOUBLE, MPI_MAX, comm );
    }
#endif

    // ----- Step 4c: compute mapped norm8wt = low-order map applied to a
    // source-norm8wt vector. For target row i this equals
    //   sum_l w_lo[i,l] * srcNorm8wt(l)
    // — equivalently, the value MCT carries in the natt+1 column of avp_o
    // after mct_sMat_avMult is applied to avp_i (whose norm8wt slot holds
    // frac post-pre-norm). When no "norm8wt" tag is available on the intx
    // side, fall back to applying the low-order map to a constant-1 vector
    // (equivalent to srcNorm8wt(l) == 1 everywhere — consistent with the
    // bounds-extraction fallback above).
    std::vector< double > mappedNorm8wt( nTargetDofs, 0.0 );
    if( hasNorm8wt )
    {
        MB_CHK_SET_ERR( loWeightMap->ApplyWeights( srcNorm8wt, mappedNorm8wt, false ),
                        "Mapped-norm8wt computation (low-order on source norm8wt) failed" );
    }
    else
    {
        std::vector< double > srcOnes( nSourceDofs, 1.0 );
        MB_CHK_SET_ERR( loWeightMap->ApplyWeights( srcOnes, mappedNorm8wt, false ),
                        "Mapped-norm8wt computation (low-order on ones) failed" );
    }

    // ----- Step 4d: scale per-row bounds by mapped norm8wt (Item 2).
    // MCT does this inside the CAAS loop:
    //     if (lnorm) then
    //        lo = lo*avp_o%rAttr(natt+1,j)
    //        hi = hi*avp_o%rAttr(natt+1,j)
    //     end if
    // Doing it once here propagates correctly into Step 5 (clipping) and
    // Step 9 (redistribution) which both use lcl_lo/lcl_hi. Note: g_lo/g_hi
    // were already snapshotted above and remain UNSCALED (matching MCT's
    // gmins/gmaxs which are the global min/max of the unscaled per-row bounds).
    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        const double w = mappedNorm8wt[i];
        lcl_lo[i] *= w;
        lcl_hi[i] *= w;
    }

    // ----- Get target areas (per matrix-row) ------------------------------
    // For BFB with MCT we MUST use the same area values MCT does. MCT uses
    // 'aream' (= area_b from the netcdf map file, loaded once when the map
    // is read). iMOAB_LoadMapFile populates the 'aream' tag on the target
    // mesh from area_b when arearead != 0 (e.g. arearead=3 for F-maps).
    //
    // The CAAS path MUST NOT recompute spherical-polygon areas from mesh
    // geometry via lHuiller (or any other re-derivation): doing so differs
    // from area_b at FP precision and silently breaks BfB with MCT. If the
    // caller has not loaded an area-bearing map and there are no online
    // areas (m_dTargetAreas) either, fail the run loudly so the caller can
    // fix their map-load configuration instead of getting silent non-BfB
    // results.
    std::vector< double > tgtAreas( nTargetDofs, 0.0 );
    {
        std::vector< moab::EntityHandle > tentVec;
        tentVec.reserve( tents.size() );
        for( moab::Range::iterator it = tents.begin(); it != tents.end(); ++it )
            tentVec.push_back( *it );

        bool got_areas = false;

        // Preferred: pull the 'aream' tag from the target MOAB mesh — this
        // is the area_b value loaded by iMOAB_LoadMapFile and is byte-identical
        // to the 'aream' field MCT uses in seq_nlmap_avNormArr.
        moab::Tag aream_tag = nullptr;
        moab::ErrorCode rval = m_interface->tag_get_handle( "aream", aream_tag );
        if( MB_SUCCESS == rval && aream_tag != nullptr && !tentVec.empty() )
        {
            const size_t nents = std::min< size_t >( tentVec.size(), nTargetDofs );
            std::vector< double > aream_vals( nents, 0.0 );
            rval = m_interface->tag_get_data( aream_tag, &tentVec[0], (int)nents, &aream_vals[0] );
            if( MB_SUCCESS == rval )
            {
                for( size_t i = 0; i < nents; i++ ) tgtAreas[i] = aream_vals[i];
                got_areas = true;
            }
        }

        // Fallback: areas were computed online and live in OfflineMap's
        // m_dTargetAreas (indexed by matrix row). This is BFB with MCT only
        // when the same online-area code path is used on both couplers; it
        // is acceptable for runs that build the map online.
        if( !got_areas )
        {
            const DataArray1D< double >& dTargetAreas = this->GetTargetAreas();
            const size_t nRows = dTargetAreas.GetRows();
            if( nRows >= nTargetDofs )
            {
                for( size_t i = 0; i < nTargetDofs; i++ )
                {
                    int r = row_dtoc_dofmap[i];
                    if( r >= 0 && (size_t)r < nRows )
                        tgtAreas[i] = dTargetAreas[r];
                }
                got_areas = true;
            }
        }

        // No fallback to lHuiller. Recomputing areas from mesh geometry is
        // not BFB with MCT and there is no safe silent default — abort.
        if( !got_areas )
        {
            MB_SET_ERR( moab::MB_FAILURE,
                        "ApplyWeightsWithDualMap: no target-cell areas available. "
                        "Neither the 'aream' tag (from iMOAB_LoadMapFile with "
                        "arearead != 0) nor OfflineMap::GetTargetAreas() (from an "
                        "online map build) provided areas. Recomputing areas from "
                        "mesh geometry is not bit-for-bit with MCT and is no longer "
                        "permitted in the CAAS path. Re-load the map file with an "
                        "area-bearing arearead setting (e.g. arearead=3 for F-maps), "
                        "or build the online map so target areas are populated." );
        }
    }

    // ----- Step 5: build per-cell CAAS weights ----------------------------
    // For BFB summation, we accumulate per-row (gid, value) pairs and reduce
    // them deterministically.
    std::vector< int >    rowGids( nTargetDofs, -1 );
    std::vector< double > massLowPerRow( nTargetDofs, 0.0 );
    std::vector< double > massHiUnclippedPerRow( nTargetDofs, 0.0 );  // y_hi BEFORE clipping
    std::vector< double > clipDefectPerRow( nTargetDofs, 0.0 );
    std::vector< double > capLowPerRow( nTargetDofs, 0.0 );
    std::vector< double > capHighPerRow( nTargetDofs, 0.0 );

    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        int r = row_dtoc_dofmap[i];
        if( r < 0 || r >= (int)row_gdofmap.size() )
            rowGids[i] = -1;  // not owned by this rank
        else
            rowGids[i] = (int)row_gdofmap[r];

        const double area = tgtAreas[i];
        const double y    = solTTagVals[i];        // y_hi (pre-clip)
        const double lo   = lcl_lo[i];
        const double hi   = lcl_hi[i];
        double yc         = y;                     // clipped value
        double dm         = 0.0;
        if( y < lo )
        {
            yc = lo;
            dm = ( y - lo ) * area;                // negative: cell exceeded below
        }
        else if( y > hi )
        {
            yc = hi;
            dm = ( y - hi ) * area;                // positive: cell exceeded above
        }
        clipDefectPerRow[i]      = dm;
        capLowPerRow[i]          = ( yc - lo ) * area;  // room to subtract
        capHighPerRow[i]         = ( hi - yc ) * area;  // room to add
        massLowPerRow[i]         = yLow[i] * area;
        // Per-row mass of the UNCLIPPED high-order projection. MCT reduces
        // exactly this quantity to obtain glbl_masses(natt+k) (M_hi_unclipped),
        // and then computes dM_total = dM_clip + (M_low - M_hi_unclipped) in
        // that 2-step order. We store the unclipped y here (rather than the
        // clipped yc as before) so MOAB's reprosum byte-matches MCT's, which
        // in turn lets the MCT-matching dM_total formula below produce the
        // same last bits as MCT.
        massHiUnclippedPerRow[i] = y * area;
        // Update solTTagVals to the clipped value for the next stage
        solTTagVals[i] = yc;
    }


    // ----- Step 6: BFB-deterministic global reductions --------------------

    // Reproducible global reductions via Worley's integer-vector algorithm
    // (moab::IntegerReprosum) — bit-identical to MCT's shr_reprosum_int
    // regardless of MPI rank count, mesh decomposition, or local iteration
    // order. This is MCT's default reprosum path (the namelist default
    // repro_sum_use_ddpdd=.false. on the MCT side). The integer-vector
    // algorithm is order-independent by construction (MPI_Allreduce with
    // MPI_SUM on int64), which eliminates the cross-rank-count ULP drift
    // that an order-sensitive reducer (e.g. Kahan or DDPDD) would otherwise
    // leak into the CAAS bounds and mass totals.
#ifdef MOAB_HAVE_MPI
    MPI_Comm reduce_comm = m_pcomm ? m_pcomm->comm() : MPI_COMM_SELF;
#else
    int reduce_comm = 0;  // serial build: comm unused but kept for API symmetry
#endif
#ifdef MOAB_HAVE_MPI
    moab::IntegerReprosum repro( reduce_comm );
#else
    moab::IntegerReprosum repro;
#endif
    // Build the ownership mask once (rowGids[i] >= 0 ↔ owned).
    const std::vector< int >& reduce_mask = rowGids;
    // Batched reduction: one MPI_Allreduce for the per-field metadata
    // (gmax_exp / gmin_exp / max_nsummands across the 5 fields) and one
    // MPI_Allreduce for the concatenated integer-vector encoding of all 5
    // fields. Bit-for-bit identical to calling sum_masked() five times
    // separately (each field still uses its own per-field metadata and
    // decode pass), but goes from 10 collective calls to 2.
    const std::vector< std::vector< double > > caasFields = {
        massLowPerRow, massHiUnclippedPerRow, clipDefectPerRow,
        capLowPerRow, capHighPerRow };
    std::vector< double > caasGsums;
    repro.sum_masked_batch( caasFields, reduce_mask, caasGsums );
    const double M_low          = caasGsums[0];
    const double M_hi_unclipped = caasGsums[1];
    const double dM_clip        = caasGsums[2];
    const double cap_low_g      = caasGsums[3];
    const double cap_high_g     = caasGsums[4];

    // ----- Step 8: total mass deficit between low-order and clipped high-order
    // The redistribution must drive the (clipped) high-order solution back to
    // the low-order mass. MCT's seq_nlmap_avNormArr (line 616) computes this
    // in EXACTLY the following 2-step form, and floating-point rounding makes
    // it FP-different from the algebraically-equivalent (M_low - M_hi_clip):
    //
    //     ! MCT (Fortran array assignment, evaluated element-wise)
    //     gwts(k) = gwts(k)          ! gwts(k) holds dM_clip after reprosum
    //             + (glbl_masses(k)  ! M_low
    //                - glbl_masses(natt+k))   ! M_hi_unclipped
    //
    // Reproducing MCT's exact bit pattern requires:
    //   (a) reducing the UNCLIPPED per-cell high-order mass directly via
    //       reprosum, NOT deriving it as M_hi_clip + dM_clip — that derivation
    //       drops 1-2 ULP because reprosum is exact only on its inputs.
    //       => see massHiUnclippedPerRow above.
    //   (b) computing dM_total in MCT's order: subtract first, then add.
    //
    // Sign: dM_total > 0 means low-order carries more mass than the clipped
    // high-order, so we need to ADD mass; dM_total < 0 means we need to REMOVE.
    //
    const double diff     = M_low - M_hi_unclipped;           // step 1: subtraction
    const double dM_total = dM_clip + diff;                   // step 2: addition (MCT order)

    // ----- Step 9: redistribute -------------------------------------------
    // For BfB with MCT seq_nlmap_avNormArr we MUST match its FP operation
    // order exactly. MCT does, per cell:
    //     y = max(lo, min(hi, nl_avp_o(k,j)))                 ! re-clip
    //     nl_avp_o(k,j) = y + ((hi - y)/tmp)*gwts(k)           ! line 648 / 662
    // i.e. divide-then-multiply, with the loop-invariant denominator
    // (cap_high_g or cap_low_g) and numerator (dM_total) NOT precomputed
    // into a single `scale = dM_total/cap_g`. Doing so introduces ULP-level
    // per-cell differences (a/b*c reorders to (c/b)*a). Similarly the
    // earlier MOAB pattern computed `room = (hi-yc)*area` and then
    // `room/area`, which doesn't algebraically cancel in FP and added two
    // extra roundings per cell. The straightforward `(hi - yc)/cap_g *
    // dM_total` form below matches MCT bit-for-bit.
    if( dM_total > 0.0 && cap_high_g > 0.0 )
    {
        for( size_t i = 0; i < nTargetDofs; i++ )
        {
            const double area = tgtAreas[i];
            const double yc   = solTTagVals[i];
            if( area > 0.0 )
                solTTagVals[i] = yc + ( ( lcl_hi[i] - yc ) / cap_high_g ) * dM_total;
        }
    }
    else if( dM_total < 0.0 && cap_low_g > 0.0 )
    {
        for( size_t i = 0; i < nTargetDofs; i++ )
        {
            const double area = tgtAreas[i];
            const double yc   = solTTagVals[i];
            if( area > 0.0 )
                solTTagVals[i] = yc + ( ( yc - lcl_lo[i] ) / cap_low_g ) * dM_total;
        }
    }

    // Final hard clip for floating-point safety, against UNSCALED global
    // extrema (Item 4). MCT's seq_nlmap_avNormArr does:
    //   if (avp_o(k,j) == 0) cycle             ! 0-mask skip
    //   nl_avp_o%rAttr(k,j) = max(gmins(k), min(gmaxs(k), nl_avp_o%rAttr(k,j)))
    // Per-row bounds (lcl_lo/lcl_hi) are now SCALED by mapped_norm8wt and so
    // would be a tighter clip than MCT applies; using global g_lo/g_hi keeps
    // the safety net loose, as the reference algorithm intends. The 0-mask
    // skip is critical: without it, target cells that were zeroed in Step 4
    // (yLow == 0) get bumped from 0 up to g_lo when g_lo > 0 (e.g.
    // positive-only fields like temperature/pressure), and the post-norm
    // divide in the driver then amplifies that wrong value by 1/wghts at
    // coastal coverage cells where wghts is tiny but nonzero.
    for( size_t i = 0; i < nTargetDofs; i++ )
    {
        if( fabs(yLow[i]) < 1E-40 ) continue;
        if( solTTagVals[i] < g_lo ) solTTagVals[i] = g_lo;
        if( solTTagVals[i] > g_hi ) solTTagVals[i] = g_hi;
    }

    // Store result back to the target tag
    MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                    "Setting target tag data failed" );

    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::DefineAnalyticalSolution( moab::Tag& solnTag,
                                                                  const std::string& solnName,
                                                                  moab::Remapper::IntersectionContext ctx,
                                                                  sample_function testFunction,
                                                                  moab::Tag* clonedSolnTag,
                                                                  std::string cloneSolnName )
{
    const bool outputEnabled = ( is_root );
    int discOrder;
    DiscretizationType discMethod;
    // moab::EntityHandle meshset;
    moab::Range entities;
    Mesh* trmesh;
    switch( ctx )
    {
        case Remapper::SourceMesh:
            // meshset    = m_remapper->m_covering_source_set;
            trmesh     = m_remapper->m_covering_source;
            entities   = ( m_remapper->point_cloud_source ? m_remapper->m_covering_source_vertices
                                                          : m_remapper->m_covering_source_entities );
            discOrder  = m_nDofsPEl_Src;
            discMethod = m_eInputType;
            break;

        case Remapper::TargetMesh:
            // meshset = m_remapper->m_target_set;
            trmesh = m_remapper->m_target;
            entities =
                ( m_remapper->point_cloud_target ? m_remapper->m_target_vertices : m_remapper->m_target_entities );
            discOrder  = m_nDofsPEl_Dest;
            discMethod = m_eOutputType;
            break;

        default:
            if( outputEnabled )
                std::cout << "Invalid context specified for defining an analytical solution tag" << std::endl;
            return moab::MB_FAILURE;
    }

    // Let us create teh solution tag with appropriate information for name, discretization order
    // (DoF space)
    MB_CHK_ERR( m_interface->tag_get_handle( solnName.c_str(), discOrder * discOrder, MB_TYPE_DOUBLE, solnTag,
                                             MB_TAG_DENSE | MB_TAG_CREAT ) );
    if( clonedSolnTag != nullptr )
    {
        if( cloneSolnName.size() == 0 )
        {
            cloneSolnName = solnName + std::string( "Cloned" );
        }
        MB_CHK_ERR( m_interface->tag_get_handle( cloneSolnName.c_str(), discOrder * discOrder, MB_TYPE_DOUBLE,
                                                 *clonedSolnTag, MB_TAG_DENSE | MB_TAG_CREAT ) );
    }

    // Triangular quadrature rule
    const int TriQuadratureOrder = 10;

    if( outputEnabled ) std::cout << "Using triangular quadrature of order " << TriQuadratureOrder << std::endl;

    TriangularQuadratureRule triquadrule( TriQuadratureOrder );

    const int TriQuadraturePoints = triquadrule.GetPoints();

    const DataArray2D< double >& TriQuadratureG = triquadrule.GetG();
    const DataArray1D< double >& TriQuadratureW = triquadrule.GetW();

    // Output data
    DataArray1D< double > dVar;
    DataArray1D< double > dVarMB;  // re-arranged local MOAB vector

    // Nodal geometric area
    DataArray1D< double > dNodeArea;

    // Calculate element areas
    // trmesh->CalculateFaceAreas(fContainsConcaveFaces);

    if( discMethod == DiscretizationType_CGLL || discMethod == DiscretizationType_DGLL )
    {
        /* Get the spectral points and sample the functionals accordingly */
        const bool fGLL          = true;
        const bool fGLLIntegrate = false;

        // Generate grid metadata
        DataArray3D< int > dataGLLNodes;
        DataArray3D< double > dataGLLJacobian;

        GenerateMetaData( *trmesh, discOrder, false, dataGLLNodes, dataGLLJacobian );

        // Number of elements
        int nElements = trmesh->faces.size();

        // Verify all elements are quadrilaterals
        for( int k = 0; k < nElements; k++ )
        {
            const Face& face = trmesh->faces[k];

            if( face.edges.size() != 4 )
            {
                _EXCEPTIONT( "Non-quadrilateral face detected; "
                             "incompatible with --gll" );
            }
        }

        // Number of unique nodes (CGLL) or total element-local DOFs (DGLL)
        const bool fDiscontinuous = ( discMethod == DiscretizationType_DGLL );
        int iMaxNode = 0;
        if( fDiscontinuous )
        {
            // DGLL: each element has independent DOFs
            iMaxNode = nElements * discOrder * discOrder;
        }
        else
        {
            // CGLL: shared nodes at element boundaries
            for( int i = 0; i < discOrder; i++ )
                for( int j = 0; j < discOrder; j++ )
                    for( int k = 0; k < nElements; k++ )
                        if( dataGLLNodes[i][j][k] > iMaxNode )
                            iMaxNode = dataGLLNodes[i][j][k];
        }

        // Get Gauss-Lobatto quadrature nodes
        DataArray1D< double > dG;
        DataArray1D< double > dW;

        GaussLobattoQuadrature::GetPoints( discOrder, 0.0, 1.0, dG, dW );

        // Get Gauss quadrature nodes
        const int nGaussP = 10;

        DataArray1D< double > dGaussG;
        DataArray1D< double > dGaussW;

        GaussQuadrature::GetPoints( nGaussP, 0.0, 1.0, dGaussG, dGaussW );

        // Allocate data
        dVar.Allocate( iMaxNode );
        dVarMB.Allocate( discOrder * discOrder * nElements );
        dNodeArea.Allocate( iMaxNode );

        // Sample data
        for( int k = 0; k < nElements; k++ )
        {
            const Face& face = trmesh->faces[k];

            // Sample data at GLL nodes
            if( fGLL )
            {
                for( int i = 0; i < discOrder; i++ )
                {
                    for( int j = 0; j < discOrder; j++ )
                    {

                        // Apply local map
                        Node node;
                        Node dDx1G;
                        Node dDx2G;

                        ApplyLocalMap( face, trmesh->nodes, dG[i], dG[j], node, dDx1G, dDx2G );

                        // Sample data at this point
                        double dNodeLon = atan2( node.y, node.x );
                        if( dNodeLon < 0.0 )
                        {
                            dNodeLon += 2.0 * M_PI;
                        }
                        double dNodeLat = asin( node.z );

                        double dSample = ( *testFunction )( dNodeLon, dNodeLat );

                        if( fDiscontinuous )
                            dVar[k * discOrder * discOrder + j * discOrder + i] = dSample;
                        else
                            dVar[dataGLLNodes[j][i][k] - 1] = dSample;
                    }
                }
                // High-order Gaussian integration over basis function
            }
            else
            {
                DataArray2D< double > dCoeff( discOrder, discOrder );

                for( int p = 0; p < nGaussP; p++ )
                {
                    for( int q = 0; q < nGaussP; q++ )
                    {

                        // Apply local map
                        Node node;
                        Node dDx1G;
                        Node dDx2G;

                        ApplyLocalMap( face, trmesh->nodes, dGaussG[p], dGaussG[q], node, dDx1G, dDx2G );

                        // Cross product gives local Jacobian
                        Node nodeCross = CrossProduct( dDx1G, dDx2G );

                        double dJacobian =
                            sqrt( nodeCross.x * nodeCross.x + nodeCross.y * nodeCross.y + nodeCross.z * nodeCross.z );

                        // Find components of quadrature point in basis
                        // of the first Face
                        SampleGLLFiniteElement( 0, discOrder, dGaussG[p], dGaussG[q], dCoeff );

                        // Sample data at this point
                        double dNodeLon = atan2( node.y, node.x );
                        if( dNodeLon < 0.0 )
                        {
                            dNodeLon += 2.0 * M_PI;
                        }
                        double dNodeLat = asin( node.z );

                        double dSample = ( *testFunction )( dNodeLon, dNodeLat );

                        // Integrate
                        for( int i = 0; i < discOrder; i++ )
                        {
                            for( int j = 0; j < discOrder; j++ )
                            {

                                double dNodalArea = dCoeff[i][j] * dGaussW[p] * dGaussW[q] * dJacobian;

                                dVar[dataGLLNodes[i][j][k] - 1] += dSample * dNodalArea;

                                dNodeArea[dataGLLNodes[i][j][k] - 1] += dNodalArea;
                            }
                        }
                    }
                }
            }
        }

        // Divide by area
        if( fGLLIntegrate )
        {
            for( size_t i = 0; i < dVar.GetRows(); i++ )
            {
                dVar[i] /= dNodeArea[i];
            }
        }

        // Let us rearrange the data based on DoF ID specification
        if( ctx == Remapper::SourceMesh )
        {
            for( unsigned j = 0; j < entities.size(); j++ )
                for( int p = 0; p < discOrder; p++ )
                    for( int q = 0; q < discOrder; q++ )
                    {
                        const int offsetDOF = j * discOrder * discOrder + p * discOrder + q;
                        dVarMB[offsetDOF]   = dVar[col_dtoc_dofmap[offsetDOF]];
                    }
        }
        else
        {
            for( unsigned j = 0; j < entities.size(); j++ )
                for( int p = 0; p < discOrder; p++ )
                    for( int q = 0; q < discOrder; q++ )
                    {
                        const int offsetDOF = j * discOrder * discOrder + p * discOrder + q;
                        dVarMB[offsetDOF]   = dVar[row_dtoc_dofmap[offsetDOF]];
                    }
        }

        // Set the tag data
        MB_CHK_ERR( m_interface->tag_set_data( solnTag, entities, &dVarMB[0] ) );
    }
    else
    {
        // assert( discOrder == 1 );
        if( discMethod == DiscretizationType_FV )
        {
            /* Compute an element-wise integral to store the sampled solution based on Quadrature
             * rules */
            // Resize the array
            dVar.Allocate( trmesh->faces.size() );

            std::vector< Node >& nodes = trmesh->nodes;

            // Loop through all Faces
            for( size_t i = 0; i < trmesh->faces.size(); i++ )
            {
                const Face& face = trmesh->faces[i];

                // Loop through all sub-triangles
                for( size_t j = 0; j < face.edges.size() - 2; j++ )
                {

                    const Node& node0 = nodes[face[0]];
                    const Node& node1 = nodes[face[j + 1]];
                    const Node& node2 = nodes[face[j + 2]];

                    // Triangle area
                    Face faceTri( 3 );
                    faceTri.SetNode( 0, face[0] );
                    faceTri.SetNode( 1, face[j + 1] );
                    faceTri.SetNode( 2, face[j + 2] );

                    double dTriangleArea = CalculateFaceArea( faceTri, nodes );

                    // Calculate the element average
                    double dTotalSample = 0.0;

                    // Loop through all quadrature points
                    for( int k = 0; k < TriQuadraturePoints; k++ )
                    {
                        Node node( TriQuadratureG[k][0] * node0.x + TriQuadratureG[k][1] * node1.x +
                                       TriQuadratureG[k][2] * node2.x,
                                   TriQuadratureG[k][0] * node0.y + TriQuadratureG[k][1] * node1.y +
                                       TriQuadratureG[k][2] * node2.y,
                                   TriQuadratureG[k][0] * node0.z + TriQuadratureG[k][1] * node1.z +
                                       TriQuadratureG[k][2] * node2.z );

                        double dMagnitude = node.Magnitude();
                        node.x /= dMagnitude;
                        node.y /= dMagnitude;
                        node.z /= dMagnitude;

                        double dLon = atan2( node.y, node.x );
                        if( dLon < 0.0 )
                        {
                            dLon += 2.0 * M_PI;
                        }
                        double dLat = asin( node.z );

                        double dSample = ( *testFunction )( dLon, dLat );

                        dTotalSample += dSample * TriQuadratureW[k] * dTriangleArea;
                    }

                    dVar[i] += dTotalSample / trmesh->vecFaceArea[i];
                }
            }
            MB_CHK_ERR( m_interface->tag_set_data( solnTag, entities, &dVar[0] ) );
        }
        else /* discMethod == DiscretizationType_PCLOUD */
        {
            /* Get the coordinates of the vertices and sample the functionals accordingly */
            std::vector< Node >& nodes = trmesh->nodes;

            // Resize the array
            dVar.Allocate( nodes.size() );

            for( size_t j = 0; j < nodes.size(); j++ )
            {
                Node& node        = nodes[j];
                double dMagnitude = node.Magnitude();
                node.x /= dMagnitude;
                node.y /= dMagnitude;
                node.z /= dMagnitude;
                double dLon = atan2( node.y, node.x );
                if( dLon < 0.0 )
                {
                    dLon += 2.0 * M_PI;
                }
                double dLat = asin( node.z );

                double dSample = ( *testFunction )( dLon, dLat );
                dVar[j]        = dSample;
            }

            MB_CHK_ERR( m_interface->tag_set_data( solnTag, entities, &dVar[0] ) );
        }
    }

    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::ComputeMetrics( moab::Remapper::IntersectionContext ctx,
                                                        moab::Tag& exactTag,
                                                        moab::Tag& approxTag,
                                                        std::map< std::string, double >& metrics,
                                                        bool verbose )
{
    const bool outputEnabled = ( is_root );
    int discOrder;
    // DiscretizationType discMethod;
    // moab::EntityHandle meshset;
    moab::Range entities;
    // Mesh* trmesh;
    switch( ctx )
    {
        case Remapper::SourceMesh:
            // meshset    = m_remapper->m_covering_source_set;
            // trmesh     = m_remapper->m_covering_source;
            entities  = ( m_remapper->point_cloud_source ? m_remapper->m_covering_source_vertices
                                                         : m_remapper->m_covering_source_entities );
            discOrder = m_nDofsPEl_Src;
            // discMethod = m_eInputType;
            break;

        case Remapper::TargetMesh:
            // meshset = m_remapper->m_target_set;
            // trmesh  = m_remapper->m_target;
            entities =
                ( m_remapper->point_cloud_target ? m_remapper->m_target_vertices : m_remapper->m_target_entities );
            discOrder = m_nDofsPEl_Dest;
            // discMethod = m_eOutputType;
            break;

        default:
            if( outputEnabled )
                std::cout << "Invalid context specified for defining an analytical solution tag" << std::endl;
            return moab::MB_FAILURE;
    }

    // Let us create teh solution tag with appropriate information for name, discretization order
    // (DoF space)
    std::string exactTagName, projTagName;
    const int ntotsize = entities.size() * discOrder * discOrder;
    std::vector< double > exactSolution( ntotsize, 0.0 ), projSolution( ntotsize, 0.0 );
    MB_CHK_ERR( m_interface->tag_get_name( exactTag, exactTagName ) );
    MB_CHK_ERR( m_interface->tag_get_data( exactTag, entities, &exactSolution[0] ) );
    MB_CHK_ERR( m_interface->tag_get_name( approxTag, projTagName ) );
    MB_CHK_ERR( m_interface->tag_get_data( approxTag, entities, &projSolution[0] ) );

    const auto& ovents = m_remapper->m_overlap_entities;

    std::vector< double > errnorms( 4, 0.0 ), globerrnorms( 4, 0.0 );  //  L1Err, L2Err, LinfErr
    double sumarea = 0.0;
    for( size_t i = 0; i < ovents.size(); ++i )
    {
        const int srcidx = m_remapper->m_overlap->vecSourceFaceIx[i];
        if( srcidx < 0 ) continue;  // Skip non-overlapping entities
        const int tgtidx = m_remapper->m_overlap->vecTargetFaceIx[i];
        if( tgtidx < 0 ) continue;  // skip ghost target faces
        const double ovarea = m_remapper->m_overlap->vecFaceArea[i];
        const double error  = fabs( exactSolution[tgtidx] - projSolution[tgtidx] );
        errnorms[0] += ovarea * error;
        errnorms[1] += ovarea * error * error;
        errnorms[3] = ( error > errnorms[3] ? error : errnorms[3] );
        sumarea += ovarea;
    }
    errnorms[2] = sumarea;
#ifdef MOAB_HAVE_MPI
    if( m_pcomm )
    {
        MPI_Reduce( &errnorms[0], &globerrnorms[0], 3, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
        MPI_Reduce( &errnorms[3], &globerrnorms[3], 1, MPI_DOUBLE, MPI_MAX, 0, m_pcomm->comm() );
    }
#else
    for( int i = 0; i < 4; ++i )
        globerrnorms[i] = errnorms[i];
#endif

    globerrnorms[0] = ( globerrnorms[0] / globerrnorms[2] );
    globerrnorms[1] = std::sqrt( globerrnorms[1] / globerrnorms[2] );

    metrics.clear();
    metrics["L1Error"]   = globerrnorms[0];
    metrics["L2Error"]   = globerrnorms[1];
    metrics["LinfError"] = globerrnorms[3];

    if( verbose && is_root )
    {
        std::cout << "Error metrics when comparing " << projTagName << " against " << exactTagName << std::endl;
        std::cout << "\t Total Intersection area = " << globerrnorms[2] << std::endl;
        std::cout << "\t L_1 error   = " << globerrnorms[0] << std::endl;
        std::cout << "\t L_2 error   = " << globerrnorms[1] << std::endl;
        std::cout << "\t L_inf error = " << globerrnorms[3] << std::endl;
    }

    return moab::MB_SUCCESS;
}
