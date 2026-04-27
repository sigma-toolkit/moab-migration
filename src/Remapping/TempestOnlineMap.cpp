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
        col_dtoc_dofmap.resize( m_remapper->m_covering_source_vertices.size(), UINT_MAX );
        src_soln_gdofs.resize( m_remapper->m_covering_source_vertices.size(), UINT_MAX );
        MB_CHK_ERR(
            m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_covering_source_vertices, &src_soln_gdofs[0] ) );
        srcTagSize = 1;
    }
    else
    {
        col_gdofmap.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX );
        col_dtoc_dofmap.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX );
        src_soln_gdofs.resize( m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX );
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
                    col_gdofmap[localDOF]      = src_soln_gdofs[offsetDOF] - 1;
                    col_dtoc_dofmap[offsetDOF] = localDOF;
                    if( vprint )
                        std::cout << "Col: " << offsetDOF << ", " << localDOF << ", " << col_gdofmap[offsetDOF] << ", "
                                  << m_nTotDofs_SrcCov << "\n";
                }
            }
        }
    }

    if( m_remapper->point_cloud_source )
    {
        assert( m_nDofsPEl_Src == 1 );
        srccol_gdofmap.resize( m_remapper->m_source_vertices.size(), UINT_MAX );
        srccol_dtoc_dofmap.resize( m_remapper->m_covering_source_vertices.size(), UINT_MAX );
        locsrc_soln_gdofs.resize( m_remapper->m_source_vertices.size(), UINT_MAX );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagSrc, m_remapper->m_source_vertices, &locsrc_soln_gdofs[0] ) );
    }
    else
    {
        srccol_gdofmap.resize( m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX );
        srccol_dtoc_dofmap.resize( m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX );
        locsrc_soln_gdofs.resize( m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX );
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
                    srccol_gdofmap[localDOF]      = locsrc_soln_gdofs[offsetDOF] - 1;
                    srccol_dtoc_dofmap[offsetDOF] = localDOF;
                }
            }
        }
    }

    int tgtTagSize = ( m_eOutputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Dest * m_nDofsPEl_Dest );
    if( m_remapper->point_cloud_target )
    {
        assert( m_nDofsPEl_Dest == 1 );
        row_gdofmap.resize( m_remapper->m_target_vertices.size(), UINT_MAX );
        row_dtoc_dofmap.resize( m_remapper->m_target_vertices.size(), UINT_MAX );
        tgt_soln_gdofs.resize( m_remapper->m_target_vertices.size(), UINT_MAX );
        MB_CHK_ERR( m_interface->tag_get_data( m_dofTagDest, m_remapper->m_target_vertices, &tgt_soln_gdofs[0] ) );
        tgtTagSize = 1;
    }
    else
    {
        row_gdofmap.resize( m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX );
        row_dtoc_dofmap.resize( m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX );
        tgt_soln_gdofs.resize( m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX );
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
                    row_gdofmap[localDOF]      = tgt_soln_gdofs[offsetDOF] - 1;
                    row_dtoc_dofmap[offsetDOF] = localDOF;
                    if( vprint )
                        std::cout << "Row: " << offsetDOF << ", " << localDOF << ", " << row_gdofmap[offsetDOF] << ", "
                                  << m_nTotDofs_Dest << "\n";
                }
            }
        }
    }

    // Let us also allocate the local representation of the sparse matrix
#if defined( MOAB_HAVE_EIGEN3 ) && defined( VERBOSE )
    if( vprint )
    {
        std::cout << "[" << rank << "]" << "DoFs: row = " << m_nTotDofs_Dest << ", " << row_gdofmap.size()
                  << ", col = " << m_nTotDofs_Src << ", " << m_nTotDofs_SrcCov << ", " << col_gdofmap.size() << "\n";
        // std::cout << "Max col_dofmap: " << maxcol << ", Min col_dofmap" << mincol << "\n";
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

        for( auto it : setMethodStrings )
        {
            // Piecewise constant monotonicity
            if( it == "mono2" )
            {
                if( nMonotoneType != 0 )
                {
                    _EXCEPTIONT( "Multiple monotonicity specifications found (--mono) or (--method \"mono#\")" );
                }
                if( ( m_eInputType == DiscretizationType_FV ) && ( m_eOutputType == DiscretizationType_FV ) )
                {
                    _EXCEPTIONT( "--method \"mono2\" is only used when remapping to/from CGLL or DGLL grids" );
                }
                nMonotoneType = 2;

                // Piecewise linear monotonicity
            }
            else if( it == "mono3" )
            {
                if( nMonotoneType != 0 )
                {
                    _EXCEPTIONT( "Multiple monotonicity specifications found (--mono) or (--method \"mono#\")" );
                }
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
            // Verify that overlap mesh is in the correct order (sanity check)
            assert( m_meshOverlap->vecSourceFaceIx.size() == m_meshOverlap->vecTargetFaceIx.size() );

            // Calculate Face areas
            if( is_root ) dbgprint.printf( 0, "Calculating overlap mesh Face areas\n" );
            local_areas[2] =
                m_meshOverlap->CalculateFaceAreas( mapOptions.fSourceConcave || mapOptions.fTargetConcave );

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
                dbgprint.printf( 0, "Overlap Mesh Recovered Area: %1.15e\n", global_areas[2] );
            }

            // Correct areas to match the areas calculated in the overlap mesh
            constexpr bool fCorrectAreas = true;
            if( fCorrectAreas )  // In MOAB-TempestRemap, we will always keep this to be true
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
                if( is_root ) dbgprint.printf( 0, "Calculating map (invdist)\n" );
                if( m_meshInputCov->faces.size() )
                    LinearRemapFVtoFVInvDist( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
            }
            else if( strMapAlgorithm == "delaunay" )
            {
                if( is_root ) dbgprint.printf( 0, "Calculating map (delaunay)\n" );
                if( m_meshInputCov->faces.size() )
                    LinearRemapTriangulation( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
            }
            else if( strMapAlgorithm == "fvintbilin" )
            {
                if( is_root ) dbgprint.printf( 0, "Calculating map (intbilin)\n" );
                if( m_meshInputCov->faces.size() )
                    LinearRemapIntegratedBilinear( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
            }
            else if( strMapAlgorithm == "fvintbilingb" )
            {
                if( is_root ) dbgprint.printf( 0, "Calculating map (intbilingb)\n" );
                if( m_meshInputCov->faces.size() )
                    LinearRemapIntegratedGeneralizedBarycentric( *m_meshInputCov, *m_meshOutput, *m_meshOverlap,
                                                                 *this );
            }
            else if( strMapAlgorithm == "fvbilin" )
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
                if( is_root ) dbgprint.printf( 0, "Calculating map (bilin)\n" );
                if( m_meshInputCov->faces.size() )
                    LinearRemapBilinear( *m_meshInputCov, *m_meshOutput, *m_meshOverlap, *this );
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
                                         mapOptions.fNoConservation );
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
    std::vector< double > dColumnSourceAreas;
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
    ierr = MPI_Gatherv( &dColumnsUnique[0], m_nTotDofs_SrcCov, MPI_INT, &dColumnIndices[0], rcount.data(),
                        displs.data(), MPI_INT, rootProc, m_pcomm->comm() );
    if( ierr != MPI_SUCCESS ) return -1;
    ierr = MPI_Gatherv( &dColumnSums[0], m_nTotDofs_SrcCov, MPI_DOUBLE, &dColumnSumsTotal[0], rcount.data(),
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

    // === Dual-map CAAS: compute bounds from low-order stencil, clip high-order result ===

    const size_t nTargetDofs = solTTagVals.size();
    const size_t nSourceDofs = solSTagVals.size();

    // Step 1: Compute absolute bounds from the low-order weight matrix stencil.
    // For each target row, find min/max of source field values over nonzero columns.
    WeightMatrix& loW = loWeightMap->GetWeightMatrix();
    std::vector< double > absLo( nTargetDofs, 1e308 );
    std::vector< double > absHi( nTargetDofs, -1e308 );

    for( size_t r = 0; r < nTargetDofs && r < (size_t)loW.outerSize(); r++ )
    {
        for( WeightMatrix::InnerIterator it( loW, r ); it; ++it )
        {
            int c = it.col();
            if( c >= 0 && c < (int)nSourceDofs )
            {
                absLo[r] = std::min( absLo[r], solSTagVals[c] );
                absHi[r] = std::max( absHi[r], solSTagVals[c] );
            }
        }
        // If the row had no nonzero entries, don't enforce bounds
        if( absLo[r] > absHi[r] )
        {
            absLo[r] = -1e308;
            absHi[r] = 1e308;
        }
    }

    // Step 2: Get target areas for conservative mass redistribution.
    // Try TempestRemap mesh face areas first; fall back to aream tag.
    std::vector< double > tgtAreas( nTargetDofs, 1.0 );

    if( m_meshOutput && (size_t)m_meshOutput->faces.size() >= nTargetDofs )
    {
        const DataArray1D< double >& faceAreas = m_meshOutput->vecFaceArea;
        for( size_t i = 0; i < nTargetDofs; i++ )
            tgtAreas[i] = faceAreas[i];
    }
    else
    {
        Tag areamTag;
        moab::ErrorCode tagRval = m_interface->tag_get_handle( "aream", areamTag );
        if( tagRval == moab::MB_SUCCESS )
            m_interface->tag_get_data( areamTag, tents, &tgtAreas[0] );
    }

    // Step 3: CAAS iteration loop — clip to bounds, redistribute mass defect conservatively
    std::string tgtSolnTagName;
    m_interface->tag_get_name( tgtSolutionTag, tgtSolnTagName );

    constexpr int nmax_caas_iterations = 10;
    constexpr double convergence_tol   = 1e-15;

    for( int iter = 0; iter < nmax_caas_iterations; iter++ )
    {
        // Clip target values to [absLo, absHi] and compute mass change from clipping
        double localMassBefore = 0.0, localMassAfter = 0.0;
        for( size_t i = 0; i < nTargetDofs; i++ )
        {
            localMassBefore += tgtAreas[i] * solTTagVals[i];
            solTTagVals[i] = std::max( absLo[i], std::min( absHi[i], solTTagVals[i] ) );
            localMassAfter += tgtAreas[i] * solTTagVals[i];
        }

        double localMassDefect  = localMassBefore - localMassAfter;
        double globalMassDefect = localMassDefect;

#ifdef MOAB_HAVE_MPI
        MPI_Allreduce( &localMassDefect, &globalMassDefect, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#endif

        if( m_remapper->verbose && is_root )
        {
            printf( "DualMap CAAS {%s}: iter %d, mass defect (global): %3.6e\n", tgtSolnTagName.c_str(), iter + 1,
                    globalMassDefect );
        }

        // Update target tag with clipped values
        MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                        "Setting target tag data failed" );

        if( fabs( globalMassDefect ) < convergence_tol ) break;

        // Redistribute mass defect proportionally within remaining room
        double localRoomUp = 0.0, localRoomDn = 0.0;
        for( size_t i = 0; i < nTargetDofs; i++ )
        {
            localRoomUp += tgtAreas[i] * ( absHi[i] - solTTagVals[i] );
            localRoomDn += tgtAreas[i] * ( solTTagVals[i] - absLo[i] );
        }

        double globalRoomUp = localRoomUp, globalRoomDn = localRoomDn;
#ifdef MOAB_HAVE_MPI
        double localRooms[2]  = { localRoomUp, localRoomDn };
        double globalRooms[2] = { 0.0, 0.0 };
        MPI_Allreduce( localRooms, globalRooms, 2, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
        globalRoomUp = globalRooms[0];
        globalRoomDn = globalRooms[1];
#endif

        double totalRoom = ( globalMassDefect > 0 ) ? globalRoomUp : globalRoomDn;

        if( fabs( totalRoom ) > 1e-20 )
        {
            for( size_t i = 0; i < nTargetDofs; i++ )
            {
                double room = ( globalMassDefect > 0 ) ? ( absHi[i] - solTTagVals[i] )
                                                       : ( solTTagVals[i] - absLo[i] );
                solTTagVals[i] += globalMassDefect * room / totalRoom;
            }

            // Update target tag after redistribution
            MB_CHK_SET_ERR( m_interface->tag_set_data( tgtSolutionTag, tents, &solTTagVals[0] ),
                            "Setting target tag data after redistribution failed" );
        }
    }

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

        // Number of unique nodes
        int iMaxNode = 0;
        for( int i = 0; i < discOrder; i++ )
        {
            for( int j = 0; j < discOrder; j++ )
            {
                for( int k = 0; k < nElements; k++ )
                {
                    if( dataGLLNodes[i][j][k] > iMaxNode )
                    {
                        iMaxNode = dataGLLNodes[i][j][k];
                    }
                }
            }
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
