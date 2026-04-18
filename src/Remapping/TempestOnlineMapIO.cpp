/*
 * =====================================================================================
 *
 *       Filename:  TempestOnlineMapIO.cpp
 *
 *    Description:  All I/O implementations related to TempestOnlineMap
 *
 *        Version:  1.0
 *        Created:  02/06/2021 02:35:41
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *        Company:  Argonne National Lab
 *
 * =====================================================================================
 */

#include "FiniteElementTools.h"
#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/TupleList.hpp"

#ifdef MOAB_HAVE_NETCDFPAR
#include "netcdfcpp_par.hpp"
#else
#include "netcdfcpp.h"
#endif

#ifdef MOAB_HAVE_PNETCDF
#include <pnetcdf.h>

#define ERR_PARNC( err )                                                              \
    if( err != NC_NOERR )                                                             \
    {                                                                                 \
        fprintf( stderr, "Error at line %d: %s\n", __LINE__, ncmpi_strerror( err ) ); \
        MPI_Abort( MPI_COMM_WORLD, 1 );                                               \
    }

#endif

#ifdef MOAB_HAVE_MPI

#define MPI_CHK_ERR( err )                                          \
    if( err )                                                       \
    {                                                               \
        std::cout << "MPI Failure. ErrorCode (" << ( err ) << ") "; \
        std::cout << "\nMPI Aborting... \n";                        \
        return moab::MB_FAILURE;                                    \
    }

#ifdef MOAB_HAVE_EIGEN3

// Function to serialize a sparse matrix to disk in text format
template < typename SparseMatrixType >
void moab::TempestOnlineMap::serializeSparseMatrix( const SparseMatrixType& mat, const std::string& filename )
{
    std::ofstream ofs( filename );
    if( !ofs.is_open() )
    {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write matrix dimensions and number of non-zeros
    int rows                             = mat.rows();
    int cols                             = mat.cols();
    typename SparseMatrixType::Index nnz = mat.nonZeros();
    ofs << rows << " " << cols << " " << nnz << "\n";

    // Iterate over non-zero elements
    for( int k = 0; k < mat.outerSize(); ++k )
    {
        for( typename SparseMatrixType::InnerIterator it( mat, k ); it; ++it )
        {
            // int row    = it.row();  // row index
            // int col    = it.col();  // col index (equals k)
            int row    = 1 + this->GetRowGlobalDoF( it.row() );  // row index
            int col    = 1 + this->GetColGlobalDoF( it.col() );  // col index
            auto value = it.value();
            ofs << row << " " << col << " " << value << "\n";
        }
    }
    ofs.close();
}

#endif

int moab::TempestOnlineMap::rearrange_arrays_by_dofs( const std::vector< unsigned int >& gdofmap,
                                                      DataArray1D< double >& vecFaceArea,
                                                      DataArray1D< double >& dCenterLon,
                                                      DataArray1D< double >& dCenterLat,
                                                      DataArray2D< double >& dVertexLon,
                                                      DataArray2D< double >& dVertexLat,
                                                      std::vector< int >& masks,
                                                      unsigned& N,  // will have the local, after
                                                      int nv,
                                                      int& maxdof )
{
    // first decide maxdof, for partitioning
    unsigned int localmax = 0;
    for( unsigned i = 0; i < N; i++ )
        if( gdofmap[i] > localmax ) localmax = gdofmap[i];

    // decide partitioning based on maxdof/size
    MPI_Allreduce( &localmax, &maxdof, 1, MPI_INT, MPI_MAX, m_pcomm->comm() );
    // maxdof is 0 based, so actual number is +1
    // maxdof
    int size_per_task = ( maxdof + 1 ) / size;  // based on this, processor to process dof x is x/size_per_task
    // so we decide to reorder by actual dof, such that task 0 has dofs from [0 to size_per_task), etc
    moab::TupleList tl;
    unsigned numr = 2 * nv + 3;         //  doubles: area, centerlon, center lat, nv (vertex lon, vertex lat)
    tl.initialize( 3, 0, 0, numr, N );  // to proc, dof, then
    tl.enableWriteAccess();

    // populate
    for( unsigned i = 0; i < N; i++ )
    {
        int gdof    = gdofmap[i];
        int to_proc = gdof / size_per_task;
        int mask    = (i >= masks.size() ? 1: masks[i]); // assume mask=1 if size is deficient (typically for SE-FV)
        if( to_proc >= size ) to_proc = size - 1;  // the last ones go to last proc
        int n                  = tl.get_n();
        tl.vi_wr[3 * n]        = to_proc;
        tl.vi_wr[3 * n + 1]    = gdof;
        tl.vi_wr[3 * n + 2]    = mask;
        tl.vr_wr[n * numr]     = vecFaceArea[i];
        tl.vr_wr[n * numr + 1] = dCenterLon[i];
        tl.vr_wr[n * numr + 2] = dCenterLat[i];
        for( int j = 0; j < nv; j++ )
        {
            tl.vr_wr[n * numr + 3 + j]      = dVertexLon[i][j];
            tl.vr_wr[n * numr + 3 + nv + j] = dVertexLat[i][j];
        }
        tl.inc_n();
    }

    // now do the heavy communication
    ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tl, 0 );

    // after communication, on each processor we should have tuples coming in
    // still need to order by global dofs; then rearrange input vectors
    moab::TupleList::buffer sort_buffer;
    sort_buffer.buffer_init( tl.get_n() );
    tl.sort( 1, &sort_buffer );
    // count how many are unique, and collapse
    int nb_unique = 1;
    for( unsigned i = 0; i < tl.get_n() - 1; i++ )
    {
        if( tl.vi_wr[3 * i + 1] != tl.vi_wr[3 * i + 4] ) nb_unique++;
    }
    vecFaceArea.Allocate( nb_unique );
    dCenterLon.Allocate( nb_unique );
    dCenterLat.Allocate( nb_unique );
    dVertexLon.Allocate( nb_unique, nv );
    dVertexLat.Allocate( nb_unique, nv );
    masks.resize( nb_unique );
    int current_size = 1;
    vecFaceArea[0]   = tl.vr_wr[0];
    dCenterLon[0]    = tl.vr_wr[1];
    dCenterLat[0]    = tl.vr_wr[2];
    masks[0]         = tl.vi_wr[2];
    for( int j = 0; j < nv; j++ )
    {
        dVertexLon[0][j] = tl.vr_wr[3 + j];
        dVertexLat[0][j] = tl.vr_wr[3 + nv + j];
    }
    for( unsigned i = 0; i < tl.get_n() - 1; i++ )
    {
        int i1 = i + 1;
        if( tl.vi_wr[3 * i + 1] != tl.vi_wr[3 * i + 4] )
        {
            vecFaceArea[current_size] = tl.vr_wr[i1 * numr];
            dCenterLon[current_size]  = tl.vr_wr[i1 * numr + 1];
            dCenterLat[current_size]  = tl.vr_wr[i1 * numr + 2];
            for( int j = 0; j < nv; j++ )
            {
                dVertexLon[current_size][j] = tl.vr_wr[i1 * numr + 3 + j];
                dVertexLat[current_size][j] = tl.vr_wr[i1 * numr + 3 + nv + j];
            }
            masks[current_size] = tl.vi_wr[3 * i1 + 2];
            current_size++;
        }
        else
        {
            vecFaceArea[current_size - 1] += tl.vr_wr[i1 * numr];  // accumulate areas; will come here only for cgll ?
        }
    }

    N = current_size;  // or nb_unique, should be the same
    return 0;
}
#endif

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::WriteParallelMap( const std::string& strFilename,
                                                          const std::map< std::string, std::string >& attrMap )
{
    size_t lastindex      = strFilename.find_last_of( "." );
    std::string extension = strFilename.substr( lastindex + 1, strFilename.size() );

    // Write the map file to disk in parallel
    if( extension == "nc" )
    {
#if !defined( MOAB_HAVE_NETCDFPAR )
        // Without parallel NetCDF, the SCRIP writer cannot handle multiple MPI ranks
        // writing to the same file. Fall back to the HDF5 format with a .h5m extension
        // and then the map can be converted to SCRIP format offline if needed.
        if( this->size > 1 )
        {
            std::string h5mFilename = strFilename.substr( 0, lastindex ) + ".h5m";
            if( !this->rank )
            {
                std::cout << "  [WriteParallelMap]: Parallel NetCDF not available; writing map to "
                          << "HDF5 format (" << h5mFilename << ") instead of SCRIP (.nc)\n";
            }
            MB_CHK_ERR( this->WriteHDF5MapFile( h5mFilename.c_str() ) );
            return moab::MB_SUCCESS;
        }
#endif
        /* Invoke the actual call to write the parallel map to disk in SCRIP format */
        MB_CHK_ERR( this->WriteSCRIPMapFile( strFilename.c_str(), attrMap ) );
    }
    else
    {
        /* Write to the parallel H5M format */
        MB_CHK_ERR( this->WriteHDF5MapFile( strFilename.c_str() ) );
    }

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::WriteSCRIPMapFile( const std::string& strFilename,
                                                           const std::map< std::string, std::string >& attrMap )
{
    NcError error( NcError::silent_nonfatal );

#ifdef MOAB_HAVE_NETCDFPAR
    bool is_independent = true;
    ParNcFile ncMap( m_pcomm->comm(), MPI_INFO_NULL, strFilename.c_str(), NcFile::Replace, NcFile::Netcdf4 );
    // ParNcFile ncMap( m_pcomm->comm(), MPI_INFO_NULL, strFilename.c_str(), NcmpiFile::replace, NcmpiFile::classic5 );
#else
    NcFile ncMap( strFilename.c_str(), NcFile::Replace );
#endif

    if( !ncMap.is_valid() )
    {
        _EXCEPTION1( "Unable to open output map file \"%s\"", strFilename.c_str() );
    }

    // Attributes
    // ncMap.add_att( "Title", "MOAB-TempestRemap Online Regridding Weight Generator" );
    auto it = attrMap.begin();
    while( it != attrMap.end() )
    {
        // set the map attributes
        ncMap.add_att( it->first.c_str(), it->second.c_str() );
        // increment iterator
        it++;
    }

    /**
     * Need to get the global maximum of number of vertices per element
     * Key issue is that when calling InitializeCoordinatesFromMeshFV, the allocation for
     *dVertexLon/dVertexLat are made based on the maximum vertices in the current process. However,
     *when writing this out, other processes may have a different size for the same array. This is
     *hence a mess to consolidate in h5mtoscrip eventually.
     **/

    /* Let us compute all relevant data for the current original source mesh on the process */
    DataArray1D< double > vecSourceFaceArea, vecTargetFaceArea;
    DataArray1D< double > dSourceCenterLon, dSourceCenterLat, dTargetCenterLon, dTargetCenterLat;
    DataArray2D< double > dSourceVertexLon, dSourceVertexLat, dTargetVertexLon, dTargetVertexLat;
    if( m_srcDiscType == DiscretizationType_FV || m_srcDiscType == DiscretizationType_PCLOUD )
    {
        this->InitializeCoordinatesFromMeshFV(
            *m_meshInput, dSourceCenterLon, dSourceCenterLat, dSourceVertexLon, dSourceVertexLat,
            ( this->m_remapper->m_source_type == moab::TempestRemapper::RLL ), /* fLatLon = false */
            m_remapper->max_source_edges );

        vecSourceFaceArea.Allocate( m_meshInput->vecFaceArea.GetRows() );
        for( unsigned i = 0; i < m_meshInput->vecFaceArea.GetRows(); ++i )
            vecSourceFaceArea[i] = m_meshInput->vecFaceArea[i];
    }
    else
    {
        DataArray3D< double > dataGLLJacobianSrc;
        this->InitializeCoordinatesFromMeshFE( *m_meshInput, m_nDofsPEl_Src, dataGLLNodesSrc, dSourceCenterLon,
                                               dSourceCenterLat, dSourceVertexLon, dSourceVertexLat );

        // Generate the continuous Jacobian for input mesh
        GenerateMetaData( *m_meshInput, m_nDofsPEl_Src, false /* fBubble */, dataGLLNodesSrc, dataGLLJacobianSrc );

        if( m_srcDiscType == DiscretizationType_CGLL )
        {
            GenerateUniqueJacobian( dataGLLNodesSrc, dataGLLJacobianSrc, vecSourceFaceArea );
        }
        else
        {
            GenerateDiscontinuousJacobian( dataGLLJacobianSrc, vecSourceFaceArea );
        }
    }

    if( m_destDiscType == DiscretizationType_FV || m_destDiscType == DiscretizationType_PCLOUD )
    {
        this->InitializeCoordinatesFromMeshFV(
            *m_meshOutput, dTargetCenterLon, dTargetCenterLat, dTargetVertexLon, dTargetVertexLat,
            ( this->m_remapper->m_target_type == moab::TempestRemapper::RLL ), /* fLatLon = false */
            m_remapper->max_target_edges );

        vecTargetFaceArea.Allocate( m_meshOutput->vecFaceArea.GetRows() );
        for( unsigned i = 0; i < m_meshOutput->vecFaceArea.GetRows(); ++i )
        {
            vecTargetFaceArea[i] = m_meshOutput->vecFaceArea[i];
        }
    }
    else
    {
        DataArray3D< double > dataGLLJacobianDest;
        this->InitializeCoordinatesFromMeshFE( *m_meshOutput, m_nDofsPEl_Dest, dataGLLNodesDest, dTargetCenterLon,
                                               dTargetCenterLat, dTargetVertexLon, dTargetVertexLat );

        // Generate the continuous Jacobian for input mesh
        GenerateMetaData( *m_meshOutput, m_nDofsPEl_Dest, false /* fBubble */, dataGLLNodesDest, dataGLLJacobianDest );

        if( m_destDiscType == DiscretizationType_CGLL )
        {
            GenerateUniqueJacobian( dataGLLNodesDest, dataGLLJacobianDest, vecTargetFaceArea );
        }
        else
        {
            GenerateDiscontinuousJacobian( dataGLLJacobianDest, vecTargetFaceArea );
        }
    }

    // Map dimensions
    unsigned nA = ( vecSourceFaceArea.GetRows() );
    unsigned nB = ( vecTargetFaceArea.GetRows() );

    std::vector< int > masksA, masksB;
    MB_CHK_SET_ERR( m_remapper->GetIMasks( moab::Remapper::SourceMesh, masksA ), "Trouble getting masks for source" );
    MB_CHK_SET_ERR( m_remapper->GetIMasks( moab::Remapper::TargetMesh, masksB ), "Trouble getting masks for target" );

    // Number of nodes per Face
    int nSourceNodesPerFace = dSourceVertexLon.GetColumns();
    int nTargetNodesPerFace = dTargetVertexLon.GetColumns();

    // if source or target cells have triangles at poles, center of those triangles need to come from
    // the original quad, not from center in 3d, converted to 2d again
    // start copy OnlineMap.cpp tempestremap
    // right now, do this only for source  mesh; copy the logic for target mesh
    for( unsigned i = 0; i < nA; i++ )
    {
        const Face& face = m_meshInput->faces[i];

        int nNodes          = face.edges.size();
        int indexNodeAtPole = -1;
        if( 3 == nNodes )  // check if one node at the poles
        {
            for( int j = 0; j < nNodes; j++ )
                if( fabs( fabs( dSourceVertexLat[i][j] ) - 90.0 ) < 1.0e-12 )
                {
                    indexNodeAtPole = j;
                    break;
                }
        }
        if( indexNodeAtPole < 0 ) continue;  // continue i loop, do nothing
        // recompute center of cell, from 3d data; add one 2 nodes at pole, and average
        int nodeAtPole = face[indexNodeAtPole];  // use the overloaded operator
        Node nodePole  = m_meshInput->nodes[nodeAtPole];
        Node newCenter = nodePole * 2;
        for( int j = 1; j < nNodes; j++ )
        {
            int indexi       = ( indexNodeAtPole + j ) % nNodes;  // nNodes is 3 !
            const Node& node = m_meshInput->nodes[face[indexi]];
            newCenter        = newCenter + node;
        }
        newCenter = newCenter * 0.25;
        newCenter = newCenter.Normalized();

#ifdef VERBOSE
        double iniLon = dSourceCenterLon[i], iniLat = dSourceCenterLat[i];
#endif
        // dSourceCenterLon, dSourceCenterLat
        XYZtoRLL_Deg( newCenter.x, newCenter.y, newCenter.z, dSourceCenterLon[i], dSourceCenterLat[i] );
#ifdef VERBOSE
        std::cout << " modify center of triangle from " << iniLon << " " << iniLat << " to " << dSourceCenterLon[i]
                  << " " << dSourceCenterLat[i] << "\n";
#endif
    }

    // first move data if in parallel
#if defined( MOAB_HAVE_MPI )
    int max_row_dof, max_col_dof;  // output; arrays will be re-distributed in chunks [maxdof/size]
    // if (size > 1)
    {
        int ierr = rearrange_arrays_by_dofs( srccol_gdofmap, vecSourceFaceArea, dSourceCenterLon, dSourceCenterLat,
                                             dSourceVertexLon, dSourceVertexLat, masksA, nA, nSourceNodesPerFace,
                                             max_col_dof );  // now nA will be close to maxdof/size
        if( ierr != 0 )
        {
            _EXCEPTION1( "Unable to arrange source data %d ", nA );
        }
        // rearrange target data: (nB)
        //
        ierr = rearrange_arrays_by_dofs( row_gdofmap, vecTargetFaceArea, dTargetCenterLon, dTargetCenterLat,
                                         dTargetVertexLon, dTargetVertexLat, masksB, nB, nTargetNodesPerFace,
                                         max_row_dof );  // now nA will be close to maxdof/size
        if( ierr != 0 )
        {
            _EXCEPTION1( "Unable to arrange target data %d ", nB );
        }
    }
#endif

    // Number of non-zeros in the remap matrix operator
    int nS = m_weightMatrix.nonZeros();

#if defined( MOAB_HAVE_MPI ) && defined( MOAB_HAVE_NETCDFPAR )
    int locbuf[5] = { (int)nA, (int)nB, nS, nSourceNodesPerFace, nTargetNodesPerFace };
    int offbuf[3] = { 0, 0, 0 };
    int globuf[5] = { 0, 0, 0, 0, 0 };
    MPI_Scan( locbuf, offbuf, 3, MPI_INT, MPI_SUM, m_pcomm->comm() );
    MPI_Allreduce( locbuf, globuf, 3, MPI_INT, MPI_SUM, m_pcomm->comm() );
    MPI_Allreduce( &locbuf[3], &globuf[3], 2, MPI_INT, MPI_MAX, m_pcomm->comm() );

    // MPI_Scan is inclusive of data in current rank; modify accordingly.
    offbuf[0] -= nA;
    offbuf[1] -= nB;
    offbuf[2] -= nS;

#else
    int offbuf[3] = { 0, 0, 0 };
    int globuf[5] = { (int)nA, (int)nB, nS, nSourceNodesPerFace, nTargetNodesPerFace };
#endif

    std::vector< std::string > srcdimNames, tgtdimNames;
    std::vector< int > srcdimSizes, tgtdimSizes;
    {
        if( m_remapper->m_source_type == moab::TempestRemapper::RLL && m_remapper->m_source_metadata.size() )
        {
            srcdimNames.push_back( "lat" );
            srcdimNames.push_back( "lon" );
            srcdimSizes.resize( 2, 0 );
            srcdimSizes[0] = m_remapper->m_source_metadata[0];
            srcdimSizes[1] = m_remapper->m_source_metadata[1];
        }
        else
        {
            srcdimNames.push_back( "num_elem" );
            srcdimSizes.push_back( globuf[0] );
        }

        if( m_remapper->m_target_type == moab::TempestRemapper::RLL && m_remapper->m_target_metadata.size() )
        {
            tgtdimNames.push_back( "lat" );
            tgtdimNames.push_back( "lon" );
            tgtdimSizes.resize( 2, 0 );
            tgtdimSizes[0] = m_remapper->m_target_metadata[0];
            tgtdimSizes[1] = m_remapper->m_target_metadata[1];
        }
        else
        {
            tgtdimNames.push_back( "num_elem" );
            tgtdimSizes.push_back( globuf[1] );
        }
    }

    // Write output dimensions entries
    unsigned nSrcGridDims = ( srcdimSizes.size() );
    unsigned nDstGridDims = ( tgtdimSizes.size() );

    NcDim* dimSrcGridRank = ncMap.add_dim( "src_grid_rank", nSrcGridDims );
    NcDim* dimDstGridRank = ncMap.add_dim( "dst_grid_rank", nDstGridDims );

    NcVar* varSrcGridDims = ncMap.add_var( "src_grid_dims", ncInt, dimSrcGridRank );
    NcVar* varDstGridDims = ncMap.add_var( "dst_grid_dims", ncInt, dimDstGridRank );

#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varSrcGridDims, is_independent );
    ncMap.enable_var_par_access( varDstGridDims, is_independent );
#endif

    // write dimension names
    {
        char szDim[64];
        for( unsigned i = 0; i < srcdimSizes.size(); i++ )
        {
            varSrcGridDims->set_cur( nSrcGridDims - i - 1 );
            varSrcGridDims->put( &( srcdimSizes[nSrcGridDims - i - 1] ), 1 );
        }

        for( unsigned i = 0; i < srcdimSizes.size(); i++ )
        {
            snprintf( szDim, 64, "name%i", i );
            varSrcGridDims->add_att( szDim, srcdimNames[nSrcGridDims - i - 1].c_str() );
        }

        for( unsigned i = 0; i < tgtdimSizes.size(); i++ )
        {
            varDstGridDims->set_cur( nDstGridDims - i - 1 );
            varDstGridDims->put( &( tgtdimSizes[nDstGridDims - i - 1] ), 1 );
        }

        for( unsigned i = 0; i < tgtdimSizes.size(); i++ )
        {
            snprintf( szDim, 64, "name%i", i );
            varDstGridDims->add_att( szDim, tgtdimNames[nDstGridDims - i - 1].c_str() );
        }
    }

    // Source and Target mesh resolutions
    NcDim* dimNA = ncMap.add_dim( "n_a", globuf[0] );
    NcDim* dimNB = ncMap.add_dim( "n_b", globuf[1] );

    // Number of nodes per Face
    NcDim* dimNVA = ncMap.add_dim( "nv_a", globuf[3] );
    NcDim* dimNVB = ncMap.add_dim( "nv_b", globuf[4] );

    // Write coordinates
    NcVar* varYCA = ncMap.add_var( "yc_a", ncDouble, dimNA );
    NcVar* varYCB = ncMap.add_var( "yc_b", ncDouble, dimNB );

    NcVar* varXCA = ncMap.add_var( "xc_a", ncDouble, dimNA );
    NcVar* varXCB = ncMap.add_var( "xc_b", ncDouble, dimNB );

    NcVar* varYVA = ncMap.add_var( "yv_a", ncDouble, dimNA, dimNVA );
    NcVar* varYVB = ncMap.add_var( "yv_b", ncDouble, dimNB, dimNVB );

    NcVar* varXVA = ncMap.add_var( "xv_a", ncDouble, dimNA, dimNVA );
    NcVar* varXVB = ncMap.add_var( "xv_b", ncDouble, dimNB, dimNVB );

    // Write masks
    NcVar* varMaskA = ncMap.add_var( "mask_a", ncInt, dimNA );
    NcVar* varMaskB = ncMap.add_var( "mask_b", ncInt, dimNB );

#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varYCA, is_independent );
    ncMap.enable_var_par_access( varYCB, is_independent );
    ncMap.enable_var_par_access( varXCA, is_independent );
    ncMap.enable_var_par_access( varXCB, is_independent );
    ncMap.enable_var_par_access( varYVA, is_independent );
    ncMap.enable_var_par_access( varYVB, is_independent );
    ncMap.enable_var_par_access( varXVA, is_independent );
    ncMap.enable_var_par_access( varXVB, is_independent );
    ncMap.enable_var_par_access( varMaskA, is_independent );
    ncMap.enable_var_par_access( varMaskB, is_independent );
#endif

    varYCA->add_att( "units", "degrees" );
    varYCB->add_att( "units", "degrees" );

    varXCA->add_att( "units", "degrees" );
    varXCB->add_att( "units", "degrees" );

    varYVA->add_att( "units", "degrees" );
    varYVB->add_att( "units", "degrees" );

    varXVA->add_att( "units", "degrees" );
    varXVB->add_att( "units", "degrees" );

    // Verify dimensionality
    if( dSourceCenterLon.GetRows() != nA )
    {
        _EXCEPTIONT( "Mismatch between dSourceCenterLon and nA" );
    }
    if( dSourceCenterLat.GetRows() != nA )
    {
        _EXCEPTIONT( "Mismatch between dSourceCenterLat and nA" );
    }
    if( dTargetCenterLon.GetRows() != nB )
    {
        _EXCEPTIONT( "Mismatch between dTargetCenterLon and nB" );
    }
    if( dTargetCenterLat.GetRows() != nB )
    {
        _EXCEPTIONT( "Mismatch between dTargetCenterLat and nB" );
    }
    if( dSourceVertexLon.GetRows() != nA )
    {
        _EXCEPTIONT( "Mismatch between dSourceVertexLon and nA" );
    }
    if( dSourceVertexLat.GetRows() != nA )
    {
        _EXCEPTIONT( "Mismatch between dSourceVertexLat and nA" );
    }
    if( dTargetVertexLon.GetRows() != nB )
    {
        _EXCEPTIONT( "Mismatch between dTargetVertexLon and nB" );
    }
    if( dTargetVertexLat.GetRows() != nB )
    {
        _EXCEPTIONT( "Mismatch between dTargetVertexLat and nB" );
    }

    varYCA->set_cur( (long)offbuf[0] );
    varYCA->put( &( dSourceCenterLat[0] ), nA );
    varYCB->set_cur( (long)offbuf[1] );
    varYCB->put( &( dTargetCenterLat[0] ), nB );

    varXCA->set_cur( (long)offbuf[0] );
    varXCA->put( &( dSourceCenterLon[0] ), nA );
    varXCB->set_cur( (long)offbuf[1] );
    varXCB->put( &( dTargetCenterLon[0] ), nB );

    varYVA->set_cur( (long)offbuf[0] );
    varYVA->put( &( dSourceVertexLat[0][0] ), nA, nSourceNodesPerFace );
    varYVB->set_cur( (long)offbuf[1] );
    varYVB->put( &( dTargetVertexLat[0][0] ), nB, nTargetNodesPerFace );

    varXVA->set_cur( (long)offbuf[0] );
    varXVA->put( &( dSourceVertexLon[0][0] ), nA, nSourceNodesPerFace );
    varXVB->set_cur( (long)offbuf[1] );
    varXVB->put( &( dTargetVertexLon[0][0] ), nB, nTargetNodesPerFace );

    varMaskA->set_cur( (long)offbuf[0] );
    varMaskA->put( &( masksA[0] ), nA );
    varMaskB->set_cur( (long)offbuf[1] );
    varMaskB->put( &( masksB[0] ), nB );

    // Write areas
    NcVar* varAreaA = ncMap.add_var( "area_a", ncDouble, dimNA );
#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varAreaA, is_independent );
#endif
    varAreaA->set_cur( (long)offbuf[0] );
    varAreaA->put( &( vecSourceFaceArea[0] ), nA );

    NcVar* varAreaB = ncMap.add_var( "area_b", ncDouble, dimNB );
#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varAreaB, is_independent );
#endif
    varAreaB->set_cur( (long)offbuf[1] );
    varAreaB->put( &( vecTargetFaceArea[0] ), nB );

    // Write SparseMatrix entries
    DataArray1D< int > vecRow( nS );
    DataArray1D< int > vecCol( nS );
    DataArray1D< double > vecS( nS );
    DataArray1D< double > dFracA( nA );
    DataArray1D< double > dFracB( nB );

    moab::TupleList tlValRow, tlValCol;
    unsigned numr = 1;  //
    // value has to be sent to processor row/nB for for fracA and col/nA for fracB
    // vecTargetArea (indexRow ) has to be sent for fracA (index col?)
    // vecTargetFaceArea will have to be sent to col index, with its index !
    tlValRow.initialize( 2, 0, 0, numr, nS );  // to proc(row),  global row , value
    tlValCol.initialize( 3, 0, 0, numr, nS );  // to proc(col),  global row / col, value
    tlValRow.enableWriteAccess();
    tlValCol.enableWriteAccess();
    /*
        dFracA[ col ] += val / vecSourceFaceArea[ col ] * vecTargetFaceArea[ row ];
        dFracB[ row ] += val ;
     */
    int offset = 0;
#if defined( MOAB_HAVE_MPI )
    int nAbase = ( max_col_dof + 1 ) / size;  // it is nA, except last rank ( == size - 1 )
    int nBbase = ( max_row_dof + 1 ) / size;  // it is nB, except last rank ( == size - 1 )
#endif
    for( int i = 0; i < m_weightMatrix.outerSize(); ++i )
    {
        for( WeightMatrix::InnerIterator it( m_weightMatrix, i ); it; ++it )
        {
            vecRow[offset] = 1 + this->GetRowGlobalDoF( it.row() );  // row index
            vecCol[offset] = 1 + this->GetColGlobalDoF( it.col() );  // col index
            vecS[offset]   = it.value();                             // value

#if defined( MOAB_HAVE_MPI )
            {
                // value M(row, col) will contribute to procRow and procCol values for fracA and fracB
                int procRow = ( vecRow[offset] - 1 ) / nBbase;
                if( procRow >= size ) procRow = size - 1;
                int procCol = ( vecCol[offset] - 1 ) / nAbase;
                if( procCol >= size ) procCol = size - 1;
                int nrInd                     = tlValRow.get_n();
                tlValRow.vi_wr[2 * nrInd]     = procRow;
                tlValRow.vi_wr[2 * nrInd + 1] = vecRow[offset] - 1;
                tlValRow.vr_wr[nrInd]         = vecS[offset];
                tlValRow.inc_n();
                int ncInd                     = tlValCol.get_n();
                tlValCol.vi_wr[3 * ncInd]     = procCol;
                tlValCol.vi_wr[3 * ncInd + 1] = vecRow[offset] - 1;
                tlValCol.vi_wr[3 * ncInd + 2] = vecCol[offset] - 1;  // this is column
                tlValCol.vr_wr[ncInd]         = vecS[offset];
                tlValCol.inc_n();
            }

#endif
            offset++;
        }
    }
#if defined( MOAB_HAVE_MPI )
    // need to send values for their row and col processors, to compute fractions there
    // now do the heavy communication
    ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tlValCol, 0 );
    ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tlValRow, 0 );

    // we have now, for example,  dFracB[ row ] += val ;
    // so we know that on current task, we received tlValRow
    // reminder dFracA[ col ] += val / vecSourceFaceArea[ col ] * vecTargetFaceArea[ row ];
    //          dFracB[ row ] += val ;
    for( unsigned i = 0; i < tlValRow.get_n(); i++ )
    {
        // int fromProc = tlValRow.vi_wr[2 * i];
        int gRowInd       = tlValRow.vi_wr[2 * i + 1];
        int localIndexRow = gRowInd - nBbase * rank;  // modulo nBbase rank is from 0 to size - 1;
        double wgt        = tlValRow.vr_wr[i];
        assert( localIndexRow >= 0 );
        assert( nB - localIndexRow > 0 );
        dFracB[localIndexRow] += wgt;
    }
    // to compute dFracA we need vecTargetFaceArea[ row ]; we know the row, and we can get the proc we need it from

    std::set< int > neededRows;
    for( unsigned i = 0; i < tlValCol.get_n(); i++ )
    {
        int rRowInd = tlValCol.vi_wr[3 * i + 1];
        neededRows.insert( rRowInd );
        // we need vecTargetFaceAreaGlobal[ rRowInd ]; this exists on proc procRow
    }
    moab::TupleList tgtAreaReq;
    tgtAreaReq.initialize( 2, 0, 0, 0, neededRows.size() );
    tgtAreaReq.enableWriteAccess();
    for( std::set< int >::iterator sit = neededRows.begin(); sit != neededRows.end(); sit++ )
    {
        int neededRow = *sit;
        int procRow   = neededRow / nBbase;
        if( procRow >= size ) procRow = size - 1;
        int nr                       = tgtAreaReq.get_n();
        tgtAreaReq.vi_wr[2 * nr]     = procRow;
        tgtAreaReq.vi_wr[2 * nr + 1] = neededRow;
        tgtAreaReq.inc_n();
    }

    ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tgtAreaReq, 0 );
    // we need to send back the tgtArea corresponding to row
    moab::TupleList tgtAreaInfo;  // load it with tgtArea at row
    tgtAreaInfo.initialize( 2, 0, 0, 1, tgtAreaReq.get_n() );
    tgtAreaInfo.enableWriteAccess();
    for( unsigned i = 0; i < tgtAreaReq.get_n(); i++ )
    {
        int from_proc     = tgtAreaReq.vi_wr[2 * i];
        int row           = tgtAreaReq.vi_wr[2 * i + 1];
        int locaIndexRow  = row - rank * nBbase;
        double areaToSend = vecTargetFaceArea[locaIndexRow];
        // int remoteIndex = tgtAreaReq.vi_wr[3*i + 2] ;

        tgtAreaInfo.vi_wr[2 * i]     = from_proc;  // send back requested info
        tgtAreaInfo.vi_wr[2 * i + 1] = row;
        tgtAreaInfo.vr_wr[i]         = areaToSend;  // this will be tgt area at row
        tgtAreaInfo.inc_n();
    }
    ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tgtAreaInfo, 0 );

    std::map< int, double > areaAtRow;
    for( unsigned i = 0; i < tgtAreaInfo.get_n(); i++ )
    {
        // we have received from proc, value for row !
        int row        = tgtAreaInfo.vi_wr[2 * i + 1];
        areaAtRow[row] = tgtAreaInfo.vr_wr[i];
    }

    // we have now for rows the
    // it is ordered by index, so:
    // now compute reminder dFracA[ col ] += val / vecSourceFaceArea[ col ] * vecTargetFaceArea[ row ];
    // tgtAreaInfo will have at index i the area we need (from row)
    // there should be an easier way :(
    for( unsigned i = 0; i < tlValCol.get_n(); i++ )
    {
        int rRowInd     = tlValCol.vi_wr[3 * i + 1];
        int colInd      = tlValCol.vi_wr[3 * i + 2];
        double val      = tlValCol.vr_wr[i];
        int localColInd = colInd - rank * nAbase;  // < local nA
        // we need vecTargetFaceAreaGlobal[ rRowInd ]; this exists on proc procRow
        auto itMap = areaAtRow.find( rRowInd );  // it should be different from end
        if( itMap != areaAtRow.end() )
        {
            double areaRow = itMap->second;  // we fished a lot for this !
            dFracA[localColInd] += val / vecSourceFaceArea[localColInd] * areaRow;
        }
    }

#endif
    // Load in data
    NcDim* dimNS = ncMap.add_dim( "n_s", globuf[2] );

    NcVar* varRow = ncMap.add_var( "row", ncInt, dimNS );
    NcVar* varCol = ncMap.add_var( "col", ncInt, dimNS );
    NcVar* varS   = ncMap.add_var( "S", ncDouble, dimNS );
#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varRow, is_independent );
    ncMap.enable_var_par_access( varCol, is_independent );
    ncMap.enable_var_par_access( varS, is_independent );
#endif

    varRow->set_cur( (long)offbuf[2] );
    varRow->put( vecRow, nS );

    varCol->set_cur( (long)offbuf[2] );
    varCol->put( vecCol, nS );

    varS->set_cur( (long)offbuf[2] );
    varS->put( &( vecS[0] ), nS );

    // Calculate and write fractional coverage arrays
    NcVar* varFracA = ncMap.add_var( "frac_a", ncDouble, dimNA );
#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varFracA, is_independent );
#endif
    varFracA->add_att( "name", "fraction of target coverage of source dof" );
    varFracA->add_att( "units", "unitless" );
    varFracA->set_cur( (long)offbuf[0] );
    varFracA->put( &( dFracA[0] ), nA );

    NcVar* varFracB = ncMap.add_var( "frac_b", ncDouble, dimNB );
#ifdef MOAB_HAVE_NETCDFPAR
    ncMap.enable_var_par_access( varFracB, is_independent );
#endif
    varFracB->add_att( "name", "fraction of source coverage of target dof" );
    varFracB->add_att( "units", "unitless" );
    varFracB->set_cur( (long)offbuf[1] );
    varFracB->put( &( dFracB[0] ), nB );

    // Add global attributes
    // std::map<std::string, std::string>::const_iterator iterAttributes =
    //     mapAttributes.begin();
    // for (; iterAttributes != mapAttributes.end(); iterAttributes++) {
    //     ncMap.add_att(
    //         iterAttributes->first.c_str(),
    //         iterAttributes->second.c_str());
    // }

    ncMap.close();

#ifdef VERBOSE
    serializeSparseMatrix( m_weightMatrix, "map_operator_" + std::to_string( rank ) + ".txt" );
#endif
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::WriteHDF5MapFile( const std::string& strOutputFile )
{
    /**
     * Need to get the global maximum of number of vertices per element
     * Key issue is that when calling InitializeCoordinatesFromMeshFV, the allocation for
     *dVertexLon/dVertexLat are made based on the maximum vertices in the current process. However,
     *when writing this out, other processes may have a different size for the same array. This is
     *hence a mess to consolidate in h5mtoscrip eventually.
     **/

    /* Let us compute all relevant data for the current original source mesh on the process */
    DataArray1D< double > vecSourceFaceArea, vecTargetFaceArea;
    DataArray1D< double > dSourceCenterLon, dSourceCenterLat, dTargetCenterLon, dTargetCenterLat;
    DataArray2D< double > dSourceVertexLon, dSourceVertexLat, dTargetVertexLon, dTargetVertexLat;
    if( m_srcDiscType == DiscretizationType_FV || m_srcDiscType == DiscretizationType_PCLOUD )
    {
        this->InitializeCoordinatesFromMeshFV(
            *m_meshInput, dSourceCenterLon, dSourceCenterLat, dSourceVertexLon, dSourceVertexLat,
            ( this->m_remapper->m_source_type == moab::TempestRemapper::RLL ) /* fLatLon = false */,
            m_remapper->max_source_edges );

        vecSourceFaceArea.Allocate( m_meshInput->vecFaceArea.GetRows() );
        for( unsigned i = 0; i < m_meshInput->vecFaceArea.GetRows(); ++i )
            vecSourceFaceArea[i] = m_meshInput->vecFaceArea[i];
    }
    else
    {
        DataArray3D< double > dataGLLJacobianSrc;
        this->InitializeCoordinatesFromMeshFE( *m_meshInput, m_nDofsPEl_Src, dataGLLNodesSrc, dSourceCenterLon,
                                               dSourceCenterLat, dSourceVertexLon, dSourceVertexLat );

        // Generate the continuous Jacobian for input mesh
        GenerateMetaData( *m_meshInput, m_nDofsPEl_Src, false /* fBubble */, dataGLLNodesSrc, dataGLLJacobianSrc );

        if( m_srcDiscType == DiscretizationType_CGLL )
        {
            GenerateUniqueJacobian( dataGLLNodesSrc, dataGLLJacobianSrc, vecSourceFaceArea );
        }
        else
        {
            GenerateDiscontinuousJacobian( dataGLLJacobianSrc, vecSourceFaceArea );
        }
    }

    if( m_destDiscType == DiscretizationType_FV || m_destDiscType == DiscretizationType_PCLOUD )
    {
        this->InitializeCoordinatesFromMeshFV(
            *m_meshOutput, dTargetCenterLon, dTargetCenterLat, dTargetVertexLon, dTargetVertexLat,
            ( this->m_remapper->m_target_type == moab::TempestRemapper::RLL ) /* fLatLon = false */,
            m_remapper->max_target_edges );

        vecTargetFaceArea.Allocate( m_meshOutput->vecFaceArea.GetRows() );
        for( unsigned i = 0; i < m_meshOutput->vecFaceArea.GetRows(); ++i )
            vecTargetFaceArea[i] = m_meshOutput->vecFaceArea[i];
    }
    else
    {
        DataArray3D< double > dataGLLJacobianDest;
        this->InitializeCoordinatesFromMeshFE( *m_meshOutput, m_nDofsPEl_Dest, dataGLLNodesDest, dTargetCenterLon,
                                               dTargetCenterLat, dTargetVertexLon, dTargetVertexLat );

        // Generate the continuous Jacobian for input mesh
        GenerateMetaData( *m_meshOutput, m_nDofsPEl_Dest, false /* fBubble */, dataGLLNodesDest, dataGLLJacobianDest );

        if( m_destDiscType == DiscretizationType_CGLL )
        {
            GenerateUniqueJacobian( dataGLLNodesDest, dataGLLJacobianDest, vecTargetFaceArea );
        }
        else
        {
            GenerateDiscontinuousJacobian( dataGLLJacobianDest, vecTargetFaceArea );
        }
    }

    moab::EntityHandle& m_meshOverlapSet = m_remapper->m_overlap_set;
    int tot_src_ents                     = m_remapper->m_source_entities.size();
    int tot_tgt_ents                     = m_remapper->m_target_entities.size();
    int tot_src_size                     = dSourceCenterLon.GetRows();
    int tot_tgt_size                     = m_dTargetCenterLon.GetRows();
    int tot_vsrc_size                    = dSourceVertexLon.GetRows() * dSourceVertexLon.GetColumns();
    int tot_vtgt_size                    = m_dTargetVertexLon.GetRows() * m_dTargetVertexLon.GetColumns();

    const int weightMatNNZ = m_weightMatrix.nonZeros();
    moab::Tag tagMapMetaData, tagMapIndexRow, tagMapIndexCol, tagMapValues, srcEleIDs, tgtEleIDs;
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SMAT_DATA", 13, moab::MB_TYPE_INTEGER, tagMapMetaData,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SMAT_ROWS", weightMatNNZ, moab::MB_TYPE_INTEGER, tagMapIndexRow,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SMAT_COLS", weightMatNNZ, moab::MB_TYPE_INTEGER, tagMapIndexCol,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SMAT_VALS", weightMatNNZ, moab::MB_TYPE_DOUBLE, tagMapValues,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceGIDS", tot_src_size, moab::MB_TYPE_INTEGER, srcEleIDs,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetGIDS", tot_tgt_size, moab::MB_TYPE_INTEGER, tgtEleIDs,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    moab::Tag srcAreaValues, tgtAreaValues;
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceAreas", tot_src_size, moab::MB_TYPE_DOUBLE, srcAreaValues,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetAreas", tot_tgt_size, moab::MB_TYPE_DOUBLE, tgtAreaValues,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    moab::Tag tagSrcCoordsCLon, tagSrcCoordsCLat, tagTgtCoordsCLon, tagTgtCoordsCLat;
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceCoordCenterLon", tot_src_size, moab::MB_TYPE_DOUBLE,
                                                 tagSrcCoordsCLon,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceCoordCenterLat", tot_src_size, moab::MB_TYPE_DOUBLE,
                                                 tagSrcCoordsCLat,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetCoordCenterLon", tot_tgt_size, moab::MB_TYPE_DOUBLE,
                                                 tagTgtCoordsCLon,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetCoordCenterLat", tot_tgt_size, moab::MB_TYPE_DOUBLE,
                                                 tagTgtCoordsCLat,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    moab::Tag tagSrcCoordsVLon, tagSrcCoordsVLat, tagTgtCoordsVLon, tagTgtCoordsVLat;
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceCoordVertexLon", tot_vsrc_size, moab::MB_TYPE_DOUBLE,
                                                 tagSrcCoordsVLon,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceCoordVertexLat", tot_vsrc_size, moab::MB_TYPE_DOUBLE,
                                                 tagSrcCoordsVLat,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetCoordVertexLon", tot_vtgt_size, moab::MB_TYPE_DOUBLE,
                                                 tagTgtCoordsVLon,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetCoordVertexLat", tot_vtgt_size, moab::MB_TYPE_DOUBLE,
                                                 tagTgtCoordsVLat,
                                                 moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                    "Retrieving tag handles failed" );
    moab::Tag srcMaskValues, tgtMaskValues;
    if( m_iSourceMask.IsAttached() )
    {
        MB_CHK_SET_ERR( m_interface->tag_get_handle( "SourceMask", m_iSourceMask.GetRows(), moab::MB_TYPE_INTEGER,
                                                     srcMaskValues,
                                                     moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                        "Retrieving tag handles failed" );
    }
    if( m_iTargetMask.IsAttached() )
    {
        MB_CHK_SET_ERR( m_interface->tag_get_handle( "TargetMask", m_iTargetMask.GetRows(), moab::MB_TYPE_INTEGER,
                                                     tgtMaskValues,
                                                     moab::MB_TAG_CREAT | moab::MB_TAG_SPARSE | moab::MB_TAG_VARLEN ),
                        "Retrieving tag handles failed" );
    }

    std::vector< int > smatrowvals( weightMatNNZ ), smatcolvals( weightMatNNZ );
    std::vector< double > smatvals( weightMatNNZ );
    // const double* smatvals = m_weightMatrix.valuePtr();
    // Loop over the matrix entries and find the max global ID for rows and columns
    for( int k = 0, offset = 0; k < m_weightMatrix.outerSize(); ++k )
    {
        for( moab::TempestOnlineMap::WeightMatrix::InnerIterator it( m_weightMatrix, k ); it; ++it, ++offset )
        {
            smatrowvals[offset] = this->GetRowGlobalDoF( it.row() );
            smatcolvals[offset] = this->GetColGlobalDoF( it.col() );
            smatvals[offset]    = it.value();
        }
    }

    /* Set the global IDs for the DoFs */
    ////
    // col_gdofmap [ col_ldofmap [ 0 : local_ndofs ] ] = GDOF
    // row_gdofmap [ row_ldofmap [ 0 : local_ndofs ] ] = GDOF
    ////
    int maxrow = 0, maxcol = 0;
    std::vector< int > src_global_dofs( tot_src_size ), tgt_global_dofs( tot_tgt_size );
    for( int i = 0; i < tot_src_size; ++i )
    {
        src_global_dofs[i] = srccol_gdofmap[i];
        maxcol             = ( src_global_dofs[i] > maxcol ) ? src_global_dofs[i] : maxcol;
    }

    for( int i = 0; i < tot_tgt_size; ++i )
    {
        tgt_global_dofs[i] = row_gdofmap[i];
        maxrow             = ( tgt_global_dofs[i] > maxrow ) ? tgt_global_dofs[i] : maxrow;
    }

    ///////////////////////////////////////////////////////////////////////////
    // The metadata in H5M file contains the following data:
    //
    //   1. n_a: Total source entities: (number of elements in source mesh)
    //   2. n_b: Total target entities: (number of elements in target mesh)
    //   3. nv_a: Max edge size of elements in source mesh
    //   4. nv_b: Max edge size of elements in target mesh
    //   5. maxrows: Number of rows in remap weight matrix
    //   6. maxcols: Number of cols in remap weight matrix
    //   7. nnz: Number of total nnz in sparse remap weight matrix
    //   8. np_a: The order of the field description on the source mesh: >= 1
    //   9. np_b: The order of the field description on the target mesh: >= 1
    //   10. method_a: The type of discretization for field on source mesh: [0 = FV, 1 = cGLL, 2 =
    //   dGLL]
    //   11. method_b: The type of discretization for field on target mesh: [0 = FV, 1 = cGLL, 2 =
    //   dGLL]
    //   12. conserved: Flag to specify whether the remap operator has conservation constraints: [0,
    //   1]
    //   13. monotonicity: Flags to specify whether the remap operator has monotonicity constraints:
    //   [0, 1, 2]
    //
    ///////////////////////////////////////////////////////////////////////////
    int map_disc_details[6];
    map_disc_details[0] = m_nDofsPEl_Src;
    map_disc_details[1] = m_nDofsPEl_Dest;
    map_disc_details[2] = ( m_srcDiscType == DiscretizationType_FV || m_srcDiscType == DiscretizationType_PCLOUD
                                ? 0
                                : ( m_srcDiscType == DiscretizationType_CGLL ? 1 : 2 ) );
    map_disc_details[3] = ( m_destDiscType == DiscretizationType_FV || m_destDiscType == DiscretizationType_PCLOUD
                                ? 0
                                : ( m_destDiscType == DiscretizationType_CGLL ? 1 : 2 ) );
    map_disc_details[4] = ( m_bConserved ? 1 : 0 );
    map_disc_details[5] = m_iMonotonicity;

#ifdef MOAB_HAVE_MPI
    int loc_smatmetadata[13] = { tot_src_ents,
                                 tot_tgt_ents,
                                 m_remapper->max_source_edges,
                                 m_remapper->max_target_edges,
                                 maxrow + 1,
                                 maxcol + 1,
                                 weightMatNNZ,
                                 map_disc_details[0],
                                 map_disc_details[1],
                                 map_disc_details[2],
                                 map_disc_details[3],
                                 map_disc_details[4],
                                 map_disc_details[5] };
    MB_CHK_SET_ERR( m_interface->tag_set_data( tagMapMetaData, &m_meshOverlapSet, 1, &loc_smatmetadata[0] ),
                    "Setting local tag data failed" );
    int glb_smatmetadata[13] = { 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 map_disc_details[0],
                                 map_disc_details[1],
                                 map_disc_details[2],
                                 map_disc_details[3],
                                 map_disc_details[4],
                                 map_disc_details[5] };
    int loc_buf[7]           = {
        tot_src_ents, tot_tgt_ents, weightMatNNZ, m_remapper->max_source_edges, m_remapper->max_target_edges,
        maxrow,       maxcol };
    int glb_buf[4] = { 0, 0, 0, 0 };
    MPI_Reduce( &loc_buf[0], &glb_buf[0], 3, MPI_INT, MPI_SUM, 0, m_pcomm->comm() );
    glb_smatmetadata[0] = glb_buf[0];
    glb_smatmetadata[1] = glb_buf[1];
    glb_smatmetadata[6] = glb_buf[2];
    MPI_Reduce( &loc_buf[3], &glb_buf[0], 4, MPI_INT, MPI_MAX, 0, m_pcomm->comm() );
    glb_smatmetadata[2] = glb_buf[0];
    glb_smatmetadata[3] = glb_buf[1];
    glb_smatmetadata[4] = glb_buf[2];
    glb_smatmetadata[5] = glb_buf[3];
#else
    int glb_smatmetadata[13] = { tot_src_ents,
                                 tot_tgt_ents,
                                 m_remapper->max_source_edges,
                                 m_remapper->max_target_edges,
                                 maxrow,
                                 maxcol,
                                 weightMatNNZ,
                                 map_disc_details[0],
                                 map_disc_details[1],
                                 map_disc_details[2],
                                 map_disc_details[3],
                                 map_disc_details[4],
                                 map_disc_details[5] };
#endif
    // These values represent number of rows and columns. So should be 1-based.
    glb_smatmetadata[4]++;
    glb_smatmetadata[5]++;

    if( this->is_root )
    {
        std::cout << "  " << this->rank << "  Writing remap weights with size [" << glb_smatmetadata[4] << " X "
                  << glb_smatmetadata[5] << "] and NNZ = " << glb_smatmetadata[6] << std::endl;
        EntityHandle root_set = 0;
        MB_CHK_SET_ERR( m_interface->tag_set_data( tagMapMetaData, &root_set, 1, &glb_smatmetadata[0] ),
                        "Setting local tag data failed" );
    }

    int dsize;
    const int numval          = weightMatNNZ;
    const void* smatrowvals_d = smatrowvals.data();
    const void* smatcolvals_d = smatcolvals.data();
    const void* smatvals_d    = smatvals.data();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagMapIndexRow, &m_meshOverlapSet, 1, &smatrowvals_d, &numval ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagMapIndexCol, &m_meshOverlapSet, 1, &smatcolvals_d, &numval ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagMapValues, &m_meshOverlapSet, 1, &smatvals_d, &numval ),
                    "Setting local tag data failed" );

    /* Set the global IDs for the DoFs */
    const void* srceleidvals_d = src_global_dofs.data();
    const void* tgteleidvals_d = tgt_global_dofs.data();
    dsize                      = src_global_dofs.size();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( srcEleIDs, &m_meshOverlapSet, 1, &srceleidvals_d, &dsize ),
                    "Setting local tag data failed" );
    dsize = tgt_global_dofs.size();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tgtEleIDs, &m_meshOverlapSet, 1, &tgteleidvals_d, &dsize ),
                    "Setting local tag data failed" );

    /* Set the source and target areas */
    const void* srcareavals_d = vecSourceFaceArea;
    const void* tgtareavals_d = vecTargetFaceArea;
    dsize                     = tot_src_size;
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( srcAreaValues, &m_meshOverlapSet, 1, &srcareavals_d, &dsize ),
                    "Setting local tag data failed" );
    dsize = tot_tgt_size;
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tgtAreaValues, &m_meshOverlapSet, 1, &tgtareavals_d, &dsize ),
                    "Setting local tag data failed" );

    /* Set the coordinates for source and target center vertices */
    const void* srccoordsclonvals_d = &dSourceCenterLon[0];
    const void* srccoordsclatvals_d = &dSourceCenterLat[0];
    dsize                           = dSourceCenterLon.GetRows();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagSrcCoordsCLon, &m_meshOverlapSet, 1, &srccoordsclonvals_d, &dsize ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagSrcCoordsCLat, &m_meshOverlapSet, 1, &srccoordsclatvals_d, &dsize ),
                    "Setting local tag data failed" );
    const void* tgtcoordsclonvals_d = &m_dTargetCenterLon[0];
    const void* tgtcoordsclatvals_d = &m_dTargetCenterLat[0];
    dsize                           = vecTargetFaceArea.GetRows();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagTgtCoordsCLon, &m_meshOverlapSet, 1, &tgtcoordsclonvals_d, &dsize ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagTgtCoordsCLat, &m_meshOverlapSet, 1, &tgtcoordsclatvals_d, &dsize ),
                    "Setting local tag data failed" );

    /* Set the coordinates for source and target element vertices */
    const void* srccoordsvlonvals_d = &( dSourceVertexLon[0][0] );
    const void* srccoordsvlatvals_d = &( dSourceVertexLat[0][0] );
    dsize                           = dSourceVertexLon.GetRows() * dSourceVertexLon.GetColumns();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagSrcCoordsVLon, &m_meshOverlapSet, 1, &srccoordsvlonvals_d, &dsize ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagSrcCoordsVLat, &m_meshOverlapSet, 1, &srccoordsvlatvals_d, &dsize ),
                    "Setting local tag data failed" );
    const void* tgtcoordsvlonvals_d = &( m_dTargetVertexLon[0][0] );
    const void* tgtcoordsvlatvals_d = &( m_dTargetVertexLat[0][0] );
    dsize                           = m_dTargetVertexLon.GetRows() * m_dTargetVertexLon.GetColumns();
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagTgtCoordsVLon, &m_meshOverlapSet, 1, &tgtcoordsvlonvals_d, &dsize ),
                    "Setting local tag data failed" );
    MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tagTgtCoordsVLat, &m_meshOverlapSet, 1, &tgtcoordsvlatvals_d, &dsize ),
                    "Setting local tag data failed" );

    /* Set the masks for source and target meshes if available */
    if( m_iSourceMask.IsAttached() )
    {
        const void* srcmaskvals_d = m_iSourceMask;
        dsize                     = m_iSourceMask.GetRows();
        MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( srcMaskValues, &m_meshOverlapSet, 1, &srcmaskvals_d, &dsize ),
                        "Setting local tag data failed" );
    }

    if( m_iTargetMask.IsAttached() )
    {
        const void* tgtmaskvals_d = m_iTargetMask;
        dsize                     = m_iTargetMask.GetRows();
        MB_CHK_SET_ERR( m_interface->tag_set_by_ptr( tgtMaskValues, &m_meshOverlapSet, 1, &tgtmaskvals_d, &dsize ),
                        "Setting local tag data failed" );
    }

#ifdef MOAB_HAVE_MPI
    const char* writeOptions = ( this->size > 1 ? "PARALLEL=WRITE_PART" : "" );
#else
    const char* writeOptions = "";
#endif

    // EntityHandle sets[3] = {m_remapper->m_source_set, m_remapper->m_target_set, m_remapper->m_overlap_set};
    EntityHandle sets[1] = { m_remapper->m_overlap_set };
    MB_CHK_ERR( m_interface->write_file( strOutputFile.c_str(), NULL, writeOptions, sets, 1 ) );

#ifdef WRITE_SCRIP_FILE
    sstr.str( "" );
    sstr << ctx.outFilename.substr( 0, lastindex ) << "_" << proc_id << ".nc";
    std::map< std::string, std::string > mapAttributes;
    mapAttributes["Creator"] = "MOAB mbtempest workflow";
    if( !ctx.proc_id ) std::cout << "Writing offline map to file: " << sstr.str() << std::endl;
    this->Write( strOutputFile.c_str(), mapAttributes, NcFile::Netcdf4 );
    sstr.str( "" );
#endif

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

void print_progress( const int barWidth, const float progress, const char* message )
{
    std::cout << message << " [";
    int pos = barWidth * progress;
    for( int i = 0; i < barWidth; ++i )
    {
        if( i < pos )
            std::cout << "=";
        else if( i == pos )
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int( progress * 100.0 ) << " %\r";
    std::cout.flush();
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// ReadParallelMap: Read a SCRIP-format map file and distribute sparse matrix
// data across MPI ranks.
//
// Strategy (adaptive, based on NNZ count):
//   1. Serial (size == 1): rank 0 reads entire file directly.
//   2. Buffered read (size > 1, nS <= NNZ threshold): rank 0 reads the file
//      using serial NcFile in fixed-size chunks, determines row ownership for
//      each entry, and scatters data to owning ranks via MPI point-to-point.
//      This avoids the need for parallel NetCDF and scales well for small-to-
//      medium maps by reducing file system contention.
//   3. Direct parallel read (size > 1, nS > NNZ threshold): all ranks read
//      their stripe of the file simultaneously using PNetCDF (preferred) or
//      parallel HDF5-backed NetCDF (NETCDFPAR). Falls back to buffered read
//      if neither is available.
//
// After the initial read/scatter, the downstream TupleList redistribution
// (for owned_dof_ids-based repartitioning) and Eigen sparse matrix assembly
// are unchanged regardless of which read strategy was used.
//
///////////////////////////////////////////////////////////////////////////////

// Tuning constants for the buffered read strategy.
// Adjust these for scalability studies on different platforms.

/// NNZ threshold: maps with nS <= this value use the buffered read strategy.
/// Maps with nS > this value use direct parallel I/O (if available).
/// Default 3M entries corresponds to ~36 MB of raw data (row+col+S).
static constexpr int BUFFERED_READ_NNZ_THRESHOLD = 3000000;

/// Buffer size in bytes for each chunk read by rank 0 in buffered mode.
/// Each sparse matrix entry is 12 bytes (2 ints + 1 double), so 64KB holds
/// ~5461 entries. Larger buffers reduce the number of read+scatter rounds
/// but increase peak memory on rank 0.
static constexpr int BUFFERED_READ_CHUNK_BYTES = 64 * 1024;

moab::ErrorCode moab::TempestOnlineMap::ReadParallelMap( const char* strSource,
                                                         const std::vector< int >& owned_dof_ids,
                                                         int arearead,
                                                         std::vector< double >& vecAreaA,
                                                         int& nA,
                                                         std::vector< double >& vecAreaB,
                                                         int& nB )
{
    NcError error( NcError::silent_nonfatal );

    const bool readAreaA = ( 1 == arearead || 3 == arearead );
    const bool readAreaB = ( 2 == arearead || 3 == arearead );
    int nS = 0;

    // =========================================================================
    // Phase 1: Read map dimensions (nA, nB, nS) and sparse matrix data.
    //
    // The read strategy is selected adaptively:
    //   - Serial or buffered read: rank 0 opens the file with serial NcFile,
    //     reads dimensions, and (for buffered mode) scatters data in chunks.
    //   - Direct parallel read: all ranks open the file with PNetCDF or
    //     NETCDFPAR and read their stripe directly.
    // =========================================================================

    std::vector< int > vecRow, vecCol;
    std::vector< double > vecS;
    int localSize = 0;  // number of sparse matrix entries on this rank after read

    // Determine which read strategy to use. For size == 1, always serial.
    // For size > 1, decide after reading dimensions (need nS).
    // We use a two-phase approach: first read dimensions on rank 0 and broadcast,
    // then select the strategy based on nS.

#ifdef MOAB_HAVE_MPI
    if( size > 1 )
    {
        // --- Multi-process path: read dimensions on rank 0 and broadcast ---
        int dims[3] = { 0, 0, 0 };  // nA, nB, nS

        if( rank == 0 )
        {
            NcFile ncDims( strSource, NcFile::ReadOnly );
            if( !ncDims.is_valid() )
            {
                _EXCEPTION1( "Unable to open input map file \"%s\" on rank 0", strSource );
            }
            NcDim* dimNA = ncDims.get_dim( "n_a" );
            NcDim* dimNB = ncDims.get_dim( "n_b" );
            NcDim* dimNS = ncDims.get_dim( "n_s" );
            if( !dimNA || !dimNB || !dimNS )
            {
                _EXCEPTION1( "Map file \"%s\" missing required dimensions (n_a, n_b, n_s)", strSource );
            }
            dims[0] = static_cast< int >( dimNA->size() );
            dims[1] = static_cast< int >( dimNB->size() );
            dims[2] = static_cast< int >( dimNS->size() );
            ncDims.close();
        }

        MPI_Bcast( dims, 3, MPI_INT, 0, m_pcomm->comm() );
        nA = dims[0];
        nB = dims[1];
        nS = dims[2];

        // Select read strategy based on NNZ count and available parallel I/O
        bool useBufferedRead = true;  // default for small maps or no parallel I/O

        if( nS > BUFFERED_READ_NNZ_THRESHOLD )
        {
            // Large map: prefer direct parallel read if available
#if defined( MOAB_HAVE_PNETCDF ) || defined( MOAB_HAVE_NETCDFPAR )
            useBufferedRead = false;
#endif
            // If neither is available, fall back to buffered read regardless of size
        }

        if( useBufferedRead )
        {
            // =================================================================
            // Buffered read: rank 0 reads in chunks and scatters to owners.
            //
            // Row ownership is determined by trivial partitioning: row i is
            // owned by rank (i / nRowPerPart), with remainder on rank 0.
            // Each chunk is read, ownership is computed per entry, and the
            // data is scattered via MPI_Scatter + MPI_Isend/MPI_Irecv.
            // =================================================================
            if( rank == 0 )
            {
                std::cout << "  [ReadParallelMap]: Using buffered read strategy for " << nS
                          << " NNZ entries (threshold=" << BUFFERED_READ_NNZ_THRESHOLD << ")\n";
            }

            const int nNNZBytes       = 2 * sizeof( int ) + sizeof( double );
            const int nMaxPerChunk    = BUFFERED_READ_CHUNK_BYTES / nNNZBytes;
            const int nBufferedReads  = static_cast< int >( std::ceil( 1.0 * nS / nMaxPerChunk ) );

            // Row ownership: trivial partitioning of nB rows across ranks
            const int nRowPerPart   = nB / size;
            const int nRowRemainder = nB % size;
            std::vector< int > rowOwnership( size );
            rowOwnership[0] = nRowPerPart + nRowRemainder;
            for( int ip = 1; ip < size; ++ip )
                rowOwnership[ip] = rowOwnership[ip - 1] + nRowPerPart;

            // File handle and variable pointers (rank 0 only)
            NcFile* ncMap        = nullptr;
            NcVar *varRowF       = nullptr, *varColF = nullptr, *varSF = nullptr;
            NcVar *varAreaAF     = nullptr, *varAreaBF = nullptr;

            if( rank == 0 )
            {
                ncMap = new NcFile( strSource, NcFile::ReadOnly );
                if( !ncMap->is_valid() )
                {
                    _EXCEPTION1( "Unable to open map file \"%s\" for buffered read", strSource );
                }
                varRowF = ncMap->get_var( "row" );
                varColF = ncMap->get_var( "col" );
                varSF   = ncMap->get_var( "S" );
                if( readAreaA ) varAreaAF = ncMap->get_var( "area_a" );
                if( readAreaB ) varAreaBF = ncMap->get_var( "area_b" );
            }

            // Accumulate received entries per rank
            std::vector< int > localRows, localCols;
            std::vector< double > localVals;
            localRows.reserve( nS / size + nS / ( size * 10 ) );  // slight overalloc
            localCols.reserve( nS / size + nS / ( size * 10 ) );
            localVals.reserve( nS / size + nS / ( size * 10 ) );

            int nEntriesRemaining = nS;
            long fileOffset       = 0;

            for( int iRead = 0; iRead < nBufferedReads; ++iRead )
            {
                // Per-chunk data and ownership (rank 0 only)
                std::vector< int > chunkRow, chunkCol;
                std::vector< double > chunkS;
                std::vector< std::vector< int > > entriesPerProc( size );
                std::vector< int > nPerProc( size, 0 );

                if( rank == 0 )
                {
                    int chunkSize = std::min( nEntriesRemaining, nMaxPerChunk );

                    chunkRow.resize( chunkSize );
                    chunkCol.resize( chunkSize );
                    chunkS.resize( chunkSize );

                    varRowF->set_cur( fileOffset );
                    varRowF->get( chunkRow.data(), chunkSize );
                    varColF->set_cur( fileOffset );
                    varColF->get( chunkCol.data(), chunkSize );
                    varSF->set_cur( fileOffset );
                    varSF->get( chunkS.data(), chunkSize );

                    // Determine ownership of each entry by its row index (1-based in file)
                    for( int ip = 0; ip < size; ++ip )
                        entriesPerProc[ip].reserve( chunkSize / size + 64 );

                    for( int i = 0; i < chunkSize; ++i )
                    {
                        int rowIdx = chunkRow[i] - 1;  // convert to 0-based
                        int owner  = 0;
                        if( rowIdx >= rowOwnership[0] )
                        {
                            // Binary search for owner
                            owner = static_cast< int >(
                                std::upper_bound( rowOwnership.begin(), rowOwnership.end(), rowIdx ) -
                                rowOwnership.begin() );
                            if( owner >= size ) owner = size - 1;
                        }
                        entriesPerProc[owner].push_back( i );
                    }

                    fileOffset += chunkSize;
                    nEntriesRemaining -= chunkSize;

                    for( int ip = 0; ip < size; ++ip )
                        nPerProc[ip] = static_cast< int >( entriesPerProc[ip].size() );
                }

                // Scatter count of entries each rank will receive in this chunk
                int nRecv = 0;
                MPI_Scatter( nPerProc.data(), 1, MPI_INT, &nRecv, 1, MPI_INT, 0, m_pcomm->comm() );

                if( rank == 0 )
                {
                    // Send data to remote ranks via non-blocking sends
                    std::vector< MPI_Request > requests;
                    requests.reserve( 2 * ( size - 1 ) );

                    // Pack and send to each remote rank
                    std::vector< std::vector< int > > sendRowCol( size );
                    std::vector< std::vector< double > > sendVals( size );

                    for( int ip = 1; ip < size; ++ip )
                    {
                        const int nDPP = nPerProc[ip];
                        if( nDPP > 0 )
                        {
                            sendRowCol[ip].resize( 2 * nDPP );
                            sendVals[ip].resize( nDPP );
                            for( int j = 0; j < nDPP; ++j )
                            {
                                int idx                  = entriesPerProc[ip][j];
                                sendRowCol[ip][2 * j]     = chunkRow[idx];
                                sendRowCol[ip][2 * j + 1] = chunkCol[idx];
                                sendVals[ip][j]           = chunkS[idx];
                            }

                            MPI_Request rqRC, rqV;
                            MPI_Isend( sendRowCol[ip].data(), 2 * nDPP, MPI_INT, ip,
                                       iRead * 1000, m_pcomm->comm(), &rqRC );
                            MPI_Isend( sendVals[ip].data(), nDPP, MPI_DOUBLE, ip,
                                       iRead * 1000 + 1, m_pcomm->comm(), &rqV );
                            requests.push_back( rqRC );
                            requests.push_back( rqV );
                        }
                    }

                    // Process rank 0's own entries while sends are in flight
                    for( int j = 0; j < nRecv; ++j )
                    {
                        int idx = entriesPerProc[0][j];
                        localRows.push_back( chunkRow[idx] );
                        localCols.push_back( chunkCol[idx] );
                        localVals.push_back( chunkS[idx] );
                    }

                    // Wait for all sends to complete
                    if( !requests.empty() )
                    {
                        std::vector< MPI_Status > stats( requests.size() );
                        MPI_Waitall( static_cast< int >( requests.size() ), requests.data(), stats.data() );
                    }
                }
                else if( nRecv > 0 )
                {
                    // Receive data from rank 0
                    std::vector< int > recvRowCol( 2 * nRecv );
                    std::vector< double > recvVals( nRecv );

                    MPI_Request rqs[2];
                    MPI_Irecv( recvRowCol.data(), 2 * nRecv, MPI_INT, 0,
                               iRead * 1000, m_pcomm->comm(), &rqs[0] );
                    MPI_Irecv( recvVals.data(), nRecv, MPI_DOUBLE, 0,
                               iRead * 1000 + 1, m_pcomm->comm(), &rqs[1] );

                    MPI_Status sts[2];
                    MPI_Waitall( 2, rqs, sts );

                    for( int j = 0; j < nRecv; ++j )
                    {
                        localRows.push_back( recvRowCol[2 * j] );
                        localCols.push_back( recvRowCol[2 * j + 1] );
                        localVals.push_back( recvVals[j] );
                    }
                }

                MPI_Barrier( m_pcomm->comm() );
            }  // end buffered read loop

            // Read area arrays on rank 0 and broadcast (small relative to sparse matrix)
            if( readAreaA )
            {
                vecAreaA.resize( nA );
                if( rank == 0 && varAreaAF )
                {
                    varAreaAF->set_cur( 0L );
                    varAreaAF->get( vecAreaA.data(), nA );
                }
                MPI_Bcast( vecAreaA.data(), nA, MPI_DOUBLE, 0, m_pcomm->comm() );
            }
            if( readAreaB )
            {
                vecAreaB.resize( nB );
                if( rank == 0 && varAreaBF )
                {
                    varAreaBF->set_cur( 0L );
                    varAreaBF->get( vecAreaB.data(), nB );
                }
                MPI_Bcast( vecAreaB.data(), nB, MPI_DOUBLE, 0, m_pcomm->comm() );
            }

            if( rank == 0 )
            {
                ncMap->close();
                delete ncMap;
            }

            // Move accumulated data into the standard vecRow/vecCol/vecS vectors
            localSize = static_cast< int >( localRows.size() );
            vecRow.swap( localRows );
            vecCol.swap( localCols );
            vecS.swap( localVals );
        }
        else
        {
            // =================================================================
            // Direct parallel read: all ranks read their stripe simultaneously.
            //
            // Strategy priority:
            //   1. PNetCDF (preferred — collective I/O, best scalability)
            //   2. NETCDFPAR (parallel HDF5-backed NetCDF)
            // =================================================================
            if( rank == 0 )
            {
                std::cout << "  [ReadParallelMap]: Using direct parallel read for " << nS
                          << " NNZ entries (threshold=" << BUFFERED_READ_NNZ_THRESHOLD << ")\n";
            }

            // Compute this rank's stripe of the sparse matrix
            localSize        = nS / size;
            long offsetRead  = rank * localSize;
            if( rank == size - 1 ) localSize += nS % size;

            vecRow.resize( localSize );
            vecCol.resize( localSize );
            vecS.resize( localSize );

            // Compute this rank's stripe of area arrays
            int localSizeA   = nA / size;
            long offsetReadA = rank * localSizeA;
            if( rank == size - 1 ) localSizeA += nA % size;

            int localSizeB   = nB / size;
            long offsetReadB = rank * localSizeB;
            if( rank == size - 1 ) localSizeB += nB % size;

            if( readAreaA ) vecAreaA.resize( localSizeA );
            if( readAreaB ) vecAreaB.resize( localSizeB );

#ifdef MOAB_HAVE_PNETCDF
            // PNetCDF path (preferred): collective parallel I/O
            int ncfile = -1;
            ERR_PARNC( ncmpi_open( m_pcomm->comm(), strSource, NC_NOWRITE, MPI_INFO_NULL, &ncfile ) );

            MPI_Offset start = static_cast< MPI_Offset >( offsetRead );
            MPI_Offset count = static_cast< MPI_Offset >( localSize );
            int varid;

            ERR_PARNC( ncmpi_inq_varid( ncfile, "S", &varid ) );
            ERR_PARNC( ncmpi_get_vara_double_all( ncfile, varid, &start, &count, vecS.data() ) );
            ERR_PARNC( ncmpi_inq_varid( ncfile, "row", &varid ) );
            ERR_PARNC( ncmpi_get_vara_int_all( ncfile, varid, &start, &count, vecRow.data() ) );
            ERR_PARNC( ncmpi_inq_varid( ncfile, "col", &varid ) );
            ERR_PARNC( ncmpi_get_vara_int_all( ncfile, varid, &start, &count, vecCol.data() ) );

            if( readAreaA )
            {
                MPI_Offset startA = static_cast< MPI_Offset >( offsetReadA );
                MPI_Offset countA = static_cast< MPI_Offset >( localSizeA );
                ERR_PARNC( ncmpi_inq_varid( ncfile, "area_a", &varid ) );
                ERR_PARNC( ncmpi_get_vara_double_all( ncfile, varid, &startA, &countA, vecAreaA.data() ) );
            }
            if( readAreaB )
            {
                MPI_Offset startB = static_cast< MPI_Offset >( offsetReadB );
                MPI_Offset countB = static_cast< MPI_Offset >( localSizeB );
                ERR_PARNC( ncmpi_inq_varid( ncfile, "area_b", &varid ) );
                ERR_PARNC( ncmpi_get_vara_double_all( ncfile, varid, &startB, &countB, vecAreaB.data() ) );
            }
            ERR_PARNC( ncmpi_close( ncfile ) );

#elif defined( MOAB_HAVE_NETCDFPAR )
            // Parallel HDF5-backed NetCDF path
            ParNcFile ncMap( m_pcomm->comm(), MPI_INFO_NULL, strSource, NcFile::ReadOnly, NcFile::Netcdf4 );
            if( !ncMap.is_valid() )
            {
                _EXCEPTION1( "Unable to open map file \"%s\" with parallel NetCDF", strSource );
            }

            NcVar* varRowP = ncMap.get_var( "row" );
            NcVar* varColP = ncMap.get_var( "col" );
            NcVar* varSP   = ncMap.get_var( "S" );
            ncMap.enable_var_par_access( varRowP, true );
            ncMap.enable_var_par_access( varColP, true );
            ncMap.enable_var_par_access( varSP, true );

            varRowP->set_cur( offsetRead );
            varRowP->get( vecRow.data(), localSize );
            varColP->set_cur( offsetRead );
            varColP->get( vecCol.data(), localSize );
            varSP->set_cur( offsetRead );
            varSP->get( vecS.data(), localSize );

            if( readAreaA )
            {
                NcVar* varAreaAP = ncMap.get_var( "area_a" );
                ncMap.enable_var_par_access( varAreaAP, true );
                varAreaAP->set_cur( offsetReadA );
                varAreaAP->get( vecAreaA.data(), localSizeA );
            }
            if( readAreaB )
            {
                NcVar* varAreaBP = ncMap.get_var( "area_b" );
                ncMap.enable_var_par_access( varAreaBP, true );
                varAreaBP->set_cur( offsetReadB );
                varAreaBP->get( vecAreaB.data(), localSizeB );
            }
            ncMap.close();
#endif
        }  // end direct parallel read
    }
    else
#endif  // MOAB_HAVE_MPI
    {
        // =================================================================
        // Serial path (size == 1): read entire file on the single process.
        // =================================================================
        NcFile ncMap( strSource, NcFile::ReadOnly );
        if( !ncMap.is_valid() )
        {
            _EXCEPTION1( "Unable to open input map file \"%s\"", strSource );
        }

        NcDim* dimNS = ncMap.get_dim( "n_s" );
        NcDim* dimNA = ncMap.get_dim( "n_a" );
        NcDim* dimNB = ncMap.get_dim( "n_b" );
        if( !dimNS || !dimNA || !dimNB )
        {
            _EXCEPTION1( "Map file \"%s\" missing required dimensions", strSource );
        }
        nS = static_cast< int >( dimNS->size() );
        nA = static_cast< int >( dimNA->size() );
        nB = static_cast< int >( dimNB->size() );

        localSize = nS;
        vecRow.resize( nS );
        vecCol.resize( nS );
        vecS.resize( nS );

        NcVar* varRowS = ncMap.get_var( "row" );
        NcVar* varColS = ncMap.get_var( "col" );
        NcVar* varSS   = ncMap.get_var( "S" );
        varRowS->get( vecRow.data(), nS );
        varColS->get( vecCol.data(), nS );
        varSS->get( vecS.data(), nS );

        if( readAreaA )
        {
            vecAreaA.resize( nA );
            NcVar* varAreaAS = ncMap.get_var( "area_a" );
            if( varAreaAS ) varAreaAS->get( vecAreaA.data(), nA );
        }
        if( readAreaB )
        {
            vecAreaB.resize( nB );
            NcVar* varAreaBS = ncMap.get_var( "area_b" );
            if( varAreaBS ) varAreaBS->get( vecAreaB.data(), nB );
        }
        ncMap.close();
    }

    // =========================================================================
    // Phase 2: Redistribute sparse matrix entries to their final owning ranks.
    //
    // After Phase 1, each rank holds a portion of the sparse matrix entries
    // (either its owned rows from the buffered read, or a stripe from the
    // direct parallel read). The rows/cols are still 1-based (SCRIP format).
    //
    // This phase uses TupleList-based crystal router communication to send
    // entries to the rank that owns each row (trivial nB/size partitioning),
    // and optionally a second redistribution based on owned_dof_ids.
    // =========================================================================

#ifdef MOAB_HAVE_EIGEN3

    typedef Eigen::Triplet< double > Triplet;
    std::vector< Triplet > tripletList;

#ifdef MOAB_HAVE_MPI
    if( size > 1 )
    {
        // Trivial row partitioning for redistribution
        const int nPerPart = nB / size;

        moab::TupleList* tl = new moab::TupleList;
        unsigned numr       = 1;
        tl->initialize( 3, 0, 0, numr, localSize );  // to_proc, row, col, value
        tl->enableWriteAccess();

        for( int i = 0; i < localSize; i++ )
        {
            int rowval  = vecRow[i] - 1;  // convert from 1-based (SCRIP) to 0-based
            int colval  = vecCol[i] - 1;
            int to_proc = rowval / nPerPart;
            if( to_proc >= size ) to_proc = size - 1;

            int n                = tl->get_n();
            tl->vi_wr[3 * n]     = to_proc;
            tl->vi_wr[3 * n + 1] = rowval;
            tl->vi_wr[3 * n + 2] = colval;
            tl->vr_wr[n]         = vecS[i];
            tl->inc_n();
        }

        // Crystal router: redistribute entries by row ownership
        ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, *tl, 0 );

        if( owned_dof_ids.size() > 0 )
        {
            // we need to send desired dof to the rendezvous point
            moab::TupleList tl_re;                                 //
            tl_re.initialize( 2, 0, 0, 0, owned_dof_ids.size() );  // to proc, value
            tl_re.enableWriteAccess();
            // send first to rendez_vous point, decided by trivial partitioning

            for( size_t i = 0; i < owned_dof_ids.size(); i++ )
            {
                int to_proc = -1;
                int dof_val = owned_dof_ids[i] - 1;  // dofs are 1 based in the file, partition from 0 ?
                to_proc     = dof_val / nPerPart;
                if( to_proc == size ) to_proc = size - 1;

                int n                  = tl_re.get_n();
                tl_re.vi_wr[2 * n]     = to_proc;
                tl_re.vi_wr[2 * n + 1] = dof_val;

                tl_re.inc_n();
            }
            ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, tl_re, 0 );
            // now we know in tl_re where do we need to send back dof_val
            moab::TupleList::buffer sort_buffer;
            sort_buffer.buffer_init( tl_re.get_n() );
            tl_re.sort( 1, &sort_buffer );  // so now we order by value

            //sort_buffer.buffer_init( tl->get_n() );

            std::map< int, int > startDofIndex, endDofIndex;  // indices in tl_re for values we want
            int dofVal = -1;
            if( tl_re.get_n() > 0 )
            {
                dofVal = tl_re.vi_rd[1];  // first dof val on this rank  tl_re.vi_rd[2 * 0 + 1];

                startDofIndex[dofVal] = 0;
                endDofIndex[dofVal]   = 0;  // start and end
                for( unsigned k = 1; k < tl_re.get_n(); k++ )
                {
                    int newDof = tl_re.vi_rd[2 * k + 1];
                    if( dofVal == newDof )
                    {
                        endDofIndex[dofVal] = k;  // increment by 1 actually
                    }
                    else
                    {
                        dofVal                = newDof;
                        startDofIndex[dofVal] = k;
                        endDofIndex[dofVal]   = k;
                    }
                }
            }
            // basically, for each value we are interested in, index in tl_re with those values are
            // tl_re.vi_rd[2*startDofIndex+1] == valDof == tl_re.vi_rd[2*endDofIndex+1]
            // so now we have ordered
            // tl_re shows to what proc do we need to send the tuple (row, col, val)
            moab::TupleList* tl_back = new moab::TupleList;
            unsigned numr            = 1;  //
            // localSize is a good guess, but maybe it should be bigger ?
            // this could be bigger for repeated dofs
            tl_back->initialize( 3, 0, 0, numr, tl->get_n() );  // to proc, row, col, value
            tl_back->enableWriteAccess();
            // now loop over tl and tl_re to see where to send
            // form the new tuple, which will contain the desired dofs per task, per row or column distribution

            for( unsigned k = 0; k < tl->get_n(); k++ )
            {
                int valDof = tl->vi_rd[3 * k + 1];  // 1 for row, 2 for column // first value, it should be
                if( startDofIndex.find( valDof ) == startDofIndex.end() ) continue;
                for( int ire = startDofIndex[valDof]; ire <= endDofIndex[valDof]; ire++ )
                {
                    int to_proc               = tl_re.vi_rd[2 * ire];
                    int n                     = tl_back->get_n();
                    tl_back->vi_wr[3 * n]     = to_proc;
                    tl_back->vi_wr[3 * n + 1] = tl->vi_rd[3 * k + 1];  // row
                    tl_back->vi_wr[3 * n + 2] = tl->vi_rd[3 * k + 2];  // col
                    tl_back->vr_wr[n]         = tl->vr_rd[k];
                    tl_back->inc_n();
                }
            }

            // now communicate to the desired tasks:
            ( m_pcomm->proc_config().crystal_router() )->gs_transfer( 1, *tl_back, 0 );

            tl_re.reset();  // clear memory, although this will go out of scope
            tl->reset();
            tl = tl_back;
        }

        // set of row and col used on this task
        std::set< int > rowSet;
        std::set< int > colSet;
        // populate the sparsematrix, using rowMap and colMap
        int n = tl->get_n();
        for( int i = 0; i < n; i++ )
        {
            const int vecRowValue = tl->vi_wr[3 * i + 1];
            const int vecColValue = tl->vi_wr[3 * i + 2];
            rowSet.insert( vecRowValue );
            colSet.insert( vecColValue );
        }
        int index = 0;
        row_gdofmap.resize( rowSet.size() );
        for( auto setIt : rowSet )
        {
            row_gdofmap[index] = setIt;
            rowMap[setIt]      = index++;
        }
        m_nTotDofs_Dest = index;
        index           = 0;
        col_gdofmap.resize( colSet.size() );
        for( auto setIt : colSet )
        {
            col_gdofmap[index] = setIt;
            colMap[setIt]      = index++;
        }
        m_nTotDofs_SrcCov = index;

        tripletList.reserve( n );
        for( int i = 0; i < n; i++ )
        {
            const int vecRowValue = tl->vi_wr[3 * i + 1];
            const int vecColValue = tl->vi_wr[3 * i + 2];
            double value          = tl->vr_wr[i];
            tripletList.emplace_back( rowMap[vecRowValue], colMap[vecColValue], value );
        }
        tl->reset();
    }
    else
#endif
    {
        // set of row and col used on this task
        std::set< int > rowSet;
        std::set< int > colSet;
        // populate the sparsematrix, using rowMap and colMap
        for( int i = 0; i < nS; i++ )
        {
            const int vecRowValue = vecRow[i] - 1;
            const int vecColValue = vecCol[i] - 1;
            rowSet.insert( vecRowValue );
            colSet.insert( vecColValue );
        }

        int index = 0;
        row_gdofmap.resize( rowSet.size() );
        for( auto setIt : rowSet )
        {
            row_gdofmap[index] = setIt;
            rowMap[setIt]      = index++;
        }
        m_nTotDofs_Dest = index;
        index           = 0;
        col_gdofmap.resize( colSet.size() );
        for( auto setIt : colSet )
        {
            col_gdofmap[index] = setIt;
            colMap[setIt]      = index++;
        }
        m_nTotDofs_SrcCov = index;

        tripletList.reserve( nS );
        for( int i = 0; i < nS; i++ )
        {
            const int vecRowValue = vecRow[i] - 1;  // the rows, cols are 1 based in the file
            const int vecColValue = vecCol[i] - 1;  // sparse matrix will be 0 based
            double value          = vecS[i];
            tripletList.emplace_back( rowMap[vecRowValue], colMap[vecColValue], value );
        }
    }

    m_weightMatrix.resize( m_nTotDofs_Dest, m_nTotDofs_SrcCov );
    m_rowVector.resize( m_nTotDofs_Dest );
    m_colVector.resize( m_nTotDofs_SrcCov );
    m_nTotDofs_Src = m_nTotDofs_SrcCov;  // do we need both?
    m_weightMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );
    // Reset the source and target data first
    m_rowVector.setZero();
    m_colVector.setZero();
#ifdef VERBOSE
    serializeSparseMatrix( m_weightMatrix, "map_operator_" + std::to_string( rank ) + ".txt" );
#endif
// #ifdef MOAB_HAVE_EIGEN3
#endif
    // TODO: make this flexible and read the order from map with help of metadata
    m_nDofsPEl_Src  = 1;  // always assume FV-FV maps are read from file
    m_nDofsPEl_Dest = 1;  // always assume FV-FV maps are read from file

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
