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

#ifdef MOAB_HAVE_NETCDF
#ifdef MOAB_HAVE_NETCDFPAR
#include "netcdfcpp_par.hpp"
#else
#include "netcdfcpp.h"
#endif
#endif

#ifdef MOAB_HAVE_PNETCDF
#include <pnetcdf.h>


#endif

#if defined( MOAB_HAVE_NETCDF ) || defined( MOAB_HAVE_PNETCDF )
// Central NC I/O: read/write SCRIP maps through the runtime dispatch layer, which picks the
// serial-NetCDF / parallel-NetCDF / PnetCDF backend from the detected on-disk format.
#include "MBNcDispatch.hpp"
#define ERR_MBNC( err, msg )                                                      \
    do                                                                            \
    {                                                                             \
        int _mbrc = ( err );                                                      \
        if( _mbrc != NC_NOERR ) { _EXCEPTION1( "MBNcDispatch error: %s", msg ); } \
    } while( 0 )
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
#if !defined( MOAB_HAVE_NETCDFPAR ) && !defined( MOAB_HAVE_PNETCDF )
        // Without a parallel SCRIP backend (parallel NetCDF or PnetCDF), the SCRIP writer
        // cannot handle multiple MPI ranks writing to the same file.
        if( this->size > 1 )
        {
#if defined( MOAB_HAVE_HDF5 )
            // Fall back to the HDF5 format with a .h5m extension; the map can be
            // converted to SCRIP format offline if needed.
            std::string h5mFilename = strFilename.substr( 0, lastindex ) + ".h5m";
            if( !this->rank )
            {
                std::cout << "  [WriteParallelMap]: Parallel NetCDF/PnetCDF not available; writing map to "
                          << "HDF5 format (" << h5mFilename << ") instead of SCRIP (.nc)\n";
            }
            MB_CHK_ERR( this->WriteHDF5MapFile( h5mFilename.c_str() ) );
            return moab::MB_SUCCESS;
#else
            MB_CHK_SET_ERR( moab::MB_FAILURE,
                            "Parallel SCRIP write requires NETCDFPAR or PnetCDF; HDF5 fallback unavailable" );
#endif
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
#if !defined( MOAB_HAVE_NETCDF ) && !defined( MOAB_HAVE_PNETCDF )
#error "Cannot enable SCRIP writing without NetCDF or PNetCDF interfaces"
#endif
    // The SCRIP map is written below through the MBNcDispatch layer (mbnc_*), which selects
    // the NetCDF / parallel-NetCDF / PnetCDF backend at runtime from the target format. The
    // file is created only after every buffer is computed (classic define-mode -> data-mode),
    // so nothing is opened here.

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

#if defined( MOAB_HAVE_MPI )
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

    // Write SparseMatrix entries (backend-agnostic: fills vecRow/vecCol/vecS and
    // computes the fractional-coverage arrays dFracA/dFracB via crystal-router comm)
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
    for( std::set< int >::iterator sit = neededRows.begin(); sit != neededRows.end(); ++sit )
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
    // ============ Write the SCRIP map through the MBNcDispatch layer (mbnc_*) ============
    // One code path for every backend (serial NetCDF, parallel NetCDF, PnetCDF). The map is
    // written as classic CDF-5, which both libnetcdf and PnetCDF read. All dimensions,
    // variables and attributes are defined first (classic define-mode), then written
    // collectively with the per-rank hyperslab offsets computed above.
    {
        const int mapFormat  = NCFMT_CLASSIC;
        NcBackend wbackend   = mbnc_choose_backend_for_write( mapFormat, (int)size );
        if( wbackend == NCB_NONE )
            _EXCEPTION1( "No NetCDF backend available to write SCRIP map \"%s\"", strFilename.c_str() );

        int ncid        = -1;
        const int cmode = NC_CLOBBER | NC_64BIT_DATA;  // CDF-5
#ifdef MOAB_HAVE_MPI
        ERR_MBNC( mbnc_create_par( wbackend, m_pcomm->comm(), MPI_INFO_NULL, strFilename.c_str(), cmode, &ncid ),
                  "create map file" );
#else
        ERR_MBNC( mbnc_create( strFilename.c_str(), cmode, &ncid ), "create map file" );
#endif

        // ---- global attributes ----
        for( std::map< std::string, std::string >::const_iterator ait = attrMap.begin(); ait != attrMap.end();
             ++ait )
            ERR_MBNC( mbnc_put_att_text( ncid, NC_GLOBAL, ait->first.c_str(), ait->second.size(), ait->second.c_str() ),
                      "global attribute" );

        // ---- dimensions ----
        int dimSrcRank, dimDstRank, dimNAp, dimNBp, dimNVAp, dimNVBp, dimNSp;
        ERR_MBNC( mbnc_def_dim( ncid, "src_grid_rank", nSrcGridDims, &dimSrcRank ), "def src_grid_rank" );
        ERR_MBNC( mbnc_def_dim( ncid, "dst_grid_rank", nDstGridDims, &dimDstRank ), "def dst_grid_rank" );
        ERR_MBNC( mbnc_def_dim( ncid, "n_a", (size_t)globuf[0], &dimNAp ), "def n_a" );
        ERR_MBNC( mbnc_def_dim( ncid, "n_b", (size_t)globuf[1], &dimNBp ), "def n_b" );
        ERR_MBNC( mbnc_def_dim( ncid, "nv_a", (size_t)globuf[3], &dimNVAp ), "def nv_a" );
        ERR_MBNC( mbnc_def_dim( ncid, "nv_b", (size_t)globuf[4], &dimNVBp ), "def nv_b" );
        ERR_MBNC( mbnc_def_dim( ncid, "n_s", (size_t)globuf[2], &dimNSp ), "def n_s" );

        // ---- variables ----
        int vSrcGridDims, vDstGridDims, vYCA, vYCB, vXCA, vXCB, vYVA, vYVB, vXVA, vXVB, vMaskA, vMaskB, vAreaA, vAreaB,
            vRow, vCol, vS, vFracA, vFracB;
        int d1[1], d2[2];
        d1[0] = dimSrcRank;
        ERR_MBNC( mbnc_def_var( ncid, "src_grid_dims", NC_INT, 1, d1, &vSrcGridDims ), "def src_grid_dims" );
        d1[0] = dimDstRank;
        ERR_MBNC( mbnc_def_var( ncid, "dst_grid_dims", NC_INT, 1, d1, &vDstGridDims ), "def dst_grid_dims" );
        d1[0] = dimNAp;
        ERR_MBNC( mbnc_def_var( ncid, "yc_a", NC_DOUBLE, 1, d1, &vYCA ), "def yc_a" );
        ERR_MBNC( mbnc_def_var( ncid, "xc_a", NC_DOUBLE, 1, d1, &vXCA ), "def xc_a" );
        ERR_MBNC( mbnc_def_var( ncid, "mask_a", NC_INT, 1, d1, &vMaskA ), "def mask_a" );
        ERR_MBNC( mbnc_def_var( ncid, "area_a", NC_DOUBLE, 1, d1, &vAreaA ), "def area_a" );
        ERR_MBNC( mbnc_def_var( ncid, "frac_a", NC_DOUBLE, 1, d1, &vFracA ), "def frac_a" );
        d1[0] = dimNBp;
        ERR_MBNC( mbnc_def_var( ncid, "yc_b", NC_DOUBLE, 1, d1, &vYCB ), "def yc_b" );
        ERR_MBNC( mbnc_def_var( ncid, "xc_b", NC_DOUBLE, 1, d1, &vXCB ), "def xc_b" );
        ERR_MBNC( mbnc_def_var( ncid, "mask_b", NC_INT, 1, d1, &vMaskB ), "def mask_b" );
        ERR_MBNC( mbnc_def_var( ncid, "area_b", NC_DOUBLE, 1, d1, &vAreaB ), "def area_b" );
        ERR_MBNC( mbnc_def_var( ncid, "frac_b", NC_DOUBLE, 1, d1, &vFracB ), "def frac_b" );
        d2[0] = dimNAp;
        d2[1] = dimNVAp;
        ERR_MBNC( mbnc_def_var( ncid, "yv_a", NC_DOUBLE, 2, d2, &vYVA ), "def yv_a" );
        ERR_MBNC( mbnc_def_var( ncid, "xv_a", NC_DOUBLE, 2, d2, &vXVA ), "def xv_a" );
        d2[0] = dimNBp;
        d2[1] = dimNVBp;
        ERR_MBNC( mbnc_def_var( ncid, "yv_b", NC_DOUBLE, 2, d2, &vYVB ), "def yv_b" );
        ERR_MBNC( mbnc_def_var( ncid, "xv_b", NC_DOUBLE, 2, d2, &vXVB ), "def xv_b" );
        d1[0] = dimNSp;
        ERR_MBNC( mbnc_def_var( ncid, "row", NC_INT, 1, d1, &vRow ), "def row" );
        ERR_MBNC( mbnc_def_var( ncid, "col", NC_INT, 1, d1, &vCol ), "def col" );
        ERR_MBNC( mbnc_def_var( ncid, "S", NC_DOUBLE, 1, d1, &vS ), "def S" );

        // ---- variable attributes (reversed grid_dims name ordering preserved) ----
        {
            char szDim[64];
            for( unsigned i = 0; i < srcdimSizes.size(); i++ )
            {
                snprintf( szDim, 64, "name%u", i );
                const std::string& nm = srcdimNames[nSrcGridDims - i - 1];
                ERR_MBNC( mbnc_put_att_text( ncid, vSrcGridDims, szDim, nm.size(), nm.c_str() ), "src_grid_dims name" );
            }
            for( unsigned i = 0; i < tgtdimSizes.size(); i++ )
            {
                snprintf( szDim, 64, "name%u", i );
                const std::string& nm = tgtdimNames[nDstGridDims - i - 1];
                ERR_MBNC( mbnc_put_att_text( ncid, vDstGridDims, szDim, nm.size(), nm.c_str() ), "dst_grid_dims name" );
            }
            const std::string deg( "degrees" );
            const int vdeg[8] = { vYCA, vYCB, vXCA, vXCB, vYVA, vYVB, vXVA, vXVB };
            for( int k = 0; k < 8; k++ )
                ERR_MBNC( mbnc_put_att_text( ncid, vdeg[k], "units", deg.size(), deg.c_str() ), "units" );
            const std::string faName( "fraction of target coverage of source dof" );
            const std::string fbName( "fraction of source coverage of target dof" );
            const std::string unitless( "unitless" );
            ERR_MBNC( mbnc_put_att_text( ncid, vFracA, "name", faName.size(), faName.c_str() ), "frac_a name" );
            ERR_MBNC( mbnc_put_att_text( ncid, vFracA, "units", unitless.size(), unitless.c_str() ), "frac_a units" );
            ERR_MBNC( mbnc_put_att_text( ncid, vFracB, "name", fbName.size(), fbName.c_str() ), "frac_b name" );
            ERR_MBNC( mbnc_put_att_text( ncid, vFracB, "units", unitless.size(), unitless.c_str() ), "frac_b units" );
        }

        ERR_MBNC( mbnc_enddef( ncid ), "enddef" );

        // ---- collective data writes (every rank participates) ----
        size_t sA = (size_t)offbuf[0], cA = (size_t)nA;
        size_t sB = (size_t)offbuf[1], cB = (size_t)nB;
        size_t sS = (size_t)offbuf[2], cS = (size_t)nS;
        double* pYCA   = ( nA > 0 ) ? &dSourceCenterLat[0] : NULL;
        double* pXCA   = ( nA > 0 ) ? &dSourceCenterLon[0] : NULL;
        double* pYCB   = ( nB > 0 ) ? &dTargetCenterLat[0] : NULL;
        double* pXCB   = ( nB > 0 ) ? &dTargetCenterLon[0] : NULL;
        int* pMaskA    = ( nA > 0 ) ? &masksA[0] : NULL;
        int* pMaskB    = ( nB > 0 ) ? &masksB[0] : NULL;
        double* pAreaA = ( nA > 0 ) ? &vecSourceFaceArea[0] : NULL;
        double* pAreaB = ( nB > 0 ) ? &vecTargetFaceArea[0] : NULL;
        double* pFracA = ( nA > 0 ) ? &dFracA[0] : NULL;
        double* pFracB = ( nB > 0 ) ? &dFracB[0] : NULL;
        int* pRow      = ( nS > 0 ) ? &vecRow[0] : NULL;
        int* pCol      = ( nS > 0 ) ? &vecCol[0] : NULL;
        double* pS     = ( nS > 0 ) ? &vecS[0] : NULL;

        ERR_MBNC( mbnc_put_vara_double( ncid, vYCA, &sA, &cA, pYCA ), "put yc_a" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vXCA, &sA, &cA, pXCA ), "put xc_a" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vYCB, &sB, &cB, pYCB ), "put yc_b" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vXCB, &sB, &cB, pXCB ), "put xc_b" );
        ERR_MBNC( mbnc_put_vara_int( ncid, vMaskA, &sA, &cA, pMaskA ), "put mask_a" );
        ERR_MBNC( mbnc_put_vara_int( ncid, vMaskB, &sB, &cB, pMaskB ), "put mask_b" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vAreaA, &sA, &cA, pAreaA ), "put area_a" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vAreaB, &sB, &cB, pAreaB ), "put area_b" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vFracA, &sA, &cA, pFracA ), "put frac_a" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vFracB, &sB, &cB, pFracB ), "put frac_b" );
        ERR_MBNC( mbnc_put_vara_int( ncid, vRow, &sS, &cS, pRow ), "put row" );
        ERR_MBNC( mbnc_put_vara_int( ncid, vCol, &sS, &cS, pCol ), "put col" );
        ERR_MBNC( mbnc_put_vara_double( ncid, vS, &sS, &cS, pS ), "put S" );
        {
            size_t s2A[2] = { (size_t)offbuf[0], 0 }, c2A[2] = { (size_t)nA, (size_t)nSourceNodesPerFace };
            double* pYVA  = ( nA > 0 ) ? &dSourceVertexLat[0][0] : NULL;
            double* pXVA  = ( nA > 0 ) ? &dSourceVertexLon[0][0] : NULL;
            ERR_MBNC( mbnc_put_vara_double( ncid, vYVA, s2A, c2A, pYVA ), "put yv_a" );
            ERR_MBNC( mbnc_put_vara_double( ncid, vXVA, s2A, c2A, pXVA ), "put xv_a" );
            size_t s2B[2] = { (size_t)offbuf[1], 0 }, c2B[2] = { (size_t)nB, (size_t)nTargetNodesPerFace };
            double* pYVB  = ( nB > 0 ) ? &dTargetVertexLat[0][0] : NULL;
            double* pXVB  = ( nB > 0 ) ? &dTargetVertexLon[0][0] : NULL;
            ERR_MBNC( mbnc_put_vara_double( ncid, vYVB, s2B, c2B, pYVB ), "put yv_b" );
            ERR_MBNC( mbnc_put_vara_double( ncid, vXVB, s2B, c2B, pXVB ), "put xv_b" );
        }
        {
            // Small global grid_dims arrays: rank 0 writes the full array, others write 0.
            size_t sg    = 0;
            size_t cgSrc = ( rank == 0 ) ? (size_t)nSrcGridDims : 0;
            size_t cgDst = ( rank == 0 ) ? (size_t)nDstGridDims : 0;
            ERR_MBNC( mbnc_put_vara_int( ncid, vSrcGridDims, &sg, &cgSrc, ( rank == 0 ) ? &srcdimSizes[0] : NULL ),
                      "put src_grid_dims" );
            ERR_MBNC( mbnc_put_vara_int( ncid, vDstGridDims, &sg, &cgDst, ( rank == 0 ) ? &tgtdimSizes[0] : NULL ),
                      "put dst_grid_dims" );
        }

        ERR_MBNC( mbnc_close( ncid ), "close" );
    }



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

///////////////////////////////////////////////////////////////////////////////
//
// ReadParallelMap: read a SCRIP-format map file and distribute the sparse matrix
// across MPI ranks. All NetCDF I/O goes through the MBNcDispatch layer (mbnc_*),
// which selects the backend (serial / parallel NetCDF / PnetCDF / buffered) from the
// detected on-disk format; each rank reads a contiguous stripe and the entries are
// then redistributed to their owners (owned_dof_ids) before Eigen assembly.
///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::ReadParallelMap( const char* strSource,
                                                         const std::vector< int >& owned_dof_ids,
                                                         int arearead,
                                                         std::vector< double >& vecAreaA,
                                                         int& nA,
                                                         std::vector< double >& vecAreaB,
                                                         int& nB )
{
#if !defined( MOAB_HAVE_NETCDF ) && !defined( MOAB_HAVE_PNETCDF )
#error "Cannot enable SCRIP reading without NetCDF or PNetCDF interfaces"
#endif

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

    // ============ Phase 1: read the map through the MBNcDispatch layer (mbnc_*) ============
    // Detect the on-disk format, let the dispatch pick the backend (serial / parallel NetCDF /
    // PnetCDF / rank-0-buffered), and read a contiguous stripe of the sparse matrix on each
    // rank. row/col/S stripes are redistributed to their owners in Phase 2; area_a/area_b are
    // read as trivial per-rank slices, which is what the downstream aream code expects.
    {
        int fileFormat = NCFMT_UNKNOWN;
#ifdef MOAB_HAVE_MPI
        if( rank == 0 ) fileFormat = mbnc_detect_format( strSource );
        MPI_Bcast( &fileFormat, 1, MPI_INT, 0, m_pcomm->comm() );
#else
        fileFormat = mbnc_detect_format( strSource );
#endif
        NcBackend rbackend = mbnc_choose_backend_for_read( fileFormat, (int)size );
        if( rbackend == NCB_NONE )
            _EXCEPTION1( "Cannot read map file \"%s\": unrecognized format, or NetCDF-4/HDF5 without libnetcdf",
                         strSource );

        int ncid = -1;
#ifdef MOAB_HAVE_MPI
        ERR_MBNC( mbnc_open_par( rbackend, m_pcomm->comm(), MPI_INFO_NULL, strSource, 0, &ncid ), "open map" );
#else
        ERR_MBNC( mbnc_open( strSource, 0, &ncid ), "open map" );
#endif

        int did     = -1;
        size_t dlen = 0;
        ERR_MBNC( mbnc_inq_dimid( ncid, "n_a", &did ), "inq n_a" );
        ERR_MBNC( mbnc_inq_dimlen( ncid, did, &dlen ), "len n_a" );
        nA = (int)dlen;
        ERR_MBNC( mbnc_inq_dimid( ncid, "n_b", &did ), "inq n_b" );
        ERR_MBNC( mbnc_inq_dimlen( ncid, did, &dlen ), "len n_b" );
        nB = (int)dlen;
        ERR_MBNC( mbnc_inq_dimid( ncid, "n_s", &did ), "inq n_s" );
        ERR_MBNC( mbnc_inq_dimlen( ncid, did, &dlen ), "len n_s" );
        nS = (int)dlen;

        // Contiguous per-rank stripes (last rank takes the remainder).
        localSize          = nS / size;
        size_t offsetRead  = (size_t)rank * (size_t)localSize;
        if( rank == size - 1 ) localSize += nS % size;
        int localSizeA     = nA / size;
        size_t offsetReadA = (size_t)rank * (size_t)localSizeA;
        if( rank == size - 1 ) localSizeA += nA % size;
        int localSizeB     = nB / size;
        size_t offsetReadB = (size_t)rank * (size_t)localSizeB;
        if( rank == size - 1 ) localSizeB += nB % size;

        vecRow.resize( localSize );
        vecCol.resize( localSize );
        vecS.resize( localSize );

        int vid   = -1;
        size_t st = offsetRead, ct = (size_t)localSize;
        ERR_MBNC( mbnc_inq_varid( ncid, "row", &vid ), "inq row" );
        ERR_MBNC( mbnc_get_vara_int( ncid, vid, &st, &ct, localSize ? vecRow.data() : NULL ), "get row" );
        ERR_MBNC( mbnc_inq_varid( ncid, "col", &vid ), "inq col" );
        ERR_MBNC( mbnc_get_vara_int( ncid, vid, &st, &ct, localSize ? vecCol.data() : NULL ), "get col" );
        ERR_MBNC( mbnc_inq_varid( ncid, "S", &vid ), "inq S" );
        ERR_MBNC( mbnc_get_vara_double( ncid, vid, &st, &ct, localSize ? vecS.data() : NULL ), "get S" );

        if( readAreaA )
        {
            vecAreaA.resize( localSizeA );
            size_t sa = offsetReadA, ca = (size_t)localSizeA;
            ERR_MBNC( mbnc_inq_varid( ncid, "area_a", &vid ), "inq area_a" );
            ERR_MBNC( mbnc_get_vara_double( ncid, vid, &sa, &ca, localSizeA ? vecAreaA.data() : NULL ), "get area_a" );
        }
        if( readAreaB )
        {
            vecAreaB.resize( localSizeB );
            size_t sb = offsetReadB, cb = (size_t)localSizeB;
            ERR_MBNC( mbnc_inq_varid( ncid, "area_b", &vid ), "inq area_b" );
            ERR_MBNC( mbnc_get_vara_double( ncid, vid, &sb, &cb, localSizeB ? vecAreaB.data() : NULL ), "get area_b" );
        }

        ERR_MBNC( mbnc_close( ncid ), "close" );
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
    // Preserve the map file's global source-DoF count (n_a) so the migration
    // can tell a masked source mesh (fewer cells than n_a -> drop is BfB-safe)
    // from a complete one (== n_a but a column missing -> real error).
    m_nTotDofs_SrcGlobal = nA;
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
