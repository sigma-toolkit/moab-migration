///////////////////////////////////////////////////////////////////////////////
///
/// \file    TempestLinearRemap.cpp
/// \author  Vijay Mahadevan
/// \version Mar 08, 2017
///

#ifdef WIN32               /* windows */
#define _USE_MATH_DEFINES  // For M_PI
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include "Announce.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
#include "FiniteVolumeTools.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MathHelper.h"
#include "SparseMatrix.h"
#include "OverlapMesh.h"

#include "DebugOutput.hpp"
#include "moab/AdaptiveKDTree.hpp"

#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/TupleList.hpp"
#include "moab/MeshTopoUtil.hpp"

#pragma GCC diagnostic pop

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <unordered_set>

// #define VERBOSE
#define USE_ComputeAdjacencyRelations

/// <summary>
///     Face index and distance metric pair.
/// </summary>
typedef std::pair< int, int > FaceDistancePair;

/// <summary>
///     Vector storing adjacent Faces.
/// </summary>
typedef std::vector< FaceDistancePair > AdjacentFaceVector;

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::LinearRemapNN_MOAB( bool use_GID_matching, bool strict_check )
{
    /* m_mapRemap size = (m_nTotDofs_Dest X m_nTotDofs_SrcCov)  */

#ifdef VVERBOSE
    {
        std::ofstream output_file( "rowcolindices.txt", std::ios::out );
        output_file << m_nTotDofs_Dest << " " << m_nTotDofs_SrcCov << " " << row_gdofmap.size() << " "
                    << row_ldofmap.size() << " " << col_gdofmap.size() << " " << col_ldofmap.size() << "\n";
        output_file << "Rows \n";
        for( unsigned iv = 0; iv < row_gdofmap.size(); iv++ )
            output_file << row_gdofmap[iv] << " " << row_dofmap[iv] << "\n";
        output_file << "Cols \n";
        for( unsigned iv = 0; iv < col_gdofmap.size(); iv++ )
            output_file << col_gdofmap[iv] << " " << col_dofmap[iv] << "\n";
        output_file.flush();  // required here
        output_file.close();
    }
#endif

    if( use_GID_matching )
    {
        std::map< unsigned, unsigned > src_gl;
        for( unsigned it = 0; it < col_gdofmap.size(); ++it )
            src_gl[col_gdofmap[it]] = it;

        std::map< unsigned, unsigned >::iterator iter;
        for( unsigned it = 0; it < row_gdofmap.size(); ++it )
        {
            unsigned row = row_gdofmap[it];
            iter         = src_gl.find( row );
            if( strict_check && iter == src_gl.end() )
            {
                std::cout << "Searching for global target DOF " << row
                          << " but could not find correspondence in source mesh.\n";
                assert( false );
            }
            else if( iter == src_gl.end() )
            {
                continue;
            }
            else
            {
                unsigned icol = src_gl[row];
                unsigned irow = it;

                // Set the permutation matrix in local space
                m_mapRemap( irow, icol ) = 1.0;
            }
        }

        return moab::MB_SUCCESS;
    }
    else
    {
        /* Create a Kd-tree to perform local queries to find nearest neighbors */

        return moab::MB_FAILURE;
    }
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapFVtoFV_Tempest_MOAB( int nOrder )
{
    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Verify ReverseNodeArray has been calculated
    if( m_meshInputCov->faces.size() > 0 && m_meshInputCov->revnodearray.size() == 0 )
    {
        _EXCEPTIONT( "ReverseNodeArray has not been calculated for m_meshInputCov" );
    }

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( TriQuadRuleOrder );

    // Number of coefficients needed at this order
#ifdef RECTANGULAR_TRUNCATION
    int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION
    int nCoefficients = nOrder * ( nOrder + 1 ) / 2;
#endif

    // Number of faces you need
    const int nRequiredFaceSetSize = nCoefficients;

    // Fit weight exponent
    const int nFitWeightsExponent = nOrder + 2;

    // Announcements
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapFVtoFV_Tempest_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Volume to Finite Volume Projection\n" );
        dbgprint.printf( 0, "Triangular quadrature rule order %i\n", TriQuadRuleOrder );
        dbgprint.printf( 0, "Number of coefficients: %i\n", nCoefficients );
        dbgprint.printf( 0, "Required adjacency set size: %i\n", nRequiredFaceSetSize );
        dbgprint.printf( 0, "Fit weights exponent: %i\n", nFitWeightsExponent );
    }

    // Current overlap face
    int ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    DataArray2D< double > dIntArray;
    DataArray1D< double > dConstraint( nCoefficients );

    // Loop through all faces on m_meshInputCov
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        // Output every 1000 elements
#ifdef VERBOSE
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Find the set of Faces that overlap faceFirst
        int ixOverlapBegin    = ixOverlap;
        unsigned ixOverlapEnd = ixOverlapBegin;

        for( ; ixOverlapEnd < m_meshOverlap->faces.size(); ixOverlapEnd++ )
        {
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapEnd] != 0 ) break;
        }

        unsigned nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

        if( nOverlapFaces == 0 ) continue;

        // Build integration array
        BuildIntegrationArray( *m_meshInputCov, *m_meshOverlap, triquadrule, ixFirst, ixOverlapBegin, ixOverlapEnd,
                               nOrder, dIntArray );

        // Set of Faces to use in building the reconstruction and associated
        // distance metric.
        AdjacentFaceVector vecAdjFaces;

        GetAdjacentFaceVectorByEdge( *m_meshInputCov, ixFirst, nRequiredFaceSetSize, vecAdjFaces );

        // Number of adjacent Faces
        int nAdjFaces = vecAdjFaces.size();

        // Determine the conservative constraint equation
        double dFirstArea = m_meshInputCov->vecFaceArea[ixFirst];
        dConstraint.Zero();
        for( int p = 0; p < nCoefficients; p++ )
        {
            for( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                dConstraint[p] += dIntArray[p][j];
            }
            dConstraint[p] /= dFirstArea;
        }

        // Build the fit array from the integration operator
        DataArray2D< double > dFitArray;
        DataArray1D< double > dFitWeights;
        DataArray2D< double > dFitArrayPlus;

        BuildFitArray( *m_meshInputCov, triquadrule, ixFirst, vecAdjFaces, nOrder, nFitWeightsExponent, dConstraint,
                       dFitArray, dFitWeights );

        // Compute the inverse fit array
        bool fSuccess = InvertFitArray_Corrected( dConstraint, dFitArray, dFitWeights, dFitArrayPlus );

        // Multiply integration array and fit array
        DataArray2D< double > dComposedArray( nAdjFaces, nOverlapFaces );
        if( fSuccess )
        {
            // Multiply integration array and inverse fit array
            for( int i = 0; i < nAdjFaces; i++ )
            {
                for( size_t j = 0; j < nOverlapFaces; j++ )
                {
                    for( int k = 0; k < nCoefficients; k++ )
                    {
                        dComposedArray( i, j ) += dIntArray( k, j ) * dFitArrayPlus( i, k );
                    }
                }
            }

            // Unable to invert fit array, drop to 1st order.  In this case
            // dFitArrayPlus(0,0) = 1 and all other entries are zero.
        }
        else
        {
            dComposedArray.Zero();
            for( size_t j = 0; j < nOverlapFaces; j++ )
            {
                dComposedArray( 0, j ) += dIntArray( 0, j );
            }
        }

        // Put composed array into map
        for( unsigned i = 0; i < vecAdjFaces.size(); i++ )
        {
            for( unsigned j = 0; j < nOverlapFaces; j++ )
            {
                int& ixFirstFaceLoc  = vecAdjFaces[i].first;
                int& ixSecondFaceLoc = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];
                // int ixFirstFaceGlob = m_remapper->GetGlobalID(moab::Remapper::SourceMesh,
                // ixFirstFaceLoc); int ixSecondFaceGlob =
                // m_remapper->GetGlobalID(moab::Remapper::TargetMesh, ixSecondFaceLoc);

                // signal to not participate, because it is a ghost target
                if( ixSecondFaceLoc < 0 ) continue;  // do not do anything

                m_mapRemap( ixSecondFaceLoc, ixFirstFaceLoc ) +=
                    dComposedArray[i][j] / m_meshOutput->vecFaceArea[ixSecondFaceLoc];
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::PrintMapStatistics()
{
    int nrows = m_weightMatrix.rows();      // Number of rows
    int ncols = m_weightMatrix.cols();      // Number of columns
    int NNZ   = m_weightMatrix.nonZeros();  // Number of non zero values
#ifdef MOAB_HAVE_MPI
    // find out min/max for NNZ, ncols, nrows
    // should work on std c++ 11
    int arr3[6] = { NNZ, nrows, ncols, -NNZ, -nrows, -ncols };
    int rarr3[6];
    MPI_Reduce( arr3, rarr3, 6, MPI_INT, MPI_MIN, 0, m_pcomm->comm() );

    int total[3];
    MPI_Reduce( arr3, total, 3, MPI_INT, MPI_SUM, 0, m_pcomm->comm() );
    if( !rank )
        std::cout << "-> Rows (min/max/sum): (" << rarr3[1] << " / " << -rarr3[4] << " / " << total[1] << "), "
                  << " Cols (min/max/sum): (" << rarr3[2] << " / " << -rarr3[5] << " / " << total[2] << "), "
                  << " NNZ (min/max/sum): (" << rarr3[0] << " / " << -rarr3[3] << " / " << total[0] << ")\n";
#else
    std::cout << "-> Rows: " << nrows << ", Cols: " << ncols << ", NNZ: " << NNZ << "\n";
#endif
}

#ifdef MOAB_HAVE_EIGEN3
void moab::TempestOnlineMap::copy_tempest_sparsemat_to_eigen3()
{
#ifndef VERBOSE
#define VERBOSE_ACTIVATED
// #define VERBOSE
#endif

    /* Should the columns be the global size of the matrix ? */
    m_weightMatrix.resize( m_nTotDofs_Dest, m_nTotDofs_SrcCov );
    m_rowVector.resize( m_weightMatrix.rows() );
    m_colVector.resize( m_weightMatrix.cols() );

#ifdef VERBOSE
    int locrows = std::max( m_mapRemap.GetRows(), m_nTotDofs_Dest );
    int loccols = std::max( m_mapRemap.GetColumns(), m_nTotDofs_SrcCov );

    std::cout << m_weightMatrix.rows() << ", " << locrows << ", " << m_weightMatrix.cols() << ", " << loccols << "\n";
    // assert(m_weightMatrix.rows() == locrows && m_weightMatrix.cols() == loccols);
#endif

    DataArray1D< int > lrows;
    DataArray1D< int > lcols;
    DataArray1D< double > lvals;
    m_mapRemap.GetEntries( lrows, lcols, lvals );
    size_t locvals = lvals.GetRows();

    // first matrix
    typedef Eigen::Triplet< double > Triplet;
    std::vector< Triplet > tripletList;
    tripletList.reserve( locvals );
    for( size_t iv = 0; iv < locvals; iv++ )
    {
        tripletList.push_back( Triplet( lrows[iv], lcols[iv], lvals[iv] ) );
    }
    m_weightMatrix.setFromTriplets( tripletList.begin(), tripletList.end() );
    m_weightMatrix.makeCompressed();

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "tempestmatrix.txt.0000" << rank;
    std::ofstream output_file( sstr.str(), std::ios::out );
    output_file << "0 " << locrows << " 0 " << loccols << "\n";
    for( unsigned iv = 0; iv < locvals; iv++ )
    {
        // output_file << lrows[iv] << " " << row_ldofmap[lrows[iv]] << " " <<
        // row_gdofmap[row_ldofmap[lrows[iv]]] << " " << col_gdofmap[col_ldofmap[lcols[iv]]] << " "
        // << lvals[iv] << "\n";
        output_file << row_gdofmap[row_ldofmap[lrows[iv]]] << " " << col_gdofmap[col_ldofmap[lcols[iv]]] << " "
                    << lvals[iv] << "\n";
    }
    output_file.flush();  // required here
    output_file.close();
#endif

#ifdef VERBOSE_ACTIVATED
#undef VERBOSE_ACTIVATED
#undef VERBOSE
#endif
    return;
}
#endif

///////////////////////////////////////////////////////////////////////////////

template < typename T >
static std::vector< size_t > sort_indexes( const std::vector< T >& v )
{
    // initialize original index locations
    std::vector< size_t > idx( v.size() );
    std::iota( idx.begin(), idx.end(), 0 );

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort( idx.begin(), idx.end(), [&v]( size_t i1, size_t i2 ) { return fabs( v[i1] ) > fabs( v[i2] ); } );

    return idx;
}

double moab::TempestOnlineMap::QLTLimiter( int caasIteration,
                                           std::vector< double >& dataCorrectedField,
                                           std::vector< double >& dataLowerBound,
                                           std::vector< double >& dataUpperBound,
                                           std::vector< double >& dMassDefect )
{
    const size_t nrows = dataCorrectedField.size();
    double dMassL      = 0.0;
    double dMassU      = 0.0;
    std::vector< double > dataCorrection( nrows );
    double dMassDiffCum                       = 0.0;
    double dLMinusU                           = fabs( dataUpperBound[0] - dataLowerBound[0] );
    const DataArray1D< double >& dTargetAreas = this->m_remapper->m_target->vecFaceArea;

    // std::vector< size_t > sortedIdx = sort_indexes( dMassDefect );
    std::vector< std::unordered_set< int > > vecAdjTargetFaces( nrows );
    constexpr bool useMOABAdjacencies = true;
#ifdef USE_ComputeAdjacencyRelations
    if( useMOABAdjacencies )
        ComputeAdjacencyRelations( vecAdjTargetFaces, caasIteration, m_remapper->m_target_entities,
                                   useMOABAdjacencies );
    else
        ComputeAdjacencyRelations( vecAdjTargetFaces, caasIteration, m_remapper->m_target_entities, useMOABAdjacencies,
                                   this->m_remapper->m_target );
#else
    moab::MeshTopoUtil mtu( m_interface );
    ;
#endif

    for( size_t i = 0; i < nrows; i++ )
    {
        // size_t index = sortedIdx[i];
        size_t index          = i;
        dataCorrection[index] = fmax( dataLowerBound[index], fmin( dataUpperBound[index], 0.0 ) );
        // dMassDiff[index] = dMassDefect[index] - dTargetAreas[index] * dataCorrection[index];
        // dMassDiff[index] = dMassDefect[index];

        dMassL += dTargetAreas[index] * dataLowerBound[index];
        dMassU += dTargetAreas[index] * dataUpperBound[index];
        dLMinusU = fmax( dLMinusU, fabs( dataUpperBound[index] - dataLowerBound[index] ) );
        dMassDiffCum += dMassDefect[index] - dTargetAreas[index] * dataCorrection[index];

#ifndef USE_ComputeAdjacencyRelations
        vecAdjTargetFaces[index].insert( index );  // add self target face first
        {
            // Compute the adjacent faces to the target face
            if( useMOABAdjacencies )
            {
                moab::Range ents;
                // ents.insert( m_remapper->m_target_entities.index( m_remapper->m_target_entities[index] ) );
                ents.insert( m_remapper->m_target_entities[index] );
                moab::Range adjEnts;
                moab::ErrorCode rval = mtu.get_bridge_adjacencies( ents, 0, 2, adjEnts, caasIteration );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
                for( moab::Range::iterator it = adjEnts.begin(); it != adjEnts.end(); ++it )
                {
                    // int adjIndex = m_interface->id_from_handle(*it)-1;
                    int adjIndex = m_remapper->m_target_entities.index( *it );
                    // printf("rank: %d, Element %lu, entity: %lu, adjIndex %d\n", rank, index, *it, adjIndex);
                    if( adjIndex >= 0 ) vecAdjTargetFaces[index].insert( adjIndex );
                }
            }
            else
            {
                AdjacentFaceVector vecAdjFaces;
                GetAdjacentFaceVectorByEdge( *this->m_remapper->m_target, index,
                                             ( m_output_order + 1 ) * ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             //  ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             //  ( 4 ) * ( m_output_order + 1 ) * ( m_output_order + 1 ),
                                             vecAdjFaces );

                // Add the adjacent faces to the target face list
                for( auto adjFace : vecAdjFaces )
                    if( adjFace.first >= 0 )
                        vecAdjTargetFaces[index].insert( adjFace.first );  // map target face to source face
            }
        }
#endif
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localDefects( 5, 0.0 ), globalDefects( 5, 0.0 );
    localDefects[0] = dMassL;
    localDefects[1] = dMassU;
    localDefects[2] = dMassDiffCum;
    localDefects[3] = dLMinusU;
    // localDefects[4] = dMassCorrectU;

    MPI_Allreduce( localDefects.data(), globalDefects.data(), 4, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );

    dMassL       = globalDefects[0];
    dMassU       = globalDefects[1];
    dMassDiffCum = globalDefects[2];
    dLMinusU     = globalDefects[3];
    // dMassCorrectU = globalDefects[4];
#endif

    //If the upper and lower bounds are too close together, just clip
    if( fabs( dMassDiffCum ) < 1e-15 || dLMinusU < 1e-15 )
    {
        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
        return dMassDiffCum;
    }
    else
    {
        if( dMassL > dMassDiffCum )
        {
            Announce( "Lower bound mass exceeds target mass by %1.15e: CAAS will need another iteration",
                      dMassL - dMassDiffCum );
            dMassDiffCum = dMassL;
            // dMass -= dMassL;
        }
        else if( dMassU < dMassDiffCum )
        {
            Announce( "Target mass exceeds upper bound mass by %1.15e: CAAS will need another iteration",
                      dMassDiffCum - dMassU );
            dMassDiffCum = dMassU;
            // dMass -= dMassU;
        }

        // TODO: optimize away dataMassVec by a simple transient double within the loop
        // DataArray1D< double > dataMassVec( nrows );  //vector of mass redistribution
        for( size_t i = 0; i < nrows; i++ )
        {
            // size_t index   = sortedIdx[i];
            size_t index                               = i;
            const std::unordered_set< int >& neighbors = vecAdjTargetFaces[index];
            if( dMassDefect[index] > 0.0 )
            {
                double dMassCorrectU = 0.0;
                for( auto it : neighbors )
                    dMassCorrectU += dTargetAreas[it] * ( dataUpperBound[it] - dataCorrection[it] );

                // double dMassDiffCumOld = dMassDefect[index];
                for( auto it : neighbors )
                    dataCorrection[it] +=
                        dMassDefect[index] * ( dataUpperBound[it] - dataCorrection[it] ) / dMassCorrectU;
            }
            else
            {
                double dMassCorrectL = 0.0;
                for( auto it : neighbors )
                    dMassCorrectL += dTargetAreas[it] * ( dataCorrection[it] - dataLowerBound[it] );

                // double dMassDiffCumOld = dMassDefect[index];
                for( auto it : neighbors )
                    dataCorrection[it] +=
                        dMassDefect[index] * ( dataCorrection[it] - dataLowerBound[it] ) / dMassCorrectL;
            }
        }

        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
    }

    return dMassDiffCum;
}

void moab::TempestOnlineMap::CAASLimiter( std::vector< double >& dataCorrectedField,
                                          std::vector< double >& dataLowerBound,
                                          std::vector< double >& dataUpperBound,
                                          double& dMass )
{
    const size_t nrows = dataCorrectedField.size();
    double dMassL      = 0.0;
    double dMassU      = 0.0;
    std::vector< double > dataCorrection( nrows );
    const DataArray1D< double >& dTargetAreas = this->m_remapper->m_target->vecFaceArea;
    double dMassDiff                          = dMass;
    double dLMinusU                           = fabs( dataUpperBound[0] - dataLowerBound[0] );
    double dMassCorrectU                      = 0.0;
    double dMassCorrectL                      = 0.0;
    for( size_t i = 0; i < nrows; i++ )
    {
        dataCorrection[i] = fmax( dataLowerBound[i], fmin( dataUpperBound[i], 0.0 ) );
        dMassL += dTargetAreas[i] * dataLowerBound[i];
        dMassU += dTargetAreas[i] * dataUpperBound[i];
        dMassDiff -= dTargetAreas[i] * dataCorrection[i];
        dLMinusU = fmax( dLMinusU, fabs( dataUpperBound[i] - dataLowerBound[i] ) );
        dMassCorrectL += dTargetAreas[i] * ( dataCorrection[i] - dataLowerBound[i] );
        dMassCorrectU += dTargetAreas[i] * ( dataUpperBound[i] - dataCorrection[i] );
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localDefects( 5, 0.0 ), globalDefects( 5, 0.0 );
    localDefects[0] = dMassL;
    localDefects[1] = dMassU;
    localDefects[2] = dMassDiff;
    localDefects[3] = dMassCorrectL;
    localDefects[4] = dMassCorrectU;

    MPI_Allreduce( localDefects.data(), globalDefects.data(), 5, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );

    dMassL        = globalDefects[0];
    dMassU        = globalDefects[1];
    dMassDiff     = globalDefects[2];
    dMassCorrectL = globalDefects[3];
    dMassCorrectU = globalDefects[4];
#endif

    //If the upper and lower bounds are too close together, just clip
    if( fabs( dMassDiff ) < 1e-15 || fabs( dLMinusU ) < 1e-15 )
    {
        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
        return;
    }
    else
    {
        if( dMassL > dMassDiff )
        {
            Announce( "%d: Lower bound mass exceeds target mass by %1.15e: CAAS will need another iteration", rank,
                      dMassL - dMassDiff );
            dMassDiff = dMassL;
            dMass -= dMassL;
        }
        else if( dMassU < dMassDiff )
        {
            Announce( "%d: Target mass exceeds upper bound mass by %1.15e: CAAS will need another iteration", rank,
                      dMassDiff - dMassU );
            dMassDiff = dMassU;
            dMass -= dMassU;
        }

        // TODO: optimize away dataMassVec by a simple transient double within the loop
        DataArray1D< double > dataMassVec( nrows );  //vector of mass redistribution
        if( dMassDiff > 0.0 )
        {
            for( size_t i = 0; i < nrows; i++ )
            {
                dataMassVec[i] = ( dataUpperBound[i] - dataCorrection[i] ) / dMassCorrectU;
                dataCorrection[i] += dMassDiff * dataMassVec[i];
            }
        }
        else
        {
            for( size_t i = 0; i < nrows; i++ )
            {
                dataMassVec[i] = ( dataCorrection[i] - dataLowerBound[i] ) / dMassCorrectL;
                dataCorrection[i] += dMassDiff * dataMassVec[i];
            }
        }

        for( size_t i = 0; i < nrows; i++ )
            dataCorrectedField[i] += dataCorrection[i];
    }

    return;
}

std::pair< double, double > moab::TempestOnlineMap::ApplyBoundsLimiting( std::vector< double >& dataInDouble,
                                                                         std::vector< double >& dataOutDouble,
                                                                         CAASType caasType,
                                                                         int caasIteration,
                                                                         double mismatch )
{
    // Currently only implemented for FV to FV remapping
    // We should generalize this to other types of remapping
    assert( !dataGLLNodesSrcCov.IsAttached() && !dataGLLNodesDest.IsAttached() );

    std::pair< double, double > massDefect( 0.0, 0.0 );

    // Check if the source and target data are of the same size
    const size_t nTargetCount                    = dataOutDouble.size();
    const DataArray1D< double >& m_dOverlapAreas = this->m_remapper->m_overlap->vecFaceArea;

    // Apply the offline map to the data
    double dMassDiff = 0.0;
    std::vector< double > x( nTargetCount );
    std::vector< double > dataLowerBound( nTargetCount );
    std::vector< double > dataUpperBound( nTargetCount );
    std::vector< double > massVector( nTargetCount );
    std::vector< std::unordered_set< int > > vecSourceOvTarget( nTargetCount );

#undef USE_ComputeAdjacencyRelations
    constexpr bool useMOABAdjacencies = true;
#ifdef USE_ComputeAdjacencyRelations
    // Compute the adjacent faces to the source face
    // However, calling MOAB to do this does not work correctly as we need ixS to be the index
    // Cannot just iterate over all entities in the source covering mesh
    if( caasType == CAAS_QLT || caasType == CAAS_LOCAL_ADJACENT )
    {
        if( useMOABAdjacencies )
        {
            moab::ErrorCode rval =
                ComputeAdjacencyRelations( vecSourceOvTarget, caasIteration, m_remapper->m_covering_source_entities,
                                           useMOABAdjacencies );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
        }
        else
        {
            moab::ErrorCode rval =
                ComputeAdjacencyRelations( vecSourceOvTarget, caasIteration, m_remapper->m_covering_source_entities,
                                           useMOABAdjacencies, m_meshInputCov );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
        }
    }
#else
    moab::MeshTopoUtil mtu( m_interface );
#endif

    // Initialize the bounds on the given source and target data
    double dSourceMin = dataInDouble[0];
    double dSourceMax = dataInDouble[0];
    double dTargetMin = dataOutDouble[0];
    double dTargetMax = dataOutDouble[0];
    for( size_t i = 0; i < m_meshOverlap->faces.size(); i++ )
    {
        const int ixS = m_meshOverlap->vecSourceFaceIx[i];
        const int ixT = m_meshOverlap->vecTargetFaceIx[i];

        if( ixT < 0 ) continue;  // skip ghost target faces

        assert( m_dOverlapAreas[i] > 0.0 );
        assert( ixS >= 0 );
        assert( ixT >= 0 );

#ifndef USE_ComputeAdjacencyRelations
        // Compute the adjacent faces to the target face
        vecSourceOvTarget[ixT].insert( ixS );  // map target face to source face
        if( ( caasType == CAAS_QLT || caasType == CAAS_LOCAL_ADJACENT ) )
        {
            if( useMOABAdjacencies )
            {
                moab::Range ents;
                ents.insert( m_remapper->m_covering_source_entities[ixS] );
                moab::Range adjEnts;
                moab::ErrorCode rval = mtu.get_bridge_adjacencies( ents, 0, 2, adjEnts, caasIteration );MB_CHK_SET_ERR_CONT( rval, "Failed to get adjacent faces" );
                for( moab::Range::iterator it = adjEnts.begin(); it != adjEnts.end(); ++it )
                {
                    int adjIndex = m_remapper->m_covering_source_entities.index( *it );
                    if( adjIndex >= 0 ) vecSourceOvTarget[ixT].insert( adjIndex );
                }
            }
            else
            {
                // Compute the adjacent faces to the target face
                AdjacentFaceVector vecAdjFaces;
                GetAdjacentFaceVectorByEdge( *m_meshInputCov, ixS,
                                             ( caasIteration ) * ( m_input_order + 1 ) * ( m_input_order + 1 ),
                                             vecAdjFaces );

                //Compute min/max over neighboring faces
                for( size_t iadj = 0; iadj < vecAdjFaces.size(); iadj++ )
                    vecSourceOvTarget[ixT].insert( vecAdjFaces[iadj].first );  // map target face to source face
            }
        }
#endif

        // Update the min and max values of the source data
        dSourceMax = fmax( dSourceMax, dataInDouble[ixS] );
        dSourceMin = fmin( dSourceMin, dataInDouble[ixS] );

        // Update the min and max values of the target data
        dTargetMin = fmin( dTargetMin, dataOutDouble[ixT] );
        dTargetMax = fmax( dTargetMax, dataOutDouble[ixT] );

        const double locMassDiff = ( dataInDouble[ixS] * m_dOverlapAreas[i] ) -  // source mass
                                   ( dataOutDouble[ixT] * m_dOverlapAreas[i] );  // target mass

        // Update the mass difference between source and target faces
        // linked to the overlap mesh element
        dMassDiff += locMassDiff;  // target mass
        massVector[ixT] += locMassDiff;
    }

#ifdef MOAB_HAVE_MPI
    std::vector< double > localMinMaxDefects( 5, 0.0 ), globalMinMaxDefects( 5, 0.0 );
    localMinMaxDefects[0] = dSourceMin;
    localMinMaxDefects[1] = dTargetMin;
    localMinMaxDefects[2] = dSourceMax;
    localMinMaxDefects[3] = dTargetMax;
    localMinMaxDefects[4] = dMassDiff;

    if( caasType == CAAS_GLOBAL )
    {
        MPI_Allreduce( localMinMaxDefects.data(), globalMinMaxDefects.data(), 2, MPI_DOUBLE, MPI_MIN, m_pcomm->comm() );
        MPI_Allreduce( localMinMaxDefects.data() + 2, globalMinMaxDefects.data() + 2, 2, MPI_DOUBLE, MPI_MAX,
                       m_pcomm->comm() );
        dSourceMin = globalMinMaxDefects[0];
        dSourceMax = globalMinMaxDefects[2];
        dTargetMin = globalMinMaxDefects[1];
        dTargetMax = globalMinMaxDefects[3];
    }
    if( caasIteration == 1 )
        MPI_Allreduce( localMinMaxDefects.data() + 4, globalMinMaxDefects.data() + 4, 1, MPI_DOUBLE, MPI_SUM,
                       m_pcomm->comm() );
    else
        globalMinMaxDefects[4] = mismatch;

    dMassDiff = localMinMaxDefects[4];
    // massDefect.first = localMinMaxDefects[4];
    massDefect.first = globalMinMaxDefects[4];
#else

    // massDefect.first = fabs( dMassDiff / ( dSourceMax - dSourceMin ) );
    massDefect.first = dMassDiff;
#endif

    // Early exit if the values are monotone already.
    // if( ( dTargetMax <= dSourceMax && dTargetMin <= dSourceMin ) || fabs( massDefect.first ) < 1e-16 )
    if( fabs( massDefect.first ) > 1e-20 )
    {
        if( caasType == CAAS_GLOBAL )
        {
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = dSourceMin - dataOutDouble[i];
                dataUpperBound[i] = dSourceMax - dataOutDouble[i];
            }
        }  // if( caasType == CAAS_GLOBAL )
        else  // caasType == CAAS_LOCAL
        {
            // Compute the local min and max values of the target data
            std::vector< double > vecLocalUpperBound( nTargetCount );
            std::vector< double > vecLocalLowerBound( nTargetCount );
            // Loop over the target faces and compute the min and max values
            // of the source data linked to the target faces
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                assert( vecSourceOvTarget[i].size() );

                double dMinI = 1E10;   // dataInDouble[vecSourceOvTarget[i][0]];
                double dMaxI = -1E10;  // dataInDouble[vecSourceOvTarget[i][0]];

                // Compute max over intersecting source faces
                for( const auto& srcElem : vecSourceOvTarget[i] )
                {
                    dMinI = fmin( dMinI, dataInDouble[srcElem] );  // min over intersecting source faces
                    dMaxI = fmax( dMaxI, dataInDouble[srcElem] );  // max over intersecting source faces
                }

                // Update the min and max values of the target data
                vecLocalLowerBound[i] = dMinI;
                vecLocalUpperBound[i] = dMaxI;
            }

            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = vecLocalLowerBound[i] - dataOutDouble[i];
                dataUpperBound[i] = vecLocalUpperBound[i] - dataOutDouble[i];
            }
        }  // caasType == CAAS_LOCAL

        // Invoke CAAS or QLT application on the map
        if( fabs( dMassDiff ) > 1e-20 )
        {
            if( caasType == CAAS_QLT )
                dMassDiff = QLTLimiter( caasIteration, dataOutDouble, dataLowerBound, dataUpperBound, massVector );
            else
                CAASLimiter( dataOutDouble, dataLowerBound, dataUpperBound, dMassDiff );
        }

        // Announce output mass
        double dMassDiffPost = 0.0;
        for( size_t i = 0; i < m_meshOverlap->faces.size(); i++ )
        {
            const int ixS = m_meshOverlap->vecSourceFaceIx[i];
            const int ixT = m_meshOverlap->vecTargetFaceIx[i];

            if( ixT < 0 ) continue;  // skip ghost target faces

            // Update the mass difference between source and target faces
            // linked to the overlap mesh element
            dMassDiffPost += ( dataInDouble[ixS] * m_dOverlapAreas[i] ) -  // source mass
                             ( dataOutDouble[ixT] * m_dOverlapAreas[i] );  // target mass
        }
        // massDefect.second = fabs( dMassDiffPost / ( dSourceMax - dSourceMin ) );
        massDefect.second = dMassDiffPost;
    }

    // Ideally should perform an AllReduce here to get the global mass difference across all processors
    // But if we satisfy the constraint on every task, essentially, the global mass difference should be zero!
    return massDefect;
}

///////////////////////////////////////////////////////////////////////////////

extern void ForceConsistencyConservation3( const DataArray1D< double >& vecSourceArea,
                                           const DataArray1D< double >& vecTargetArea,
                                           DataArray2D< double >& dCoeff,
                                           bool fMonotone,
                                           bool fSparseConstraints = false );

///////////////////////////////////////////////////////////////////////////////

extern void ForceIntArrayConsistencyConservation( const DataArray1D< double >& vecSourceArea,
                                                  const DataArray1D< double >& vecTargetArea,
                                                  DataArray2D< double >& dCoeff,
                                                  bool fMonotone );

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapSE4_Tempest_MOAB( const DataArray3D< int >& dataGLLNodes,
                                                          const DataArray3D< double >& dataGLLJacobian,
                                                          int nMonotoneType,
                                                          bool fContinuousIn,
                                                          bool fNoConservation )
{
    // Order of the polynomial interpolant
    int nP = dataGLLNodes.GetRows();

    // Order of triangular quadrature rule
    const int TriQuadRuleOrder = 4;

    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( TriQuadRuleOrder );

    int TriQuadraturePoints = triquadrule.GetPoints();

    const DataArray2D< double >& TriQuadratureG = triquadrule.GetG();

    const DataArray1D< double >& TriQuadratureW = triquadrule.GetW();

    // Sample coefficients
    DataArray2D< double > dSampleCoeff( nP, nP );

    // GLL Quadrature nodes on quadrilateral elements
    DataArray1D< double > dG;
    DataArray1D< double > dW;
    GaussLobattoQuadrature::GetPoints( nP, 0.0, 1.0, dG, dW );

    // Announcements
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapSE4_Tempest_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Volume Projection\n" );
        dbgprint.printf( 0, "Triangular quadrature rule order %i\n", TriQuadRuleOrder );
        dbgprint.printf( 0, "Order of the FE polynomial interpolant: %i\n", nP );
    }

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // NodeVector from m_meshOverlap
    const NodeVector& nodesOverlap = m_meshOverlap->nodes;
    const NodeVector& nodesFirst   = m_meshInputCov->nodes;

    // Vector of source areas
    DataArray1D< double > vecSourceArea( nP * nP );

    DataArray1D< double > vecTargetArea;
    DataArray2D< double > dCoeff;

#ifdef VERBOSE
    std::stringstream sstr;
    sstr << "remapdata_" << rank << ".txt";
    std::ofstream output_file( sstr.str() );
#endif

    // Current Overlap Face
    int ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    // generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

    // Loop over all input Faces
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        if( faceFirst.edges.size() != 4 )
        {
            _EXCEPTIONT( "Only quadrilateral elements allowed for SE remapping" );
        }
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Need to re-number the overlap elements such that vecSourceFaceIx[a:b] = 0, then 1 and so
        // on wrt the input mesh data Then the overlap_end and overlap_begin will be correct.
        // However, the relation with MOAB and Tempest will go out of the roof

        // Determine how many overlap Faces and triangles are present
        int nOverlapFaces    = 0;
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // if( m_meshOverlap->vecTargetFaceIx[ixOverlapTemp] < 0 ) continue;  // skip ghost target faces
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 ) break;

            nOverlapFaces++;
        }

        // No overlaps
        if( nOverlapFaces == 0 ) continue;

        // Allocate remap coefficients array for meshFirst Face
        DataArray3D< double > dRemapCoeff( nP, nP, nOverlapFaces );

        // Find the local remap coefficients
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            const Face& faceOverlap = m_meshOverlap->faces[ixOverlap + j];
            if( m_meshOverlap->vecFaceArea[ixOverlap + j] < 1.e-16 )  // machine precision
            {
                Announce( "Very small overlap at index %i area polygon: (%1.10e )", ixOverlap + j,
                          m_meshOverlap->vecFaceArea[ixOverlap + j] );
                int n = faceOverlap.edges.size();
                Announce( "Number nodes: %d", n );
                for( int k = 0; k < n; k++ )
                {
                    Node nd = nodesOverlap[faceOverlap[k]];
                    Announce( "Node %d  %d  : %1.10e  %1.10e %1.10e ", k, faceOverlap[k], nd.x, nd.y, nd.z );
                }
                continue;
            }

            // #ifdef VERBOSE
            // if ( is_root )
            //     Announce ( "\tLocal ID: %i/%i = %i, areas = %2.8e", j + ixOverlap, nOverlapFaces,
            //     m_remapper->lid_to_gid_covsrc[m_meshOverlap->vecSourceFaceIx[ixOverlap + j]],
            //     m_meshOverlap->vecFaceArea[ixOverlap + j] );
            // #endif

            int nbEdges           = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center;  // not used if nbEdges == 3
            if( nbEdges > 3 )
            {  // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for( int k = 0; k < nbEdges; k++ )
                {
                    const Node& node = nodesOverlap[faceOverlap[k]];
                    center           = center + node;
                }
                center = center / nbEdges;
                center = center.Normalized();  // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriangleArea;

            // Loop over all sub-triangles of this Overlap Face
            for( int k = 0; k < nOverlapTriangles; k++ )
            {
                if( nbEdges == 3 )  // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0         = nodesOverlap[faceOverlap[0]];
                    node1         = nodesOverlap[faceOverlap[1]];
                    node2         = nodesOverlap[faceOverlap[2]];
                    dTriangleArea = CalculateFaceArea( faceOverlap, nodesOverlap );
                }
                else  // decompose polygon in triangles around the center
                {
                    node0         = center;
                    node1         = nodesOverlap[faceOverlap[k]];
                    int k1        = ( k + 1 ) % nbEdges;
                    node2         = nodesOverlap[faceOverlap[k1]];
                    nodes[0]      = center;
                    nodes[1]      = node1;
                    nodes[2]      = node2;
                    dTriangleArea = CalculateFaceArea( faceTri, nodes );
                }
                // Coordinates of quadrature Node
                for( int l = 0; l < TriQuadraturePoints; l++ )
                {
                    Node nodeQuadrature;
                    nodeQuadrature.x = TriQuadratureG[l][0] * node0.x + TriQuadratureG[l][1] * node1.x +
                                       TriQuadratureG[l][2] * node2.x;

                    nodeQuadrature.y = TriQuadratureG[l][0] * node0.y + TriQuadratureG[l][1] * node1.y +
                                       TriQuadratureG[l][2] * node2.y;

                    nodeQuadrature.z = TriQuadratureG[l][0] * node0.z + TriQuadratureG[l][1] * node1.z +
                                       TriQuadratureG[l][2] * node2.z;

                    nodeQuadrature = nodeQuadrature.Normalized();

                    // Find components of quadrature point in basis
                    // of the first Face
                    double dAlpha;
                    double dBeta;

                    ApplyInverseMap( faceFirst, nodesFirst, nodeQuadrature, dAlpha, dBeta );

                    // Check inverse map value
                    if( ( dAlpha < -1.0e-13 ) || ( dAlpha > 1.0 + 1.0e-13 ) || ( dBeta < -1.0e-13 ) ||
                        ( dBeta > 1.0 + 1.0e-13 ) )
                    {
                        _EXCEPTION4( "Inverse Map for element %d and subtriangle %d out of range "
                                     "(%1.5e %1.5e)",
                                     j, l, dAlpha, dBeta );
                    }

                    // Sample the finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nP, dAlpha, dBeta, dSampleCoeff );

                    // Add sample coefficients to the map if m_meshOverlap->vecFaceArea[ixOverlap + j] > 0
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dRemapCoeff[p][q][j] += TriQuadratureW[l] * dTriangleArea * dSampleCoeff[p][q] /
                                                    m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
            }
        }

#ifdef VERBOSE
        output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    output_file << dRemapCoeff[p][q][j] << " ";
                }
            }
        }
        output_file << std::endl;
#endif

        // Force consistency and conservation
        if( !fNoConservation )
        {
            double dTargetArea = 0.0;
            for( int j = 0; j < nOverlapFaces; j++ )
            {
                dTargetArea += m_meshOverlap->vecFaceArea[ixOverlap + j];
            }

            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    vecSourceArea[p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                }
            }

            const double areaTolerance = 1e-10;
            // Source elements are completely covered by target volumes
            if( fabs( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea ) <= areaTolerance )
            {
                vecTargetArea.Allocate( nOverlapFaces );
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }

                dCoeff.Allocate( nOverlapFaces, nP * nP );

                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }

                // Target volumes only partially cover source elements
            }
            else if( m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea > areaTolerance )
            {
                double dExtraneousArea = m_meshInputCov->vecFaceArea[ixFirst] - dTargetArea;

                vecTargetArea.Allocate( nOverlapFaces + 1 );
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    vecTargetArea[j] = m_meshOverlap->vecFaceArea[ixOverlap + j];
                }
                vecTargetArea[nOverlapFaces] = dExtraneousArea;

#ifdef VERBOSE
                Announce( "Partial volume: %i (%1.10e / %1.10e)", ixFirst, dTargetArea,
                          m_meshInputCov->vecFaceArea[ixFirst] );
#endif
                if( dTargetArea > m_meshInputCov->vecFaceArea[ixFirst] )
                {
                    _EXCEPTIONT( "Partial element area exceeds total element area" );
                }

                dCoeff.Allocate( nOverlapFaces + 1, nP * nP );

                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
                        }
                    }
                }
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] = dataGLLJacobian[p][q][ixFirst];
                    }
                }
                for( int j = 0; j < nOverlapFaces; j++ )
                {
                    for( int p = 0; p < nP; p++ )
                    {
                        for( int q = 0; q < nP; q++ )
                        {
                            dCoeff[nOverlapFaces][p * nP + q] -=
                                dRemapCoeff[p][q][j] * m_meshOverlap->vecFaceArea[ixOverlap + j];
                        }
                    }
                }
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dCoeff[nOverlapFaces][p * nP + q] /= dExtraneousArea;
                    }
                }

                // Source elements only partially cover target volumes
            }
            else
            {
                Announce( "Coverage area: %1.10e, and target element area: %1.10e)", ixFirst,
                          m_meshInputCov->vecFaceArea[ixFirst], dTargetArea );
                _EXCEPTIONT( "Target grid must be a subset of source grid" );
            }

            ForceConsistencyConservation3( vecSourceArea, vecTargetArea, dCoeff, ( nMonotoneType > 0 )
                                           /*, m_remapper->lid_to_gid_covsrc[ixFirst]*/ );

            for( int j = 0; j < nOverlapFaces; j++ )
            {
                for( int p = 0; p < nP; p++ )
                {
                    for( int q = 0; q < nP; q++ )
                    {
                        dRemapCoeff[p][q][j] = dCoeff[j][p * nP + q];
                    }
                }
            }
        }

#ifdef VERBOSE
        // output_file << "[" << m_remapper->lid_to_gid_covsrc[ixFirst] << "] \t";
        // for ( int j = 0; j < nOverlapFaces; j++ )
        // {
        //     for ( int p = 0; p < nP; p++ )
        //     {
        //         for ( int q = 0; q < nP; q++ )
        //         {
        //             output_file << dRemapCoeff[p][q][j] << " ";
        //         }
        //     }
        // }
        // output_file << std::endl;
#endif

        // Put these remap coefficients into the SparseMatrix map
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];

            // signal to not participate, because it is a ghost target
            if( ixSecondFace < 0 ) continue;  // do not do anything

            for( int p = 0; p < nP; p++ )
            {
                for( int q = 0; q < nP; q++ )
                {
                    if( fContinuousIn )
                    {
                        int ixFirstNode = dataGLLNodes[p][q][ixFirst] - 1;

                        smatMap( ixSecondFace, ixFirstNode ) += dRemapCoeff[p][q][j] *
                                                                m_meshOverlap->vecFaceArea[ixOverlap + j] /
                                                                m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                    else
                    {
                        int ixFirstNode = ixFirst * nP * nP + p * nP + q;

                        smatMap( ixSecondFace, ixFirstNode ) += dRemapCoeff[p][q][j] *
                                                                m_meshOverlap->vecFaceArea[ixOverlap + j] /
                                                                m_meshOutput->vecFaceArea[ixSecondFace];
                    }
                }
            }
        }
        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }
#ifdef VERBOSE
    output_file.flush();  // required here
    output_file.close();
#endif

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_MOAB( const DataArray3D< int >& dataGLLNodesIn,
                                                        const DataArray3D< double >& dataGLLJacobianIn,
                                                        const DataArray3D< int >& dataGLLNodesOut,
                                                        const DataArray3D< double >& dataGLLJacobianOut,
                                                        const DataArray1D< double >& dataNodalAreaOut,
                                                        int nPin,
                                                        int nPout,
                                                        int nMonotoneType,
                                                        bool fContinuousIn,
                                                        bool fContinuousOut,
                                                        bool fNoConservation )
{
    // Triangular quadrature rule
    TriangularQuadratureRule triquadrule( 8 );

    const DataArray2D< double >& dG = triquadrule.GetG();
    const DataArray1D< double >& dW = triquadrule.GetW();

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D< double > dSampleCoeffIn( nPin, nPin );
    DataArray2D< double > dSampleCoeffOut( nPout, nPout );

    // Announcemnets
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapGLLtoGLL2_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Element Projection\n" );
        dbgprint.printf( 0, "Order of the input FE polynomial interpolant: %i\n", nPin );
        dbgprint.printf( 0, "Order of the output FE polynomial interpolant: %i\n", nPout );
    }

    // Build the integration array for each element on m_meshOverlap
    DataArray3D< double > dGlobalIntArray( nPin * nPin, m_meshOverlap->faces.size(), nPout * nPout );

    // Number of overlap Faces per source Face
    DataArray1D< int > nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        // Determine how many overlap Faces and triangles are present
        int nOverlapFaces    = 0;
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];
            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 )
            {
                break;
            }

            nOverlapFaces++;
        }

        nAllOverlapFaces[ixFirst] = nOverlapFaces;

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Geometric area of each output node
    DataArray2D< double > dGeometricOutputArea( m_meshOutput->faces.size(), nPout * nPout );

    // Area of each overlap element in the output basis
    DataArray2D< double > dOverlapOutputArea( m_meshOverlap->faces.size(), nPout * nPout );

    // Loop through all faces on m_meshInputCov
    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    if( is_root ) dbgprint.printf( 0, "Building conservative distribution maps\n" );

    // generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Quantities from the First Mesh
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector& nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        if( !nOverlapFaces ) continue;

        // // Calculate total element Jacobian
        // double dTotalJacobian = 0.0;
        // for (int s = 0; s < nPin; s++) {
        //     for (int t = 0; t < nPin; t++) {
        //         dTotalJacobian += dataGLLJacobianIn[s][t][ixFirst];
        //     }
        // }

        // Loop through all Overlap Faces
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // Quantities from the overlap Mesh
            const Face& faceOverlap = m_meshOverlap->faces[ixOverlap + i];

            const NodeVector& nodesOverlap = m_meshOverlap->nodes;

            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            // signal to not participate, because it is a ghost target
            if( ixSecond < 0 ) continue;  // do not do anything

            const NodeVector& nodesSecond = m_meshOutput->nodes;

            const Face& faceSecond = m_meshOutput->faces[ixSecond];

            int nbEdges           = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center;  // not used if nbEdges == 3
            if( nbEdges > 3 )
            {  // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for( int k = 0; k < nbEdges; k++ )
                {
                    const Node& node = nodesOverlap[faceOverlap[k]];
                    center           = center + node;
                }
                center = center / nbEdges;
                center = center.Normalized();  // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriArea;

            // Loop over all sub-triangles of this Overlap Face
            for( int j = 0; j < nOverlapTriangles; j++ )
            {
                if( nbEdges == 3 )  // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0    = nodesOverlap[faceOverlap[0]];
                    node1    = nodesOverlap[faceOverlap[1]];
                    node2    = nodesOverlap[faceOverlap[2]];
                    dTriArea = CalculateFaceArea( faceOverlap, nodesOverlap );
                }
                else  // decompose polygon in triangles around the center
                {
                    node0    = center;
                    node1    = nodesOverlap[faceOverlap[j]];
                    int j1   = ( j + 1 ) % nbEdges;
                    node2    = nodesOverlap[faceOverlap[j1]];
                    nodes[0] = center;
                    nodes[1] = node1;
                    nodes[2] = node2;
                    dTriArea = CalculateFaceArea( faceTri, nodes );
                }

                for( int k = 0; k < triquadrule.GetPoints(); k++ )
                {
                    // Get the nodal location of this point
                    double dX[3];

                    dX[0] = dG( k, 0 ) * node0.x + dG( k, 1 ) * node1.x + dG( k, 2 ) * node2.x;
                    dX[1] = dG( k, 0 ) * node0.y + dG( k, 1 ) * node1.y + dG( k, 2 ) * node2.y;
                    dX[2] = dG( k, 0 ) * node0.z + dG( k, 1 ) * node1.z + dG( k, 2 ) * node2.z;

                    double dMag = sqrt( dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2] );

                    dX[0] /= dMag;
                    dX[1] /= dMag;
                    dX[2] /= dMag;

                    Node nodeQuadrature( dX[0], dX[1], dX[2] );

                    // Find the components of this quadrature point in the basis
                    // of the first Face.
                    double dAlphaIn;
                    double dBetaIn;

                    ApplyInverseMap( faceFirst, nodesFirst, nodeQuadrature, dAlphaIn, dBetaIn );

                    // Find the components of this quadrature point in the basis
                    // of the second Face.
                    double dAlphaOut;
                    double dBetaOut;

                    ApplyInverseMap( faceSecond, nodesSecond, nodeQuadrature, dAlphaOut, dBetaOut );

                    /*
                                        // Check inverse map value
                                        if ((dAlphaIn < 0.0) || (dAlphaIn > 1.0) ||
                                            (dBetaIn  < 0.0) || (dBetaIn  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaIn, dBetaIn);
                                        }

                                        // Check inverse map value
                                        if ((dAlphaOut < 0.0) || (dAlphaOut > 1.0) ||
                                            (dBetaOut  < 0.0) || (dBetaOut  > 1.0)
                                        ) {
                                            _EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
                                                dAlphaOut, dBetaOut);
                                        }
                    */
                    // Sample the First finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPin, dAlphaIn, dBetaIn, dSampleCoeffIn );

                    // Sample the Second finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPout, dAlphaOut, dBetaOut, dSampleCoeffOut );

                    // Overlap output area
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            double dNodeArea = dSampleCoeffOut[s][t] * dW[k] * dTriArea;

                            dOverlapOutputArea[ixOverlap + i][s * nPout + t] += dNodeArea;

                            dGeometricOutputArea[ixSecond][s * nPout + t] += dNodeArea;
                        }
                    }

                    // Compute overlap integral
                    int ixp = 0;
                    for( int p = 0; p < nPin; p++ )
                    {
                        for( int q = 0; q < nPin; q++ )
                        {
                            int ixs = 0;
                            for( int s = 0; s < nPout; s++ )
                            {
                                for( int t = 0; t < nPout; t++ )
                                {
                                    // Sample the Second finite element at this point
                                    dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
                                        dSampleCoeffOut[s][t] * dSampleCoeffIn[p][q] * dW[k] * dTriArea;

                                    ixs++;
                                }
                            }

                            ixp++;
                        }
                    }
                }
            }
        }

        // Coefficients
        DataArray2D< double > dCoeff( nOverlapFaces * nPout * nPout, nPin * nPin );

        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            dCoeff[i * nPout * nPout + ixs][ixp] = dGlobalIntArray[ixp][ixOverlap + i][ixs] /
                                                                   dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Source areas
        DataArray1D< double > vecSourceArea( nPin * nPin );

        for( int p = 0; p < nPin; p++ )
        {
            for( int q = 0; q < nPin; q++ )
            {
                vecSourceArea[p * nPin + q] = dataGLLJacobianIn[p][q][ixFirst];
            }
        }

        // Target areas
        DataArray1D< double > vecTargetArea( nOverlapFaces * nPout * nPout );

        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];
            int ixs = 0;
            for( int s = 0; s < nPout; s++ )
            {
                for( int t = 0; t < nPout; t++ )
                {
                    vecTargetArea[i * nPout * nPout + ixs] = dOverlapOutputArea[ixOverlap + i][nPout * s + t];

                    ixs++;
                }
            }
        }

        // Force consistency and conservation
        if( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation( vecSourceArea, vecTargetArea, dCoeff, ( nMonotoneType != 0 ) );
        }

        // Update global coefficients
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            dGlobalIntArray[ixp][ixOverlap + i][ixs] =
                                dCoeff[i * nPout * nPout + ixs][ixp] * dOverlapOutputArea[ixOverlap + i][s * nPout + t];

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

#ifdef VVERBOSE
        // Check column sums (conservation)
        for( int i = 0; i < nPin * nPin; i++ )
        {
            double dColSum = 0.0;
            for( int j = 0; j < nOverlapFaces * nPout * nPout; j++ )
            {
                dColSum += dCoeff[j][i] * vecTargetArea[j];
            }
            printf( "Col %i: %1.15e\n", i, dColSum / vecSourceArea[i] );
        }

        // Check row sums (consistency)
        for( int j = 0; j < nOverlapFaces * nPout * nPout; j++ )
        {
            double dRowSum = 0.0;
            for( int i = 0; i < nPin * nPin; i++ )
            {
                dRowSum += dCoeff[j][i];
            }
            printf( "Row %i: %1.15e\n", j, dRowSum );
        }
#endif

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Build redistribution map within target element
    if( is_root ) dbgprint.printf( 0, "Building redistribution maps on target mesh\n" );
    DataArray1D< double > dRedistSourceArea( nPout * nPout );
    DataArray1D< double > dRedistTargetArea( nPout * nPout );
    std::vector< DataArray2D< double > > dRedistributionMaps;
    dRedistributionMaps.resize( m_meshOutput->faces.size() );

    for( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        dRedistributionMaps[ixSecond].Allocate( nPout * nPout, nPout * nPout );

        for( int i = 0; i < nPout * nPout; i++ )
        {
            dRedistributionMaps[ixSecond][i][i] = 1.0;
        }

        for( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistSourceArea[s] = dGeometricOutputArea[ixSecond][s];
        }

        for( int s = 0; s < nPout * nPout; s++ )
        {
            dRedistTargetArea[s] = dataGLLJacobianOut[s / nPout][s % nPout][ixSecond];
        }

        if( !fNoConservation )
        {
            ForceIntArrayConsistencyConservation( dRedistSourceArea, dRedistTargetArea, dRedistributionMaps[ixSecond],
                                                  ( nMonotoneType != 0 ) );

            for( int s = 0; s < nPout * nPout; s++ )
            {
                for( int t = 0; t < nPout * nPout; t++ )
                {
                    dRedistributionMaps[ixSecond][s][t] *= dRedistTargetArea[s] / dRedistSourceArea[t];
                }
            }
        }
    }

    // Construct the total geometric area
    DataArray1D< double > dTotalGeometricArea( dataNodalAreaOut.GetRows() );
    for( size_t ixSecond = 0; ixSecond < m_meshOutput->faces.size(); ixSecond++ )
    {
        for( int s = 0; s < nPout; s++ )
        {
            for( int t = 0; t < nPout; t++ )
            {
                dTotalGeometricArea[dataGLLNodesOut[s][t][ixSecond] - 1] +=
                    dGeometricOutputArea[ixSecond][s * nPout + t];
            }
        }
    }

    // Compose the integration operator with the output map
    ixOverlap = 0;

    if( is_root ) dbgprint.printf( 0, "Assembling map\n" );

    // Map from source DOFs to target DOFs with redistribution applied
    DataArray2D< double > dRedistributedOp( nPin * nPin, nPout * nPout );

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        if( !nOverlapFaces ) continue;

        // Put composed array into map
        for( int j = 0; j < nOverlapFaces; j++ )
        {
            int ixSecondFace = m_meshOverlap->vecTargetFaceIx[ixOverlap + j];

            // signal to not participate, because it is a ghost target
            if( ixSecondFace < 0 ) continue;  // do not do anything

            dRedistributedOp.Zero();
            for( int p = 0; p < nPin * nPin; p++ )
            {
                for( int s = 0; s < nPout * nPout; s++ )
                {
                    for( int t = 0; t < nPout * nPout; t++ )
                    {
                        dRedistributedOp[p][s] +=
                            dRedistributionMaps[ixSecondFace][s][t] * dGlobalIntArray[p][ixOverlap + j][t];
                    }
                }
            }

            int ixp = 0;
            for( int p = 0; p < nPin; p++ )
            {
                for( int q = 0; q < nPin; q++ )
                {
                    int ixFirstNode;
                    if( fContinuousIn )
                    {
                        ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                    }
                    else
                    {
                        ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
                    }

                    int ixs = 0;
                    for( int s = 0; s < nPout; s++ )
                    {
                        for( int t = 0; t < nPout; t++ )
                        {
                            int ixSecondNode;
                            if( fContinuousOut )
                            {
                                ixSecondNode = dataGLLNodesOut[s][t][ixSecondFace] - 1;

                                if( !fNoConservation )
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dataNodalAreaOut[ixSecondNode];
                                }
                                else
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dTotalGeometricArea[ixSecondNode];
                                }
                            }
                            else
                            {
                                ixSecondNode = ixSecondFace * nPout * nPout + s * nPout + t;

                                if( !fNoConservation )
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dataGLLJacobianOut[s][t][ixSecondFace];
                                }
                                else
                                {
                                    smatMap( ixSecondNode, ixFirstNode ) +=
                                        dRedistributedOp[ixp][ixs] / dGeometricOutputArea[ixSecondFace][s * nPout + t];
                                }
                            }

                            ixs++;
                        }
                    }

                    ixp++;
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void moab::TempestOnlineMap::LinearRemapGLLtoGLL2_Pointwise_MOAB( const DataArray3D< int >& dataGLLNodesIn,
                                                                  const DataArray3D< double >& /*dataGLLJacobianIn*/,
                                                                  const DataArray3D< int >& dataGLLNodesOut,
                                                                  const DataArray3D< double >& /*dataGLLJacobianOut*/,
                                                                  const DataArray1D< double >& dataNodalAreaOut,
                                                                  int nPin,
                                                                  int nPout,
                                                                  int nMonotoneType,
                                                                  bool fContinuousIn,
                                                                  bool fContinuousOut )
{
    // Gauss-Lobatto quadrature within Faces
    DataArray1D< double > dGL;
    DataArray1D< double > dWL;

    GaussLobattoQuadrature::GetPoints( nPout, 0.0, 1.0, dGL, dWL );

    // Get SparseMatrix represntation of the OfflineMap
    SparseMatrix< double >& smatMap = this->GetSparseMatrix();

    // Sample coefficients
    DataArray2D< double > dSampleCoeffIn( nPin, nPin );

    // Announcemnets
    moab::DebugOutput dbgprint( std::cout, this->rank, 0 );
    dbgprint.set_prefix( "[LinearRemapGLLtoGLL2_Pointwise_MOAB]: " );
    if( is_root )
    {
        dbgprint.printf( 0, "Finite Element to Finite Element (Pointwise) Projection\n" );
        dbgprint.printf( 0, "Order of the input FE polynomial interpolant: %i\n", nPin );
        dbgprint.printf( 0, "Order of the output FE polynomial interpolant: %i\n", nPout );
    }

    // Number of overlap Faces per source Face
    DataArray1D< int > nAllOverlapFaces( m_meshInputCov->faces.size() );

    int ixOverlap = 0;

    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
        size_t ixOverlapTemp = ixOverlap;
        for( ; ixOverlapTemp < m_meshOverlap->faces.size(); ixOverlapTemp++ )
        {
            // const Face & faceOverlap = m_meshOverlap->faces[ixOverlapTemp];

            if( ixFirst - m_meshOverlap->vecSourceFaceIx[ixOverlapTemp] != 0 ) break;

            nAllOverlapFaces[ixFirst]++;
        }

        // Increment the current overlap index
        ixOverlap += nAllOverlapFaces[ixFirst];
    }

    // Number of times this point was found
    DataArray1D< bool > fSecondNodeFound( dataNodalAreaOut.GetRows() );

    ixOverlap = 0;
#ifdef VERBOSE
    const unsigned outputFrequency = ( m_meshInputCov->faces.size() / 10 ) + 1;
#endif
    // Loop through all faces on m_meshInputCov
    for( size_t ixFirst = 0; ixFirst < m_meshInputCov->faces.size(); ixFirst++ )
    {
#ifdef VERBOSE
        // Announce computation progress
        if( ixFirst % outputFrequency == 0 && is_root )
        {
            dbgprint.printf( 0, "Element %zu/%lu\n", ixFirst, m_meshInputCov->faces.size() );
        }
#endif
        // Quantities from the First Mesh
        const Face& faceFirst = m_meshInputCov->faces[ixFirst];

        const NodeVector& nodesFirst = m_meshInputCov->nodes;

        // Number of overlapping Faces and triangles
        int nOverlapFaces = nAllOverlapFaces[ixFirst];

        // Loop through all Overlap Faces
        for( int i = 0; i < nOverlapFaces; i++ )
        {
            // Quantities from the Second Mesh
            int ixSecond = m_meshOverlap->vecTargetFaceIx[ixOverlap + i];

            // signal to not participate, because it is a ghost target
            if( ixSecond < 0 ) continue;  // do not do anything

            const NodeVector& nodesSecond = m_meshOutput->nodes;
            const Face& faceSecond        = m_meshOutput->faces[ixSecond];

            // Loop through all nodes on the second face
            for( int s = 0; s < nPout; s++ )
            {
                for( int t = 0; t < nPout; t++ )
                {
                    size_t ixSecondNode;
                    if( fContinuousOut )
                    {
                        ixSecondNode = dataGLLNodesOut[s][t][ixSecond] - 1;
                    }
                    else
                    {
                        ixSecondNode = ixSecond * nPout * nPout + s * nPout + t;
                    }

                    if( ixSecondNode >= fSecondNodeFound.GetRows() ) _EXCEPTIONT( "Logic error" );

                    // Check if this node has been found already
                    if( fSecondNodeFound[ixSecondNode] ) continue;

                    // Check this node
                    Node node;
                    Node dDx1G;
                    Node dDx2G;

                    ApplyLocalMap( faceSecond, nodesSecond, dGL[t], dGL[s], node, dDx1G, dDx2G );

                    // Find the components of this quadrature point in the basis
                    // of the first Face.
                    double dAlphaIn;
                    double dBetaIn;

                    ApplyInverseMap( faceFirst, nodesFirst, node, dAlphaIn, dBetaIn );

                    // Check if this node is within the first Face
                    if( ( dAlphaIn < -1.0e-10 ) || ( dAlphaIn > 1.0 + 1.0e-10 ) || ( dBetaIn < -1.0e-10 ) ||
                        ( dBetaIn > 1.0 + 1.0e-10 ) )
                        continue;

                    // Node is within the overlap region, mark as found
                    fSecondNodeFound[ixSecondNode] = true;

                    // Sample the First finite element at this point
                    SampleGLLFiniteElement( nMonotoneType, nPin, dAlphaIn, dBetaIn, dSampleCoeffIn );

                    // Add to map
                    for( int p = 0; p < nPin; p++ )
                    {
                        for( int q = 0; q < nPin; q++ )
                        {
                            int ixFirstNode;
                            if( fContinuousIn )
                            {
                                ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
                            }
                            else
                            {
                                ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
                            }

                            smatMap( ixSecondNode, ixFirstNode ) += dSampleCoeffIn[p][q];
                        }
                    }
                }
            }
        }

        // Increment the current overlap index
        ixOverlap += nOverlapFaces;
    }

    // Check for missing samples
    for( size_t i = 0; i < fSecondNodeFound.GetRows(); i++ )
    {
        if( !fSecondNodeFound[i] )
        {
            _EXCEPTION1( "Can't sample point %i", i );
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////
