/*
 * =====================================================================================
 *
 *       Filename:  ApplyWeights.cpp
 *
 *    Description:  Kernels to apply the sparse-matrix onto a dense vector using various
 *                  optimized implementation algorithms. Internally, we store the maps
 *                  using Eigen3 datastructures, but for bit-for-bit reproducibility,
 *                  we also have some alternatives to experiment.
 *
 *       Revision:  none
 *
 * =====================================================================================
 */

#include "moab/Remapping/TempestOnlineMap.hpp"

#include <algorithm>
#include <utility>
#include <vector>

// ** Kahan Summation Algorithm for improved numerical accuracy **
struct KahanSum
{
    double sum        = 0.0;
    double correction = 0.0;

    void add( double value )
    {
        double y   = value - correction;  // Correct the input
        double t   = sum + y;             // Perform the sum
        correction = ( t - sum ) - y;     // Update correction
        sum        = t;                   // Store the new sum
    }

    double result() const
    {
        return sum;
    }
};

///////////////////////////////////////////////////////////////////////////////


// Pairwise summation helper function
inline double pairwiseSum( const std::set< double >& sorted )
{
    if( sorted.empty() ) return 0.0;
    if( sorted.size() == 1 ) return *sorted.begin();

    // Accumulate pairwise to minimize rounding error
    double sum = 0.0;
    for( double val : sorted )
        sum += val;
    return sum;
}

// Pairwise summation helper function
inline double pairwiseKahanSum( const std::set< double >& sorted )
{
    if( sorted.empty() ) return 0.0;
    if( sorted.size() == 1 ) return *sorted.begin();

    // Accumulate pairwise to minimize rounding error
    // Apply pairwise summation with Kahan correction
    KahanSum kahan;
    for( double val : sorted )
        kahan.add( val );
    return kahan.result();
}

// Sparse matrix-vector multiplication using pairwise summation
inline void deterministicSparseMatVecMul( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                          const typename moab::TempestOnlineMap::WeightColVector& x,
                                          typename moab::TempestOnlineMap::WeightRowVector& result )
{
    constexpr bool useKahanSum    = false;
    constexpr bool usePairwiseSum = false;

    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        std::set< double > accumulators;
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            // accumulators contains the sorted values of the product: A(row, col) * x(col)
            accumulators.insert( it.value() * x( it.col() ) );
        }
        if( usePairwiseSum ) result( row ) = pairwiseSum( accumulators );
        if( useKahanSum ) result( row ) = pairwiseKahanSum( accumulators );

        if( !usePairwiseSum && !useKahanSum )
        {
            double sum = 0.0;
            for( double val : accumulators )
                sum += val;
            result( row ) = sum;
        }
    }
}

//
// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatVecMulKahan( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                               const typename moab::TempestOnlineMap::WeightColVector& x,
                                               typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        KahanSum kahan;
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            double product = it.value() * x( it.col() );  // Compute product
            kahan.add( product );
        }

        result( row ) = kahan.result();
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatVecMulClean( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                               const typename moab::TempestOnlineMap::WeightColVector& x,
                                               typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate row-wise to enforce a fixed summation order
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            const double product = it.value() * x( it.col() );  // Compute product
            result( it.row() ) += product;
        }
    }
}

inline void deterministicSparseMatVecMulNative( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                const typename moab::TempestOnlineMap::WeightColVector& x,
                                                typename moab::TempestOnlineMap::WeightRowVector& result )
{
    result = A * x;  // Perform the matrix-vector multiplication using Eigen3
}

// Sparse matrix-vector multiplication with a per-row sum order keyed by the
// GLOBAL source-DOF index. Eigen's CSR InnerIterator visits stored entries in
// ascending LOCAL-column order, but the local-column index for the same
// global source DOF can change with the MPI rank count and with the order in
// which coverage entities were inserted. Sorting each row's contributions by
// their global column DOF before summing makes the per-row dot product
// invariant to local matrix layout, hence bit-for-bit reproducible across
// rank-count and decomposition changes.
//
// Cost: O(nnz log K) where K is the average row stencil size. For typical
// FV-FV (~9-25 nnz/row) and FV-SE (~16-100 nnz/row) coupling maps this is
// dominated by the SpMV itself; the sort is small per-row work and avoids
// any global communication.
inline void deterministicSparseMatVecMulSorted( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                const typename moab::TempestOnlineMap::WeightColVector& x,
                                                typename moab::TempestOnlineMap::WeightRowVector& result,
                                                const std::vector< unsigned >& col_gdofmap )
{
    result.setZero();

    std::vector< std::pair< unsigned, double > > terms;
    for( int row = 0; row < A.outerSize(); ++row )
    {
        terms.clear();
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            const int c          = static_cast< int >( it.col() );
            const unsigned gcol  = ( c >= 0 && c < static_cast< int >( col_gdofmap.size() ) )
                                       ? col_gdofmap[c]
                                       : static_cast< unsigned >( c );
            const double product = it.value() * x( c );
            terms.emplace_back( gcol, product );
        }
        std::sort( terms.begin(), terms.end(),
                   []( const std::pair< unsigned, double >& a,
                       const std::pair< unsigned, double >& b ) { return a.first < b.first; } );
        double s = 0.0;
        for( const auto& p : terms )
            s += p.second;
        result( row ) = s;
    }
}

// Transpose variant: y = A^T x with per-output-column accumulation in
// global-row-DOF order. Uses one accumulator vector per local matrix column;
// each accumulator sorts its (global_row, value) pairs before summing so the
// final entry is independent of the row-traversal order chosen by Eigen.
inline void deterministicSparseMatTransposeVecMulSorted( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                         const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                         typename moab::TempestOnlineMap::WeightColVector& result,
                                                         const std::vector< unsigned >& row_gdofmap )
{
    result.setZero();

    std::vector< std::vector< std::pair< unsigned, double > > > accumulators( A.cols() );

    for( int row = 0; row < A.outerSize(); ++row )
    {
        const unsigned grow = ( row >= 0 && row < static_cast< int >( row_gdofmap.size() ) )
                                  ? row_gdofmap[row]
                                  : static_cast< unsigned >( row );
        const double xr     = x( row );
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            accumulators[it.col()].emplace_back( grow, it.value() * xr );
        }
    }

    for( int col = 0; col < A.cols(); ++col )
    {
        auto& bin = accumulators[col];
        std::sort( bin.begin(), bin.end(),
                   []( const std::pair< unsigned, double >& a,
                       const std::pair< unsigned, double >& b ) { return a.first < b.first; } );
        double s = 0.0;
        for( const auto& p : bin )
            s += p.second;
        result( col ) = s;
    }
}

// Deterministic sparse matrix-vector multiplication with A^T * x using pairwise summation
inline void deterministicSparseMatTransposeVecMul( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                   const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                   typename moab::TempestOnlineMap::WeightColVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Temporary storage for pairwise summation
    std::vector< std::set< double > > accumulators( A.cols() );

    // Iterate over A row-wise, but accumulate into result as if computing A^T * x
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            accumulators[it.col()].insert( it.value() * x( row ) );
        }
    }

    // Compute final sum using pairwise summation for each entry
    for( int col = 0; col < A.cols(); ++col )
    {
        // result( col ) = pairwiseSum( accumulators[col] );
        result( col ) = pairwiseKahanSum( accumulators[col] );
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatTransposeVecMulClean( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                        const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                        typename moab::TempestOnlineMap::WeightColVector& result )
{
    result.setZero();  // Ensure no uninitialized memory issues

    // Iterate over A row-wise, but accumulate into result as if computing A^T * x
    for( int row = 0; row < A.outerSize(); ++row )
    {
        for( typename moab::TempestOnlineMap::WeightMatrix::InnerIterator it( A, row ); it; ++it )
        {
            const double product = it.value() * x( row );  // Compute product
            result( it.col() ) += product;                 // Accumulate contributions to the corresponding row in A^T
        }
    }
}

// Perform a deterministic sparse matrix-vector multiplication
inline void deterministicSparseMatTransposeVecMulNative( const typename moab::TempestOnlineMap::WeightMatrix& A,
                                                         const typename moab::TempestOnlineMap::WeightRowVector& x,
                                                         typename moab::TempestOnlineMap::WeightColVector& result )
{
    result = A.adjoint() * x;  // Perform the adjoint.matrix-vector multiplication using Eigen3
}
///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::ApplyWeights( std::vector< double >& srcVals,
                                                      std::vector< double >& tgtVals,
                                                      bool transpose )
{
    // Reset the source and target data first
    m_rowVector.setZero();
    m_colVector.setZero();

#ifdef VERBOSE
    std::stringstream sstr;
    static int callId = 0;
    callId++;
    sstr << "projection_id_" << callId << "_s_" << size << "_rk_" << rank << ".txt";
    std::ofstream output_file( sstr.str() );
#endif
    // Perform the actual projection of weights: application of weight matrix onto the source
    // solution vector

    if( transpose )
    {
        // Permute the source data first
        for( unsigned i = 0; i < srcVals.size(); ++i )
        {
            if( row_dtoc_dofmap[i] >= 0 )
                m_rowVector( row_dtoc_dofmap[i] ) = srcVals[i];  // permute and set the row (source) vector properly
        }

        // Now apply the adjoint operator: m_colVector = m_weightMatrix.adjoint() * m_rowVector;
        // Master baseline: Eigen's natural CSR per-column accumulation.
        // The local-column index assignment in col_dtoc_dofmap is already
        // partition-invariant on master (low-order SpMV is BFB across rank
        // counts under this kernel), so reordering by global row DOF —
        // attempted in deterministicSparseMatTransposeVecMulSorted — gives
        // a different (also deterministic) summation order and breaks the
        // BFB property the master baseline relies on. Keep the *Sorted
        // variant defined for future experimentation but do not call it
        // from this code path.
        deterministicSparseMatTransposeVecMulClean( m_weightMatrix, m_rowVector, m_colVector );
        // deterministicSparseMatTransposeVecMulSorted( m_weightMatrix, m_rowVector, m_colVector, row_gdofmap );
        // deterministicSparseMatTransposeVecMul( m_weightMatrix, m_rowVector, m_colVector );
        // deterministicSparseMatTransposeVecMulNative( m_weightMatrix, m_rowVector, m_colVector );

        // Permute the resulting target data back
        for( unsigned i = 0; i < tgtVals.size(); ++i )
        {
            if( col_dtoc_dofmap[i] >= 0 )
                tgtVals[i] = m_colVector( col_dtoc_dofmap[i] );  // permute and set the row (source) vector properly
        }
    }
    else
    {
#ifdef VERBOSE
        output_file << "ColVector: " << m_colVector.size() << ", SrcVals: " << srcVals.size()
                    << ", Sizes: " << m_nTotDofs_SrcCov << ", " << col_dtoc_dofmap.size() << "\n";
#endif
        for( unsigned i = 0; i < srcVals.size(); ++i )
        {
            if( col_dtoc_dofmap[i] >= 0 )
                m_colVector( col_dtoc_dofmap[i] ) = srcVals[i];  // permute and set the row (source) vector properly
#ifdef VERBOSE
            output_file << i << " " << col_gdofmap[col_dtoc_dofmap[i]] + 1 << "  " << srcVals[i] << "\n";
#endif
        }
        
        // Now apply the operator: m_rowVector = m_weightMatrix * m_colVector;
        // Master baseline: Eigen's natural CSR per-row accumulation. The
        // local-column index assignment in col_dtoc_dofmap is partition-
        // invariant on master (validated by imoab_read_compute_map's
        // baseline check up to 16 tasks for the low-order map), so the
        // per-row sum order is already BFB across rank counts. Sorting
        // by global col DOF — attempted in deterministicSparseMatVecMulSorted
        // — produces a different (also deterministic) summation order and
        // breaks the BFB property the master baseline relies on. Keep the
        // *Sorted variant defined for future experimentation but do not
        // call it from this code path.
        deterministicSparseMatVecMulClean( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMulSorted( m_weightMatrix, m_colVector, m_rowVector, col_gdofmap );
        // deterministicSparseMatVecMul( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMulNative( m_weightMatrix, m_colVector, m_rowVector );
        // deterministicSparseMatVecMulKahan( m_weightMatrix, m_colVector, m_rowVector );

        // Permute the resulting target data back
#ifdef VERBOSE
        output_file << "RowVector: " << m_rowVector.size() << ", TgtVals:" << tgtVals.size()
                    << ", Sizes: " << m_nTotDofs_Dest << ", " << row_gdofmap.size() << "\n";
#endif
        for( unsigned i = 0; i < tgtVals.size(); ++i )
        {
            if( row_dtoc_dofmap[i] >= 0 )
            {
                tgtVals[i] = m_rowVector( row_dtoc_dofmap[i] );  // permute and set the row (source) vector properly
#ifdef VERBOSE
                output_file << i << " " << row_gdofmap[row_dtoc_dofmap[i]] + 1 << "  " << tgtVals[i] << "\n";
#endif
            }
        }
    }

    // if( caasType != CAAS_NONE )
    // {
    //     constexpr int nmax_caas_iterations = 5;
    //     double mismatch                    = 1.0;
    //     int caasIteration                  = 0;
    //     while( mismatch > 1e-15 &&
    //            caasIteration++ < nmax_caas_iterations )  // iterate until convergence or a maximum of 5 iterations
    //     {
    //         std::pair< double, double > mDefect = this->ApplyCAASLimiting( srcVals, tgtVals, caasType );
    //         if( m_remapper->verbose )
    //             printf( "Rank %d: -- Iteration: %d, Net original mass defect: %3.4e, mass defect post-CAAS: %3.4e\n",
    //                     m_remapper->rank, caasIteration, mDefect.first, mDefect.second );
    //         mismatch = mDefect.second;
    //     }
    // }

#ifdef VERBOSE
    output_file.flush();  // required here
    output_file.close();
#endif

    // All done with matvec application
    return moab::MB_SUCCESS;
}


