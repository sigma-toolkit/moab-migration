/*
 * IntegerReprosum.hpp
 *
 *   Bit-reproducible global sum via Worley's integer-vector algorithm.
 *
 *   This is a C++ port of E3SM's `shr_reprosum_int` (Fortran) found in
 *   share/util/shr_reprosum_mod.F90. Each FP summand is represented as an
 *   integer vector keyed by exponent levels; the integer vector is then
 *   reduced via MPI_Allreduce on int64_t (genuinely commutative and
 *   associative), and the result is converted back to floating point with
 *   a deterministic reconstruction. The final sum is bit-identical
 *   regardless of:
 *     - MPI rank count
 *     - local iteration order on each rank
 *     - mesh decomposition
 *
 *   For BFB matching with MCT-coupler runs, set `reprosum_use_ddpdd=.false.`
 *   in `user_nl_cpl` so MCT also uses the integer-vector path (its default).
 *
 *   References:
 *     P. Worley, "Reproducibility and Performance of MPI-Reduce in MPAS-O,"
 *     Oak Ridge National Laboratory.
 *     E3SM share/util/shr_reprosum_mod.F90 :: shr_reprosum_int
 *
 *   Author: ported for MOAB CAAS dual-map BFB work.
 */

#ifndef MOAB_INTEGER_REPROSUM_HPP
#define MOAB_INTEGER_REPROSUM_HPP

#include "moab/MOABConfig.h"

#include <vector>
#include <cstdint>

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

namespace moab
{

/**
 * Reproducible global summation using the integer-vector algorithm
 * (Worley). One instance per MPI communicator; the sum() / sum_masked()
 * entry points can be called repeatedly on different value vectors.
 */
class IntegerReprosum
{
  public:
#ifdef MOAB_HAVE_MPI
    explicit IntegerReprosum( MPI_Comm comm );
#else
    /// Serial constructor; the comm parameter is unused.
    explicit IntegerReprosum( int comm = 0 );
#endif

    ~IntegerReprosum() = default;

    /**
     * Compute the global sum of every entry in `vals` across all ranks
     * on the constructor's MPI communicator. Bit-identical regardless of
     * how `vals` is partitioned across ranks or iterated locally.
     */
    double sum( const std::vector< double >& vals ) const;

    /**
     * Compute the global sum of `vals[i]` only where `mask[i] >= 0`.
     * Useful for excluding halo / not-owned entries when summing per-rank
     * partial vectors. `mask` must be the same length as `vals`.
     */
    double sum_masked( const std::vector< double >& vals,
                       const std::vector< int >& mask ) const;

    /**
     * Convenience batch entry: compute one global sum per input vector,
     * sharing the metadata (gmax/gmin exponent) reduction across all
     * fields in the batch. Equivalent to calling sum_masked() N times
     * but with one fewer MPI_Allreduce per field for the metadata.
     * All input vectors must have the same length and use the same mask.
     */
    void sum_masked_batch( const std::vector< std::vector< double > >& fields,
                           const std::vector< int >& mask,
                           std::vector< double >& gsums ) const;

  private:
#ifdef MOAB_HAVE_MPI
    MPI_Comm m_comm;
#endif

    // === implementation helpers ===

    struct Metadata
    {
        int gmax_exp;       ///< global max exponent of any non-zero summand
        int gmin_exp;       ///< global min exponent of any non-zero summand
        int max_nsummands;  ///< MPI_MAX of local count of non-zero summands
        int arr_max_shift;  ///< per-level integer shift; chosen so the level
                            ///< sum cannot overflow int64
        int max_levels;     ///< number of integer-vector levels for this field
        int extra_levels;   ///< padding levels at the high end to absorb
                            ///< overflow during cross-rank summation
    };

    /// Reduce per-field local exponent extrema and the local count of
    /// non-zero summands across the comm and derive arr_max_shift,
    /// max_levels, extra_levels.
    Metadata compute_metadata( const std::vector< double >& vals,
                               const std::vector< int >* mask ) const;

    /// Encode local summands into the integer-vector representation.
    /// `iv` is sized `max_levels + extra_levels`, indexed so that
    /// iv[idx] corresponds to algorithmic level `idx - (extra_levels - 1)`.
    void encode_local( const std::vector< double >& vals,
                       const std::vector< int >* mask,
                       const Metadata& md,
                       std::vector< int64_t >& iv ) const;

    /// MPI_Allreduce(MPI_SUM, MPI_INT64_T) on the integer vector.
    void reduce_global( std::vector< int64_t >& iv ) const;

    /// Reconstruct double from the global integer vector. Faithful port
    /// of the decode section of shr_reprosum_int (preprocess for non-
    /// overlap and same-sign, truncate at FP-representable boundary,
    /// sum the resulting r8 components smallest-to-largest).
    double decode_global( const std::vector< int64_t >& iv,
                          const Metadata& md ) const;

    // === constants ===
    static constexpr int kI8Digits = 63;  ///< mantissa bits in int64 (sign excluded)
    static constexpr int kR8Digits = 53;  ///< mantissa bits in double
    static constexpr int kRadix    = 2;
};

}  // namespace moab

#endif  // MOAB_INTEGER_REPROSUM_HPP
