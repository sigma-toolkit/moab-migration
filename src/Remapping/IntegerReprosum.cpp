/*
 * IntegerReprosum.cpp
 *
 *   C++ port of E3SM's `shr_reprosum_int` (Worley's integer-vector
 *   reproducible sum). See IntegerReprosum.hpp for design notes.
 *
 *   Algorithm summary (per call):
 *     1. Reduce gmax_exp / gmin_exp / max_nsummands across the comm.
 *     2. Derive arr_max_shift, max_levels, extra_levels.
 *     3. Encode each local summand into integer-vector levels.
 *     4. Postprocess locally to absorb overflow within each rank.
 *     5. MPI_Allreduce(MPI_SUM, MPI_INT64_T) on the level vector.
 *     6. Reconstruct double: preprocess for non-overlap and same-sign,
 *        truncate at FP-representable boundaries, sum components
 *        smallest-to-largest, restore sign.
 *
 *   Maps directly onto the Fortran reference at lines 1144–1905 of
 *   E3SM/share/util/shr_reprosum_mod.F90.
 */

#include "moab/Remapping/IntegerReprosum.hpp"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <cassert>

namespace moab
{

#ifdef MOAB_HAVE_MPI
IntegerReprosum::IntegerReprosum( MPI_Comm comm ) : m_comm( comm ) {}
#else
IntegerReprosum::IntegerReprosum( int /*comm*/ ) {}
#endif

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace
{

/// Fortran-equivalent fraction(x) for radix-2 doubles: returns f in [0.5, 1)
/// (or in (-1, -0.5] for negative x) and stores the exponent so that
///   x = f * 2^e
/// Same as std::frexp.
inline double frac_and_exp( double x, int& exp_out )
{
    return std::frexp( x, &exp_out );
}

/// Fortran-equivalent scale(x, n) for radix-2 doubles: returns x * 2^n.
/// Same as std::ldexp.
inline double scale2( double x, int n )
{
    return std::ldexp( x, n );
}

/// Fortran-equivalent set_exponent(x, e): returns a value with the same
/// fraction as x but with exponent e. Same as
///   ldexp(frexp(x, &dummy), e)
inline double set_exp( double x, int e )
{
    int dummy;
    double f = std::frexp( x, &dummy );
    return std::ldexp( f, e );
}

/// MPI_INT64_T equivalent for our int64 type. Use MPI_LONG_LONG_INT for
/// portability across MPI implementations that don't yet expose
/// MPI_INT64_T as a primary datatype.
#ifdef MOAB_HAVE_MPI
inline MPI_Datatype mpi_int64()
{
#ifdef MPI_INT64_T
    return MPI_INT64_T;
#else
    return MPI_LONG_LONG_INT;
#endif
}
#endif

/// Integer power 2^k (for small k <= 62).
inline int64_t i2pow( int k )
{
    return ( static_cast< int64_t >( 1 ) << k );
}

}  // namespace

// ---------------------------------------------------------------------------
// Phase 1: metadata reduction
// ---------------------------------------------------------------------------

IntegerReprosum::Metadata IntegerReprosum::compute_metadata( const std::vector< double >& vals,
                                                              const std::vector< int >* mask ) const
{
    // Local extrema (over non-zero values only) and total count of summands
    // (including zeros — matches Fortran shr_reprosum_int's `nsummands`,
    // which is the caller-passed array size, not a non-zero count).
    int local_max_exp = std::numeric_limits< int >::min();
    int local_min_exp = std::numeric_limits< int >::max();
    int local_count   = 0;

    for( size_t i = 0; i < vals.size(); ++i )
    {
        if( mask && ( *mask )[i] < 0 ) continue;
        ++local_count;  // total owned summands, including zeros
        const double v = vals[i];
        if( v == 0.0 ) continue;
        int e;
        std::frexp( v, &e );
        if( e > local_max_exp ) local_max_exp = e;
        if( e < local_min_exp ) local_min_exp = e;
    }

    int gmax_exp = local_max_exp;
    int gmin_exp = local_min_exp;
    int gcount   = local_count;

#ifdef MOAB_HAVE_MPI
    // One Allreduce, three quantities. Use negation trick on local_min_exp
    // so we can take a single MPI_MAX reduction; same trick the Fortran
    // does for arr_lextremes.
    int local_arr[3];
    int global_arr[3];
    local_arr[0] = local_count;       // MPI_MAX over count -> max_nsummands
    local_arr[1] = local_max_exp;     // MPI_MAX over max -> gmax_exp
    local_arr[2] = -local_min_exp;    // MPI_MAX over (-min) -> -gmin_exp
    MPI_Allreduce( local_arr, global_arr, 3, MPI_INT, MPI_MAX, m_comm );
    gcount   = global_arr[0];
    gmax_exp = global_arr[1];
    gmin_exp = -global_arr[2];
#endif

    // Mirror MCT's all-zero / no-non-zero-summand-anywhere fixup
    // (shr_reprosum_mod.F90 lines 939-941):
    //   arr_gmin_exp = min(arr_gmax_exp, arr_gmin_exp)
    // This collapses the sentinel pair (gmax = INT_MIN, gmin = INT_MAX)
    // to (INT_MIN, INT_MIN) so subsequent arithmetic on (gmax - gmin)
    // doesn't blow up.
    if( gmin_exp > gmax_exp ) gmin_exp = gmax_exp;

    Metadata md;
    md.max_nsummands = gcount;
    md.gmax_exp      = gmax_exp;
    md.gmin_exp      = gmin_exp;

    // Edge case: all summands zero (or no owned summands anywhere).
    // Use arbitrary safe defaults; sum() will skip work and return 0.
    if( md.max_nsummands == 0 )
    {
        md.arr_max_shift = kI8Digits / 4;
        md.max_levels    = 2;
        md.extra_levels  = ( kI8Digits - 1 ) / md.arr_max_shift;
        md.gmax_exp      = 0;
        md.gmin_exp      = 0;
        return md;
    }

    // Conservative bound: account for thread-then-task summation per
    // shr_reprosum_calc lines 800-802 (we have one thread, but keep the
    // identical formula for byte-equivalence with MCT's thread=1 case).
    const int omp_nthreads_local = 1;
    int max_n = ( md.max_nsummands / omp_nthreads_local ) + 1;
#ifdef MOAB_HAVE_MPI
    int nproc = 1;
    MPI_Comm_size( m_comm, &nproc );
    if( max_n < nproc * omp_nthreads_local ) max_n = nproc * omp_nthreads_local;
#endif

    // arr_max_shift = digits(int64) - (exponent(real(max_n)) + 1)
    //
    // Fortran's exponent(x) for x = f * 2^E with f in [0.5, 1) returns E.
    // C++ frexp(x, &e) gives the same e. The +1 here accounts for the
    // upper bound: max_n < 2^(e+1), so summing max_n integers each less
    // than 2^arr_max_shift in absolute value gives a sum less than
    // 2^(arr_max_shift + e + 1). For this to fit in int64 (without
    // overflow during MPI_Allreduce(MPI_SUM)) we need
    //   arr_max_shift + e + 1 <= digits(int64)
    // hence arr_max_shift = digits(int64) - (e + 1).
    int e_of_max_n;
    std::frexp( static_cast< double >( max_n ), &e_of_max_n );
    md.arr_max_shift = kI8Digits - ( e_of_max_n + 1 );

    if( md.arr_max_shift < 2 )
    {
        // Too many summands. The Fortran aborts here. We do too — caller
        // shouldn't call us with > ~2^60 summands per rank. Returning a
        // bogus value is worse than aborting clearly.
        std::abort();
    }

    // max_levels = 2 + (digits(double) + (gmax_exp - gmin_exp)) / arr_max_shift
    md.max_levels = 2 + ( kR8Digits + ( md.gmax_exp - md.gmin_exp ) ) / md.arr_max_shift;
    if( md.max_levels < 2 ) md.max_levels = 2;

    // extra_levels = (digits(int64) - 1) / arr_max_shift
    md.extra_levels = ( kI8Digits - 1 ) / md.arr_max_shift;
    if( md.extra_levels < 1 ) md.extra_levels = 1;

    return md;
}

// ---------------------------------------------------------------------------
// Phase 3: encode local summands into integer-vector representation
// ---------------------------------------------------------------------------
//
// Layout: iv has size (max_levels + extra_levels). The Fortran indexes
// from -(extra_levels-1) up to +max_levels, where positive levels carry
// the most-significant digits. We use a 0-based index where
//   iv_index(level) = level + (extra_levels - 1)
// so:
//   level = -(extra_levels - 1) -> iv_index = 0          (least significant)
//   level =  +1                 -> iv_index = extra_levels
//   level =  +max_levels        -> iv_index = max_levels + extra_levels - 1
// ---------------------------------------------------------------------------

void IntegerReprosum::encode_local( const std::vector< double >& vals,
                                     const std::vector< int >* mask,
                                     const Metadata& md,
                                     std::vector< int64_t >& iv ) const
{
    const int max_levels    = md.max_levels;
    const int extra_levels  = md.extra_levels;
    const int arr_max_shift = md.arr_max_shift;
    const int gmax_exp      = md.gmax_exp;

    iv.assign( max_levels + extra_levels, 0 );
    if( md.max_nsummands == 0 ) return;

    auto iv_index = [extra_levels]( int level ) -> int { return level + ( extra_levels - 1 ); };

    for( size_t i = 0; i < vals.size(); ++i )
    {
        if( mask && ( *mask )[i] < 0 ) continue;
        const double x = vals[i];
        if( x == 0.0 ) continue;

        int arr_exp;
        const double arr_frac = std::frexp( x, &arr_exp );

        // gmax_exp is supposed to be a global upper bound; if a summand
        // exceeds it, the metadata reduction was wrong. Caller bug —
        // skip silently rather than overflow the integer vector.
        if( arr_exp > gmax_exp ) continue;

        int arr_shift = arr_max_shift - ( gmax_exp - arr_exp );
        int ilevel;

        if( arr_shift < 1 )
        {
            ilevel = ( 1 + ( gmax_exp - arr_exp ) ) / arr_max_shift;
            arr_shift = ilevel * arr_max_shift - ( gmax_exp - arr_exp );
            while( arr_shift < 1 )
            {
                arr_shift += arr_max_shift;
                ilevel    += 1;
            }
        }
        else
        {
            ilevel = 1;
        }

        if( ilevel > max_levels ) continue;  // smaller than smallest representable

        // First shift / truncate / accumulate.
        double remainder = scale2( arr_frac, arr_shift );
        int64_t i_part   = static_cast< int64_t >( remainder );  // truncates toward 0
        iv[iv_index( ilevel )] += i_part;
        remainder -= static_cast< double >( i_part );

        // Continue while remainder is non-zero and we still have levels.
        while( remainder != 0.0 && ilevel < max_levels )
        {
            ++ilevel;
            remainder = scale2( remainder, arr_max_shift );
            i_part    = static_cast< int64_t >( remainder );
            iv[iv_index( ilevel )] += i_part;
            remainder -= static_cast< double >( i_part );
        }
    }

    // Postprocess: walk levels high-to-low (most-significant to least),
    // moving overflow into the lower-significance level. Same as the
    // Fortran "(a)" comment block at lines 1410–1432 of
    // shr_reprosum_mod.F90. Required so the integer vector cannot
    // overflow during the subsequent MPI_Allreduce(SUM).
    const int64_t shift_factor = i2pow( arr_max_shift );
    const int min_level        = -( extra_levels - 1 );
    for( int level = max_levels; level >= min_level + 1; --level )
    {
        const int idx = iv_index( level );
        if( std::llabs( iv[idx] ) >= shift_factor )
        {
            const int64_t carry = iv[idx] / shift_factor;
            iv[iv_index( level - 1 )] += carry;
            iv[idx] -= carry * shift_factor;
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 4: MPI_Allreduce on the integer vector
// ---------------------------------------------------------------------------

void IntegerReprosum::reduce_global( std::vector< int64_t >& iv ) const
{
#ifdef MOAB_HAVE_MPI
    if( iv.empty() ) return;
    std::vector< int64_t > out( iv.size() );
    MPI_Allreduce( iv.data(), out.data(), static_cast< int >( iv.size() ), mpi_int64(), MPI_SUM,
                   m_comm );
    iv.swap( out );
#else
    (void)iv;
#endif
}

// ---------------------------------------------------------------------------
// Phase 5: reconstruct double from the global integer vector
// ---------------------------------------------------------------------------
//
// Faithful port of the decode section of shr_reprosum_int
// (Fortran lines 1593–1867). Steps mirror the comments in the reference.
// ---------------------------------------------------------------------------

double IntegerReprosum::decode_global( const std::vector< int64_t >& iv_in,
                                        const Metadata& md ) const
{
    if( md.max_nsummands == 0 ) return 0.0;
    const int max_levels    = md.max_levels;
    const int extra_levels  = md.extra_levels;
    const int arr_max_shift = md.arr_max_shift;
    const int gmax_exp      = md.gmax_exp;
    const int min_level     = -( extra_levels - 1 );

    auto iv_index = [extra_levels]( int level ) -> int { return level + ( extra_levels - 1 ); };

    // Working copy.
    std::vector< int64_t > iv = iv_in;
    const int64_t shift_factor = i2pow( arr_max_shift );

    // ---- (a)(i) propagate carries to non-overlap (high to low) ------------
    for( int level = max_levels; level >= min_level + 1; --level )
    {
        const int idx = iv_index( level );
        if( std::llabs( iv[idx] ) >= shift_factor )
        {
            const int64_t carry = iv[idx] / shift_factor;
            iv[iv_index( level - 1 )] += carry;
            iv[idx] -= carry * shift_factor;
        }
    }

    // Find the first non-zero level (low to high).
    int first_level = max_levels;
    for( int level = min_level; level <= max_levels; ++level )
    {
        if( iv[iv_index( level )] != 0 )
        {
            first_level = level;
            break;
        }
    }
    if( first_level == max_levels && iv[iv_index( max_levels )] == 0 )
    {
        // Sum is exactly zero.
        return 0.0;
    }

    // Determine sign of the sum (sign of first non-zero level).
    int64_t sign = ( iv[iv_index( first_level )] < 0 ) ? -1 : 1;

    // ---- (a)(ii) make all components have the same sign ------------------
    if( first_level < max_levels )
    {
        for( int j = first_level; j <= max_levels - 1; ++j )
        {
            const int j_idx  = iv_index( j );
            const int j1_idx = iv_index( j + 1 );
            const int64_t s_here = ( iv[j_idx] < 0 ) ? -1 : ( iv[j_idx] > 0 ? 1 : 0 );
            const int64_t s_next = ( iv[j1_idx] < 0 ) ? -1 : ( iv[j1_idx] > 0 ? 1 : 0 );
            // Treat 0 at the next level as "different sign so always
            // borrow", matching the Fortran condition
            //   sign(jlevel) /= sign(jlevel+1) .or. iv(jlevel+1) == 0
            if( s_here != s_next || iv[j1_idx] == 0 )
            {
                iv[j_idx]  -= sign;
                iv[j1_idx] += sign * shift_factor;
            }
        }
    }

    // ---- (a)(iii) flip to positive temporarily ---------------------------
    if( sign < 0 )
    {
        for( int level = first_level; level <= max_levels; ++level )
            iv[iv_index( level )] = -iv[iv_index( level )];
    }

    // ---- (a)(iv) re-impose non-overlap (carries) -------------------------
    for( int level = max_levels; level >= min_level + 1; --level )
    {
        const int idx = iv_index( level );
        if( std::llabs( iv[idx] ) >= shift_factor )
        {
            const int64_t carry = iv[idx] / shift_factor;
            iv[iv_index( level - 1 )] += carry;
            iv[idx] -= carry * shift_factor;
        }
    }

    // ---- (b)(c)(d) iterate: truncate at FP-representable digit, convert
    //                          to FP, append to summand_vector ------------
    std::vector< double > summand_vector;
    summand_vector.reserve( static_cast< size_t >( max_levels + extra_levels ) );

    bool first_iteration = true;
    int arr_shift_curr   = gmax_exp - min_level * arr_max_shift;
    int digit_count      = 0;
    int begin_level      = min_level;

    while( begin_level <= max_levels )
    {
        // Determine the level at which the cumulative number of integer
        // digits equals or exceeds digits(double) = 53. That's where
        // truncation needs to happen.
        int trunc_loc   = 0;
        int trunc_level = max_levels;

        for( int level = begin_level; level <= max_levels; ++level )
        {
            int LX;
            if( first_iteration )
            {
                if( digit_count == 0 )
                {
                    if( iv[iv_index( level )] != 0 )
                    {
                        const double Xf = static_cast< double >( iv[iv_index( level )] );
                        int e_of_X;
                        std::frexp( Xf, &e_of_X );
                        LX = e_of_X;
                    }
                    else
                    {
                        LX = 0;
                    }
                }
                else
                {
                    LX = arr_max_shift;
                }
            }
            else
            {
                if( level == begin_level && digit_count != 0 )
                    LX = 0;
                else
                    LX = arr_max_shift;
            }

            if( digit_count + LX >= kR8Digits )
            {
                trunc_level = level;
                trunc_loc   = ( digit_count + LX ) - kR8Digits;
                break;
            }
            else
            {
                digit_count += LX;
            }
        }
        first_iteration = false;

        // Compute the truncated value at trunc_level and the remainder
        // (the bits that didn't fit in digits(double)).
        int64_t trunc_level_rem = 0;
        if( trunc_loc != 0 )
        {
            const int64_t pow_trunc = i2pow( trunc_loc );
            const int64_t kept      = iv[iv_index( trunc_level )] / pow_trunc;
            const int64_t kept_full = kept * pow_trunc;
            trunc_level_rem         = iv[iv_index( trunc_level )] - kept_full;
            iv[iv_index( trunc_level )] = kept_full;
        }

        // Convert truncated integer-vector segment [begin_level..trunc_level]
        // to FP and accumulate into a fresh summand_vector entry.
        double seg_sum = 0.0;
        for( int level = begin_level; level <= trunc_level; ++level )
        {
            const int64_t v = iv[iv_index( level )];
            if( v != 0 )
            {
                const double Xf = static_cast< double >( v );
                int e_of_X;
                std::frexp( Xf, &e_of_X );
                const int curr_exp = e_of_X + arr_shift_curr;
                const int min_exp  = std::numeric_limits< double >::min_exponent;
                if( curr_exp >= min_exp )
                {
                    seg_sum += set_exp( Xf, curr_exp );
                }
                else
                {
                    // Subnormal-region scaling: split into two ldexp's so
                    // intermediate stays representable. Mirrors the Fortran
                    // set_exponent + scale combo.
                    const double rxv = set_exp( Xf, curr_exp - min_exp );
                    seg_sum += scale2( rxv, min_exp );
                }
            }

            // Step the arr_shift down by arr_max_shift unless we're
            // staying at the same trunc_level for the next iteration
            // (which happens when trunc_loc > 0).
            if( level < trunc_level || trunc_loc == 0 )
            {
                arr_shift_curr -= arr_max_shift;
            }
        }

        summand_vector.push_back( seg_sum );

        if( trunc_loc == 0 )
        {
            digit_count = 0;
            begin_level = trunc_level + 1;
        }
        else
        {
            digit_count = trunc_loc;
            begin_level = trunc_level;
            // The remainder at trunc_level becomes the new starting value
            // for the next iteration.
            iv[iv_index( trunc_level )] = trunc_level_rem;
        }
    }

    // ---- (e) sum smallest to largest -------------------------------------
    double result = 0.0;
    for( auto it = summand_vector.rbegin(); it != summand_vector.rend(); ++it )
        result += *it;

    // ---- (f) restore sign ------------------------------------------------
    if( sign < 0 ) result = -result;

    return result;
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

double IntegerReprosum::sum( const std::vector< double >& vals ) const
{
    return sum_masked( vals, std::vector< int >() );
}

double IntegerReprosum::sum_masked( const std::vector< double >& vals,
                                     const std::vector< int >& mask ) const
{
    const std::vector< int >* mask_ptr = mask.empty() ? nullptr : &mask;

    Metadata md = compute_metadata( vals, mask_ptr );
    if( md.max_nsummands == 0 ) return 0.0;

    std::vector< int64_t > iv;
    encode_local( vals, mask_ptr, md, iv );
    reduce_global( iv );
    return decode_global( iv, md );
}

void IntegerReprosum::sum_masked_batch( const std::vector< std::vector< double > >& fields,
                                         const std::vector< int >& mask,
                                         std::vector< double >& gsums ) const
{
    gsums.assign( fields.size(), 0.0 );
    // For now this is just a loop; the metadata Allreduce could be batched
    // across fields for a small saving but each field still needs its own
    // integer-vector Allreduce. Leave as-is until profiling says otherwise.
    for( size_t f = 0; f < fields.size(); ++f )
        gsums[f] = sum_masked( fields[f], mask );
}

}  // namespace moab
