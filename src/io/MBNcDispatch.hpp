//-------------------------------------------------------------------------
// Filename      : MBNcDispatch.hpp
//
// Purpose       : Runtime dispatch layer for NetCDF / Parallel-NetCDF I/O.
//
//   Before this layer existed, src/io/{ReadNC,WriteNC,NCHelper*} chose
//   between the standard NetCDF C API (nc_*) and Parallel-NetCDF
//   (ncmpi_*) at *compile time* via the NCFUNC* macros in ReadNC.hpp /
//   WriteNC.hpp. The compile-time choice is wrong: PNetCDF cannot read
//   NetCDF-4/HDF5 files at all, so any build configured with
//   MOAB_HAVE_PNETCDF (and not MOAB_HAVE_NETCDFPAR) silently fails to
//   load SCRIP / ESMF / MPAS grid files that happen to be NetCDF-4.
//
//   This header replaces the compile-time choice with runtime dispatch
//   based on the *detected* file format. Both pnetcdf.h and netcdf.h
//   are included whenever the corresponding feature flag is on, so the
//   wrapper functions can route each call to the backend that actually
//   handles the file.
//
//   Selection logic (per ReadParallelMap pattern in
//   src/Remapping/TempestOnlineMapIO.cpp):
//     - Classic (CDF-1/2/5)   : prefer PNetCDF; fall back to nc_*_par;
//                               last resort, buffered serial+scatter.
//     - NetCDF-4 (HDF5)       : nc_*_par; last resort, buffered serial.
//     - Serial (mpi_size==1)  : plain nc_* (handles either format).
//     - Unknown / unreadable  : error.
//
// File-id tagging
//   Library file IDs returned by nc_open() and ncmpi_open() come from
//   independent namespaces and CAN collide. Rather than maintain a
//   registry keyed by the library id, this layer tags the backend
//   identity into the high nibble of the fileId returned to callers:
//
//        bit 31         bit 28           bit 0
//        +-------+---------------------------+
//        | tag   |     library file id       |
//        +-------+---------------------------+
//
//   Every wrapper function untags the input id and dispatches on the
//   tag. Callers store the *tagged* id wherever they currently store
//   fileId; nothing else changes.
//
// Type bridging
//   ncmpi_* takes MPI_Offset arrays for starts/counts and for some
//   inquiry outputs (inq_dim length, inq_att length). nc_* takes
//   size_t. Wrappers standardize on size_t at the API boundary and
//   convert to MPI_Offset on the stack when routing to PNetCDF.
//
// Nonblocking calls
//   PNetCDF provides ncmpi_iget_* + ncmpi_wait_all for request
//   aggregation. Standard NetCDF has no equivalent. On non-PNetCDF
//   backends, mbnc_iget_* degrades to an immediate blocking collective
//   call and sets the request handle to a sentinel that mbnc_wait_all
//   recognizes and ignores. Existing #ifdef MOAB_HAVE_PNETCDF blocks
//   in callers therefore remain correct after the macro swap.
//
//   Creator       : Vijay Mahadevan, 2026-06-11
//-------------------------------------------------------------------------

#ifndef MB_NC_DISPATCH_HPP
#define MB_NC_DISPATCH_HPP

#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

// Always include the standard NetCDF C API when MOAB has any NetCDF support.
// libnetcdf is a transitive dependency of libpnetcdf, so its headers are
// reachable whenever PNetCDF is configured. The wrapper layer needs both
// pnetcdf.h (for ncmpi_*) AND netcdf.h (for nc_*) compiled in simultaneously.
#ifdef MOAB_HAVE_NETCFF
#include "netcdf.h"

#ifdef MOAB_HAVE_NETCDFPAR
#include "netcdf_par.h"
#endif

#endif // MOAB_HAVE_NETCFF

#ifdef MOAB_HAVE_PNETCDF
#include "pnetcdf.h"
#endif

#include <cstddef>  // size_t, ptrdiff_t

namespace moab
{

// ============================================================================
// Backend tags
// ============================================================================

/// Identifies which underlying NetCDF library will service a given file.
/// Encoded into the high nibble of every tagged file id this layer hands
/// back, so dispatch is branch-free and stateless.
enum NcBackend
{
    NCB_NONE          = 0,  ///< Invalid / closed handle
    NCB_NETCDF_SERIAL = 1,  ///< Plain nc_*, single rank or rank-0 reads
    NCB_PNETCDF       = 2,  ///< ncmpi_*, classic CDF-1/2/5 parallel
    NCB_NETCDF_PAR    = 3,  ///< nc_*_par, parallel HDF5 / parallel CDF (when libnetcdf has PNetCDF backend)
    NCB_BUFFERED      = 4   ///< nc_* on rank 0, MPI scatter/broadcast on every read
};

/// File format detected from the on-disk magic bytes. Independent of
/// which library will end up servicing the file (a classic file can be
/// opened by either PNetCDF or libnetcdf-parallel, for instance).
enum NcFileFormat
{
    NCFMT_UNKNOWN = 0,
    NCFMT_CLASSIC = 1,  ///< NetCDF-3: CDF-1, CDF-2, CDF-5
    NCFMT_NETCDF4 = 2   ///< NetCDF-4 / HDF5
};

// ============================================================================
// Tagged file id encoding
// ============================================================================

/// Bits reserved for the backend tag in the top of a tagged file id.
/// 4 bits leaves 28 bits (268M) for the underlying library file id,
/// which is far more than any NetCDF library will hand out in practice.
constexpr int MBNC_TAG_BITS  = 4;
constexpr int MBNC_TAG_SHIFT = 32 - MBNC_TAG_BITS;
constexpr int MBNC_ID_MASK   = ( 1 << MBNC_TAG_SHIFT ) - 1;

inline int mbnc_make_tagged( int libId, NcBackend backend )
{
    return ( static_cast< int >( backend ) << MBNC_TAG_SHIFT ) | ( libId & MBNC_ID_MASK );
}

inline int mbnc_lib_id( int taggedId )
{
    return taggedId & MBNC_ID_MASK;
}

inline NcBackend mbnc_backend_of( int taggedId )
{
    // Use unsigned shift to avoid sign-extension on negative-looking ids.
    return static_cast< NcBackend >( static_cast< unsigned >( taggedId ) >> MBNC_TAG_SHIFT );
}

// ============================================================================
// Format probe + backend chooser (pure functions — no I/O state)
// ============================================================================

/// Inspect the file header to classify NetCDF format. Reads only the
/// first 8 bytes — relies on the well-known "CDF\xNN" and HDF5 magic
/// signatures. Returns NCFMT_UNKNOWN if the file cannot be opened or
/// the signature does not match. Safe to call from any single rank.
int mbnc_detect_format( const char* path );

/// Pick the best available backend for reading a file of the given
/// format with the given MPI size. Returns NCB_NONE if no compatible
/// backend is configured (caller should error out with a precise
/// message naming the missing capability).
///
/// Decision matrix at runtime, given compile-time MOAB_HAVE_PNETCDF
/// and MOAB_HAVE_NETCDFPAR flags (visible to this function via #ifdef):
///
///   format=classic, mpi_size=1                 -> NCB_NETCDF_SERIAL
///   format=netcdf4, mpi_size=1                 -> NCB_NETCDF_SERIAL
///   format=classic, parallel, have PNetCDF     -> NCB_PNETCDF
///   format=classic, parallel, have NETCDFPAR   -> NCB_NETCDF_PAR
///   format=classic, parallel, neither          -> NCB_BUFFERED
///   format=netcdf4, parallel, have NETCDFPAR   -> NCB_NETCDF_PAR
///   format=netcdf4, parallel, only PNetCDF     -> NCB_BUFFERED
///   format=unknown                             -> NCB_NONE
NcBackend mbnc_choose_backend_for_read( int format, int mpi_size );

/// Pick a backend for writing. The format is provided by the caller
/// (typically derived from a user option / file extension), not probed.
NcBackend mbnc_choose_backend_for_write( int requested_format, int mpi_size );

// ============================================================================
// Buffered-fallback context registry
//
// Only NCB_BUFFERED files need extra metadata (the MPI communicator and
// rank/size) so the wrappers can scatter / broadcast on each call. The
// registry is keyed by the *tagged* file id and consulted only by the
// buffered code path.
// ============================================================================

#ifdef MOAB_HAVE_MPI
void mbnc_register_buffered( int taggedFileId, MPI_Comm comm );
void mbnc_unregister_buffered( int taggedFileId );
MPI_Comm mbnc_buffered_comm( int taggedFileId );
int      mbnc_buffered_rank( int taggedFileId );
int      mbnc_buffered_size( int taggedFileId );
#endif

// ============================================================================
// File open / close / create
//
// These return TAGGED file ids; pass them unmodified to every other
// mbnc_* function. The caller must close via mbnc_close() (never
// nc_close / ncmpi_close directly).
// ============================================================================

#ifdef MOAB_HAVE_MPI
/// Parallel open with explicit backend (chosen via mbnc_choose_backend_for_read).
/// Returns NC_NOERR on success, NetCDF error code otherwise; *taggedFileId
/// is set only on success.
int mbnc_open_par( NcBackend backend, MPI_Comm comm, MPI_Info info, const char* path, int omode, int* taggedFileId );

/// Parallel create with explicit backend (chosen via mbnc_choose_backend_for_write).
int mbnc_create_par( NcBackend backend, MPI_Comm comm, MPI_Info info, const char* path, int cmode,
                     int* taggedFileId );
#endif

/// Serial open. Tags the returned id as NCB_NETCDF_SERIAL.
int mbnc_open( const char* path, int omode, int* taggedFileId );

/// Serial create.
int mbnc_create( const char* path, int cmode, int* taggedFileId );

/// Close. Untags + dispatches.
int mbnc_close( int taggedFileId );

// ============================================================================
// Define mode (writes only)
// ============================================================================

int mbnc_redef( int taggedFileId );
int mbnc_enddef( int taggedFileId );
int mbnc_def_dim( int taggedFileId, const char* name, size_t len, int* dimid );
int mbnc_def_var( int taggedFileId, const char* name, nc_type xtype, int ndims, const int* dimids, int* varid );

// ============================================================================
// PNetCDF independent / collective mode toggles
//
// These are PNetCDF-specific knobs that switch the library between
// collective and independent I/O modes. No-ops on non-PNetCDF backends.
// ============================================================================

int mbnc_begin_indep_data( int taggedFileId );
int mbnc_end_indep_data( int taggedFileId );

// ============================================================================
// Inquiry — file / dimension / variable / attribute metadata
// ============================================================================

int mbnc_inq_natts( int taggedFileId, int* nattsp );
int mbnc_inq_ndims( int taggedFileId, int* ndimsp );
int mbnc_inq_nvars( int taggedFileId, int* nvarsp );

int mbnc_inq_dimid( int taggedFileId, const char* name, int* dimidp );
int mbnc_inq_dim( int taggedFileId, int dimid, char* name, size_t* lenp );
int mbnc_inq_dimlen( int taggedFileId, int dimid, size_t* lenp );

int mbnc_inq_varid( int taggedFileId, const char* name, int* varidp );
int mbnc_inq_varname( int taggedFileId, int varid, char* name );
int mbnc_inq_vartype( int taggedFileId, int varid, nc_type* xtypep );
int mbnc_inq_varndims( int taggedFileId, int varid, int* ndimsp );
int mbnc_inq_vardimid( int taggedFileId, int varid, int* dimids );
int mbnc_inq_varnatts( int taggedFileId, int varid, int* nattsp );

int mbnc_inq_attname( int taggedFileId, int varid, int attnum, char* name );
int mbnc_inq_att( int taggedFileId, int varid, const char* name, nc_type* xtypep, size_t* lenp );

// ============================================================================
// Attribute get / put
// ============================================================================

int mbnc_get_att_text( int taggedFileId, int varid, const char* name, char* value );
int mbnc_get_att_int( int taggedFileId, int varid, const char* name, int* value );
int mbnc_get_att_short( int taggedFileId, int varid, const char* name, short* value );
int mbnc_get_att_long( int taggedFileId, int varid, const char* name, long* value );
int mbnc_get_att_float( int taggedFileId, int varid, const char* name, float* value );
int mbnc_get_att_double( int taggedFileId, int varid, const char* name, double* value );

int mbnc_put_att_text( int taggedFileId, int varid, const char* name, size_t len, const char* value );
int mbnc_put_att_int( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const int* value );
int mbnc_put_att_short( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const short* value );
int mbnc_put_att_float( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const float* value );
int mbnc_put_att_double( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len,
                         const double* value );

// ============================================================================
// Variable get_vara / put_vara — collective by default on parallel backends
//
// All array arguments are size_t. PNetCDF wrappers convert to MPI_Offset
// on the stack (NetCDF variables rarely have more than ~6 dimensions, so
// the conversion is a few-cycle loop).
// ============================================================================

int mbnc_get_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data );
int mbnc_get_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data );
int mbnc_get_vara_long( int taggedFileId, int varid, const size_t* start, const size_t* count, long* data );
int mbnc_get_vara_text( int taggedFileId, int varid, const size_t* start, const size_t* count, char* data );

int mbnc_get_vars_double( int taggedFileId, int varid, const size_t* start, const size_t* count,
                          const ptrdiff_t* stride, double* data );

int mbnc_put_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, const double* data );
int mbnc_put_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data );
int mbnc_put_vara_text( int taggedFileId, int varid, const size_t* start, const size_t* count, const char* data );

// Independent-mode get/put variants — used inside begin_indep_data /
// end_indep_data brackets on PNetCDF. On non-PNetCDF backends the
// collective/independent distinction is per-variable (set externally
// via nc_var_par_access on NETCDF_PAR; meaningless for SERIAL), so
// these route to the same nc_get_vara_* / nc_put_vara_* call as the
// collective wrappers and the caller is responsible for any per-var
// access-mode toggling on the libnetcdf side.
int mbnc_get_vara_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data );
int mbnc_get_vara_int_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data );
int mbnc_get_vara_long_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, long* data );
int mbnc_get_vara_text_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, char* data );
int mbnc_get_vars_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count,
                                const ptrdiff_t* stride, double* data );

int mbnc_put_vara_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count,
                                const double* data );
int mbnc_put_vara_int_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data );
int mbnc_put_vara_text_indep( int taggedFileId, int varid, const size_t* start, const size_t* count,
                              const char* data );

// ============================================================================
// Nonblocking variable get + wait (PNetCDF request aggregation)
//
// On non-PNetCDF backends:
//   - mbnc_iget_* performs an immediate blocking collective call and
//     sets *req = MBNC_REQ_NULL.
//   - mbnc_wait_all is a no-op that returns NC_NOERR and writes
//     NC_NOERR into every entry of statuses[] whose corresponding
//     request equals MBNC_REQ_NULL.
//
// Net effect: existing #ifdef MOAB_HAVE_PNETCDF blocks in callers
// continue to compile and execute correctly even when the runtime
// backend isn't PNetCDF.
// ============================================================================

constexpr int MBNC_REQ_NULL = -1;  ///< sentinel for "no real request pending"

int mbnc_iget_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data,
                           int* req );
int mbnc_iget_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data, int* req );

// Nonblocking puts (symmetric to iget_*). On non-PNetCDF backends: blocking
// collective put; *req = MBNC_REQ_NULL. mbnc_wait_all handles both iget and
// iput requests on PNetCDF (the underlying ncmpi_wait_all is direction-agnostic).
int mbnc_iput_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, const double* data,
                           int* req );
int mbnc_iput_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data,
                        int* req );

int mbnc_wait_all( int taggedFileId, int nreq, int* requests, int* statuses );

}  // namespace moab

#endif  // MB_NC_DISPATCH_HPP
