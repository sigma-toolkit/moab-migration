//-------------------------------------------------------------------------
// Filename      : MBNcDispatch.cpp
//
// Purpose       : Runtime dispatch implementations. See MBNcDispatch.hpp
//                 for design notes.
//
//   Creator     : Vijay Mahadevan, 2026-06-11
//-------------------------------------------------------------------------

#include "MBNcDispatch.hpp"

#include <cstdio>   // std::fopen / std::fread / std::fclose for format probe
#include <cstring>  // std::memset

#ifdef MOAB_HAVE_MPI
#include <map>
#endif

// Phase-1 limitation: this dispatch layer always needs the standard NetCDF
// C API for the serial / parallel-NetCDF / buffered backends. Pure-PNetCDF
// builds (HAVE_PNETCDF without HAVE_NETCDF) are rare and currently fall
// outside this refactor's scope — they continue to work via the legacy
// NCFUNC macro path until that case is wired up here.
#ifndef MOAB_HAVE_NETCDF
#error "MBNcDispatch currently requires MOAB_HAVE_NETCDF. PNetCDF-only builds are not yet supported by this layer."
#endif

namespace moab
{

// ============================================================================
// Internal helpers
// ============================================================================

namespace
{

/// Maximum NetCDF variable rank we expect to encounter. Plenty for any
/// climate / mesh file (typical max is 4–5). Used for stack-allocated
/// size_t ↔ MPI_Offset conversion buffers in PNetCDF dispatches.
constexpr int kMaxDims = 16;

#ifdef MOAB_HAVE_PNETCDF
/// Inquire variable rank and convert size_t start/count arrays to MPI_Offset
/// in caller-supplied stack buffers. Returns NC_NOERR on success or the
/// underlying ncmpi error code.
inline int to_mpi_offset_pair( int libId, int varid, const size_t* start, const size_t* count, MPI_Offset* outStart,
                               MPI_Offset* outCount, int* outNdims )
{
    int ndims;
    int rc = ncmpi_inq_varndims( libId, varid, &ndims );
    if( rc != NC_NOERR ) return rc;
    if( ndims > kMaxDims ) return NC_EMAXDIMS;
    for( int i = 0; i < ndims; ++i )
    {
        outStart[i] = static_cast< MPI_Offset >( start[i] );
        outCount[i] = static_cast< MPI_Offset >( count[i] );
    }
    *outNdims = ndims;
    return NC_NOERR;
}
#endif

#ifdef MOAB_HAVE_MPI
/// Buffered-context registry. Only NCB_BUFFERED files appear here; the
/// other backends carry all needed state internally. Map is keyed by the
/// tagged file id (so PNetCDF / NetCDF id collisions are impossible).
struct BufferedCtx
{
    MPI_Comm comm;
    int      rank;
    int      size;
};

std::map< int, BufferedCtx >& buffered_registry()
{
    static std::map< int, BufferedCtx > reg;
    return reg;
}
#endif

}  // namespace

// ============================================================================
// Format probe
// ============================================================================

int mbnc_detect_format( const char* path )
{
    std::FILE* fp = std::fopen( path, "rb" );
    if( !fp ) return NCFMT_UNKNOWN;

    unsigned char magic[8];
    std::memset( magic, 0, sizeof( magic ) );
    const size_t nread = std::fread( magic, 1, sizeof( magic ), fp );
    std::fclose( fp );

    if( nread < 4 ) return NCFMT_UNKNOWN;

    // Classic NetCDF families: "CDF" + version byte (0x01, 0x02, or 0x05)
    if( magic[0] == 'C' && magic[1] == 'D' && magic[2] == 'F' )
    {
        if( magic[3] == 0x01 || magic[3] == 0x02 || magic[3] == 0x05 ) return NCFMT_CLASSIC;
    }

    // HDF5 superblock signature (used by NetCDF-4)
    if( nread >= 8 && magic[0] == 0x89 && magic[1] == 'H' && magic[2] == 'D' && magic[3] == 'F' && magic[4] == 0x0D &&
        magic[5] == 0x0A && magic[6] == 0x1A && magic[7] == 0x0A )
    {
        return NCFMT_NETCDF4;
    }

    return NCFMT_UNKNOWN;
}

// ============================================================================
// Backend chooser
// ============================================================================

NcBackend mbnc_choose_backend_for_read( int format, int mpi_size )
{
    if( format != NCFMT_CLASSIC && format != NCFMT_NETCDF4 ) return NCB_NONE;

    // Serial: plain nc_* always handles both formats (libnetcdf v4+).
    if( mpi_size <= 1 ) return NCB_NETCDF_SERIAL;

    if( format == NCFMT_CLASSIC )
    {
#ifdef MOAB_HAVE_PNETCDF
        return NCB_PNETCDF;  // best fit for classic in parallel
#elif defined( MOAB_HAVE_NETCDFPAR )
        return NCB_NETCDF_PAR;  // works if libnetcdf was built with PNetCDF backend
#else
        return NCB_NONE;  // buffered fallback comes in Phase 17 (task #17)
#endif
    }

    // format == NCFMT_NETCDF4
#ifdef MOAB_HAVE_NETCDFPAR
    return NCB_NETCDF_PAR;
#else
    return NCB_NONE;  // PNetCDF alone cannot read NetCDF-4. Buffered fallback in Phase 17.
#endif
}

NcBackend mbnc_choose_backend_for_write( int requested_format, int mpi_size )
{
    // For writes, format is the user's intent — same matrix.
    return mbnc_choose_backend_for_read( requested_format, mpi_size );
}

// ============================================================================
// Buffered-context registry
// ============================================================================

#ifdef MOAB_HAVE_MPI
void mbnc_register_buffered( int taggedFileId, MPI_Comm comm )
{
    BufferedCtx ctx;
    ctx.comm = comm;
    MPI_Comm_rank( comm, &ctx.rank );
    MPI_Comm_size( comm, &ctx.size );
    buffered_registry()[taggedFileId] = ctx;
}

void mbnc_unregister_buffered( int taggedFileId )
{
    buffered_registry().erase( taggedFileId );
}

MPI_Comm mbnc_buffered_comm( int taggedFileId )
{
    auto it = buffered_registry().find( taggedFileId );
    return ( it == buffered_registry().end() ) ? MPI_COMM_NULL : it->second.comm;
}

int mbnc_buffered_rank( int taggedFileId )
{
    auto it = buffered_registry().find( taggedFileId );
    return ( it == buffered_registry().end() ) ? -1 : it->second.rank;
}

int mbnc_buffered_size( int taggedFileId )
{
    auto it = buffered_registry().find( taggedFileId );
    return ( it == buffered_registry().end() ) ? 0 : it->second.size;
}
#endif

// ============================================================================
// File open / close / create
// ============================================================================

#ifdef MOAB_HAVE_MPI
int mbnc_open_par( NcBackend backend, MPI_Comm comm, MPI_Info info, const char* path, int omode, int* taggedFileId )
{
    int libId = -1;
    int rc    = NC_EBADID;

    switch( backend )
    {
#ifdef MOAB_HAVE_PNETCDF
        case NCB_PNETCDF:
            rc = ncmpi_open( comm, path, omode, info, &libId );
            break;
#endif
#ifdef MOAB_HAVE_NETCDFPAR
        case NCB_NETCDF_PAR:
            rc = nc_open_par( path, omode | NC_MPIIO, comm, info, &libId );
            break;
#endif
        case NCB_BUFFERED: {
            // Rank 0 opens; other ranks defer to rank 0 for every subsequent call.
            int rank = 0;
            MPI_Comm_rank( comm, &rank );
            if( rank == 0 ) rc = nc_open( path, omode, &libId );
            MPI_Bcast( &rc, 1, MPI_INT, 0, comm );
            // libId is meaningful only on rank 0, but every rank still gets a tagged
            // handle — the tag carries the backend, and the registry carries the comm.
            if( rc == NC_NOERR )
            {
                *taggedFileId = mbnc_make_tagged( libId, NCB_BUFFERED );
                mbnc_register_buffered( *taggedFileId, comm );
            }
            return rc;
        }
        case NCB_NETCDF_SERIAL:
            // Caller asked for serial open via the parallel entry point — honor it.
            rc = nc_open( path, omode, &libId );
            break;
        default:
            return NC_EBADID;
    }

    if( rc == NC_NOERR ) *taggedFileId = mbnc_make_tagged( libId, backend );
    return rc;
}

int mbnc_create_par( NcBackend backend, MPI_Comm comm, MPI_Info info, const char* path, int cmode, int* taggedFileId )
{
    int libId = -1;
    int rc    = NC_EBADID;

    switch( backend )
    {
#ifdef MOAB_HAVE_PNETCDF
        case NCB_PNETCDF:
            rc = ncmpi_create( comm, path, cmode, info, &libId );
            break;
#endif
#ifdef MOAB_HAVE_NETCDFPAR
        case NCB_NETCDF_PAR:
            rc = nc_create_par( path, cmode | NC_MPIIO, comm, info, &libId );
            break;
#endif
        case NCB_BUFFERED: {
            int rank = 0;
            MPI_Comm_rank( comm, &rank );
            if( rank == 0 ) rc = nc_create( path, cmode, &libId );
            MPI_Bcast( &rc, 1, MPI_INT, 0, comm );
            if( rc == NC_NOERR )
            {
                *taggedFileId = mbnc_make_tagged( libId, NCB_BUFFERED );
                mbnc_register_buffered( *taggedFileId, comm );
            }
            return rc;
        }
        case NCB_NETCDF_SERIAL:
            rc = nc_create( path, cmode, &libId );
            break;
        default:
            return NC_EBADID;
    }

    if( rc == NC_NOERR ) *taggedFileId = mbnc_make_tagged( libId, backend );
    return rc;
}
#endif  // MOAB_HAVE_MPI

int mbnc_open( const char* path, int omode, int* taggedFileId )
{
    int libId = -1;
    int rc    = nc_open( path, omode, &libId );
    if( rc == NC_NOERR ) *taggedFileId = mbnc_make_tagged( libId, NCB_NETCDF_SERIAL );
    return rc;
}

int mbnc_create( const char* path, int cmode, int* taggedFileId )
{
    int libId = -1;
    int rc    = nc_create( path, cmode, &libId );
    if( rc == NC_NOERR ) *taggedFileId = mbnc_make_tagged( libId, NCB_NETCDF_SERIAL );
    return rc;
}

int mbnc_close( int taggedFileId )
{
    const int libId           = mbnc_lib_id( taggedFileId );
    const NcBackend backend   = mbnc_backend_of( taggedFileId );
    int rc                    = NC_EBADID;

    switch( backend )
    {
#ifdef MOAB_HAVE_PNETCDF
        case NCB_PNETCDF:
            rc = ncmpi_close( libId );
            break;
#endif
#ifdef MOAB_HAVE_NETCDFPAR
        case NCB_NETCDF_PAR:
            rc = nc_close( libId );  // parallel handle uses nc_close
            break;
#endif
        case NCB_NETCDF_SERIAL:
            rc = nc_close( libId );
            break;
#ifdef MOAB_HAVE_MPI
        case NCB_BUFFERED: {
            // Only rank 0 actually has an open file.
            const int rank = mbnc_buffered_rank( taggedFileId );
            if( rank == 0 ) rc = nc_close( libId );
            MPI_Bcast( &rc, 1, MPI_INT, 0, mbnc_buffered_comm( taggedFileId ) );
            mbnc_unregister_buffered( taggedFileId );
            return rc;
        }
#endif
        default:
            break;
    }

    return rc;
}

// ============================================================================
// Define mode
// ============================================================================

int mbnc_redef( int taggedFileId )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_redef( libId );
#endif
    (void)backend;
    return nc_redef( libId );
}

int mbnc_enddef( int taggedFileId )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_enddef( libId );
#endif
    (void)backend;
    return nc_enddef( libId );
}

int mbnc_def_dim( int taggedFileId, const char* name, size_t len, int* dimid )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_def_dim( libId, name, static_cast< MPI_Offset >( len ), dimid );
#endif
    (void)backend;
    return nc_def_dim( libId, name, len, dimid );
}

int mbnc_def_var( int taggedFileId, const char* name, nc_type xtype, int ndims, const int* dimids, int* varid )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_def_var( libId, name, xtype, ndims, dimids, varid );
#endif
    (void)backend;
    return nc_def_var( libId, name, xtype, ndims, dimids, varid );
}

// ============================================================================
// PNetCDF independent / collective mode toggles
//
// On non-PNetCDF backends these are no-ops returning NC_NOERR — standard
// NetCDF does not separate independent and collective I/O modes the way
// PNetCDF does.
// ============================================================================

int mbnc_begin_indep_data( int taggedFileId )
{
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_begin_indep_data( mbnc_lib_id( taggedFileId ) );
#endif
    (void)backend;
    return NC_NOERR;
}

int mbnc_end_indep_data( int taggedFileId )
{
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_end_indep_data( mbnc_lib_id( taggedFileId ) );
#endif
    (void)backend;
    return NC_NOERR;
}

// ============================================================================
// Inquiry wrappers
//
// All take size_t pointers for length-style outputs. The PNetCDF branch
// receives into a local MPI_Offset and converts on return.
// ============================================================================

int mbnc_inq_natts( int taggedFileId, int* nattsp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_natts( libId, nattsp );
#endif
    (void)backend;
    return nc_inq_natts( libId, nattsp );
}

int mbnc_inq_ndims( int taggedFileId, int* ndimsp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_ndims( libId, ndimsp );
#endif
    (void)backend;
    return nc_inq_ndims( libId, ndimsp );
}

int mbnc_inq_nvars( int taggedFileId, int* nvarsp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_nvars( libId, nvarsp );
#endif
    (void)backend;
    return nc_inq_nvars( libId, nvarsp );
}

int mbnc_inq_dimid( int taggedFileId, const char* name, int* dimidp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_dimid( libId, name, dimidp );
#endif
    (void)backend;
    return nc_inq_dimid( libId, name, dimidp );
}

int mbnc_inq_dim( int taggedFileId, int dimid, char* name, size_t* lenp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset tmp = 0;
        int rc         = ncmpi_inq_dim( libId, dimid, name, &tmp );
        if( lenp ) *lenp = static_cast< size_t >( tmp );
        return rc;
    }
#endif
    (void)backend;
    return nc_inq_dim( libId, dimid, name, lenp );
}

int mbnc_inq_dimlen( int taggedFileId, int dimid, size_t* lenp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset tmp = 0;
        int rc         = ncmpi_inq_dimlen( libId, dimid, &tmp );
        if( lenp ) *lenp = static_cast< size_t >( tmp );
        return rc;
    }
#endif
    (void)backend;
    return nc_inq_dimlen( libId, dimid, lenp );
}

int mbnc_inq_varid( int taggedFileId, const char* name, int* varidp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_varid( libId, name, varidp );
#endif
    (void)backend;
    return nc_inq_varid( libId, name, varidp );
}

int mbnc_inq_varname( int taggedFileId, int varid, char* name )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_varname( libId, varid, name );
#endif
    (void)backend;
    return nc_inq_varname( libId, varid, name );
}

int mbnc_inq_vartype( int taggedFileId, int varid, nc_type* xtypep )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_vartype( libId, varid, xtypep );
#endif
    (void)backend;
    return nc_inq_vartype( libId, varid, xtypep );
}

int mbnc_inq_varndims( int taggedFileId, int varid, int* ndimsp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_varndims( libId, varid, ndimsp );
#endif
    (void)backend;
    return nc_inq_varndims( libId, varid, ndimsp );
}

int mbnc_inq_vardimid( int taggedFileId, int varid, int* dimids )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_vardimid( libId, varid, dimids );
#endif
    (void)backend;
    return nc_inq_vardimid( libId, varid, dimids );
}

int mbnc_inq_varnatts( int taggedFileId, int varid, int* nattsp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_varnatts( libId, varid, nattsp );
#endif
    (void)backend;
    return nc_inq_varnatts( libId, varid, nattsp );
}

int mbnc_inq_attname( int taggedFileId, int varid, int attnum, char* name )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_inq_attname( libId, varid, attnum, name );
#endif
    (void)backend;
    return nc_inq_attname( libId, varid, attnum, name );
}

int mbnc_inq_att( int taggedFileId, int varid, const char* name, nc_type* xtypep, size_t* lenp )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset tmp = 0;
        int rc         = ncmpi_inq_att( libId, varid, name, xtypep, &tmp );
        if( lenp ) *lenp = static_cast< size_t >( tmp );
        return rc;
    }
#endif
    (void)backend;
    return nc_inq_att( libId, varid, name, xtypep, lenp );
}

// ============================================================================
// Attribute get
// ============================================================================

int mbnc_get_att_text( int taggedFileId, int varid, const char* name, char* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_text( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_text( libId, varid, name, value );
}

int mbnc_get_att_int( int taggedFileId, int varid, const char* name, int* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_int( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_int( libId, varid, name, value );
}

int mbnc_get_att_short( int taggedFileId, int varid, const char* name, short* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_short( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_short( libId, varid, name, value );
}

int mbnc_get_att_long( int taggedFileId, int varid, const char* name, long* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_long( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_long( libId, varid, name, value );
}

int mbnc_get_att_float( int taggedFileId, int varid, const char* name, float* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_float( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_float( libId, varid, name, value );
}

int mbnc_get_att_double( int taggedFileId, int varid, const char* name, double* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_get_att_double( libId, varid, name, value );
#endif
    (void)backend;
    return nc_get_att_double( libId, varid, name, value );
}

// ============================================================================
// Attribute put
// ============================================================================

int mbnc_put_att_text( int taggedFileId, int varid, const char* name, size_t len, const char* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
        return ncmpi_put_att_text( libId, varid, name, static_cast< MPI_Offset >( len ), value );
#endif
    (void)backend;
    return nc_put_att_text( libId, varid, name, len, value );
}

int mbnc_put_att_int( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const int* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
        return ncmpi_put_att_int( libId, varid, name, xtype, static_cast< MPI_Offset >( len ), value );
#endif
    (void)backend;
    return nc_put_att_int( libId, varid, name, xtype, len, value );
}

int mbnc_put_att_short( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const short* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
        return ncmpi_put_att_short( libId, varid, name, xtype, static_cast< MPI_Offset >( len ), value );
#endif
    (void)backend;
    return nc_put_att_short( libId, varid, name, xtype, len, value );
}

int mbnc_put_att_float( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const float* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
        return ncmpi_put_att_float( libId, varid, name, xtype, static_cast< MPI_Offset >( len ), value );
#endif
    (void)backend;
    return nc_put_att_float( libId, varid, name, xtype, len, value );
}

int mbnc_put_att_double( int taggedFileId, int varid, const char* name, nc_type xtype, size_t len, const double* value )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
        return ncmpi_put_att_double( libId, varid, name, xtype, static_cast< MPI_Offset >( len ), value );
#endif
    (void)backend;
    return nc_put_att_double( libId, varid, name, xtype, len, value );
}

// ============================================================================
// Variable get_vara — collective by default on parallel PNetCDF
// ============================================================================

int mbnc_get_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_double_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_double( libId, varid, start, count, data );
}

int mbnc_get_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_int_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_int( libId, varid, start, count, data );
}

int mbnc_get_vara_long( int taggedFileId, int varid, const size_t* start, const size_t* count, long* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_long_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_long( libId, varid, start, count, data );
}

int mbnc_get_vara_text( int taggedFileId, int varid, const size_t* start, const size_t* count, char* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_text_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_text( libId, varid, start, count, data );
}

int mbnc_get_vars_double( int taggedFileId, int varid, const size_t* start, const size_t* count,
                          const ptrdiff_t* stride, double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims], st[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        for( int i = 0; i < ndims; ++i )
            st[i] = static_cast< MPI_Offset >( stride[i] );
        return ncmpi_get_vars_double_all( libId, varid, s, c, st, data );
    }
#endif
    (void)backend;
    return nc_get_vars_double( libId, varid, start, count, stride, data );
}

// ============================================================================
// Variable put_vara
// ============================================================================

int mbnc_put_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, const double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_double_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_double( libId, varid, start, count, data );
}

int mbnc_put_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_int_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_int( libId, varid, start, count, data );
}

int mbnc_put_vara_text( int taggedFileId, int varid, const size_t* start, const size_t* count, const char* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_text_all( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_text( libId, varid, start, count, data );
}

// ============================================================================
// Nonblocking get + wait
//
// PNetCDF: real ncmpi_iget_* request aggregation.
// Other backends: immediate blocking collective; *req = MBNC_REQ_NULL.
// mbnc_wait_all is a no-op for non-PNetCDF (returns NC_NOERR, writes
// NC_NOERR into statuses[] for any entries whose request is MBNC_REQ_NULL).
// ============================================================================

int mbnc_iget_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data,
                           int* req )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_iget_vara_double( libId, varid, s, c, data, req );
    }
#endif
    (void)backend;
    if( req ) *req = MBNC_REQ_NULL;
    return mbnc_get_vara_double( taggedFileId, varid, start, count, data );
}

int mbnc_iget_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data, int* req )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_iget_vara_int( libId, varid, s, c, data, req );
    }
#endif
    (void)backend;
    if( req ) *req = MBNC_REQ_NULL;
    return mbnc_get_vara_int( taggedFileId, varid, start, count, data );
}

// ============================================================================
// Independent-mode get/put — used inside begin_indep_data brackets on PNetCDF.
// On non-PNetCDF backends these route to the same plain nc_*_vara_* calls
// (per-var access mode is the caller's responsibility outside PNetCDF).
// ============================================================================

int mbnc_get_vara_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_double( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_double( libId, varid, start, count, data );
}

int mbnc_get_vara_int_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, int* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_int( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_int( libId, varid, start, count, data );
}

int mbnc_get_vara_long_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, long* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_long( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_long( libId, varid, start, count, data );
}

int mbnc_get_vara_text_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, char* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_get_vara_text( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_get_vara_text( libId, varid, start, count, data );
}

int mbnc_get_vars_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count,
                                const ptrdiff_t* stride, double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims], st[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        for( int i = 0; i < ndims; ++i )
            st[i] = static_cast< MPI_Offset >( stride[i] );
        return ncmpi_get_vars_double( libId, varid, s, c, st, data );
    }
#endif
    (void)backend;
    return nc_get_vars_double( libId, varid, start, count, stride, data );
}

int mbnc_put_vara_double_indep( int taggedFileId, int varid, const size_t* start, const size_t* count,
                                const double* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_double( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_double( libId, varid, start, count, data );
}

int mbnc_put_vara_int_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_int( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_int( libId, varid, start, count, data );
}

int mbnc_put_vara_text_indep( int taggedFileId, int varid, const size_t* start, const size_t* count, const char* data )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_put_vara_text( libId, varid, s, c, data );
    }
#endif
    (void)backend;
    return nc_put_vara_text( libId, varid, start, count, data );
}

int mbnc_iput_vara_double( int taggedFileId, int varid, const size_t* start, const size_t* count, const double* data,
                           int* req )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_iput_vara_double( libId, varid, s, c, data, req );
    }
#endif
    (void)backend;
    if( req ) *req = MBNC_REQ_NULL;
    return mbnc_put_vara_double( taggedFileId, varid, start, count, data );
}

int mbnc_iput_vara_int( int taggedFileId, int varid, const size_t* start, const size_t* count, const int* data,
                        int* req )
{
    const int libId         = mbnc_lib_id( taggedFileId );
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF )
    {
        MPI_Offset s[kMaxDims], c[kMaxDims];
        int ndims = 0;
        int rc    = to_mpi_offset_pair( libId, varid, start, count, s, c, &ndims );
        if( rc != NC_NOERR ) return rc;
        return ncmpi_iput_vara_int( libId, varid, s, c, data, req );
    }
#endif
    (void)backend;
    if( req ) *req = MBNC_REQ_NULL;
    return mbnc_put_vara_int( taggedFileId, varid, start, count, data );
}

int mbnc_wait_all( int taggedFileId, int nreq, int* requests, int* statuses )
{
    const NcBackend backend = mbnc_backend_of( taggedFileId );
#ifdef MOAB_HAVE_PNETCDF
    if( backend == NCB_PNETCDF ) return ncmpi_wait_all( mbnc_lib_id( taggedFileId ), nreq, requests, statuses );
#endif
    (void)backend;
    // No-op on non-PNetCDF backends: the corresponding iget_* already
    // performed the blocking call. Mark statuses as OK for sentinel reqs.
    if( statuses )
    {
        for( int i = 0; i < nreq; ++i )
        {
            if( !requests || requests[i] == MBNC_REQ_NULL ) statuses[i] = NC_NOERR;
        }
    }
    return NC_NOERR;
}

}  // namespace moab
