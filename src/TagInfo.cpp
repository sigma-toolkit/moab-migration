#include "TagInfo.hpp"
#include "moab/Error.hpp"
#include "moab/ErrorHandler.hpp"
#include <cstring> /* memcpy */
#include <cstdlib> /* realloc & free */
#include <cassert>

namespace moab
{

TagInfo::TagInfo( const char* name, int size, DataType type, const void* default_value, int default_value_size )
    : mDefaultValue( NULL ), mMeshValue( NULL ), mDefaultValueSize( default_value_size ), mMeshValueSize( 0 ),
      mDataSize( size ), dataType( type )
{
    if( default_value )
    {
        mDefaultValue = malloc( mDefaultValueSize );
        memcpy( mDefaultValue, default_value, mDefaultValueSize );
    }
    if( name ) mTagName = name;
}

TagInfo::~TagInfo()
{
    free( mDefaultValue );
    mDefaultValue     = 0;
    mDefaultValueSize = 0;
}

// enum MOAB_EXPORT DataType
// {
//     MB_TYPE_OPAQUE   = 0, /**< byte array */
//     MB_TYPE_INTEGER  = 1, /**< native 'int' type */
//     MB_TYPE_UNSIGNED_INTEGER  = 2, /**< native 'int' type */
//     MB_TYPE_LONG = 3, /**< native 'long' type */
//     MB_TYPE_UNSIGNED_LONG = 4, /**< native 'unsigned long' type */
//     MB_TYPE_UNSIGNED_LONG_LONG = 5, /**< native 'unsigned long' type */
//     MB_TYPE_FLOAT    = 6, /**< native 'float' type */
//     MB_TYPE_DOUBLE   = 7, /**< native 'double' type */
//     MB_TYPE_BIT      = 8, /**< mandatory type for tags with MB_TAG_BIT storage */
//     MB_TYPE_HANDLE   = 8, /**< EntityHandle */
//     MB_MAX_DATA_TYPE = MB_TYPE_HANDLE
// };
int TagInfo::size_from_data_type( DataType t )
{
    static const int sizes[] = { 1,  // MB_TYPE_OPAQUE
                                sizeof( int ), // MB_TYPE_INTEGER
                                sizeof( unsigned int ), // MB_TYPE_UNSIGNED_INTEGER
                                sizeof( long ), // MB_TYPE_LONG
                                sizeof( unsigned long ), // MB_TYPE_UNSIGNED_LONG
                                sizeof( unsigned long long ), // MB_TYPE_UNSIGNED_LONG_LONG
                                sizeof( float ), // MB_TYPE_FLOAT
                                sizeof( double ), // MB_TYPE_DOUBLE
                                1, // MB_TYPE_BIT
                                sizeof( EntityHandle ), // MB_TYPE_HANDLE
                                0 };
    return sizes[t];
}

bool TagInfo::equals_default_value( const void* data, int size ) const
{
    if( !get_default_value() ) return false;

    if( variable_length() && size != get_default_value_size() ) return false;

    if( !variable_length() && size >= 0 && size != get_size() ) return false;

    if( get_data_type() == MB_TYPE_BIT )
    {
        assert( get_size() <= 8 && get_default_value_size() == 1 );
        unsigned char byte1 = *reinterpret_cast< const unsigned char* >( data );
        unsigned char byte2 = *reinterpret_cast< const unsigned char* >( get_default_value() );
        unsigned char mask  = (unsigned char)( ( 1u << get_size() ) - 1 );
        return ( byte1 & mask ) == ( byte2 & mask );
    }
    else
    {
        return !memcmp( data, get_default_value(), get_default_value_size() );
    }
}

// Check that all lengths are valid multiples of the type size.
// Returns true if all lengths are valid, false otherwise.
bool TagInfo::check_valid_sizes( const int* sizes, int num_sizes ) const
{
    const unsigned size = size_from_data_type( get_data_type() );
    if( 1 == size ) return true;

    unsigned sum = 0;
    for( int i = 0; i < num_sizes; ++i )
        sum |= ( (unsigned)sizes[i] ) % size;

    return ( 0 == sum );
}

ErrorCode TagInfo::validate_lengths( Error* /* error_handler */, const int* lengths, size_t num_lengths ) const
{
    int bits = 0;
    if( variable_length() )
    {
        if( !lengths )
        {
            MB_SET_ERR( MB_VARIABLE_DATA_LENGTH, "No size specified for variable-length tag" );
        }
        const unsigned type_size = size_from_data_type( get_data_type() );
        if( type_size == 1 ) return MB_SUCCESS;
        for( size_t i = 0; i < num_lengths; ++i )
            bits |= lengths[i] % type_size;
    }
    else if( lengths )
    {
        for( size_t i = 0; i < num_lengths; ++i )
            bits |= lengths[i] - get_size();
    }
    if( 0 == bits ) return MB_SUCCESS;

    MB_SET_ERR( MB_INVALID_SIZE, "Tag data with invalid size" );
}

}  // namespace moab
