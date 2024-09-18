/*
 * ModifyMap.cpp
 *
 *  read an existing map, modify col variable using old global ids from source file
 */

#include "netcdf.h"
#include <iostream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"

#include "moab/ProgOptions.hpp"

using namespace moab;

int ncFile1;

#define ERR_NC( e )                                \
    {                                              \
        printf( "Error: %s\n", nc_strerror( e ) ); \
        exit( 2 );                                 \
    }

#define GET_DIM1( ncdim, name, val )                            \
    {                                                           \
        int gdfail = nc_inq_dimid( ncFile1, name, &( ncdim ) ); \
        if( NC_NOERR == gdfail )                                \
        {                                                       \
            size_t tmp_val;                                     \
            gdfail = nc_inq_dimlen( ncFile1, ncdim, &tmp_val ); \
            if( NC_NOERR != gdfail )                            \
            {                                                   \
                ERR_NC( gdfail )                                \
            }                                                   \
            else                                                \
                ( val ) = tmp_val;                              \
        }                                                       \
        else                                                    \
            ( val ) = 0;                                        \
    }

#define GET_VAR1( name, id, dims )                                     \
    {                                                                  \
        ( id )     = -1;                                               \
        int gvfail = nc_inq_varid( ncFile1, name, &( id ) );           \
        if( NC_NOERR == gvfail )                                       \
        {                                                              \
            int ndims;                                                 \
            gvfail = nc_inq_varndims( ncFile1, id, &ndims );           \
            if( NC_NOERR == gvfail )                                   \
            {                                                          \
                ( dims ).resize( ndims );                              \
                gvfail = nc_inq_vardimid( ncFile1, id, &( dims )[0] ); \
                if( NC_NOERR != gvfail )                               \
                {                                                      \
                    ERR_NC( gvfail )                                   \
                }                                                      \
            }                                                          \
        }                                                              \
    }

#define GET_1D_INT_VAR1( name, id, vals )                                               \
    {                                                                                   \
        GET_VAR1( name, id, vals );                                                     \
        if( -1 != ( id ) )                                                              \
        {                                                                               \
            size_t ntmp;                                                                \
            int ivfail = nc_inq_dimlen( ncFile1, ( vals )[0], &ntmp );                  \
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
            ( vals ).resize( ntmp );                                                    \
            size_t ntmp1 = 0;                                                           \
            ivfail       = nc_get_vara_int( ncFile1, id, &ntmp1, &ntmp, &( vals )[0] ); \
            if( NC_NOERR != ivfail )                                                    \
            {                                                                           \
                ERR_NC( ivfail )                                                        \
            }                                                                           \
        }                                                                               \
    }



int main( int argc, char* argv[] )
{

    ProgOptions opts;

    std::string inputfile1, inputfile2;
    opts.addOpt< std::string >( "existingMap,i", "input map file", &inputfile1 );
    opts.addOpt< std::string >( "meshFile,m", "mesh file with ids", &inputfile2 );

    opts.parseCommandLine( argc, argv );

    // Open netcdf/exodus file to read/write
    int fail = nc_open( inputfile1.c_str(), NC_NOWRITE, &ncFile1 );
    if( NC_NOWRITE != fail )
    {
        ERR_NC( fail )
    }
    Core moab;
    Interface* mb = &moab;
	EntityHandle sf;
	ErrorCode rval = mb->create_meshset( MESHSET_SET, sf );MB_CHK_ERR( rval );

	rval = mb->load_file( inputfile2.c_str(), &sf );MB_CHK_ERR( rval );
	// get all 2d cells, and reset global id
	Range cells;
	rval = mb->get_entities_by_dimension(0, 2, cells);MB_CHK_ERR( rval );
	Tag global_id_tag = mb->globalId_tag();

	std::vector<int> globalIds(cells.size());
	rval = mb->tag_get_data(global_id_tag, cells, &globalIds[0]);MB_CHK_SET_ERR( rval, "Can't get global id vals" );
	Tag oldIdTag;
    rval = mb->tag_get_handle( "OLD_GLOBAL_ID", oldIdTag );MB_CHK_SET_ERR( rval, "Can't get old id tag" );
    std::vector<int> oldIds(cells.size());

	rval = mb->tag_get_data(oldIdTag, cells, &oldIds[0]);MB_CHK_SET_ERR( rval, "Can't get old global id vals" );

	// construct a map between global ids and old global ids
	// will modify col variable accordingly
	std::map<int, int> idMap;
	for (size_t i=0; i<globalIds.size(); i++ )
	{
		idMap[globalIds[i]] = oldIds[i];
	}
	int temp_dim;
	int ns1, idcol1;
	GET_DIM1( temp_dim, "n_s", ns1 );
	std::vector< int > col1( ns1 );
	GET_1D_INT_VAR1( "col", idcol1, col1 );

	for (int i=0; i<ns1; i++)
		col1[i] = idMap[col1[i]];
	// write now the new col in the netcdf file
	size_t start=0;
	size_t count =ns1;
	int ierr = nc_put_vara_int(ncFile1, idcol1, &start, &count,  &col1[0]);
	std::cout << " to put new val for col:  ierr:" << ierr << "\n";
	// now rewrite the col variable, such that col will have old values
	ierr = nc_close(ncFile1);
	std::cout << " close the file  ierr:" << ierr << "\n";

}
