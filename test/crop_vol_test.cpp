/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/*
 * crop_vol_test.cpp
 * \brief Driver to create a volume from 2 surfaces (top and bottom), a user defined polyline
 * (contour) and a direction
 *
 */

#include "moab/Core.hpp"
#include <iostream>
#include <fstream>
#include <cstring>
#include "moab/FBEngine.hpp"
#include "moab/GeomTopoTool.hpp"
#include "TestUtil.hpp"

std::string filename_top;
std::string filename_bot;
std::string polygon_file_name;
double min_dot = 0.8;
std::string vol_file;

int number_tests_successful = 0;
int number_tests_failed     = 0;

using namespace moab;

ErrorCode volume_test( FBEngine* pFacet );

void handle_error_code( ErrorCode rv, int& number_failed, int& number_successful )
{
    if( rv == MB_SUCCESS )
    {
        std::cout << "Success";
        number_successful++;
    }
    else
    {
        std::cout << "Failure";
        number_failed++;
    }
}

int main( int argc, char* argv[] )
{
    filename_bot      = TestDir + "unittest/BedCrop2.h5m";
    filename_top      = TestDir + "unittest/SECrop2.h5m";  //"/polyPB.txt";
    polygon_file_name = TestDir + "unittest/poly14.txt";
    vol_file          = "volIce.h5m";  // output

    number_tests_successful = 0;
    bool remove_output      = true;

    if( argc == 1 )
    {
        std::cout << "Using default input files: " << filename_bot << " " << filename_top << " " << polygon_file_name
                  << " " << min_dot << " " << vol_file << std::endl;
        std::cout << "    default output file: " << vol_file << " will be deleted \n";
    }
    else if( argc == 6 )
    {
        remove_output     = false;
        filename_bot      = argv[1];
        filename_top      = argv[2];
        polygon_file_name = argv[3];
        min_dot           = atof( argv[4] );
        vol_file          = argv[5];
    }
    else
    {
        std::cout << " wrong number of arguments\n";
        return 1;  // fail
    }
    Core mbcore;
    Interface* mb = &mbcore;

    MB_CHK_SET_ERR( mb->load_file( filename_bot.c_str() ), "failed to load bed file" );

    // load the second face (we know we have one face in each file)
    MB_CHK_SET_ERR( mb->load_file( filename_top.c_str() ), "failed to load top file" );

    FBEngine* pFacet = new FBEngine( mb, NULL, true );  // smooth facetting, no OBB tree passed

    if( !pFacet ) return 1;  // error

    // should the init be part of constructor or not?
    // this is where the obb tree is constructed, and smooth faceting initialized, too.
    MB_CHK_SET_ERR( pFacet->Init(), "failed to initialize smoothing" );

    std::cout << "volume creation test: ";
    handle_error_code( volume_test( pFacet ), number_tests_failed, number_tests_successful );
    std::cout << "\n";

    if( remove_output )
    {
        remove( vol_file.c_str() );
    }
    return number_tests_failed;
}

ErrorCode volume_test( FBEngine* pFacet )
{
    // there should be 2 faces in the model
    std::ifstream datafile( polygon_file_name.c_str(), std::ifstream::in );
    if( !datafile )
    {
        std::cout << "can't read file\n";
        return MB_FAILURE;
    }
    //
    char temp[100];
    double direction[3];  // normalized
    double gridSize;
    datafile.getline( temp, 100 );  // first line

    // get direction and mesh size along polygon segments, from file
    sscanf( temp, " %lf %lf %lf %lf ", direction, direction + 1, direction + 2, &gridSize );
    // NORMALIZE(direction);// just to be sure

    std::vector< double > xyz;
    while( !datafile.eof() )
    {
        datafile.getline( temp, 100 );
        // int id = 0;
        double x, y, z;
        int nr = sscanf( temp, "%lf %lf %lf", &x, &y, &z );
        if( nr == 3 )
        {
            xyz.push_back( x );
            xyz.push_back( y );
            xyz.push_back( z );
        }
    }
    int sizePolygon = (int)xyz.size() / 3;
    if( sizePolygon < 3 )
    {
        std::cerr << " Not enough points in the polygon" << std::endl;
        return MB_FAILURE;
    }
    EntityHandle root_set;
    MB_CHK_SET_ERR( pFacet->getRootSet( &root_set ), "ERROR : getRootSet failed!" );

    int top = 2;  //  iBase_FACE;

    Range faces;
    MB_CHK_SET_ERR( pFacet->getEntities( root_set, top, faces ), "Failed to get faces in volume_test." );

    if( 2 != faces.size() )
    {
        std::cerr << "expected 2 faces, top and bottom\n";
        return MB_FAILURE;
    }
    EntityHandle newFace1;  // first test is with closed surface
    MB_CHK_SET_ERR( pFacet->split_surface_with_direction( faces[0], xyz, direction, /*closed*/ 1, min_dot, newFace1 ),
                    "Failed to crop first face." );

    EntityHandle newFace2;  // first test is with closed surface
    MB_CHK_SET_ERR( pFacet->split_surface_with_direction( faces[1], xyz, direction, /*closed*/ 1, min_dot, newFace2 ),
                    "Failed to crop second face." );

    // here we make the assumption that the edges, vertices are matching fine...
    EntityHandle volume;
    MB_CHK_SET_ERR( pFacet->create_volume_with_direction( newFace1, newFace2, direction, volume ),
                    "Failed to create volume." );

    Interface* mb = pFacet->moab_instance();
    pFacet->delete_smooth_tags();
    // duplicate -- extract a new geom topo tool
    GeomTopoTool* gtt       = pFacet->get_gtt();
    GeomTopoTool* duplicate = NULL;
    std::vector< EntityHandle > gents;
    gents.push_back( volume );
    MB_CHK_SET_ERR( gtt->duplicate_model( duplicate, &gents ), "Failed to extract volume." );
    EntityHandle newRootSet = duplicate->get_root_model_set();
    delete pFacet;
    pFacet = NULL;  // try not to write the obb tree
    delete duplicate;

    MB_CHK_ERR( mb->write_file( vol_file.c_str(), NULL, NULL, &newRootSet, 1 ) );

    return MB_SUCCESS;
}
