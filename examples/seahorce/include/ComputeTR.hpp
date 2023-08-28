#ifndef __compute_tempestremap_hpp__
#define __compute_tempestremap_hpp__

#include "RemapMPASROMS.hpp"
#include "TempestRemapAPI.h"
#include "FiniteVolumeTools.h"

const char* template_map_output_filename = "mpas_roms_map_";

moab::ErrorCode CloneToTRMesh( moab::Interface* m_interface, Mesh& mesh, moab::EntityHandle mesh_set )
{
    using namespace moab;

    moab::ErrorCode rval;
    moab::Range elems, verts;

    NodeVector& nodes = mesh.nodes;
    FaceVector& faces = mesh.faces;

    rval = m_interface->get_entities_by_dimension( mesh_set, 2, elems );MB_CHK_ERR( rval );

    // resize the number of elements in Tempest mesh
    faces.resize( elems.size() );

    // let us now get the vertices from all the elements
    rval = m_interface->get_connectivity( elems, verts );MB_CHK_ERR( rval );
    if( verts.size() == 0 )
    {
        rval = m_interface->get_entities_by_dimension( mesh_set, 0, verts );MB_CHK_ERR( rval );
    }
    // assert(verts.size() > 0); // If not, this may be an invalid mesh ! possible for unbalanced
    // loads

    std::map< EntityHandle, int > indxMap;
    {
        int j = 0;
        for( Range::iterator it = verts.begin(); it != verts.end(); it++ )
            indxMap[*it] = j++;
    }

    for( size_t iface = 0; iface < elems.size(); ++iface )
    {
        Face& face           = faces[iface];
        EntityHandle ehandle = elems[iface];

        // get the connectivity for each edge
        const EntityHandle* connectface;
        int nnodesf;
        rval = m_interface->get_connectivity( ehandle, connectface, nnodesf );MB_CHK_ERR( rval );

        face.edges.resize( nnodesf );
        for( int iverts = 0; iverts < nnodesf; ++iverts )
        {
            int indx = indxMap[connectface[iverts]];
            assert( indx >= 0 );
            face.SetNode( iverts, indx );
        }
    }

    size_t nnodes = verts.size();
    nodes.resize( nnodes );

    // Set the data for the vertices
    std::vector< double > coordx( nnodes ), coordy( nnodes ), coordz( nnodes );
    rval = m_interface->get_coords( verts, &coordx[0], &coordy[0], &coordz[0] );MB_CHK_ERR( rval );
    for( size_t inode = 0; inode < nnodes; ++inode )
    {
        Node& node = nodes[inode];
        node.x     = coordx[inode];
        node.y     = coordy[inode];
        node.z     = coordz[inode];
    }
    coordx.clear();
    coordy.clear();
    coordz.clear();

    mesh.RemoveZeroEdges();
    mesh.RemoveCoincidentNodes();

    // Generate reverse node array and edge map
    mesh.ConstructEdgeMap( false );
    // mesh.ConstructReverseNodeArray();

    // mesh.Validate();

    return MB_SUCCESS;
}

moab::ErrorCode LoadTempestRemapWeights( RuntimeContext& context,
                                         moab::EntityHandle src_set,
                                         moab::EntityHandle tgt_set,
                                         std::string strMethod )
{
    CloneToTRMesh( context.moab_interface, context.meshInput, src_set );
    CloneToTRMesh( context.moab_interface, context.meshOutput, tgt_set );

    context.meshInput.ConstructEdgeMap();
    context.meshOutput.ConstructEdgeMap();

    // load the 2D intersection mesh from disk
    // context.meshOverlap = Mesh( "mesh_intersection.g" );

    // next read the map file
    NcError ncerror( NcError::silent_nonfatal );
    std::string map_output_filename =
        std::string( template_map_output_filename ) + ( strMethod.size() ? strMethod : "fv" ) + ".nc";
    std::cout << "Reading TempestRemap map file: " << map_output_filename << "\n";
    context.weightMap.Read( map_output_filename );

    return moab::MB_SUCCESS;
}

moab::ErrorCode ComputeTempestRemapWeights( RuntimeContext& context,
                                            moab::EntityHandle src_set,
                                            moab::EntityHandle tgt_set,
                                            std::string strMethod,
                                            bool ensureMonotonicity )
{
    // err = remapper.ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( err );
    // err = remapper.ConvertMeshToTempest( moab::Remapper::TargetMesh );MB_CHK_ERR( err );
    CloneToTRMesh( context.moab_interface, context.meshInput, src_set );
    CloneToTRMesh( context.moab_interface, context.meshOutput, tgt_set );

    context.meshInput.ConstructEdgeMap();
    context.meshOutput.ConstructEdgeMap();

    // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
    if( context.proc_id == 0 )
        std::cout << "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes" << std::endl;
    // err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );
    bool concaveMeshA = false, concaveMeshB = false, allowNoOverlap = true, verbose = false;
    int ierr =
        GenerateOverlapWithMeshes( context.meshInput, context.meshOutput, context.meshOverlap, "" /*outFilename*/,
                                   "Netcdf4", "exact", concaveMeshA, concaveMeshB, allowNoOverlap, verbose );
    if( ierr )
    {
        MB_CHK_SET_ERR( moab::MB_FAILURE, "TempestRemap: Can't compute the intersection of meshes on the sphere" );
    }

    if( context.proc_id == 0 ) std::cout << "\nSetup computation of weights" << std::endl;

    // Call to generate the remapping weights with the tempest meshes
    std::string map_output_filename =
        std::string( template_map_output_filename ) + ( strMethod.size() ? strMethod : "fv" ) + ".nc";
    GenerateOfflineMapAlgorithmOptions mapOptions;
    mapOptions.nPin             = 1;
    mapOptions.nPout            = 1;
    mapOptions.fSourceConcave   = false;
    mapOptions.fTargetConcave   = false;
    mapOptions.strMethod        = strMethod;  // invdist, bilin, intbilin, delaunay
    mapOptions.fMonotone        = ensureMonotonicity;
    mapOptions.fNoCorrectAreas  = false;
    mapOptions.fNoCheck         = true;
    mapOptions.strOutputMapFile = map_output_filename;  // ask TR to write it out
    mapOptions.strOutputFormat  = "Netcdf4";

    if( context.proc_id == 0 ) std::cout << "Compute weights with TempestRemap" << std::endl;
    ierr = GenerateOfflineMapWithMeshes( context.meshInput,    // Mesh inputMesh
                                         context.meshOutput,   // Mesh outputMesh,
                                         context.meshOverlap,  // Mesh overlapMesh,
                                         "fv",                 // std::string inputDiscretization,
                                         "fv",                 // std::string outputDiscretization,
                                         mapOptions,           // const GenerateOfflineMapAlgorithmOptions& options
                                         context.weightMap );

    // check the generated weights and output information
    {
        const double dNormalTolerance = 1.0E-8;
        const double dStrictTolerance = 1.0E-12;
        context.weightMap.CheckMap( true, true, ensureMonotonicity, dNormalTolerance, dStrictTolerance );
    }

    // Write the map to disk
#ifdef WRITE_MAP_FILE
    {
        // First write out the overlap mesh to disk
        context.meshOverlap.Write( "mesh_intersection.g" );

        // Next prepare set of attributes to add to the map NC file
        typedef std::map< std::string, std::string > AttributeMap;
        typedef AttributeMap::value_type AttributePair;

        AttributeMap mapAttributes;
        mapAttributes.insert( AttributePair( "grid_file_src", "mpas_grid.h5m" ) );
        mapAttributes.insert( AttributePair( "grid_file_dst", "roms_grid.h5m" ) );
        mapAttributes.insert( AttributePair( "grid_file_ovr", "mesh_intersection.g" ) );
        mapAttributes.insert(
            AttributePair( "concave_src", ( mapOptions.fSourceConcave ) ? ( "true" ) : ( "false" ) ) );
        mapAttributes.insert(
            AttributePair( "concave_dst", ( mapOptions.fTargetConcave ) ? ( "true" ) : ( "false" ) ) );
        if( mapOptions.strSourceMeta != "" )
        {
            mapAttributes.insert( AttributePair( "meta_src", mapOptions.strSourceMeta ) );
        }
        if( mapOptions.strTargetMeta != "" )
        {
            mapAttributes.insert( AttributePair( "meta_dst", mapOptions.strTargetMeta ) );
        }
        mapAttributes.insert( AttributePair( "type_src", "fv" ) );
        mapAttributes.insert( AttributePair( "type_dst", "fv" ) );
        mapAttributes.insert( AttributePair( "np_src", std::to_string( (long long)mapOptions.nPin ) ) );
        mapAttributes.insert( AttributePair( "np_dst", std::to_string( (long long)mapOptions.nPout ) ) );
        mapAttributes.insert( AttributePair( "mono", ( mapOptions.fMonotone ) ? ( "true" ) : ( "false" ) ) );
        mapAttributes.insert( AttributePair( "nobubble", "false" ) );
        mapAttributes.insert( AttributePair( "nocorrectareas", "false" ) );
        mapAttributes.insert( AttributePair( "noconserve", "false" ) );
        mapAttributes.insert( AttributePair( "sparse_constraints", "false" ) );
        mapAttributes.insert( AttributePair( "method", mapOptions.strMethod ) );
        mapAttributes.insert( AttributePair( "version", "RemapMPASROMS v0.1" ) );

        if( context.proc_id == 0 ) std::cout << "\nWrite the weights to " << map_output_filename << std::endl;
        context.weightMap.Write( mapOptions.strOutputMapFile, mapAttributes, NcFile::Netcdf4Classic );

        // // Write the map file to disk in parallel using either HDF5 or SCRIP interface
        // err = weightMap.WriteParallelMap( output_filename.c_str() );MB_CHK_ERR( err );
    }
#endif
    return moab::MB_SUCCESS;
}

double ApplyCAASLimiting( OfflineMap& mapOperator,
                          Mesh& meshInput,
                          Mesh& meshOverlap,
                          const int nPin,
                          DataArray1D< double >& dataInDouble,
                          DataArray1D< double >& dataOutDouble,
                          bool useCAASLocal )
{
    const size_t nSourceCount                   = dataInDouble.GetRows();
    const size_t nTargetCount                   = dataOutDouble.GetRows();
    const DataArray1D< double >& m_dSourceAreas = mapOperator.GetSourceAreas();
    const DataArray1D< double >& m_dTargetAreas = mapOperator.GetTargetAreas();

    // Announce input mass
    double dSourceMass = 0.0;
    double dSourceMin  = dataInDouble[0];
    double dSourceMax  = dataInDouble[0];
    for( size_t i = 0; i < nSourceCount; i++ )
    {
        dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
        dSourceMax = fmax( dSourceMax, dataInDouble[i] );
        dSourceMin = fmin( dSourceMin, dataInDouble[i] );
    }

    // Apply the offline map to the data
    {
        DataArray1D< double > x( nTargetCount );
        DataArray1D< double > dataLowerBound( nTargetCount );
        DataArray1D< double > dataUpperBound( nTargetCount );

        double dMassDiff = dSourceMass;

        double dTargetMin = dataOutDouble[0];
        double dTargetMax = dataOutDouble[0];
        for( size_t i = 0; i < nTargetCount; i++ )
        {
            dMassDiff -= dataOutDouble[i] * m_dTargetAreas[i];
            dTargetMax = fmax( dTargetMax, dataOutDouble[i] );
            dTargetMin = fmin( dTargetMin, dataOutDouble[i] );
        }

        // Early exit if the values are monotone already.
        if( dTargetMax <= dSourceMax && dTargetMin <= dSourceMin ) return 0.0;

        if( useCAASLocal )
        {
            double dMinI;
            double dMaxI;

            int nTargetFaces = nTargetCount;
            std::vector< double > vecLocalUpperBound( nTargetCount );
            std::vector< double > vecLocalLowerBound( nTargetCount );

            std::vector< std::vector< int > > vecSourceOvTarget( nTargetFaces );
            for( size_t i = 0; i < meshOverlap.faces.size(); i++ )
            {

                int ixT = meshOverlap.vecTargetFaceIx[i];
                int ixS = meshOverlap.vecSourceFaceIx[i];
                vecSourceOvTarget[ixT].push_back( ixS );
            }

            //FV to FV
            {
                for( size_t i = 0; i < nTargetCount; i++ )
                {
                    if( !vecSourceOvTarget[i].size() ) continue;
                    dMaxI = dataInDouble[vecSourceOvTarget[i][0]];
                    dMinI = dataInDouble[vecSourceOvTarget[i][0]];

                    //Compute max over interstecting source faces

                    for( size_t j = 0; j < vecSourceOvTarget[i].size(); j++ )
                    {
                        int k = vecSourceOvTarget[i][j];
                        dMaxI = fmax( dMaxI, dataInDouble[k] );
                        dMinI = fmin( dMinI, dataInDouble[k] );
                    }

                    if( useCAASLocal )
                    {
                        double dMaxIAdj = dMaxI;
                        double dMinIAdj = dMinI;

                        AdjacentFaceVector vecAdjFaces;

                        GetAdjacentFaceVectorByEdge( meshInput, vecSourceOvTarget[i][0], ( nPin + 1 ) * ( nPin + 1 ),
                                                     vecAdjFaces );

                        //Compute max over neighboring faces
                        for( size_t j = 0; j < vecAdjFaces.size(); j++ )
                        {
                            int k = vecAdjFaces[j].first;

                            dMaxIAdj = fmax( dMaxIAdj, dataInDouble[k] );
                            dMinIAdj = fmin( dMinIAdj, dataInDouble[k] );
                        }

                        vecLocalLowerBound[i] = dMinIAdj;
                        vecLocalUpperBound[i] = dMaxIAdj;
                    }
                    else
                    {
                        vecLocalLowerBound[i] = dMinI;
                        vecLocalUpperBound[i] = dMaxI;
                    }
                }
            }

            for( size_t i = 0; i < dataLowerBound.GetRows(); i++ )
            {
                dataLowerBound[i] = vecLocalLowerBound[i] - dataOutDouble[i];
                dataUpperBound[i] = vecLocalUpperBound[i] - dataOutDouble[i];
            }

        }     // if( useCAASLocal )
        else  // useCAASGlobal
        {
            for( size_t i = 0; i < nTargetCount; i++ )
            {
                dataLowerBound[i] = dSourceMin - dataLowerBound[i];
                dataUpperBound[i] = dSourceMax - dataUpperBound[i];
            }
        }

        // Invoke CAAS application on the offline map
        mapOperator.CAAS( dataOutDouble, dataLowerBound, dataUpperBound, dMassDiff );
    }

    // Announce output mass
    double dTargetMass = 0.0;
    double dTargetMin  = dataOutDouble[0];
    double dTargetMax  = dataOutDouble[0];
    for( size_t i = 0; i < nTargetCount; i++ )
    {
        dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
        if( dataOutDouble[i] < dTargetMin )
        {
            dTargetMin = dataOutDouble[i];
        }
        if( dataOutDouble[i] > dTargetMax )
        {
            dTargetMax = dataOutDouble[i];
        }
    }

    return ( dTargetMass - dSourceMass );
}

double ApplyCAASLimiting_ABC( OfflineMap& mapOperator,
                              Mesh& meshInput,
                              Mesh& meshOverlap,
                              const int nPin,
                              DataArray1D< double >& dataInDouble,
                              DataArray1D< double >& dataOutDouble,
                              bool useCAAS,
                              bool useCAASLocal )
{

    const int nSourceCount                      = dataInDouble.GetRows();
    const int nTargetCount                      = dataOutDouble.GetRows();
    const DataArray1D< double >& m_dSourceAreas = mapOperator.GetSourceAreas();
    const DataArray1D< double >& m_dTargetAreas = mapOperator.GetTargetAreas();

    // Announce input mass
    double dSourceMass = 0.0;
    double dSourceMin  = dataInDouble[0];
    double dSourceMax  = dataInDouble[0];
    for( int i = 0; i < nSourceCount; i++ )
    {
        dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
        if( dataInDouble[i] < dSourceMin )
        {
            dSourceMin = dataInDouble[i];
        }
        if( dataInDouble[i] > dSourceMax )
        {
            dSourceMax = dataInDouble[i];
        }
    }

    // Apply the offline map to the data

    if( useCAASLocal || useCAAS )
    {
        DataArray1D< double > l = dataOutDouble;
        DataArray1D< double > u = dataOutDouble;
        DataArray1D< double > x( nTargetCount );
        double b = dSourceMass;

        for( size_t i = 0; i < l.GetRows(); i++ )
        {
            b -= dataOutDouble[i] * m_dTargetAreas[i];
        }

        if( useCAASLocal )
        {
            // int GLLSizeIn  = 0;  // FV
            // int GLLSizeOut = 0;  // FV
            // int pOut       = dataGLLNodesOut.GetSize( 0 );
            // int qOut       = dataGLLNodesOut.GetSize( 1 );
            // int pIn        = dataGLLNodesIn.GetSize( 0 );
            // int qIn        = dataGLLNodesIn.GetSize( 1 );
            // int pIn          = nPin;
            double f_maxI    = 0.0;
            double f_minI    = 0.0;
            int nTargetFaces = nTargetCount;

            std::vector< std::vector< int > > SourceOvTarget( nTargetFaces );

            for( size_t i = 0; i < meshOverlap.faces.size(); i++ )
            {
                int ixT = meshOverlap.vecTargetFaceIx[i];
                int ixS = meshOverlap.vecSourceFaceIx[i];
                SourceOvTarget[ixT].push_back( ixS );
            }

            std::vector< double > local_UB( nTargetCount );
            std::vector< double > local_LB( nTargetCount );

            for( int i = 0; i < nTargetCount; i++ )
            {
                AdjacentFaceVector vecAdjFaces;

                if( !SourceOvTarget[i].size() ) continue;

                GetAdjacentFaceVectorByEdge( meshInput, SourceOvTarget[i][0], ( nPin + 1 ) * ( nPin + 1 ),
                                             vecAdjFaces );
                f_maxI = dataInDouble[vecAdjFaces[0].first];
                f_minI = dataInDouble[vecAdjFaces[0].first];
                for( size_t j = 0; j < vecAdjFaces.size(); j++ )
                {
                    int k  = vecAdjFaces[j].first;
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                for( size_t j = 0; j < SourceOvTarget[i].size(); j++ )
                {

                    int k  = SourceOvTarget[i][j];
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                // f_minI=fmax(f_minI,0.0);

                local_UB[i] = f_maxI;
                local_LB[i] = f_minI;
            }

            // double mt = 0.0;
            // for( int i = 0; i < nTargetCount; i++ )
            // {
            //     mt += m_dTargetAreas[i] * ( local_LB[i] - dataOutDouble[i] );
            // }

            for( size_t i = 0; i < l.GetRows(); i++ )
            {
                l[i] = local_LB[i] - l[i];
                u[i] = local_UB[i] - u[i];
            }

            // Adjust mass of lower bound if greater than b
            double mL = 0.0;

            for( int i = 0; i < nTargetCount; i++ )
            {
                mL += m_dTargetAreas[i] * l[i];
            }
            if( mL > b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mL   = mL - m_dTargetAreas[i] * l[i] + m_dTargetAreas[i] * ( dSourceMin - dataOutDouble[i] );
                    l[i] = dSourceMin - dataOutDouble[i];
                    if( mL < b )
                    {
                        break;
                    }
                }
            }

            // Adjust mass of upper bound if less than b
            double mU = 0.0;

            if( mU < b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mU   = mU - m_dTargetAreas[i] * u[i] + m_dTargetAreas[i] * ( dSourceMax - dataOutDouble[i] );
                    u[i] = dSourceMax - dataOutDouble[i];
                    if( mU > b )
                    {
                        break;
                    }
                }
            }
        }  // if( useCAASLocal )
        else if( useCAAS )
        {
            for( size_t i = 0; i < l.GetRows(); i++ )
            {
                l[i] = dSourceMin - l[i];
                u[i] = dSourceMax - u[i];
            }
        }

        // Invoke CAAS application on the offline map
        mapOperator.CAAS( x, l, u, b );

        // Add correction
        for( size_t i = 0; i < l.GetRows(); i++ )
        {
            dataOutDouble[i] += x[i];
        }
    }

    // Announce output mass
    double dTargetMass = 0.0;
    double dTargetMin  = dataOutDouble[0];
    double dTargetMax  = dataOutDouble[0];
    for( int i = 0; i < nTargetCount; i++ )
    {
        dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
        if( dataOutDouble[i] < dTargetMin )
        {
            dTargetMin = dataOutDouble[i];
        }
        if( dataOutDouble[i] > dTargetMax )
        {
            dTargetMax = dataOutDouble[i];
        }
    }

    return ( dTargetMass - dSourceMass );
}

#endif  // __compute_tempestremap_hpp__
