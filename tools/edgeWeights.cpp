/*
 * compareMaps.cpp
 * this tool will take 2 existing map files in nc format, and compare their sparse matrices
 * we can compare xc, yc, areas, fractions with ncdiff from nco
 * maybe there is another utility in nco, need to ask Charlie Zender
 *
 * example of usage:
 * ./mbcmpmaps -i map1.nc -j map2.nc
 * will look for row, col, S entries, and use eigen3 sparse matrix constructor
 *
 * can be built only if netcdf and eigen3 are available
 *
 *
 */
#include "moab/MOABConfig.h"

#ifndef MOAB_HAVE_HDF5
#error edgeWeights tool requires HDF5
#endif

// MOAB includes
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "DebugOutput.hpp"
#include "moab/MeshTopoUtil.hpp"

#include "nanoflann.hpp"

#ifdef MOAB_HAVE_MPI
// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

    bool
    closestPointDistanceToLine();

// And this is the "dataset to kd-tree" adaptor class:
template < typename Derived >
struct PointCloudAdaptor
{
    const Derived& obj;  //!< A const ref to the data set origin

    /// The constructor that sets the data set source
    PointCloudAdaptor( const Derived& obj_ ) : obj( obj_ ) {}

    /// CRTP helper method
    inline const Derived& derived() const
    {
        return obj;
    }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return derived().size();
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline double kdtree_get_pt( const size_t idx, const size_t dim ) const
    {
        const size_t offset = idx * 3;
        return derived()[offset + dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template < class BBOX >
    bool kdtree_get_bbox( BBOX& /*bb*/ ) const
    {
        return false;
    }

};  // end of PointCloudAdaptor

int main( int argc, char* argv[] )
{
    using namespace moab;
    ProgOptions opts;

    closestPointDistanceToLine();
    // return 0;

#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
#endif

    std::string sourceMeshFile, intersectionMeshFile, outputMeshFile( "output.h5m" );
    opts.addOpt< std::string >( "source", "source input mesh", &sourceMeshFile );
    opts.addOpt< std::string >( "intersection", "intersection input mesh", &intersectionMeshFile );
    opts.addOpt< std::string >( "output", "output mesh with edge weights", &outputMeshFile );

    opts.parseCommandLine( argc, argv );

    if( sourceMeshFile.empty() || intersectionMeshFile.empty() )
    {
        opts.printHelp();
        exit( 1 );
    }

    moab::Interface* mbCore = new( std::nothrow ) moab::Core;

    if( NULL == mbCore )
    {
        return 1;
    }

    ErrorCode rval;
    moab::EntityHandle sourceSet, intersectionSet;

    rval = mbCore->create_meshset( 0, sourceSet );MB_CHK_SET_ERR( rval, "Couldn't create master/slave set" );
    rval = mbCore->create_meshset( 0, intersectionSet );MB_CHK_SET_ERR( rval, "Couldn't create master/slave set" );

#ifdef MOAB_HAVE_MPI
    ParallelComm* pc = new ParallelComm( mbCore, MPI_COMM_WORLD );
    // if( debug ) options << "DEBUG_IO=1;CPUTIME;";
    const std::string options =
        ( pc->size() > 1 ? "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;" : "" );
#else
    const std::string options = "";
#endif

    // Load intersection mesh file
    std::cout << "Opening " << intersectionMeshFile << "\n";
    rval = mbCore->load_file( intersectionMeshFile.c_str(), &intersectionSet, options.c_str() );MB_CHK_ERR( rval );
    // rval = mbCore->load_file( intersectionMeshFile.c_str(), 0, options.c_str() );MB_CHK_ERR( rval );

    Range intersectionVertices, intersectionEdges, intersectionElems;
    rval = mbCore->get_entities_by_dimension( intersectionSet, 2, intersectionElems, true );MB_CHK_ERR( rval );
    // rval = mbCore->get_entities_by_dimension( intersectionSet, 0, intersectionVertices, true );MB_CHK_ERR( rval );
    rval = mbCore->get_connectivity( intersectionElems, intersectionVertices );MB_CHK_ERR( rval );
    rval =
        mbCore->get_adjacencies( intersectionVertices, 1 /*edges*/, true, intersectionEdges, moab::Interface::UNION );MB_CHK_ERR( rval );
    rval = moab::IntxUtils::ScaleToRadius( mbCore, intersectionSet, 1.0 );MB_CHK_ERR( rval );

    rval = mbCore->add_entities( intersectionSet, intersectionEdges );MB_CHK_ERR( rval );
    rval = pc->assign_global_ids( intersectionSet, 1, 1 );MB_CHK_ERR( rval );

    std::cout << "Intersection set contains: " << intersectionEdges.size() << " edges and "
              << intersectionVertices.size() << " vertices\n";

    // Load source mesh file
    std::cout << "Opening " << sourceMeshFile << "\n";
    rval = mbCore->load_file( sourceMeshFile.c_str(), &sourceSet, options.c_str() );MB_CHK_ERR( rval );
    rval = moab::IntxUtils::ScaleToRadius( mbCore, sourceSet, 1.0 );MB_CHK_ERR( rval );

    Tag edgeWeightTag;
    double defaultDouble = 0.0;
    rval = mbCore->tag_get_handle( "EdgeWeight", 1, MB_TYPE_DOUBLE, edgeWeightTag, MB_TAG_DENSE | MB_TAG_CREAT,
                                   &defaultDouble );MB_CHK_ERR( rval );
    // rval = mbCore->tag_get_data( edgeWeightTag, &rootset, 1, smat_metadata_glb );MB_CHK_ERR( rval );

    Tag sourceParentTag;
    rval = mbCore->tag_get_handle( "SourceParent", sourceParentTag );MB_CHK_ERR( rval );

    // rval = mbCore->get_adjacencies( intersectionElems, 1 /*edges*/, true, intersectionEdges, moab::Interface::UNION );
    // MB_CHK_ERR( rval );
    // MeshTopoUtil( mbCore ).get_bridge_adjacencies( intersectionVertices, 1, 1, intersectionEdges, 10 );
    // MB_CHK_SET_ERR( rval, "Failed to get bridge adjacencies" );
    // rval = mbCore->get_adjacencies( intersectionEdges, 2 /*elements*/, true, intersectionElems, moab::Interface::UNION );
    // MB_CHK_ERR( rval );

    Range sourceVertices, sourceEdges, sourceElems;
    rval = mbCore->get_entities_by_dimension( sourceSet, 0, sourceVertices );
    MB_CHK_ERR( rval );
    rval = mbCore->get_entities_by_dimension( sourceSet, 1, sourceEdges );
    MB_CHK_ERR( rval );
    rval = mbCore->get_entities_by_dimension( sourceSet, 2, sourceElems );
    MB_CHK_ERR( rval );

    std::cout << "Source set contains: " << sourceEdges.size() << " edges and " << sourceVertices.size()
              << " vertices\n";

    using PCKD = PointCloudAdaptor< std::vector< double > >;
    // std::vector< double > cloud( sourceVertices.size() * 3 );
    // rval = mbCore->get_coords( sourceVertices, cloud.data() );MB_CHK_ERR( rval );
    // const PCKD pckd( cloud );  // The adaptor

    // construct a kd-tree index:
    constexpr int dimension = 3;
    using sourceKdTree =
        nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor< double, PCKD >, PCKD, dimension >;
    // sourceKdTree index1( dimension, pckd, { 10 /* max leaf */ } );

    Tag gidTag = mbCore->globalId_tag();
    std::map< int, EntityHandle > gidSrcMap;
    int gidElem;
    for( auto it = sourceElems.begin(); it != sourceElems.end(); ++it )
    {
        EntityHandle elem = *it;
        rval              = mbCore->tag_get_data( gidTag, &elem, 1, &gidElem );MB_CHK_ERR( rval );
        gidSrcMap[gidElem] = elem;
    }

    auto edgeVolume = []( const double a[3], const double b[3] ) -> double {
        // Compute the distance between the two edges
        CartVect ab( b[0] - a[0], b[1] - a[1], b[2] - a[2] );
        return ab.length();
    };

    auto do_knn_search = []( const sourceKdTree& index, double queryPt[3] ) {
        // do a knn search
        const size_t num_results = 1;
        size_t ret_index;
        double out_dist_sqr;
        nanoflann::KNNResultSet< double > resultSet( num_results );
        // double query_pt[3] = { 0.5, 0.5, 0.5 };

        resultSet.init( &ret_index, &out_dist_sqr );
        index.findNeighbors( resultSet, queryPt );

        std::cout << "knnSearch(nn=" << num_results << "): \n";
        std::cout << "ret_index=" << ret_index << " out_dist_sqr=" << out_dist_sqr << std::endl;
    };

    // A function to calculate the distance from a point to a line defined by two endpoints in 3D
    auto pointToLineDistanceA = []( const double A[3], const double B[3],
                                   const double P[3] ) -> double {
        // Vector from A B
        std::array< double, 3 > AB = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };

        // Vector from A to P
        std::array< double, 3 > AP = { P[0] - A[0], P[1] - A[1], P[2] - A[2] };

        // Calculate the magnitude of AB
        double ABMagnitude = std::sqrt( AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2] );

        // If AB is zero vector (A and B are the same point), return distance from P to A
        if( ABMagnitude == 0 )
        {
            return std::sqrt( AP[0] * AP[0] + AP[1] * AP[1] + AP[2] * AP[2] );
        }

        // Normalize the direction vector AB
        std::array< double, 3 > ABNormalized = { AB[0] / ABMagnitude, AB[1] / ABMagnitude, AB[2] / ABMagnitude };

        // Calculate the projection of AP onto ABNormalized
        double t = ( AP[0] * ABNormalized[0] + AP[1] * ABNormalized[1] + AP[2] * ABNormalized[2] );

        // Clamp t to the range [0, 1] to find the closest point on the segment
        t = std::max( 0.0, std::min( 1.0, t ) );

        // Calculate the closest point Q on the line segment
        std::array< double, 3 > Q = { A[0] + t * ABNormalized[0], A[1] + t * ABNormalized[1],
                                      A[2] + t * ABNormalized[2] };

        // Calculate the distance from P to Q
        double distance = std::sqrt( ( P[0] - Q[0] ) * ( P[0] - Q[0] ) + ( P[1] - Q[1] ) * ( P[1] - Q[1] ) +
                                     ( P[2] - Q[2] ) * ( P[2] - Q[2] ) );

        return distance;
    };

    // A function to calculate the distance from a point to a line defined by two endpoints in 3D
    auto pointToLineDistance = []( const double A[3], const double B[3], const double P[3],
                                   std::pair< double, double >& distances ) -> bool {
        // Vector from A B: AB = B - A
        // std::array< double, 3 > AB = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
        moab::CartVect AB( B[0] - A[0], B[1] - A[1], B[2] - A[2] );

        // Vector from A to P: AP = P - A
        // std::array< double, 3 > AP = { P[0] - A[0], P[1] - A[1], P[2] - A[2] };
        moab::CartVect AP( P[0] - A[0], P[1] - A[1], P[2] - A[2] );
        // AP.normalize();

        // Calculate the magnitude of AB
        distances.first          = AP.length();
        distances.second         = distances.first;
        const double ABMagnitude = AB.length();

        // If AB is zero vector (A and B are the same point), return distance from P to A
        if( ABMagnitude == 0 ) return true;

        // Normalize the direction vector AB to get a unit vector
        AB.normalize();
        AP.normalize();

        // Calculate the projection of AP onto AB: t = AP (.) AB
        double t = ( AP[0] * AB[0] + AP[1] * AB[1] + AP[2] * AB[2] );

        // check if the cosine angle between the two vectors is 0 or 180 degrees
        if( fabs( t ) < 1.0 - 1e-13 ) return false;

        moab::CartVect BP( P[0] - B[0], P[1] - B[1], P[2] - B[2] );
        distances.second = BP.length();

        // point is on the line; so return the distance from the point to the line
        return true;
    };

    // Loop over intersection edges and see where they lie on the source mesh.
    // This is a completely independent process and embarassingly parallel
    // We can parallelize this loop with OpenMP as well.
    std::pair< double, double > distances;
    for( auto it = intersectionEdges.begin(); it != intersectionEdges.end(); ++it )
    {
        EntityHandle edge = *it;
        // Get the vertices of the edge
        const EntityHandle* edgeVertices;
        int nnodes;
        rval = mbCore->get_connectivity( edge, edgeVertices, nnodes, true );MB_CHK_ERR( rval );

        assert( nnodes == 2 );

        // Get the coordinates of the vertices
        double coords[6], edgeCenterCoords[3];
        rval = mbCore->get_coords( edgeVertices, 2, coords );MB_CHK_ERR( rval );
        rval = mbCore->get_coords( &edge, 1, edgeCenterCoords );MB_CHK_ERR( rval );

        std::vector< EntityHandle > intxAdjs;
        rval = mbCore->get_adjacencies( &edge, 1, 2, true, intxAdjs, moab::Interface::UNION );MB_CHK_ERR( rval );

        // assert( intxAdjs.size() <= 2 );

        std::vector< int > srcParents( intxAdjs.size() );
        rval = mbCore->tag_get_data( sourceParentTag, &intxAdjs[0], intxAdjs.size(), &srcParents[0] );MB_CHK_ERR( rval );

        // std::cout << "Number of edge adjacencies: " << intxAdjs.size() << ";" << srcParents[0] << std::endl;
        // Compute the edge weight
        double edgeWeight = 0.0;
        // if ( sourceVertices.index( edgeVertices[0] ) >= 0 ) // first vertex part of original source mesh
        // {
        //     // compute edge distance to the closest vertex on the source mesh and set the edgeWeight tag
        //     // continue;
        //     std::vector< EntityHandle > sourceAdjs;
        //     rval = mbCore->get_adjacencies( &edgeVertices[0], 1, 1 /*edges*/, true, sourceAdjs, moab::Interface::UNION );MB_CHK_ERR( rval );

        //     std::cout << "Found vertex A on source" << std::endl;

        //     // Loop over the source edges and find the closest edge that contains the other vertex
        //     for ( auto jit = sourceAdjs.begin(); jit != sourceAdjs.end(); ++jit )
        //     {
        //         EntityHandle sourceEdge = *jit;

        //         // if the adjacent edge is the same as the current edge, skip it
        //         if( sourceEdge == edge ) continue;

        //         const EntityHandle* sourceEdgeVertices;
        //         int nsourceEdgeNodes;
        //         rval = mbCore->get_connectivity( sourceEdge, sourceEdgeVertices, nsourceEdgeNodes, true );MB_CHK_ERR( rval );

        //         assert( nsourceEdgeNodes == 2 );

        //         // Get the coordinates of the vertices
        //         double sourceCoords[6];
        //         rval = mbCore->get_coords( sourceEdgeVertices, 2, sourceCoords );MB_CHK_ERR( rval );

        //         if( pointToLineDistance( sourceCoords, sourceCoords + 3, coords + 3, distances ) )
        //         {
        //             edgeWeight = edgeVolume( coords + 3, sourceCoords ) / edgeVolume( sourceCoords, sourceCoords + 3 );
        //             std::cout << "A: Found intersection edge coincident on source edge with distance: "
        //                       << distances.second << " and weight: " << edgeWeight << std::endl;
        //             break;
        //         }
        //     }
        // }
        // else if( sourceVertices.index( edgeVertices[1] ) >= 0 )  // second vertex part of original source mesh
        // {
        //     // compute edge distance to the closest vertex on the source mesh and set the edgeWeight tag
        //     // continue;
        //     std::vector< EntityHandle > sourceAdjs;
        //     rval = mbCore->get_adjacencies( &edgeVertices[1], 1, 1 /*edges*/, true, sourceAdjs, moab::Interface::UNION );MB_CHK_ERR( rval );

        //     std::cout << "Found vertex B on source" << std::endl;
        //     // Loop over the source edges and find the closest edge that contains the other vertex
        //     for ( auto jit = sourceAdjs.begin(); jit != sourceAdjs.end(); ++jit )
        //     {
        //         EntityHandle sourceEdge = *jit;

        //         // if the adjacent edge is the same as the current edge, skip it
        //         if( sourceEdge == edge ) continue;

        //         const EntityHandle* sourceEdgeVertices;
        //         int nsourceEdgeNodes;
        //         rval = mbCore->get_connectivity( sourceEdge, sourceEdgeVertices, nsourceEdgeNodes, true );MB_CHK_ERR( rval );

        //         assert( nsourceEdgeNodes == 2 );

        //         // Get the coordinates of the vertices
        //         double sourceCoords[6];
        //         rval = mbCore->get_coords( sourceEdgeVertices, 2, sourceCoords );MB_CHK_ERR( rval );

        //         if( pointToLineDistance( sourceCoords, sourceCoords + 3, coords, distances ) )
        //         {
        //             edgeWeight = edgeVolume( coords, sourceCoords + 3 ) / edgeVolume( sourceCoords, sourceCoords + 3 );
        //             std::cout << "B: Found intersection edge coincident on source edge with distance: "
        //                       << distances.second << " and weight: " << edgeWeight << std::endl;
        //         }
        //     }
        // }
        // else
        {
            // // none of the edges are part of the original source mesh itself
            // // but they could be part of the edges in the source mesh
            // // a more complicated case;

            for (auto isrcParent : srcParents)
            // EntityHandle isrcParent = srcParents[0];
            {
                auto srcPIter = gidSrcMap.find(isrcParent);

                if (srcPIter == gidSrcMap.end())
                {
                    std::cerr << "Parent " << isrcParent << " not found for edge: " << edge << std::endl;
                    continue;
                }

                EntityHandle srcEnt = srcPIter->second;
                const EntityHandle* sourceElemVertices;
                int nsourceElemNodes;
                rval = mbCore->get_connectivity( srcEnt, sourceElemVertices, nsourceElemNodes, true );MB_CHK_ERR( rval );

                std::vector< EntityHandle > sourceAdjs;
                rval = mbCore->get_adjacencies( sourceElemVertices, nsourceElemNodes, 1 /*edges*/, true, sourceAdjs, moab::Interface::UNION );
                MB_CHK_ERR( rval );

                // std::cout << "Parent: " << srcEnt << ", edges: " << sourceAdjs.size() << std::endl;

                // Loop over the source edges and find the closest edge that contains the other vertex
                for( auto jit = sourceAdjs.begin(); jit != sourceAdjs.end(); ++jit )
                {
                    EntityHandle sourceEdge = *jit;

                    // if the adjacent edge is the same as the current edge, skip it
                    if( sourceEdge == edge )
                    {
                        std::cout << "Full weight for edge: " << edge << std::endl;
                        edgeWeight = 1.0;
                        // break;
                    }
                    else
                    {
                        const EntityHandle* sourceEdgeVertices;
                        int nsourceEdgeNodes;
                        rval = mbCore->get_connectivity( sourceEdge, sourceEdgeVertices, nsourceEdgeNodes, true );MB_CHK_ERR( rval );

                        assert( nsourceEdgeNodes == 2 );

                        // Get the coordinates of the vertices
                        double sourceCoords[6];
                        rval = mbCore->get_coords( sourceEdgeVertices, 2, sourceCoords );MB_CHK_ERR( rval );

                        if( // pointToLineDistance( sourceCoords, sourceCoords + 3, edgeCenterCoords, distances ) ||
                            ( pointToLineDistance( sourceCoords, sourceCoords + 3, coords, distances ) ||
                              pointToLineDistance( sourceCoords, sourceCoords + 3, coords + 3, distances ) ) )
                        {
                            // edgeWeight =
                            //     edgeVolume( coords, coords + 3 ) / edgeVolume( sourceCoords, sourceCoords + 3 );
                            // std::cout << "A. Found intersection edge coincident on source edge with distance: "
                            //           << distances.first << " and weight: " << edgeWeight << std::endl;
                            edgeWeight = 1.0;
                            // break;
                        }
                        // std::cout << "Didn't find edge midpoint in source: " << distances.first << std::endl;
                        // if( pointToLineDistance( sourceCoords, sourceCoords + 3, coords, distances ) &&
                        //     pointToLineDistance( sourceCoords, sourceCoords + 3, coords + 3, distances ) )
                        // {
                        //     edgeWeight =
                        //         edgeVolume( coords, coords + 3 ) / edgeVolume( sourceCoords, sourceCoords + 3 );
                        //     std::cout << "A. Found intersection edge coincident on source edge with distance: "
                        //             << distances.first << " and weight: " << edgeWeight << std::endl;
                        //     break;
                        // }
                        // else if( pointToLineDistance( sourceCoords, sourceCoords + 3, coords, distances )  )
                        // {
                        //     edgeWeight = distances.first / edgeVolume( sourceCoords, sourceCoords + 3 );
                        //     std::cout << "B. Found intersection edge coincident on source edge with distance: "
                        //                 << distances.first << " and weight: " << edgeWeight << std::endl;
                        //     break;
                        // }
                        // else if( pointToLineDistance( sourceCoords, sourceCoords + 3, coords + 3, distances ) )
                        // {
                        //     edgeWeight = distances.first / edgeVolume( sourceCoords, sourceCoords + 3 );
                        //     std::cout << "C. Found intersection edge coincident on source edge with distance: "
                        //             << distances.first << " and weight: " << edgeWeight << std::endl;
                        //     break;
                        // }
                        // else continue;
                    }
                }
            }
        }

        // edgeWeight = 1.0;
        // Compute the edge weight
        rval = mbCore->tag_set_data( edgeWeightTag, &edge, 1, &edgeWeight );MB_CHK_ERR( rval );
    }

    moab::Range sourceEdgesWithWeights = moab::intersect( sourceEdges, intersectionEdges );
    std::cout << "Source edges: " << sourceEdges.size() << "; Intersection edges: " << intersectionEdges.size()
              << "; Source edges intersected: " << sourceEdgesWithWeights.size() << std::endl;
    // std::vector< double > edgeWeights( intersectionEdges.size(), 1.0 );
    // rval = mbCore->tag_set_data( edgeWeightTag, intersectionEdges, edgeWeights.data() );MB_CHK_ERR( rval );

    rval = mbCore->write_file( outputMeshFile.c_str(), 0, 0, &intersectionSet, 1 );MB_CHK_ERR( rval );

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}

bool closestPointDistanceToLine( )
{

    // A function to calculate the distance from a point to a line defined by two endpoints in 3D
    auto pointToLineDistance = []( const double A[3], const double B[3], const double P[3] ) -> double {
        // Vector from A B: AB = B - A
        // std::array< double, 3 > AB = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
        moab::CartVect AB( B[0] - A[0], B[1] - A[1], B[2] - A[2] );

        // Vector from A to P: AP = P - A
        // std::array< double, 3 > AP = { P[0] - A[0], P[1] - A[1], P[2] - A[2] };
        moab::CartVect AP( P[0] - A[0], P[1] - A[1], P[2] - A[2] );
        // AP.normalize();

        // Calculate the magnitude of AB
        const double ABMagnitude = AB.length();
        const double APMagnitude = AP.length();

        // If AB is zero vector (A and B are the same point), return distance from P to A
        if( ABMagnitude == 0 ) return APMagnitude;

        // Normalize the direction vector AB to get a unit vector
        AB.normalize();
        AP.normalize();

        // Calculate the projection of AP onto AB: t = AP (.) AB
        double t = ( AP[0] * AB[0] + AP[1] * AB[1] + AP[2] * AB[2] );

        std::cout.precision( 15 );
        std::cout << "t: " << t << std::endl;
        if ( fabs(t) < 1-1e-15 )
            return false;

        const double q = t * t;

        // if( q <= 1.0 ) point p( x, y ) is on the line
        // else p( x, y ) is not on line
        std::cout << "Distance from point to first point in line segment: " << APMagnitude << std::endl;
        // std::cout << "q: " << q << std::endl;
        return true;

        // // Clamp t to the range [0, 1] to find the closest point on the segment
        // // This should give the barycentric coordinates of the projection of P onto the line segment
        // t = std::max( 0.0, std::min( 1.0, t ) );

        // // Calculate the closest point Q on the line segment
        // CartVect Q = { A[0] + t * AB[0], A[1] + t * AB[1],
        //                               A[2] + t * AB[2] };

        // // Calculate the distance from P to Q
        // double distance = std::sqrt( ( P[0] - Q[0] ) * ( P[0] - Q[0] ) + ( P[1] - Q[1] ) * ( P[1] - Q[1] ) +
        //                              ( P[2] - Q[2] ) * ( P[2] - Q[2] ) );

        // return distance;

    };

    //  usage
    double xP[3] = { 2.0, 2.0, 2.0 };  // Point
    double xQ[3] = { -2.5, -2.5, -2.5 };  // Point
    double xA[3] = { 1.0, 1.0, 1.0 };  // First endpoint of line
    double xB[3] = { 3.0, 3.0, 3.0 };  // Second endpoint of the line

    if( pointToLineDistance( xA, xB, xP ) )
        std::cout << "Point on the line!" << std::endl;
    else
        std::cout << "Point not on the line!" << std::endl;

    if( pointToLineDistance( xA, xB, xQ ) )
        std::cout << "Point on the line!" << std::endl;
    else
        std::cout << "Point not on the line!" << std::endl;
}
