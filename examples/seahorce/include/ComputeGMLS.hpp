#ifndef __computeGMLS_hpp__
#define __computeGMLS_hpp__

#include "RemapMPASROMS.hpp"

// Include compadre headers
#include <Compadre_Config.h>
#include <Compadre_GMLS.hpp>
#include <Compadre_Evaluator.hpp>
#include <Compadre_PointCloudSearch.hpp>

moab::ErrorCode ComputeGMLSInterpolant( RuntimeContext&,
                                        const std::vector< double >& src_xyz,
                                        const std::vector< double >& src_tdata,
                                        const std::vector< double >& dst_xyz,
                                        std::vector< double >& dst_tdata )
{
    int argc    = 0;
    char** argv = nullptr;
    using namespace Compadre;
    const int order                   = 2;
    const int dimension               = 3;
    const size_t number_source_coords = src_tdata.size();
    const size_t number_target_coords = dst_tdata.size();
    const std::string constraint_name = "conserve";
    const std::string solver_name     = "GMLS";
    const std::string problem_name    = "GMLS-MPAS-ROMS";
    const int min_neighbors           = 1 * Compadre::GMLS::getNP( order );
    const int max_neighbors           = 1 * Compadre::GMLS::getNP( order ) * 1.25;

    Kokkos::initialize( argc, argv );
    {

        // Kokkos::View< int**, Kokkos::DefaultExecutionSpace > neighbor_lists( "neighbor lists", number_target_coords,
        //                                                          max_neighbors + 1 );  // first column is # of neighbors
        Kokkos::View< double**, Kokkos::DefaultExecutionSpace > source_coords( "neighbor coordinates", number_source_coords,
                                                                   dimension );
        Kokkos::View< double*, Kokkos::DefaultExecutionSpace > sampling_data( "samples of true solution",
                                                                  source_coords.extent( 0 ) );

        Kokkos::View< double*, Kokkos::DefaultExecutionSpace > epsilon( "h supports", number_target_coords );
        Kokkos::View< double**, Kokkos::DefaultExecutionSpace > target_coords( "target coordinates", number_target_coords,
                                                                   dimension );

        for( size_t i = 0; i < number_source_coords; i++ )
        {
            const size_t offset   = i * 3;
            source_coords( i, 0 ) = src_xyz[offset + 0];
            source_coords( i, 1 ) = src_xyz[offset + 1];
            source_coords( i, 2 ) = src_xyz[offset + 2];
            sampling_data( i )    = src_tdata[i];
        }
        for( size_t i = 0; i < number_target_coords; i++ )
        {
            const size_t offset   = i * 3;
            target_coords( i, 0 ) = dst_xyz[offset + 0];
            target_coords( i, 1 ) = dst_xyz[offset + 1];
            target_coords( i, 2 ) = dst_xyz[offset + 2];
            epsilon( i ) = 1e-3;
        }

        GMLS my_GMLS( order, dimension, solver_name.c_str(), problem_name.c_str(), constraint_name.c_str(),
                      2 /*manifold order*/ );

        // number_of_neighbors_list must be the same size as the number of target sites so that it can be populated
        // with the number of neighbors for each target site.
        Kokkos::View< int* > number_of_neighbors_list( "number of neighbor lists",
                                                          number_target_coords );  // first column is # of neighbors
        // Point cloud construction for neighbor search
        // CreatePointCloudSearch constructs an object of type PointCloudSearch, but deduces the templates for you
        auto point_cloud_search( CreatePointCloudSearch( source_coords, dimension ) );

        // each row is a neighbor list for a target site, with the first column of each row containing
        // the number of neighbors for that rows corresponding target site
        double epsilon_multiplier = 1.5;
        int estimated_upper_bound_number_neighbors =
            point_cloud_search.getEstimatedNumberNeighborsUpperBound( min_neighbors, dimension, epsilon_multiplier );

        std::cout << "Estimated number of maximum neighbors = " << estimated_upper_bound_number_neighbors << std::endl;

        // neighbor_lists will contain all neighbor lists (for each target site) in a compressed row format
        // Initially, we do a dry-run to calculate neighborhood sizes before actually storing the result. This is
        // why we can start with a neighbor_lists size of 0.
        Kokkos::View< int* > neighbor_lists( "neighbor lists",
                                                0 );  // first column is # of neighbors

        size_t storage_size =
            point_cloud_search.generateCRNeighborListsFromKNNSearch( true /*dry run*/, target_coords, neighbor_lists,
                                                                     number_of_neighbors_list, epsilon, min_neighbors,
                                                                     epsilon_multiplier );

        // resize neighbor_lists so as to be large enough to contain all neighborhoods
        Kokkos::resize( neighbor_lists, storage_size );

        std::cout << "Storage size = " << storage_size << std::endl;

        int maxnn = number_of_neighbors_list( 0 ), minnn = number_of_neighbors_list( 0 );
        for( size_t i = 0; i < number_target_coords; ++i )
        {
            maxnn = ( number_of_neighbors_list( i ) > maxnn ? number_of_neighbors_list( i ) : maxnn );
            minnn = ( number_of_neighbors_list( i ) < minnn ? number_of_neighbors_list( i ) : minnn );
        }
        std::cout << "number_neighbors: max = " << maxnn << ", min = " << minnn << std::endl;

        // query the point cloud a second time, but this time storing results into neighbor_lists
        point_cloud_search.generateCRNeighborListsFromKNNSearch( false /*not dry run*/, target_coords, neighbor_lists,
                                                                 number_of_neighbors_list, epsilon, min_neighbors,
                                                                 epsilon_multiplier );

        std::cout << "finished computing neighbor lists " << std::endl;
        // my_GMLS.setProblemData( neighbor_lists, source_coords, target_coords, epsilon );
        my_GMLS.setProblemData( neighbor_lists, number_of_neighbors_list, source_coords, target_coords, epsilon );
        std::cout << "finished setting problem data " << std::endl;

        my_GMLS.setWeightingParameter( 10 );

        std::vector< TargetOperation > lro( 1 );
        lro[0] = ScalarPointEvaluation;

        // and then pass them to the GMLS class
        my_GMLS.addTargets( lro );

        // sets the weighting kernel function from WeightingFunctionType
        my_GMLS.setWeightingType( WeightingFunctionType::Power );

        // power to use in that weighting kernel function
        my_GMLS.setWeightingParameter( 2 );
        std::cout << "finished setting all parameters " << std::endl;

        // generate the alphas that to be combined with data for each target operation requested in lro
        my_GMLS.generateAlphas( 16,
                                false /* keep polynomial coefficients, only needed for a test later in this program */ );
        std::cout << "finished generating alphas " << std::endl;

        Evaluator gmls_evaluator( &my_GMLS );
        auto output_value = gmls_evaluator.applyAlphasToDataAllComponentsAllTargetSites< double*, Kokkos::DefaultExecutionSpace >(
            sampling_data, ScalarPointEvaluation );
        for( size_t i = 0; i < number_target_coords; i++ )
        {
            // dst_tdata[i] = gmls_evaluator.applyAlphasToDataSingleComponentSingleTargetSite(
            //     sampling_data, 0, ScalarPointEvaluation, i, 0, 0, 0, 0, 0 );
            dst_tdata[i] = output_value( i );
        }
    }

    Kokkos::finalize();

    return moab::MB_SUCCESS;
}

#endif  // __computeGMLS_hpp__