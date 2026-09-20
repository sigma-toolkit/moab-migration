/**
 * @example ParallelRemapTemplate.cpp
 * @brief Parallel spherical scalar remapping with iMOAB and TempestRemap.
 *
 * Direct requirements: an MPI-enabled MOAB build with TempestRemap enabled.
 * iMOAB is part of libMOAB, not a separate dependency. NetCDF is optional and
 * is used only by --map-file.
 *
 * This example runs source, target, and intersection applications on one MPI
 * communicator. It writes an analytic field to the intersection application's
 * covering mesh, which is the source representation consumed by
 * iMOAB_ApplyScalarProjectionWeights. Production component-to-coupler
 * applications should instead establish a communication graph and migrate
 * their component field to that coverage mesh.
 */

#include "moab/MOABConfig.h"

#ifndef MOAB_HAVE_MPI
#error "ParallelRemapTemplate requires MOAB to be configured with MPI"
#endif

#ifndef MOAB_HAVE_TEMPESTREMAP
#error "ParallelRemapTemplate requires MOAB to be configured with TempestRemap"
#endif

#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "moab/ProgOptions.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

std::string to_lower( std::string value )
{
    std::transform( value.begin(), value.end(), value.begin(),
                    []( unsigned char c ) { return static_cast< char >( std::tolower( c ) ); } );
    return value;
}

[[noreturn]] void throw_iMOAB_error( ErrCode error, const std::string& operation, MPI_Comm comm )
{
    int rank = -1;
    MPI_Comm_rank( comm, &rank );
    std::ostringstream message;
    message << "Rank " << rank << ": " << operation << " failed with iMOAB error code " << error;
    throw std::runtime_error( message.str() );
}

void check_iMOAB( ErrCode error, const std::string& operation, MPI_Comm comm )
{
    if( error != moab::MB_SUCCESS ) throw_iMOAB_error( error, operation, comm );
}

void check_mpi( int error, const std::string& operation )
{
    if( error == MPI_SUCCESS ) return;
    char error_string[MPI_MAX_ERROR_STRING] = { 0 };
    int error_length                        = 0;
    MPI_Error_string( error, error_string, &error_length );
    throw std::runtime_error( operation + " failed: " + std::string( error_string, error_length ) );
}

int checked_storage_length( size_t size, const std::string& description )
{
    if( size > static_cast< size_t >( std::numeric_limits< int >::max() ) )
        throw std::runtime_error( description + " exceeds the iMOAB integer storage limit" );
    return static_cast< int >( size );
}

class MpiSession
{
  public:
    MpiSession( int* argc, char*** argv ) : owns_mpi_( false )
    {
        int initialized = 0;
        check_mpi( MPI_Initialized( &initialized ), "MPI_Initialized" );
        if( !initialized )
        {
            check_mpi( MPI_Init( argc, argv ), "MPI_Init" );
            owns_mpi_ = true;
        }
    }

    ~MpiSession() noexcept
    {
        if( !owns_mpi_ ) return;
        int finalized = 0;
        if( MPI_Finalized( &finalized ) == MPI_SUCCESS && !finalized ) MPI_Finalize();
    }

    MpiSession( const MpiSession& )            = delete;
    MpiSession& operator=( const MpiSession& ) = delete;

    MPI_Comm communicator() const { return MPI_COMM_WORLD; }

    int rank() const
    {
        int result = -1;
        check_mpi( MPI_Comm_rank( communicator(), &result ), "MPI_Comm_rank" );
        return result;
    }

    int size() const
    {
        int result = 0;
        check_mpi( MPI_Comm_size( communicator(), &result ), "MPI_Comm_size" );
        return result;
    }

  private:
    bool owns_mpi_;
};

class iMOABRuntime
{
  public:
    iMOABRuntime( int argc, char** argv, MPI_Comm comm ) : initialized_( false )
    {
        check_iMOAB( iMOAB_Initialize( argc, argv ), "iMOAB_Initialize", comm );
        initialized_ = true;
    }

    ~iMOABRuntime() noexcept
    {
        if( initialized_ ) iMOAB_Finalize();
    }

    iMOABRuntime( const iMOABRuntime& )            = delete;
    iMOABRuntime& operator=( const iMOABRuntime& ) = delete;

  private:
    bool initialized_;
};

enum class DiscretizationMethod : std::uint8_t  // NOLINT(performance-enum-size)
{
    FiniteVolume,
    ContinuousGll,
    DiscontinuousGll
};

class Discretization
{
  public:
    Discretization( const std::string& name, int order ) : method_( parse( name ) ), order_( order )
    {
        if( order_ < 1 ) throw std::runtime_error( "Discretization order must be positive" );
        if( method_ != DiscretizationMethod::FiniteVolume && order_ > std::numeric_limits< int >::max() / order_ )
            throw std::runtime_error( "Discretization order is too large" );
    }

    const char* name() const
    {
        switch( method_ )
        {
            case DiscretizationMethod::FiniteVolume: return "fv";
            case DiscretizationMethod::ContinuousGll: return "cgll";
            case DiscretizationMethod::DiscontinuousGll: return "dgll";
        }
        return "";
    }

    const char* dof_tag() const { return method_ == DiscretizationMethod::FiniteVolume ? "GLOBAL_ID" : "GLOBAL_DOFS"; }

    int order() const { return order_; }

    int components_per_entity() const { return method_ == DiscretizationMethod::FiniteVolume ? 1 : order_ * order_; }

  private:
    static DiscretizationMethod parse( const std::string& value )
    {
        const std::string normalized = to_lower( value );
        if( normalized == "fv" ) return DiscretizationMethod::FiniteVolume;
        if( normalized == "cgll" ) return DiscretizationMethod::ContinuousGll;
        if( normalized == "dgll" ) return DiscretizationMethod::DiscontinuousGll;
        throw std::runtime_error( "Unsupported discretization '" + value + "'. Use fv, cgll, or dgll." );
    }

    DiscretizationMethod method_;
    int order_;
};

struct Field
{
    Field( std::string field_name, const Discretization& field_discretization )
        : name( std::move( field_name ) ), components( field_discretization.components_per_entity() )
    {
    }

    std::string name;
    int components;
};

class RemapConfiguration
{
  public:
    RemapConfiguration( std::string source_mesh, std::string target_mesh, std::optional< std::string > output_mesh,
                         std::optional< std::string > map_file,
                         std::string read_options, std::string weights_id, Discretization source_discretization_value,
                         Discretization target_discretization_value, std::string source_field, std::string target_field,
                         int ghost_layers, int filter_type, bool write_local_intersection )
        : source_mesh( std::move( source_mesh ) ), target_mesh( std::move( target_mesh ) ), output_mesh( std::move( output_mesh ) ),
          map_file( std::move( map_file ) ), read_options( std::move( read_options ) ), weights_id( std::move( weights_id ) ),
          source_discretization( source_discretization_value ), target_discretization( target_discretization_value ),
          source_field( std::move( source_field ), source_discretization ), target_field( std::move( target_field ), target_discretization ),
          ghost_layers( ghost_layers ), filter_type( filter_type ), write_local_intersection( write_local_intersection )
    {
        if( this->source_mesh.empty() || this->target_mesh.empty() )
            throw std::runtime_error( "Both --source and --target are required" );
        if( ghost_layers < 0 ) throw std::runtime_error( "ghost-layers must be non-negative" );
        if( filter_type < 0 || filter_type > 3 ) throw std::runtime_error( "filter-type must be in [0, 3]" );
    }

    std::string source_mesh;
    std::string target_mesh;
    std::optional< std::string > output_mesh;
    std::optional< std::string > map_file;
    std::string read_options;
    std::string weights_id;
    Discretization source_discretization;
    Discretization target_discretization;
    Field source_field;
    Field target_field;
    int ghost_layers;
    int filter_type;
    bool write_local_intersection;
};

class iMOABApplication
{
  public:
    iMOABApplication( const char* name, int component_id, MPI_Comm comm )
        : comm_( comm ), component_id_( component_id ), application_id_( -1 ), registered_( false )
    {
        check_iMOAB( iMOAB_RegisterApplication( name, &comm_, &component_id_, id() ),
                     std::string( "iMOAB_RegisterApplication(" ) + name + ")", comm_ );
        registered_ = true;
    }

    virtual ~iMOABApplication() noexcept
    {
        if( registered_ ) iMOAB_DeregisterApplication( id() );
    }

    iMOABApplication( const iMOABApplication& )            = delete;
    iMOABApplication& operator=( const iMOABApplication& ) = delete;

    void load( const std::string& mesh_file, const std::string& options, int ghost_layers )
    {
        check_iMOAB( iMOAB_LoadMesh( id(), mesh_file.c_str(), options.c_str(), &ghost_layers ),
                     "iMOAB_LoadMesh(" + mesh_file + ")", comm_ );
    }

    void define( const Field& field )
    {
        int tag_type  = DENSE_DOUBLE;
        int tag_index = 0;
        int components = field.components;
        check_iMOAB( iMOAB_DefineTagStorage( id(), field.name.c_str(), &tag_type, &components, &tag_index ),
                     "iMOAB_DefineTagStorage(" + field.name + ")", comm_ );
    }

    int visible_element_count() const
    {
        int vertices[3] = { 0, 0, 0 };
        int elements[3] = { 0, 0, 0 };
        check_iMOAB( iMOAB_GetMeshInfo( id(), vertices, elements, nullptr, nullptr, nullptr ), "iMOAB_GetMeshInfo", comm_ );
        return elements[2];
    }

    std::vector< double > field_values( const Field& field, int entities ) const
    {
        std::vector< double > values( static_cast< size_t >( entities ) * field.components, 0.0 );
        if( values.empty() ) return values;
        int storage_length = checked_storage_length( values.size(), field.name );
        int entity_type    = 1;
        check_iMOAB( iMOAB_GetDoubleTagStorage( id(), field.name.c_str(), &storage_length, &entity_type, values.data() ),
                     "iMOAB_GetDoubleTagStorage(" + field.name + ")", comm_ );
        return values;
    }

    void write( const std::string& mesh_file, const MpiSession& mpi ) const
    {
        const char* options = mpi.size() > 1 ? "PARALLEL=WRITE_PART" : "";
        check_iMOAB( iMOAB_WriteMesh( id(), mesh_file.c_str(), options ), "iMOAB_WriteMesh(" + mesh_file + ")", comm_ );
    }

    // The C API needs this opaque application handle to compose applications.
    iMOAB_AppID native_id() const { return id(); }

  protected:
    iMOAB_AppID id() const { return const_cast< iMOAB_AppID >( &application_id_ ); }
    MPI_Comm communicator() const { return comm_; }

  private:
    MPI_Comm comm_;
    int component_id_;
    int application_id_;
    bool registered_;
};

class IntersectionApplication : public iMOABApplication
{
  public:
    IntersectionApplication( const char* name, int component_id, MPI_Comm comm ) : iMOABApplication( name, component_id, comm ) {}

    void compute( const iMOABApplication& source, const iMOABApplication& target )
    {
        check_iMOAB( iMOAB_ComputeMeshIntersectionOnSphere( source.native_id(), target.native_id(), id() ),
                     "iMOAB_ComputeMeshIntersectionOnSphere", communicator() );
    }

    void seed_coverage_field( const Field& field )
    {
        // ApplyScalarProjectionWeights reads the source field from this app's coverage mesh.
        define( field );
        int coverage_entities = 0;
        check_iMOAB( iMOAB_GetCoverageMeshInfo( id(), &coverage_entities, nullptr, nullptr ), "iMOAB_GetCoverageMeshInfo",
                     communicator() );
        require_nonempty( coverage_entities, "Source coverage mesh" );
        if( coverage_entities == 0 ) return;

        std::vector< int > global_ids( static_cast< size_t >( coverage_entities ) );
        check_iMOAB( iMOAB_GetCoverageMeshInfo( id(), &coverage_entities, global_ids.data(), nullptr ),
                     "iMOAB_GetCoverageMeshInfo(ids)", communicator() );
        std::vector< double > values = analytic_values( global_ids, field.components );
        int storage_length           = checked_storage_length( values.size(), field.name );
        check_iMOAB( iMOAB_SetDoubleTagStorageOnCoverage( id(), field.name.c_str(), &storage_length, values.data() ),
                     "iMOAB_SetDoubleTagStorageOnCoverage(" + field.name + ")", communicator() );
    }

    void compute_weights( const RemapConfiguration& configuration )
    {
        int no_bubble        = 1;
        int monotone         = 0;
        int volumetric       = 0;
        int inverse_distance = 0;
        int no_conservation  = 0;
        int validate         = 0;
        int source_order     = configuration.source_discretization.order();
        int target_order     = configuration.target_discretization.order();
        check_iMOAB( iMOAB_ComputeScalarProjectionWeights(
                         id(), configuration.weights_id.c_str(), configuration.source_discretization.name(), &source_order,
                         configuration.target_discretization.name(), &target_order, nullptr, &no_bubble, &monotone, &volumetric,
                         &inverse_distance, &no_conservation, &validate, configuration.source_discretization.dof_tag(),
                         configuration.target_discretization.dof_tag() ),
                     "iMOAB_ComputeScalarProjectionWeights", communicator() );
    }

    void apply_weights( const RemapConfiguration& configuration )
    {
        int filter_type = configuration.filter_type;
        check_iMOAB( iMOAB_ApplyScalarProjectionWeights( id(), &filter_type, configuration.weights_id.c_str(),
                                                         configuration.source_field.name.c_str(), configuration.target_field.name.c_str() ),
                     "iMOAB_ApplyScalarProjectionWeights", communicator() );
    }

    void write_map( const std::string& weights_id, const std::string& map_file )
    {
        check_iMOAB( iMOAB_WriteMapFile( id(), weights_id.c_str(), map_file.c_str() ), "iMOAB_WriteMapFile(" + map_file + ")",
                     communicator() );
    }

    void write_local_mesh( char* prefix )
    {
        check_iMOAB( iMOAB_WriteLocalMesh( id(), prefix ), "iMOAB_WriteLocalMesh", communicator() );
    }

  private:
    void require_nonempty( int local_count, const std::string& mesh_name ) const
    {
        int global_count = 0;
        check_mpi( MPI_Allreduce( &local_count, &global_count, 1, MPI_INT, MPI_SUM, communicator() ), "MPI_Allreduce" );
        if( global_count == 0 ) throw std::runtime_error( mesh_name + " contains no visible 2D entities" );
    }

    static std::vector< double > analytic_values( const std::vector< int >& global_ids, int components )
    {
        std::vector< double > values( global_ids.size() * static_cast< size_t >( components ), 0.0 );
        for( size_t entity = 0; entity < global_ids.size(); ++entity )
            for( int component = 0; component < components; ++component )
                values[entity * static_cast< size_t >( components ) + component] =
                    static_cast< double >( global_ids[entity] ) + 0.01 * component;
        return values;
    }
};

class RemapWorkflow
{
  public:
    RemapWorkflow( const RemapConfiguration& configuration, const MpiSession& mpi, int argc, char** argv )
        : configuration_( configuration ), mpi_( mpi ), runtime_( argc, argv, mpi.communicator() ),
          source_( "PAR_REMAP_SOURCE", 101, mpi.communicator() ), target_( "PAR_REMAP_TARGET", 202, mpi.communicator() ),
          intersection_( "PAR_REMAP_INTERSECTION", 303, mpi.communicator() )
    {
    }

    void run()
    {
        load_meshes();
        define_fields();
        intersect_meshes();
        populate_coverage_field();
        generate_weights();
        project_field();
        report_target_field();
        write_outputs();
    }

  private:
    void load_meshes()
    {
        source_.load( configuration_.source_mesh, configuration_.read_options, configuration_.ghost_layers );
        target_.load( configuration_.target_mesh, configuration_.read_options, configuration_.ghost_layers );
        source_entities_ = source_.visible_element_count();
        target_entities_ = target_.visible_element_count();
        require_globally_nonempty( source_entities_, "Source mesh" );
        require_globally_nonempty( target_entities_, "Target mesh" );
    }

    void define_fields()
    {
        source_.define( configuration_.source_field );
        target_.define( configuration_.target_field );
    }

    void intersect_meshes()
    {
        if( mpi_.rank() == 0 )
            std::cout << "Computing mesh intersection for " << configuration_.source_mesh << " -> " << configuration_.target_mesh << '\n';
        intersection_.compute( source_, target_ );
    }

    void populate_coverage_field() { intersection_.seed_coverage_field( configuration_.source_field ); }

    void generate_weights() { intersection_.compute_weights( configuration_ ); }

    void project_field() { intersection_.apply_weights( configuration_ ); }

    void report_target_field() const
    {
        const std::vector< double > values = target_.field_values( configuration_.target_field, target_entities_ );
        const double local_sum = std::accumulate( values.begin(), values.end(), 0.0 );
        const double local_min = values.empty() ? std::numeric_limits< double >::infinity()
                                                : *std::min_element( values.begin(), values.end() );
        const double local_max = values.empty() ? -std::numeric_limits< double >::infinity()
                                                : *std::max_element( values.begin(), values.end() );
        const long long local_count = static_cast< long long >( values.size() );

        double global_sum = 0.0;
        double global_min = 0.0;
        double global_max = 0.0;
        long long global_count = 0;
        check_mpi( MPI_Reduce( &local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_.communicator() ), "MPI_Reduce(sum)" );
        check_mpi( MPI_Reduce( &local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, mpi_.communicator() ), "MPI_Reduce(min)" );
        check_mpi( MPI_Reduce( &local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, mpi_.communicator() ), "MPI_Reduce(max)" );
        check_mpi( MPI_Reduce( &local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, 0, mpi_.communicator() ),
                   "MPI_Reduce(count)" );
        if( mpi_.rank() == 0 )
        {
            const double mean = global_count == 0 ? 0.0 : global_sum / static_cast< double >( global_count );
            std::cout << configuration_.target_field.name << ": count=" << global_count << " min=" << global_min
                      << " max=" << global_max << " mean=" << mean << '\n';
        }
    }

    void write_outputs()
    {
#ifdef MOAB_HAVE_NETCDF
        if( configuration_.map_file ) intersection_.write_map( configuration_.weights_id, *configuration_.map_file );
#else
        if( mpi_.rank() == 0 && configuration_.map_file )
            std::cout << "Skipping map-file write: this MOAB build has no NetCDF support\n";
#endif
        if( configuration_.write_local_intersection )
        {
            char prefix[] = "parallel_remap_intx";
            intersection_.write_local_mesh( prefix );
        }
        if( configuration_.output_mesh ) target_.write( *configuration_.output_mesh, mpi_ );
    }

    void require_globally_nonempty( int local_count, const std::string& mesh_name ) const
    {
        int global_count = 0;
        check_mpi( MPI_Allreduce( &local_count, &global_count, 1, MPI_INT, MPI_SUM, mpi_.communicator() ), "MPI_Allreduce" );
        if( global_count == 0 ) throw std::runtime_error( mesh_name + " contains no visible 2D entities" );
    }

    const RemapConfiguration& configuration_;
    const MpiSession& mpi_;
    iMOABRuntime runtime_;
    iMOABApplication source_;
    iMOABApplication target_;
    IntersectionApplication intersection_;
    int source_entities_ = 0;
    int target_entities_ = 0;
};

RemapConfiguration parse_command_line( int argc, char* argv[] )
{
    std::string source_mesh;
    std::string target_mesh;
    std::string output_mesh = "remap_target_out.h5m";
    std::string map_file;
    std::string read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
    std::string weights_id   = "scalar";
    std::string source_disc  = "fv";
    std::string target_disc  = "fv";
    std::string source_field = "SRC_FIELD";
    std::string target_field = "TGT_FIELD";
    int source_order         = 1;
    int target_order         = 1;
    int ghost_layers         = 0;
    int filter_type          = 0;
    bool write_local_intx    = false;

    ProgOptions options;
    options.addOpt< std::string >( "source,s", "Source spherical mesh file", &source_mesh );
    options.addOpt< std::string >( "target,t", "Target spherical mesh file", &target_mesh );
    options.addOpt< std::string >( "output,o", "Output target mesh with projected field", &output_mesh );
    options.addOpt< std::string >( "map-file", "Optional NetCDF map file", &map_file );
    options.addOpt< std::string >( "read-options", "Parallel iMOAB_LoadMesh options", &read_options );
    options.addOpt< std::string >( "weights-id", "Computed weight-map identifier", &weights_id );
    options.addOpt< std::string >( "source-disc", "Source discretization: fv, cgll, or dgll", &source_disc );
    options.addOpt< std::string >( "target-disc", "Target discretization: fv, cgll, or dgll", &target_disc );
    options.addOpt< int >( "source-order", "Source discretization order", &source_order );
    options.addOpt< int >( "target-order", "Target discretization order", &target_order );
    options.addOpt< std::string >( "source-field", "Source scalar tag", &source_field );
    options.addOpt< std::string >( "target-field", "Target scalar tag", &target_field );
    options.addOpt< int >( "ghost-layers", "Requested mesh ghost layers", &ghost_layers );
    options.addOpt< int >( "filter-type", "0 none, 1 global, 2 local, 3 patch", &filter_type );
    options.addOpt< void >( "write-local-intx", "Write local intersection meshes", &write_local_intx );
    options.parseCommandLine( argc, argv );

    const std::optional< std::string > requested_output = output_mesh.empty() ? std::nullopt : std::optional< std::string >{ output_mesh };
    const std::optional< std::string > requested_map    = map_file.empty() ? std::nullopt : std::optional< std::string >{ map_file };
    return RemapConfiguration( source_mesh, target_mesh, requested_output, requested_map, read_options, weights_id,
                               Discretization( source_disc, source_order ), Discretization( target_disc, target_order ),
                               source_field, target_field, ghost_layers, filter_type, write_local_intx );
}

}  // namespace

int main( int argc, char* argv[] )
try
{
    MpiSession mpi( &argc, &argv );
    try
    {
        const RemapConfiguration configuration = parse_command_line( argc, argv );
        RemapWorkflow( configuration, mpi, argc, argv ).run();
    }
    catch( const std::exception& error )
    {
        std::cerr << "Rank " << mpi.rank() << ": " << error.what() << '\n';
        MPI_Abort( mpi.communicator(), EXIT_FAILURE );
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
catch( const std::exception& error )
{
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
}
