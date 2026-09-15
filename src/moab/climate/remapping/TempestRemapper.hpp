/*
 * =====================================================================================
 *
 *       Filename:  TempestRemapper.hpp
 *
 *    Description:  Interface to the TempestRemap library to enable intersection and
 *                  high-order conservative remapping of climate solution from
 *                  arbitrary resolution of source and target grids on the sphere.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#ifndef MB_TEMPESTREMAPPER_HPP
#define MB_TEMPESTREMAPPER_HPP

#include "moab/climate/remapping/Remapper.hpp"
#include "moab/climate/intx_mesh/Intx2MeshOnSphere.hpp"
#include "moab/climate/intx_mesh/IntxUtils.hpp"

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

namespace moab
{

// Forward declare our friend, the mapper
class TempestOnlineMap;

class TempestRemapper : public Remapper
{
  public:
#ifdef MOAB_HAVE_MPI
    TempestRemapper( moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL, bool offlineMode = false )
        : Remapper( mbInt, pcomm ),
#else
    TempestRemapper( moab::Interface* mbInt, bool offlineMode = false )
        : Remapper( mbInt ),
#endif
          offlineWorkflow( offlineMode ), meshValidate( false ), constructEdgeMap( false ), m_source_type( DEFAULT ),
          m_target_type( DEFAULT )
    {
    }

    virtual ~TempestRemapper();

    // Mesh type with a correspondence to Tempest/Climate formats
    enum TempestMeshType
    {
        DEFAULT        = -1,
        CS             = 0,
        RLL            = 1,
        ICO            = 2,
        ICOD           = 3,
        OVERLAP_FILES  = 4,
        OVERLAP_MEMORY = 5,
        OVERLAP_MOAB   = 6
    };

    friend class TempestOnlineMap;

  public:
    /**
     * @brief Initialize the TempestRemapper object internal data structures including the mesh sets
     *        and TempestRemap mesh references.
     *
     * @param initialize_fsets Flag to initialize the mesh sets (default: true)
     * @return ErrorCode indicating the status of the initialization
     */
    virtual ErrorCode initialize( bool initialize_fsets = true );

    /**
     * @brief Deallocate and clear any memory initialized in the TempestRemapper object
     *
     * @return ErrorCode indicating the status of the clear operation
     */
    virtual ErrorCode clear();

    /**
     * @brief Generate a mesh in memory of given type (CS/RLL/ICO/MPAS(structured)) and store it
     *        under the context specified by the user.
     *
     * @param ctx Intersection context
     * @param type Type of mesh to generate
     * @return ErrorCode indicating the status of the mesh generation
     */
    moab::ErrorCode GenerateMesh( Remapper::IntersectionContext ctx, TempestMeshType type );

    /**
     * @brief Load a mesh from disk of given type and store it under the context specified by the user.
     *
     * @param ctx Intersection context
     * @param inputFilename File name of the mesh to load
     * @param type Type of mesh to load
     * @return ErrorCode indicating the status of the mesh loading
     */
    moab::ErrorCode LoadMesh( Remapper::IntersectionContext ctx, std::string inputFilename, TempestMeshType type );

    /**
     * @brief Construct a source covering mesh such that it completely encompasses the target grid in
     *        parallel. This operation is critical to ensure that the parallel advancing-front
     *        intersection algorithm can compute the intersection mesh only locally without any process
     *        communication.
     *
     * @param tolerance Tolerance for the covering mesh construction (default: 1e-8)
     * @param radius_src Radius of the source mesh (default: 1.0)
     * @param radius_tgt Radius of the target mesh (default: 1.0)
     * @param boxeps Box epsilon value (default: 0.1)
     * @param regional_mesh Flag to indicate if the mesh is regional (default: false)
     * @param gnomonic Flag to indicate if the mesh is gnomonic (default: true)
     * @param nb_ghost_layers Number of ghost layers (default: 0)
     * @return ErrorCode indicating the status of the covering mesh construction
     */
    moab::ErrorCode ConstructCoveringSet( double tolerance    = 1e-8,
                                          double radius_src   = 1.0,
                                          double radius_tgt   = 1.0,
                                          double boxeps       = 0.1,
                                          bool regional_mesh  = false,
                                          bool gnomonic       = true,
                                          int nb_ghost_layers = 0 );

    /**
     * @brief Compute the intersection mesh between the source and target grids that have been
     *        instantiated in the Remapper. This function invokes the parallel advancing-front
     *        intersection algorithm internally for spherical meshes and can handle arbitrary
     *        unstructured grids (CS, RLL, ICO, MPAS) with and without holes.
     *
     * @param kdtree_search Flag to enable k-d tree search (default: true)
     * @param use_tempest Flag to use TempestRemap (default: false)
     * @return ErrorCode indicating the status of the intersection mesh computation
     */
    moab::ErrorCode ComputeOverlapMesh( bool kdtree_search = true, bool use_tempest = false );

    /**
     * @brief Convert the TempestRemap mesh object to a corresponding MOAB mesh representation
     *        according to the intersection context.
     *
     * @param ctx Intersection context
     * @return ErrorCode indicating the status of the mesh conversion
     */
    moab::ErrorCode ConvertTempestMesh( Remapper::IntersectionContext ctx );

    /**
     * @brief Convert the MOAB mesh representation to a corresponding TempestRemap mesh object
     *        according to the intersection context.
     *
     * @param ctx Intersection context
     * @return ErrorCode indicating the status of the mesh conversion
     */
    moab::ErrorCode ConvertMeshToTempest( Remapper::IntersectionContext ctx );

    /**
     * @brief Get the TempestRemap mesh object according to the intersection context.
     *
     * @param ctx Intersection context
     * @return Pointer to the TempestRemap mesh object
     */
    Mesh* GetMesh( Remapper::IntersectionContext ctx );

    /**
     * @brief Set the TempestRemap mesh object according to the intersection context.
     *
     * @param ctx Intersection context
     * @param mesh Pointer to the TempestRemap mesh object
     * @param overwrite Flag to overwrite the existing mesh (default: true)
     */
    void SetMesh( Remapper::IntersectionContext ctx, Mesh* mesh, bool overwrite = true );

    /**
     * @brief Set the mesh set according to the intersection context.
     *
     * @param ctx Intersection context
     * @param mset MOAB mesh set handle
     * @param entities MOAB range of entities (optional)
     */
    void SetMeshSet( Remapper::IntersectionContext ctx, moab::EntityHandle mset, moab::Range* entities = nullptr );

    /**
     * @brief Get the covering mesh (TempestRemap) object.
     *
     * @return Pointer to the covering mesh object
     */
    Mesh* GetCoveringMesh();

    /**
     * @brief Get the MOAB mesh set corresponding to the intersection context.
     *
     * @param ctx Intersection context
     * @return MOAB mesh set handle
     */
    moab::EntityHandle& GetMeshSet( Remapper::IntersectionContext ctx );

    /**
     * @brief Const overload. Get the MOAB mesh set corresponding to the intersection context.
     *
     * @param ctx Intersection context
     * @return MOAB mesh set handle
     */
    moab::EntityHandle GetMeshSet( Remapper::IntersectionContext ctx ) const;

    /**
     * @brief Get the mesh element entities corresponding to the intersection context.
     *
     * @param ctx Intersection context
     * @return MOAB range of mesh element entities
     */
    moab::Range& GetMeshEntities( Remapper::IntersectionContext ctx );

    /**
     * @brief Const overload. Get the mesh element entities corresponding to the intersection context.
     *
     * @param ctx Intersection context
     * @return MOAB range of mesh element entities
     */
    const moab::Range& GetMeshEntities( Remapper::IntersectionContext ctx ) const;

    /**
     * @brief Get the mesh vertices corresponding to the intersection context. Useful for point-cloud
     *        meshes.
     *
     * @param ctx Intersection context
     * @return MOAB range of mesh vertices
     */
    moab::Range& GetMeshVertices( Remapper::IntersectionContext ctx );

    /**
     * @brief Const overload. Get the mesh vertices corresponding to the intersection context. Useful
     *        for point-cloud meshes.
     *
     * @param ctx Intersection context
     * @return MOAB range of mesh vertices
     */
    const moab::Range& GetMeshVertices( Remapper::IntersectionContext ctx ) const;

    /**
     * @brief Get access to the underlying source covering set if available. Else return the source
     *        set.
     *
     * @return MOAB entity handle of the covering set
     */
    moab::EntityHandle& GetCoveringSet();

    /**
     * @brief Set the mesh type corresponding to the intersection context
     *
     * @param ctx Intersection context
     * @param metadata Vector of mesh type metadata
     */
    void SetMeshType( Remapper::IntersectionContext ctx, const std::vector< int >& metadata );

    /**
     * @brief Reconstruct mesh, used now only for IO; need a better solution maybe
     *
     * @param ctx Intersection context
     * @param meshSet MOAB mesh set handle
     */
    void ResetMeshSet( Remapper::IntersectionContext ctx, moab::EntityHandle meshSet );

    /**
     * @brief Get the mesh type corresponding to the intersection context
     *
     * @param ctx Intersection context
     * @return Mesh type
     */
    TempestMeshType GetMeshType( Remapper::IntersectionContext ctx ) const;

    /**
     * @brief Gather the overlap mesh and associated source/target data and write it out to disk
     *        using the TempestRemap output interface. This information can then be used with the
     *        "GenerateOfflineMap" tool in TempestRemap as needed.
     *
     * @param strOutputFileName Output file name
     * @param fAllParallel Flag to write all parallel data (default: false)
     * @param fInputConcave Flag to indicate if the input mesh is concave (default: false)
     * @param fOutputConcave Flag to indicate if the output mesh is concave (default: false)
     * @return ErrorCode indicating the status of the write operation
     */
    moab::ErrorCode WriteTempestIntersectionMesh( std::string strOutputFileName,
                                                  const bool fAllParallel,
                                                  const bool fInputConcave,
                                                  const bool fOutputConcave );

    /**
     * @brief Generate the necessary metadata and specifically the GLL node numbering for DoFs for a
     *        CS mesh. This negates the need for running external code like HOMME to output the
     *        numbering needed for computing maps. The functionality is used through the `mbconvert`
     *        tool to compute processor-invariant Global DoF IDs at GLL nodes.
     *
     * @param ntot_elements Total number of elements
     * @param entities MOAB range of entities
     * @param secondary_entities MOAB range of secondary entities (optional)
     * @param dofTagName Name of the DoF tag
     * @param nP Number of points
     * @return ErrorCode indicating the status of the metadata generation
     */
    moab::ErrorCode GenerateCSMeshMetadata( const int ntot_elements,
                                            moab::Range& entities,
                                            moab::Range* secondary_entities,
                                            const std::string& dofTagName,
                                            int nP );

    /**
     * @brief Generate the necessary metadata for DoF node numbering in a given mesh.
     *        Currently, only the functionality to generate numbering on CS grids is supported.
     *
     * @param mesh Mesh object
     * @param ntot_elements Total number of elements
     * @param entities MOAB range of entities
     * @param secondary_entities MOAB range of secondary entities (optional)
     * @param dofTagName Name of the DoF tag
     * @param nP Number of points
     * @return ErrorCode indicating the status of the metadata generation
     */
    moab::ErrorCode GenerateMeshMetadata( Mesh& mesh,
                                          const int ntot_elements,
                                          moab::Range& entities,
                                          moab::Range* secondary_entities,
                                          const std::string dofTagName,
                                          int nP );

    /**
     * @brief Get all the ghosted overlap entities that were accumulated to enable conservation in
     *        parallel
     *
     * @param sharedGhostEntities MOAB range of ghosted overlap entities
     * @return ErrorCode indicating the status of the get operation
     */
    moab::ErrorCode GetOverlapAugmentedEntities( moab::Range& sharedGhostEntities );

#ifndef MOAB_HAVE_MPI
    /**
     * @brief Internal method to assign vertex and element global IDs if one does not exist already
     *
     * @param idtag Tag for the global IDs
     * @param this_set MOAB entity handle of the mesh set
     * @param dimension Dimension of the mesh (default: 2)
     * @param start_id Starting ID (default: 1)
     * @return ErrorCode indicating the status of the assignment
     */
    moab::ErrorCode assign_vertex_element_IDs( Tag idtag,
                                               EntityHandle this_set,
                                               const int dimension = 2,
                                               const int start_id  = 1 );
#endif

    /**
     * @brief Get the masks that could have been defined
     *
     * @param ctx Intersection context
     * @param masks Vector of masks
     * @return ErrorCode indicating the status of the get operation
     */
    ErrorCode GetIMasks( Remapper::IntersectionContext ctx, std::vector< int >& masks );

  public:  // public members
    /**
     * @brief Flag indicating whether the workflow is in offline mode.
     *
     * This flag is used to determine the context of the workflow, specifically
     * whether it is running in an offline mode (mbtempest).
     */
    const bool offlineWorkflow;

    /**
     * @brief Flag to enable mesh validation after loading from file.
     *
     * If set to true, the mesh will be validated after it is loaded from a file.
     */
    bool meshValidate;

    /**
     * @brief Flag to construct the edge map within the TempestRemap data structures.
     *
     * If set to true, the edge map will be constructed within the TempestRemap
     * data structures.
     */
    bool constructEdgeMap;

    /**
     * @brief Global verbosity flag.
     *
     * This flag controls the verbosity of the output. If set to true, more
     * detailed output will be generated.
     */
    static const bool verbose = true;

  private:
    /**
     * @brief Convert all MOAB meshes to TempestRemap format.
     *
     * Utility method primarily used in online workflows to convert all meshes
     * to TempestRemap format.
     *
     * @return ErrorCode Status of the conversion
     */
    ErrorCode ConvertAllMeshesToTempest();

    /**
     * @brief Convert overlap mesh to source-ordered format.
     *
     * Transforms the overlap mesh into a source-ordered format compatible
     * with TempestRemap.
     *
     * @return moab::ErrorCode Status of the conversion
     */
    moab::ErrorCode ConvertOverlapMeshSourceOrdered();

    // Private methods

    /**
     * @brief Load a mesh from disk in TempestRemap format.
     *
     * @param inputFilename Path to the input mesh file
     * @param tempest_mesh Pointer to store the loaded TempestRemap mesh
     * @return moab::ErrorCode Status of the load operation
     */
    moab::ErrorCode load_tempest_mesh_private( std::string inputFilename, Mesh** tempest_mesh );

    /**
     * @brief Convert a MOAB mesh to TempestRemap format.
     *
     * @param mesh Target TempestRemap mesh pointer
     * @param meshset MOAB mesh set to convert
     * @param entities Range of entities to convert
     * @param pverts Optional pointer to vertex range
     * @return moab::ErrorCode Status of the conversion
     */
    moab::ErrorCode convert_mesh_to_tempest_private( Mesh* mesh,
                                                     moab::EntityHandle meshset,
                                                     moab::Range& entities,
                                                     moab::Range* pverts );

    /**
     * @brief Convert a TempestRemap mesh to MOAB format.
     *
     * @param type Type of the TempestRemap mesh
     * @param mesh Source TempestRemap mesh
     * @param meshset Target MOAB mesh set
     * @param entities Range to store converted entities
     * @param vertices Optional range to store vertices
     * @return moab::ErrorCode Status of the conversion
     */
    moab::ErrorCode convert_tempest_mesh_private( TempestMeshType type,
                                                  Mesh* mesh,
                                                  moab::EntityHandle& meshset,
                                                  moab::Range& entities,
                                                  moab::Range* vertices );

    /**
     * @brief Augment overlap mesh with ghosted entities.
     *
     * Adds ghosted entities to the overlap mesh to ensure conservation
     * in parallel operations.
     *
     * @return moab::ErrorCode Status of the augmentation
     */
    moab::ErrorCode AugmentOverlapSet();

    /* Source meshset, mesh and entity references */
    Mesh* m_source;
    TempestMeshType m_source_type;
    moab::Range m_source_entities;
    moab::Range m_source_vertices;
    moab::EntityHandle m_source_set;
    int max_source_edges;
    bool point_cloud_source;
    std::vector< int > m_source_metadata;

    /* Target meshset, mesh and entity references */
    Mesh* m_target;
    TempestMeshType m_target_type;
    moab::Range m_target_entities;
    moab::Range m_target_vertices;
    moab::EntityHandle m_target_set;
    int max_target_edges;
    bool point_cloud_target;
    std::vector< int > m_target_metadata;

    /* Overlap meshset, mesh and entity references */
    Mesh* m_overlap;
    TempestMeshType m_overlap_type;
    moab::Range m_overlap_entities;
    moab::EntityHandle m_overlap_set;
    std::vector< std::pair< int, int > > m_sorted_overlap_order;

    /* Intersection context on a sphere */
    moab::Intx2MeshOnSphere* mbintx;

    /* Parallel - migrated mesh that is in the local view */
    Mesh* m_covering_source;
    moab::EntityHandle m_covering_source_set;
    moab::Range m_covering_source_entities;
    moab::Range m_covering_source_vertices;

    /* local to glboal and global to local ID maps */
    // std::map< int, int > gid_to_lid_src, gid_to_lid_covsrc, gid_to_lid_tgt;
    // std::map< int, int > lid_to_gid_src, lid_to_gid_covsrc, lid_to_gid_tgt;

    IntxAreaUtils::AreaMethod m_area_method;

    bool rrmgrids;
    bool is_parallel, is_root;
    int rank, size;
};

// Inline functions
inline Mesh* TempestRemapper::GetMesh( Remapper::IntersectionContext ctx )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source;
        case Remapper::TargetMesh:
            return m_target;
        case Remapper::OverlapMesh:
            return m_overlap;
        case Remapper::CoveringMesh:
            return m_covering_source;
        case Remapper::DEFAULT:
        default:
            return NULL;
    }
}

inline void TempestRemapper::SetMesh( Remapper::IntersectionContext ctx, Mesh* mesh, bool overwrite )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            if( !overwrite && m_source ) return;
            if( overwrite && m_source ) delete m_source;
            m_source = mesh;
            break;
        case Remapper::TargetMesh:
            if( !overwrite && m_target ) return;
            if( overwrite && m_target ) delete m_target;
            m_target = mesh;
            break;
        case Remapper::OverlapMesh:
            if( !overwrite && m_overlap ) return;
            if( overwrite && m_overlap ) delete m_overlap;
            m_overlap = mesh;
            break;
        case Remapper::CoveringMesh:
            if( !overwrite && m_covering_source ) return;
            if( overwrite && m_covering_source ) delete m_covering_source;
            m_covering_source = mesh;
            break;
        case Remapper::DEFAULT:
        default:
            break;
    }
}

inline void TempestRemapper::ResetMeshSet( Remapper::IntersectionContext ctx, moab::EntityHandle meshSet )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            delete m_source;
            m_source     = new Mesh;
            m_source_set = meshSet;
            convert_mesh_to_tempest_private( m_source, m_source_set, m_source_entities, &m_source_vertices );
            m_source->CalculateFaceAreas( false );  // fInputConcave is false ?
            break;
        case Remapper::TargetMesh:
            // not needed yet
            break;
        case Remapper::OverlapMesh:
            // not needed yet
            break;
        case Remapper::CoveringMesh:
            // not needed yet
            break;
        case Remapper::DEFAULT:
        default:
            break;
    }
}

inline moab::EntityHandle& TempestRemapper::GetMeshSet( Remapper::IntersectionContext ctx )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_set;
        case Remapper::TargetMesh:
            return m_target_set;
        case Remapper::OverlapMesh:
            return m_overlap_set;
        case Remapper::CoveringMesh:
            return m_covering_source_set;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_overlap_set );
    }
}

inline moab::EntityHandle TempestRemapper::GetMeshSet( Remapper::IntersectionContext ctx ) const
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_set;
        case Remapper::TargetMesh:
            return m_target_set;
        case Remapper::OverlapMesh:
            return m_overlap_set;
        case Remapper::CoveringMesh:
            return m_covering_source_set;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_overlap_set );
    }
}

inline moab::Range& TempestRemapper::GetMeshEntities( Remapper::IntersectionContext ctx )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_entities;
        case Remapper::TargetMesh:
            return m_target_entities;
        case Remapper::OverlapMesh:
            return m_overlap_entities;
        case Remapper::CoveringMesh:
            return m_covering_source_entities;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_overlap_entities );
    }
}

inline const moab::Range& TempestRemapper::GetMeshEntities( Remapper::IntersectionContext ctx ) const
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_entities;
        case Remapper::TargetMesh:
            return m_target_entities;
        case Remapper::OverlapMesh:
            return m_overlap_entities;
        case Remapper::CoveringMesh:
            return m_covering_source_entities;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_overlap_entities );
    }
}

inline moab::Range& TempestRemapper::GetMeshVertices( Remapper::IntersectionContext ctx )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_vertices;
        case Remapper::TargetMesh:
            return m_target_vertices;
        case Remapper::CoveringMesh:
            return m_covering_source_vertices;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_source_vertices );
    }
}

inline const moab::Range& TempestRemapper::GetMeshVertices( Remapper::IntersectionContext ctx ) const
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_vertices;
        case Remapper::TargetMesh:
            return m_target_vertices;
        case Remapper::CoveringMesh:
            return m_covering_source_vertices;
        case Remapper::DEFAULT:
        default:
            MB_SET_ERR_RET_VAL( "Invalid context passed to GetMeshSet", m_source_vertices );
    }
}

inline void TempestRemapper::SetMeshType( Remapper::IntersectionContext ctx, const std::vector< int >& metadata )
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            m_source_type = static_cast< moab::TempestRemapper::TempestMeshType >( metadata[0] );
            if( metadata[0] == 1 )  // RLL mesh
            {
                m_source_metadata.resize( 2 );
                m_source_metadata[0] = metadata[1];
                m_source_metadata[1] = metadata[2];
            }
            else
            {
                m_source_metadata.resize( 1 );
                m_source_metadata[0] = metadata[1];
            }
            break;
        case Remapper::TargetMesh:
            m_target_type = static_cast< moab::TempestRemapper::TempestMeshType >( metadata[0] );
            if( metadata[0] == 1 )  // RLL mesh
            {
                m_target_metadata.resize( 2 );
                m_target_metadata[0] = metadata[1];
                m_target_metadata[1] = metadata[2];
            }
            else
            {
                m_target_metadata.resize( 1 );
                m_target_metadata[0] = metadata[1];
            }
            break;
        case Remapper::OverlapMesh:
            m_overlap_type = moab::TempestRemapper::OVERLAP_FILES;
        default:
            break;
    }
}

inline TempestRemapper::TempestMeshType TempestRemapper::GetMeshType( Remapper::IntersectionContext ctx ) const
{
    switch( ctx )
    {
        case Remapper::SourceMesh:
            return m_source_type;
        case Remapper::TargetMesh:
            return m_target_type;
        case Remapper::OverlapMesh:
            return m_overlap_type;
        case Remapper::DEFAULT:
        default:
            return TempestRemapper::DEFAULT;
    }
}

inline Mesh* TempestRemapper::GetCoveringMesh()
{
    return m_covering_source;
}

inline moab::EntityHandle& TempestRemapper::GetCoveringSet()
{
    return m_covering_source_set;
}

// inline int TempestRemapper::GetGlobalID( Remapper::IntersectionContext ctx, int localID )
// {
//     switch( ctx )
//     {
//         case Remapper::SourceMesh:
//             return lid_to_gid_src[localID];
//         case Remapper::TargetMesh:
//             return lid_to_gid_tgt[localID];
//         case Remapper::CoveringMesh:
//             return lid_to_gid_covsrc[localID];
//         case Remapper::OverlapMesh:
//         case Remapper::DEFAULT:
//         default:
//             return -1;
//     }
// }

// inline int TempestRemapper::GetLocalID( Remapper::IntersectionContext ctx, int globalID )
// {
//     switch( ctx )
//     {
//         case Remapper::SourceMesh:
//             return gid_to_lid_src[globalID];
//         case Remapper::TargetMesh:
//             return gid_to_lid_tgt[globalID];
//         case Remapper::CoveringMesh:
//             return gid_to_lid_covsrc[globalID];
//         case Remapper::DEFAULT:
//         case Remapper::OverlapMesh:
//         default:
//             return -1;
//     }
// }

}  // namespace moab

#endif  // MB_TEMPESTREMAPPER_HPP
