/*! \page userguide User's Guide

  \tableofcontents

  \section introduction 1. Introduction

In scientific computing, systems of partial differential equations (PDEs) are solved on computers.  One of the most widely used methods to solve PDEs numerically is to solve over discrete neighborhoods or "elements" of the domain.  Popular discretization methods include Finite Difference (FD), Finite Element (FE), and Finite Volume (FV).  These methods require the decomposition of the domain into a discretized representation, which is referred to as a "mesh".  The mesh is one of the fundamental types of data linking the various tools in the analysis process (mesh generation, analysis, visualization, etc.).  Thus, the representation of mesh data and operations on those data play a very important role in PDE-based simulations.

MOAB is a component for representing and evaluating mesh data. It is part of the Scalable Interfaces for Geometry and Mesh based Applications (SIGMA) toolkit [1]. MOAB can store structured and unstructured mesh, consisting of elements in the finite element "zoo", along with polygons and polyhedra.  The functional interface to MOAB is simple, consisting of only four fundamental data types.  This data is quite powerful, allowing the representation of most types of metadata commonly found on the mesh.  Internally, MOAB uses array-based storage for fine-grained data, which in many cases provides more efficient access, especially for large portions of mesh and associated data.  MOAB is optimized for efficiency in space and time, based on access to mesh in chunks rather than through individual entities, while also versatile enough to support individual entity access. MOAB is also unique in that it maintains a completely parallel view of the mesh database so that queries to traverse the cells and manipulation of the grid in memory is scalable in parallel.

The MOAB data model consists of the following four fundamental types: mesh interface instance, mesh entities (vertex, edge, tri, etc.), sets, and tags.  Entities are addressed through handles rather than pointers, to allow the underlying representation of an entity to change without changing the handle to that entity.  Sets are arbitrary groupings of mesh entities and other sets.  Sets also support parent/child relationships as a relation distinct from sets containing other sets.  The directed graph provided by set parent/child relationships is useful for embedding graphs whose nodes include collections of mesh entities; this approach has been used to represent a wide variety of application-specific data, including geometric model topology, processor partitions, and various types of search trees.  Tags are named data which can be assigned to the mesh as a whole, individual entities, or sets.  Tags are a mechanism for attaching data to individual entities, and sets are a mechanism for describing relations between entities; the combination of these two mechanisms is a powerful yet simple interface for representing metadata or application-specific data.

Various mesh-related tools are provided with MOAB or can be used directly with MOAB.  These tools can be used for mesh format translation (mbconvert), mesh skinning (Skinner class), solution transfer between meshes (MBCoupler tool), ray tracing and other geometric searches (OrientedBoxTreeTool, AdaptiveKDTree), visualization (vtkMOABReader tool), mesh optimization/smoothing (Mesquite interfaces), and relation between mesh and geometric models (the separately-packed Lasso tool).  These tools are described later in this document.

MOAB is written in the C++ programming language, with applications interacting with MOAB mostly through its moab::Interface class.  All of the MOAB functions and classes are isolated in and accessed through the moab namespace<sup>1</sup>. The remainder of this report gives class and function names without the "moab::" namespace qualification; unless otherwise noted, the namespace qualifier should be added to all class and function names referenced here.  MOAB also implements the iMOAB interface, which is a language-agnostic, lightweight interface that can be called directly from C, C++, Fortran, and Python.  Almost all of the functionality in MOAB can be accessed through the iMOAB interface.  MOAB is developed and supported on the Linux and MacOS operating systems, as well as various HPC operating systems.  MOAB can be used on parallel computing systems as well, including laptops, workstations, clusters and high-end petascale/exascale parallel systems (e.g., NERSC-Perlmutter, OLCF-Frontier, ALCF-Aurora).  MOAB is released under a standard LGPL open source software license.

MOAB is used in several ways in various applications.  MOAB serves as the underlying mesh data representation in several scientific computing applications [2], [3], [4], [5], [6].  MOAB can also be used as a mesh format translator, using readers and writers included in MOAB.  MOAB has also been used as a bridge to couple results in multi-physics analysis and to link these applications with other mesh services [7] in order to consistently solve problems at scale.

The remainder of this guide is organized as follows.  Section 2, "Getting Started", provides a few simple examples of using MOAB to perform simple tasks on a mesh.  Section 3 discusses the MOAB data model in more detail, including some aspects of the implementation.  Section 4 summarizes the MOAB function API and related mesh services.  Section 5 describes some of the tools included with MOAB, and the implementation of mesh readers/writers for MOAB.  Section 6 describes how to build MOAB-based applications.  Section 7 contains a detailed description of the iMOAB mesh interface and its capabilities.  Sections 8 and 9 discuss MOAB's representations of structured and spectral element meshes, respectively.  Section 10 gives helpful hints for accessing MOAB in an efficient manner from applications.  Section 11 discusses parallel computing capabilities in MOAB.  Section 12 covers advanced features including mesh quality assessment, remapping, and coupling.  Section 13 gives a conclusion and future plans for MOAB development.  Section 14 gives references cited in this report.

Several other sources of information about MOAB may also be of interest to readers.  Meta-data conventions define how sets and /or tags are used together to represent various commonly-used simulation constructs in MOAB [8]. MOAB also has an active mailing list [9] to provide guidance and encourage discussions about usage and development of algorithms using MOAB. A separate mailing list [10] for MOAB-related release announcements and future guidance may also be subscribed by users.  Potential users are encouraged to interact with the MOAB team using these mailing lists, and watch the SIGMA website [1] for updates.

<sup>1</sup> Non-namespaced names are also provided for backward compatibility, with the "MB" prefix added to the class or variable name.

 \ref contents

  \section interface 2. MOAB Data Model
The MOAB data model describes the basic types used in MOAB and the language used to communicate that data to applications.  This chapter describes that data model, along with some of the reasons for some of the design choices in MOAB.

 \ref contents

 \subsection twoone 2.1. MOAB Interface

MOAB is written in C++.  The primary interface with applications is through member functions of the abstract base class Interface.  The MOAB library is created by instantiating Core, which implements the Interface API.  Multiple instances of MOAB can exist concurrently in the same application; mesh entities are not shared between these instances<sup>2</sup>.  MOAB is most easily viewed as a database of mesh objects accessed through the instance.  No other assumptions explicitly made about the nature of the mesh stored there; for example, there is no fundamental requirement that elements fill space or do not overlap each other geometrically.

<sup>2</sup> One exception to this statement is when the parallel interface to MOAB is used; in this case, entity sharing between instances is handled explicitly using message passing.  This is described in more detail in Section 11 of this document.

 \ref contents

 \subsection twotwo 2.2. Mesh Entities
MOAB represents the following topological mesh entities: vertex, edge, triangle, quadrilateral, polygon, tetrahedron, pyramid, prism, knife, hexahedron, polyhedron.  MOAB uses the EntityType enumeration to refer to these entity types (see Table 1).  This enumeration has several special characteristics, chosen intentionally: the types begin with vertex, entity types are grouped by topological dimension, with lower-dimensional entities appearing before higher dimensions; the enumeration includes an entity type for sets (described in the next section); and MBMAXTYPE is included at the end of this enumeration, and can be used to terminate loops over type.  In addition to these defined values, the an increment operator (++) is defined such that variables of type EntityType can be used as iterators in loops.
MOAB refers to entities using "handles".  Handles are implemented as long integer data types, with the four highest-order bits used to store the entity type (mesh vertex, edge, tri, etc.) and the remaining bits storing the entity id.  This scheme is convenient for applications because:
- Handles sort lexicographically by type and dimension; this can be useful for grouping and iterating over entities by type.
- The type of an entity is indicated by the handle itself, without needing to call a function.
- Entities allocated in sequence will typically have contiguous handles; this characteristic can be used to efficiently store and operate on large lists of handles.
.

This handle implementation is exposed to applications intentionally, because of optimizations that it enables, and is unlikely to change in future versions.

  \subsection tableone Table 1: Values defined for the EntityType enumerated type.
<table border="1">
<tr>
<td>MBVERTEX = 0</td>
<td>MBPRISM</td>
</tr>
<tr>
<td>MBEDGE</td>
<td>MBKNIFE</td>
</tr>
<tr>
<td>MBTRI</td>
<td>MBHEX</td>
</tr>
<tr>
<td>MBQUAD</td>
<td>MBPOLYHEDRON</td>
</tr>
<tr>
<td>MBPOLYGON</td>
<td>MBENTITYSET</td>
</tr>
<tr>
<td>MBTET</td>
<td>MBMAXTYPE</td>
</tr>
<tr>
<td>MBPYRAMID</td>
<td></td>
</tr>
</table>

MOAB defines a special class for storing lists of entity handles, named Range.  This class stores handles as a series of (start_handle, end_handle) subrange tuples.  If a list of handles has large contiguous ranges, it can be represented in almost constant size using Range.  Since entities are typically created in groups, e.g. during mesh generation or file import, a high degree of contiguity in handle space is typical.  Range provides an interface similar to C++ STL containers like std::vector, containing iterator data types and functions for initializing and iterating over entity handles stored in the range.  Range also provides functions for efficient Boolean operations like subtraction and intersection.  Most API functions in MOAB come in both range-based and vector-based variants.  By definition, a list of entities stored in an Range is always sorted, and can contain a given entity handle only once.  Range cannot store the handle 0 (zero).

Typical usage of an Range object would look like:

\code
using namespace moab;
   int my_function(Range &from_range) {
          int num_in_range = from_range.size();
          Range to_range;
          Range::iterator rit;
    for (rit = from_range.begin(); rit != from_range.end(); ++rit) {
            EntityHandle this_ent = *rit;
            to_range.insert(this_ent);
          }
        }
\endcode

Here, the range is iterated similar to how std::vector is iterated.

  \subsection adjacencies 2.3. Adjacencies & AEntities

The term adjacencies is used to refer to those entities topologically connected to a given entity, e.g. the faces bounded by a given edge or the vertices bounding a given region.  The same term is used for both higher-dimensional (or bounded) and lower-dimensional (or bounding) adjacent entities.  MOAB provides functions for querying adjacent entities by target dimension, using the same functions for higher- and lower-dimension adjacencies.  By default, MOAB stores the minimum data necessary to recover adjacencies between entities.  When a mesh is initially loaded into MOAB, only entity-vertex (i.e. "downward") adjacencies are stored, in the form of entity connectivity.  When "upward" adjacencies are requested for the first time, e.g. from vertices to regions, MOAB stores all vertex-entity adjacencies explicitly, for all entities in the mesh.  Non-vertex entity to entity adjacencies are never stored, unless explicitly requested by the application.

In its most fundamental form, a mesh need only be represented by its vertices and the entities of maximal topological dimension.  For example, a hexahedral mesh can be represented as the connectivity of the hex elements and the vertices forming the hexes.  Edges and faces in a 3D mesh need not be explicitly represented.  We refer to such entities as "AEntities", where 'A' refers to "Auxiliary", "Ancillary", and a number of other words mostly beginning with 'A'.  Individual AEntities are created only when requested by applications, either using mesh modification functions or by requesting adjacencies with a special "create if missing" flag passed as "true".  This reduces the overall memory usage when representing large meshes.  Note entities must be explicitly represented before they can be assigned tag values or added to entity sets (described in following Sections).

\ref contents

 \subsection twothree 2.4. Entity Sets
Entity sets are also known as "mesh sets", or when the context is clear, not to be confused with std::set, just "sets". Entity sets are used to store arbitrary collections of entities and other sets.  Sets are used for a variety of things in mesh-based applications, from the set of entities discretizing a given geometric model entity to the entities partitioned to a specific processor in a parallel finite element application.  MOAB entity sets can also store parent/child relations with other entity sets, with these relations distinct from contains relations.  Parent/child relations are useful for building directed graphs with graph nodes representing collections of mesh entities; this construct can be used, for example, to represent an interface of mesh faces shared by two distinct collections of mesh regions.  MOAB also defines one special set, the "root set" or the interface itself; all entities are part of this set by definition.  Defining a root set allows the use of a single set of MOAB API functions to query entities in the overall mesh as well as its subsets.

MOAB entity sets can be one of two distinct types: list-type entity sets preserve the order in which entities are added to the set, and can store a given entity handle multiple times in the same set; set-type sets are always ordered by handle, regardless of the order of addition to the set, and can store a given entity handle only once.  This characteristic is assigned when the set is created, and cannot be changed during the set's lifetime.

MOAB provides the option to track or not track entities in a set.  When entities (and sets) are deleted by other operations in MOAB, they will also be removed from containing sets for which tracking has been enabled.  This behavior is assigned when the set is created, and cannot be changed during the set's lifetime.  The cost of turning tracking on for a given set is sizeof(EntityHandle) for each entity added to the set; MOAB stores containing sets in the same list which stores adjacencies to other entities.

Using an entity set looks like the following:
\code
using namespace moab;
// load a file using MOAB, putting the loaded mesh into a file set
EntityHandle file_set;
MB_CHK_ERR( moab->create_meshset(MESHSET_SET, file_set) );
MB_CHK_SET_ERR( moab->load_file("fname.vtk", &file_set), "Unable to load fname.vtk" );
Range set_ents;
// get all the 3D entities in the set
MB_CHK_ERR( moab->get_entities_by_dimension(file_set, 3, set_ents) );
\endcode

Entity sets are often used in conjunction with tags (described in the next section), and provide a powerful mechanism to store a variety of meta-data with meshes.

\ref contents

 \subsection twofour 2.5. Tags

Applications of a mesh database often need to attach data to mesh entities.  The types of attached data are often not known at compile time, and can vary across individual entities and entity types.  MOAB refers to this attached data as a "tag".  Tags can be thought of loosely as a variable, which can be given a distinct value for individual entities, entity sets, or for the interface itself.  A tag is referenced using a handle, similarly to how entities are referenced in MOAB.  Each MOAB tag has the following characteristics, which can be queried through the MOAB interface:
- Name
- Size (in bytes)
- Storage type
- Data type (integer, double, opaque, entity handle)
- Handle
.

The storage type determines how tag values are stored on entities.

- Dense: Dense tag values are stored in arrays which match arrays of contiguous entity handles.  Dense tags are more efficient in both storage and memory if large numbers of entities are assigned the same tag.  Storage for a given dense tag is not allocated until a tag value is set on an entity; memory for a given dense tag is allocated for all entities in a given sequence at the same time.
- Sparse: Sparse tags are stored as a list of (entity handle, tag value) tuples, one list per sparse tag, sorted by entity handle.
- Bit: Bit tags are stored similarly to dense tags, but with special handling to allow allocation in bit-size amounts per entity.
.

MOAB also supports variable-length tags, which can have a different length for each entity they are assigned to.  Variable length tags are stored similarly to sparse tags.

The data type of a tag can either be one understood at compile time (integer, double, entity handle), in which case the tag value can be saved and restored properly to/from files and between computers of different architecture (MOAB provides a native HDF5-based save/restore format for this purpose; see Section 4.6).  The opaque data type is used for character strings, or for allocating "raw memory" for use by applications (e.g. for storage application-defined structures or other abstract data types).  These tags are saved and restored as raw memory, with no special handling for endian or precision differences.

An application would use the following code to attach a double-precision tag to vertices in a mesh, e.g. to assign a temperature field to those vertices:

\code
using namespace moab;
// load a file using MOAB and get the vertices
MB_CHK_ERR( moab->load_file("fname.vtk") );
Range verts;
rval = moab->get_entities_by_dimension(0, 0, verts);
// create a tag called "TEMPERATURE"
Tag temperature;
double def_val = -1.0d-300, new_val = 273.0;
MB_CHK_ERR( moab->tag_create("TEMPERATURE", sizeof(double), MB_TAG_DENSE,
                        MB_TYPE_DOUBLE, temperature, &def_val) );
// assign a value to vertices
for (Range::iterator vit = verts.begin();
     vit != verts.end(); ++vit)
  MB_CHK_ERR( moab->tag_set_data(temperature, &(*rit), 1, &new_val) );

\endcode

The semantic meaning of a tag is determined by applications using it.  However, to promote interoperability between applications, there are a number of tag names reserved by MOAB which are intended to be used by convention.  Mesh readers and writers in MOAB use these tag conventions, and applications can use them as well to access the same data. Ref. [1] maintains an up-to-date list of conventions for meta-data usage in MOAB.

  \ref contents

  \section api 3. MOAB API Design Philosophy and Summary

This section describes the design philosophy behind MOAB, and summarizes the functions, data types and enumerated variables in the MOAB API.  A complete description of the MOAB API is available in online documentation in the MOAB distribution [11].

MOAB is designed to operate efficiently on collections of entities.  Entities are often created or referenced in groups (e.g. the mesh faces discretizing a given geometric face, the 3D elements read from a file), with those groups having some form of temporal or spatial locality.  The interface provides special mechanisms for reading data directly into the native storage used in MOAB, and for writing large collections of entities directly from that storage, to avoid data copies.  MOAB applications structured to take advantage of that locality will typically operate more efficiently.

MOAB has been designed to maximize the flexibility of mesh data which can be represented.  There is no explicit constraint on the geometric structure of meshes represented in MOAB, or on the connectivity between elements.  In particular, MOAB allows the representation of multiple entities with the same exact connectivity; however, in these cases, explicit adjacencies must be used to distinguish adjacencies with AEntities bounding such entities.

The number of vertices used to represent a given topological entity can vary, depending on analysis needs; this is often the case in FEA.  For example, applications often use "quadratic" or 10-vertex tetrahedral, with vertices at edge midpoints as well as corners.  MOAB does not distinguish these variants by entity type, referring to all variants as "tetrahedra".  The number of vertices for a given entity is used to distinguish the variants, with canonical numbering conventions used to determine placement of the vertices [12].  This is similar to how such variations are represented in the Exodus [13] and Patran [14] file formats.  In practice, we find that this simplifies coding in applications, since in many cases the handling of entities depends only on the number of corner vertices in the element.  Some MOAB API functions provide a flag which determines whether corner or all vertices are requested.

The MOAB API is designed to balance complexity and ease of use.  This balance is evident in the following general design characteristics:

- Entity lists: Lists of entities are passed to and from MOAB in a variety of forms.  Lists output from MOAB are passed as either STL vector or Range data types.  Either of these constructs may be more efficient in both time and memory, depending on the semantics of the data being requested.  Input lists are passed as either Range's, or as a pointer to EntityHandle and a size.  The latter allows the same function to be used when passing individual entities, without requiring construction of an otherwise unneeded STL vector.
- Entity sets: Most query functions accept an entity set as input.  Applications can pass zero to indicate a request for the whole interface.  Note that this convention applies only to query functions; attempts to add or subtract entities to/from the interface using set-based modification functions, or to add parents or children to the interface set, will fail.  Allowing specification of the interface set in this manner avoids the need for a separate set of API functions to query the database as a whole.
- Implicit Booleans in output lists: A number of query functions in MOAB allow specification of a Boolean operation (Interface::INTERSECT or Interface::UNION).  This operation is applied to the results of the query, often eliminating the need for code the application would need to otherwise implement.  For example, to find the set of vertices shared by a collection of quadrilaterals, the application would pass that list of quadrilaterals to a request for vertex adjacencies, with Interface::INTERSECT passed for the Boolean flag.  The list of vertices returned would be the same as if the application called that function for each individual entity, and computed the intersection of the results over all the quadrilaterals.  Applications may also input non-empty lists to store the results, in which case the intersection is also performed with entities already in the list.  In many cases, this allows optimizations in both time and memory inside the MOAB implementation.

Since these objectives are at odds with each other, tradeoffs had to be made between them.  Some specific issues that came up are:

- Using ranges: Where possible, entities can be referenced using either ranges (which allow efficient storage of long lists) or STL vectors (which allow list order to be preserved), in both input and output arguments.
- Entities in sets: Accessing the entities in a set is done using the same functions which access entities in the entire mesh.  The whole mesh is referenced by specifying a set handle of zero<sup>3</sup>.
- Entity vectors on input: Functions which could normally take a single entity as input are specified to take a vector of handles instead.  Single entities are specified by taking the address of that entity handle and specifying a list length of one.  This minimizes the number of functions, while preserving the ability to input single entities.<sup>4</sup>

<sup>3</sup> This convention applies only to query functions; attempts to add or subtract entities to/from the interface using set-based modification functions, or to add parents or children to the interface set, will fail.

<sup>4</sup> This design choice was made to minimize the number of functions in the API, while preserving the ability to input single entities.  The tradeoff is that single entities must be passed as pointers, which can be confusing to some users.

  \ref contents

  \section apiservices 4. MOAB Function API and Related Mesh Services

This section provides a comprehensive overview of the MOAB function API and the related mesh services that are available with MOAB. The MOAB API is designed to provide efficient access to mesh data while maintaining flexibility for various application needs.

\subsection api_overview 4.1. API Overview

The MOAB API is organized into several functional categories:

- **Entity Management**: Functions for creating, querying, and modifying mesh entities
- **Set Operations**: Functions for working with entity sets and their relationships
- **Tag Operations**: Functions for attaching and querying data on entities
- **Adjacency Queries**: Functions for traversing mesh connectivity
- **File I/O**: Functions for reading and writing mesh files
- **Parallel Operations**: Functions for distributed mesh operations

\subsection core_functions 4.2. Core API Functions

The core MOAB API provides fundamental operations on mesh entities:

\subsubsection entity_ops Entity Operations
- \c get_entities_by_type(): Retrieve entities of a specific type
- \c get_entities_by_dimension(): Retrieve entities of a specific dimension
- \c get_connectivity(): Get vertex connectivity for entities
- \c get_coords(): Get coordinates for vertices
- \c create_vertices(): Create new vertices
- \c create_elements(): Create new elements

\subsubsection set_ops Set Operations
- \c create_meshset(): Create a new entity set
- \c add_entities(): Add entities to a set
- \c remove_entities(): Remove entities from a set
- \c get_entities_by_handle(): Get entities in a set
- \c add_parent_child(): Establish parent-child relationships

\subsubsection tag_ops Tag Operations
- \c tag_get_handle(): Get or create a tag
- \c tag_set_data(): Set tag values on entities
- \c tag_get_data(): Get tag values from entities
- \c tag_delete(): Delete a tag
- \c tag_get_length(): Get tag size information

\subsubsection adjacency_ops Adjacency Operations
- \c get_adjacencies(): Get adjacent entities
- \c get_entities_by_handle(): Get entities in a set
- \c get_entities_by_type_and_tag(): Get entities by type and tag values

\subsection mesh_services 4.3. Related Mesh Services

MOAB provides several mesh-related services and tools that extend its core functionality:

\subsubsection format_conversion Format Conversion
MOAB includes comprehensive support for mesh format conversion through the \c mbconvert tool and built-in readers/writers:

- **Input Formats**: Exodus, VTK, Gmsh, Cubit, Nastran, and many others
- **Output Formats**: HDF5 (native), Exodus, VTK, Gmsh, and others
- **Batch Processing**: Convert multiple files with different formats
- **Metadata Preservation**: Maintain tags and sets across format conversions

\subsubsection mesh_analysis Mesh Analysis Tools
MOAB provides tools for analyzing mesh properties and quality:

- **Topology Analysis**: Analyze mesh connectivity and topology
- **Geometry Analysis**: Compute geometric properties (areas, volumes, etc.)
- **Quality Assessment**: Evaluate mesh quality metrics
- **Statistics Generation**: Generate comprehensive mesh statistics

\subsubsection mesh_modification Mesh Modification
Tools for modifying and manipulating meshes:

- **Skinner Class**: Extract mesh boundaries and skins
- **Mesh Refinement**: Tools for mesh refinement and coarsening
- **Mesh Smoothing**: Improve mesh quality through vertex smoothing
- **Entity Creation**: Create new entities from existing ones

\subsubsection spatial_queries Spatial Query Tools
Efficient spatial search and query capabilities:

- **OrientedBoxTreeTool**: Fast spatial queries using oriented bounding boxes
- **AdaptiveKDTree**: Efficient k-d tree implementation for spatial searches
- **Point Location**: Find containing elements for points
- **Ray Tracing**: Perform ray-mesh intersection queries

\subsubsection coupling_tools Coupling and Transfer Tools
Tools for coupling different mesh-based applications:

- **MBCoupler**: Primary tool for coupling different mesh-based applications
- **Solution Transfer**: Transfer solution data between different meshes
- **Interface Management**: Manage interfaces between different physics domains
- **Parallel Coupling**: Support for parallel coupling operations

\subsubsection visualization_tools Visualization Support
Tools for visualizing mesh data:

- **vtkMOABReader**: VTK-based visualization support
- **HDF5 Visualization**: Tools for visualizing HDF5-based mesh data
- **Quality Visualization**: Visualize mesh quality metrics
- **Statistics Visualization**: Generate visual reports of mesh statistics

\subsubsection quality_tools Quality Assessment and Improvement
Integration with the Mesquite library for mesh quality:

- **Quality Metrics**: Comprehensive mesh quality assessment
- **Quality Improvement**: Tools for improving mesh quality
- **Smoothing Algorithms**: Various smoothing and optimization algorithms
- **Quality Reports**: Generate detailed quality assessment reports

\subsection api_examples 4.4. API Usage Examples

Here are some common examples of using the MOAB API:

\subsubsection basic_usage Basic Usage Example
\code
using namespace moab;
Core moab;
Interface* mb = &moab;

// Load a mesh file
ErrorCode rval = mb->load_file("mesh.h5m");
MB_CHK_ERR(rval);

// Get all vertices
Range vertices;
rval = mb->get_entities_by_dimension(0, 0, vertices);
MB_CHK_ERR(rval);

// Get coordinates
std::vector<double> coords;
rval = mb->get_coords(vertices, coords);
MB_CHK_ERR(rval);
\endcode

\subsubsection tag_usage Tag Usage Example
\code
// Create a temperature tag
Tag temp_tag;
double default_temp = 273.15;
rval = mb->tag_get_handle("TEMPERATURE", 1, MB_TYPE_DOUBLE,
                         temp_tag, MB_TAG_DENSE, &default_temp);
MB_CHK_ERR(rval);

// Set temperature values on vertices
std::vector<double> temps(vertices.size(), 300.0);
rval = mb->tag_set_data(temp_tag, vertices, temps.data());
MB_CHK_ERR(rval);
\endcode

\subsubsection set_usage Set Usage Example
\code
// Create a set for boundary faces
EntityHandle boundary_set;
rval = mb->create_meshset(MESHSET_SET, boundary_set);
MB_CHK_ERR(rval);

// Get all faces
Range faces;
rval = mb->get_entities_by_dimension(0, 2, faces);
MB_CHK_ERR(rval);

// Add faces to boundary set
rval = mb->add_entities(boundary_set, faces);
MB_CHK_ERR(rval);
\endcode

\subsection performance_considerations 4.5. Performance Considerations

When using the MOAB API, consider the following performance guidelines:

- **Bulk Operations**: Use bulk operations when possible instead of individual entity operations
- **Range Usage**: Use Range objects for large entity collections
- **Tag Storage**: Choose appropriate tag storage types (dense vs. sparse)
- **Adjacency Caching**: Be aware of when adjacencies are created and cached
- **Memory Management**: Use appropriate data structures for your use case

\subsection error_handling 4.6. Error Handling

MOAB uses a comprehensive error handling system:

- **ErrorCode**: All functions return an ErrorCode indicating success or failure
- **Error Macros**: Use MB_CHK_ERR and related macros for consistent error handling
- **Error Messages**: Provide descriptive error messages for debugging
- **Recovery**: Design applications to handle and recover from errors gracefully

  \ref contents

  \section tools 5. MOAB Tools and Utilities

MOAB comes with a comprehensive set of tools and utilities for mesh processing and analysis.

\subsection conversion 5.1. Mesh Format Conversion

The mbconvert tool provides format conversion capabilities:

- **Input Formats**: Support for reading various mesh formats (Exodus, VTK, Gmsh, etc.)
- **Output Formats**: Support for writing to various mesh formats
- **Batch Processing**: Support for converting multiple files
- **Parallel Conversion**: Parallel I/O for large mesh files

\subsection visualization 5.2. Visualization Support

MOAB provides several visualization tools:

- **VTK Integration**: Native VTK reader for ParaView integration
- **HDF5 Visualization**: Tools for visualizing HDF5-based mesh data
- **Mesh Statistics**: Tools for generating mesh statistics and reports
- **Quality Visualization**: Tools for visualizing mesh quality metrics

\subsection analysis 5.3. Mesh Analysis Tools

MOAB provides various mesh analysis capabilities:

- **Mesh Statistics**: Comprehensive mesh statistics and reporting
- **Topology Analysis**: Tools for analyzing mesh topology
- **Geometry Analysis**: Tools for analyzing geometric properties
- **Performance Analysis**: Tools for analyzing mesh operation performance

  \ref contents

  \section building 6. Building MOAB-Based Applications

This section describes how to build applications that use MOAB, including compilation, linking, and deployment considerations.

\subsection compilation 6.1. Compilation Requirements

MOAB applications require the following for compilation:

- **C++ Compiler**: MOAB requires a C++11 (or later) compatible compiler
- **CMake**: MOAB uses CMake as its build system
- **Dependencies**: Various optional dependencies depending on features used
- **Platform Support**: Linux, macOS, and various HPC platforms

\subsection cmake_integration 6.2. CMake Integration

The recommended way to build MOAB applications is using CMake:

\code
cmake_minimum_required(VERSION 3.10)
project(my_moab_app)

# Find MOAB
find_package(MOAB REQUIRED)

# Add executable
add_executable(my_app main.cpp)

# Link with MOAB
target_link_libraries(my_app MOAB::MOAB)

# Include directories
target_include_directories(my_app PRIVATE ${MOAB_INCLUDE_DIRS})
\endcode

\subsection dependencies 6.3. Dependencies

MOAB has several optional dependencies:

- **HDF5**: Required for native H5M format support
- **NetCDF**: Required for Exodus format support
- **MPI**: Required for parallel computing features
- **VTK**: Optional for VTK format support
- **CGNS**: Optional for CGNS format support

\subsection linking 6.4. Linking Considerations

When linking MOAB applications:

- **Static vs. Dynamic**: MOAB can be linked statically or dynamically
- **Runtime Libraries**: Ensure all runtime libraries are available
- **MPI**: Use the same MPI implementation for parallel applications
- **Thread Safety**: MOAB is thread-safe for read operations

\subsection deployment 6.5. Deployment

Considerations for deploying MOAB applications:

- **Library Paths**: Ensure MOAB libraries are in the library path
- **Runtime Dependencies**: Include all required runtime libraries
- **Platform Compatibility**: Test on target platforms
- **Performance Tuning**: Optimize for target hardware

  \ref contents

  \section imoab 7. iMOAB Interface

The iMOAB interface provides a lightweight, language-agnostic C/Fortran interface to MOAB that exposes the underlying MOAB Interface object so that queries can be made seamlessly. This interface is designed to be simple, efficient, and portable across different programming languages and computing platforms. The iMOAB interface is the recommended approach for applications that need to integrate MOAB functionality across multiple programming languages or require a simplified interface to MOAB's capabilities.

\subsection imoab_design 7.1. Design Philosophy

The iMOAB interface follows several key design principles:

- **POD Types**: All data in the interface are exposed via Plain Old Data (POD) types to ensure maximum compatibility across languages.
- **Reference Passing**: Everything is passed by reference to avoid language-specific calling conventions.
- **User-Allocated Arrays**: Arrays are allocated by the user code, eliminating concerns about de-allocation of data by the interface.
- **Explicit Sizing**: Always pass the pointer to the start of array along with the total allocated size for the array.
- **Direct MOAB Access**: iMOAB provides direct access to the underlying MOAB Interface object.

\subsection imoab_features 7.2. Key Features

The iMOAB interface provides comprehensive access to MOAB's capabilities:

- **Entity Management**: Create, query, and modify mesh entities
- **Set Operations**: Work with entity sets and their relationships
- **Tag Operations**: Attach and query data on entities
- **File I/O**: Read and write mesh files in various formats
- **Parallel Operations**: Support for distributed mesh operations
- **Spatial Queries**: Efficient spatial search capabilities
- **Quality Assessment**: Mesh quality evaluation and improvement

\subsection imoab_languages 7.3. Language Bindings

iMOAB provides native bindings for C and Fortran programming languages:

- **C/C++**: Direct C interface with C++ wrapper classes
- **Fortran**: Support for Fortran 77, 90, 2003, and 2008
- **Cross-Language**: Seamless integration between C and Fortran bindings

Note: MOAB itself provides Python bindings through the pymoab package, which is separate from the iMOAB interface.

\subsection imoab_examples 7.4. Usage Examples

Here are examples of using iMOAB in C and Fortran:

**C Example:**
\code
#include "iMOAB.h"

// Initialize iMOAB
int appID;
iMOAB_Initialize(&appID);

// Load a mesh file
int fileID;
iMOAB_LoadMesh(appID, "mesh.h5m", &fileID);

// Get mesh information
int num_vertices, num_elements;
iMOAB_GetMeshInfo(appID, fileID, &num_vertices, &num_elements);

// Get vertex coordinates
double* coords = malloc(3 * num_vertices * sizeof(double));
iMOAB_GetVertexCoords(appID, fileID, coords, num_vertices);

// Clean up
free(coords);
iMOAB_Finalize(appID);
\endcode

**Fortran Example:**
\code
program imoab_example
  use imoab
  implicit none

  integer :: appID, fileID, num_vertices, num_elements
  real(8), allocatable :: coords(:)

  ! Initialize iMOAB
  call imoab_initialize(appID)

  ! Load a mesh file
  call imoab_load_mesh(appID, "mesh.h5m", fileID)

  ! Get mesh information
  call imoab_get_mesh_info(appID, fileID, num_vertices, num_elements)

  ! Get vertex coordinates
  allocate(coords(3 * num_vertices))
  call imoab_get_vertex_coords(appID, fileID, coords, num_vertices)

  ! Clean up
  deallocate(coords)
  call imoab_finalize(appID)
end program imoab_example
\endcode

\subsection imoab_parallel 7.5. Parallel Computing Support

iMOAB provides comprehensive support for parallel computing:

- **Parallel I/O**: Efficient parallel reading and writing of mesh files
- **Ghost Layer Management**: Automatic management of ghost entities
- **Load Balancing**: Tools for redistributing mesh data across processes
- **Inter-Process Communication**: Efficient communication patterns for mesh operations
- **Parallel Coupling**: Support for coupling different mesh-based applications in parallel

\subsection imoab_advanced 7.6. Advanced Features

iMOAB supports advanced MOAB features:

- **Spectral Elements**: Support for high-order finite elements
- **Structured Meshes**: Support for structured mesh operations
- **Quality Assessment**: Integration with mesh quality tools
- **Spatial Queries**: Efficient spatial search capabilities
- **Remapping**: Tools for transferring data between different meshes

  \ref contents

  \section structured 8. Structured Mesh Support

MOAB provides comprehensive support for structured meshes through the ScdInterface class. This section describes the structured mesh capabilities and usage.

\subsection structured_representation 8.1. Structured Mesh Representation

MOAB represents structured meshes using a parametric (i,j,k) coordinate system:

- **Parametric Space**: Each structured block has a parametric space defined by (imin,imax) × (jmin,jmax) × (kmin,kmax)
- **Entity Allocation**: Entities are allocated with i varying fastest, then j, then k
- **Periodic Support**: Support for periodic meshes in i and/or j directions
- **Parallel Support**: Full parallel support for structured meshes

\subsection structured_creation 8.2. Creating Structured Meshes

Structured meshes can be created using the ScdInterface:

\code
using namespace moab;

// Create a structured mesh interface
ScdInterface* scdi;
ErrorCode rval = moab->query_interface(scdi);

// Create a structured box
ScdBox* box;
rval = scdi->construct_box(HomCoord(0,0,0), HomCoord(10,10,10),
                          HomCoord(0,0,0), HomCoord(10,10,10), box);

// Get the structured mesh data
EntityHandle start_vertex, start_element;
rval = box->get_vertex_handle(start_vertex);
rval = box->get_element_handle(start_element);
\endcode

\subsection structured_operations 8.3. Structured Mesh Operations

MOAB provides various operations on structured meshes:

- **Parametric Access**: Access entities using parametric coordinates
- **Adjacency Queries**: Efficient adjacency queries for structured meshes
- **Mesh Modification**: Tools for modifying structured meshes
- **Parallel Operations**: Full parallel support for structured mesh operations

  \ref contents

  \section spectral 9. Spectral Element Support

MOAB provides specialized support for spectral element methods, which are commonly used in high-order finite element applications.

\subsection spectral_elements 9.1. Spectral Elements

Spectral elements are high-order finite elements that use specific quadrature points:

- **Gauss-Lobatto Points**: Standard quadrature points for spectral elements
- **Higher-Order Support**: Support for elements with multiple nodes per edge/face
- **Spectral Interpolation**: Tools for spectral interpolation and projection
- **Quality Assessment**: Specialized quality metrics for spectral elements

\subsection spectral_tools 9.2. Spectral Element Tools

MOAB provides several tools for working with spectral elements:

- **SpectralMeshTool**: Primary tool for spectral element operations
- **Higher-Order Factory**: Tools for creating higher-order elements
- **Spectral Interpolation**: Tools for interpolating between different spectral orders
- **Quality Metrics**: Specialized quality metrics for spectral elements

  \ref contents

  \section optimization 10. Performance Optimization

This section provides guidelines for optimizing MOAB performance in your applications.

\subsection memory 10.1. Memory Management

- **Range Usage**: Use Range objects for large entity lists to reduce memory usage
- **Tag Storage**: Choose appropriate tag storage types (dense vs. sparse) based on usage patterns
- **Entity Creation**: Create entities in batches rather than individually
- **Memory Cleanup**: Properly clean up unused entities and tags

\subsection algorithm_optimization 10.2. Algorithm Optimization

- **Batch Operations**: Use batch operations whenever possible
- **Spatial Locators**: Use spatial locators for repeated point location queries
- **Adjacency Caching**: Cache adjacency information for repeated queries
- **Parallel Efficiency**: Use parallel operations efficiently

\subsection io 10.3. I/O Optimization

- **Parallel I/O**: Use parallel I/O for large mesh files
- **File Format Selection**: Choose appropriate file formats for your use case
- **Compression**: Use compression for large mesh files
- **Incremental Loading**: Load mesh data incrementally when possible

  \ref contents

  \section parallel 11. Parallel Computing in MOAB

MOAB provides comprehensive support for parallel computing, enabling scalable mesh operations on distributed memory systems. This section describes the parallel computing capabilities and best practices for using MOAB in parallel applications.

\subsection parallel_architecture 11.1. Parallel Architecture

MOAB's parallel architecture is based on the following principles:

- **Distributed Memory Model**: Each process owns a portion of the mesh
- **Ghost Layer Exchange**: Processes maintain copies of boundary entities from neighboring processes
- **Global Consistency**: MOAB maintains global consistency of mesh data across all processes
- **Scalable Communication**: Efficient communication patterns for large-scale parallel systems

\subsection parallel_comm 11.2. Parallel Communication

The ParallelComm class provides the core parallel communication functionality:

- **Entity Sharing**: Automatic detection and management of shared entities
- **Ghost Layer Management**: Exchange of ghost entities between processes
- **Load Balancing**: Redistribution of mesh data across processes
- **Inter-Process Communication**: Efficient communication patterns for mesh operations

\subsection parallel_io 11.3. Parallel I/O

MOAB supports efficient parallel I/O operations:

- **Parallel File Reading**: Multiple processes can read different portions of a mesh file simultaneously
- **Parallel File Writing**: Distributed mesh data can be written efficiently to a single file
- **HDF5 Integration**: Native support for HDF5-based parallel I/O
- **Format Support**: Parallel I/O for various mesh formats including Exodus, VTK, and MOAB's native H5M format

\subsection parallel_tools 11.4. Parallel Tools and Utilities

MOAB provides several tools for parallel mesh operations:

- **MBCoupler**: Parallel coupling tool for transferring data between different meshes
- **Parallel Mesh Refinement**: Tools for parallel mesh refinement and coarsening
- **Load Balancing**: Utilities for redistributing mesh data to improve load balance
- **Parallel Quality Assessment**: Tools for assessing mesh quality in parallel

  \ref contents

  \section advanced 12. Advanced Features

MOAB provides several advanced features that extend its capabilities beyond basic mesh operations.

\subsection quality 12.1. Mesh Quality Assessment

MOAB integrates with the Mesquite library to provide comprehensive mesh quality assessment and improvement:

- **Quality Metrics**: Various metrics for assessing element quality (aspect ratio, Jacobian, etc.)
- **Quality Improvement**: Tools for improving mesh quality through vertex smoothing and optimization
- **Quality Monitoring**: Real-time monitoring of mesh quality during mesh generation or modification

\subsection remapping 12.2. Remapping and Interpolation

MOAB provides advanced remapping capabilities, particularly for climate science applications:

- **TempestRemap Integration**: Support for high-order conservative remapping
- **Bilinear Interpolation**: Standard bilinear interpolation for structured meshes
- **Conservative Remapping**: Mass-conserving interpolation schemes
- **Online and Offline Maps**: Support for both pre-computed and on-the-fly remapping

\subsection coupling 12.3. Multi-Physics Coupling

MOAB provides tools for coupling multiple physics applications:

- **MBCoupler**: Primary tool for coupling different mesh-based applications
- **Data Transfer**: Efficient transfer of solution data between different meshes
- **Interface Management**: Tools for managing interfaces between different physics domains
- **Load Balancing**: Support for load balancing across coupled applications

\subsection spatial 12.4. Spatial Queries and Geometric Operations

MOAB provides efficient spatial query capabilities:

- **AdaptiveKDTree**: Efficient k-d tree implementation for spatial queries
- **OrientedBoxTreeTool**: Tool for oriented bounding box queries
- **Ray Tracing**: Efficient ray tracing for geometric queries
- **Point Location**: Fast point location in unstructured meshes

\subsection spectral_advanced 12.5. Spectral Element Support

MOAB provides specialized support for spectral element methods:

- **Higher-Order Elements**: Support for high-order finite elements
- **Spectral Mesh Tools**: Specialized tools for spectral element meshes
- **Gauss-Lobatto Points**: Support for Gauss-Lobatto quadrature points
- **Spectral Interpolation**: Tools for spectral interpolation and projection

  \ref contents

  \section conclusion 13. Conclusion and Future Plans

MOAB continues to evolve to meet the needs of the scientific computing community. Current development efforts focus on:

- **Enhanced Parallel Performance**: Continued improvements in parallel scalability
- **Advanced Coupling**: Enhanced tools for multi-physics coupling
- **New File Formats**: Support for emerging mesh file formats
- **Cloud Computing**: Support for cloud-based mesh operations
- **Machine Learning Integration**: Integration with machine learning frameworks for mesh generation and optimization

MOAB's open-source development model encourages community contributions and ensures that the library remains responsive to user needs. The active development community and comprehensive documentation make MOAB an excellent choice for mesh-based scientific computing applications.

  \ref contents

  \section references 14. References

[1] Mahadevan, Vijay S., Iulian Grindeanu, Rajeev Jain, Patrick Shriwise, Navamita Ray, Paul Wilson, Tautges, Timothy J., "SIGMA -- MOAB.", URL: http://sigma.mcs.anl.gov/

[2] Mahadevan, Vijay S., Iulian Grindeanu, Robert Jacob, and Jason Sarich. "Improving climate model coupling through a complete mesh representation: a case study with E3SM (v1) and MOAB (v5. x).", Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2018-280, in review, 2018.

[3] Mahadevan, Vijay S., Elia Merzari, Timothy Tautges, Rajeev Jain, Aleksandr Obabko, Michael Smith, and Paul Fischer. "High-resolution coupled physics solvers for analysing fine-scale nuclear reactor design problems." Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences 372, no. 2021 (2014): 20130381.

[4] Yan, Mi, Kirk Jordan, Dinesh Kaushik, Michael Perrone, Vipin Sachdeva, Timothy J. Tautges, and John Magerlein. "Coupling a basin modeling and a seismic code using MOAB." Procedia Computer Science 9 (2012): 986-993.

[5] Jacob, Robert, Jayesh Krishna, Xiabing Xu, Tim Tautges, Iulian Grindeanu, Rob Latham, Kara Peterson et al. "ParNCL and ParGAL: Data-parallel tools for postprocessing of large-scale Earth science data." Procedia Computer Science 18 (2013): 1245-1254.

[6] Bohm, Tim D., Mohamed E. Sawan, Steve T. Jackson, and Paul PH Wilson. "Detailed nuclear analysis of ITER ELM coils." Fusion Engineering and Design 87, no. 5-6 (2012): 657-661.

[7] Tautges, Timothy J., and Alvaro Caceres. "Scalable parallel solution coupling for multiphysics reactor simulation." In Journal of Physics: Conference Series, vol. 180, no. 1, p. 012017. IOP Publishing, 2009.

[8] Tautges, Timothy J. "MOAB-SD: integrated structured and unstructured mesh representation." Engineering With Computers 20, no. 3 (2004): 286-293.

[9] MOAB Users and Developers Email List., moab-dev@mcs.anl.gov.

[10] MOAB Announcement Email List., moab-announce@mcs.anl.gov.

[11] MOAB online documentation., http://ftp.mcs.anl.gov/pub/fathom/moab-docs/index.html

[12] T.J. Tautges, “Canonical numbering systems for finite-element codes,” Communications in Numerical Methods in Engineering,  vol. Online, Mar. 2009.

[13] L.A. Schoof and V.R. Yarberry, EXODUS II: A Finite Element Data Model,  Albuquerque, NM: Sandia National Laboratories, 1994.

[14] M. PATRAN, “PATRAN User’s Manual,” 2005.

[15] VisIt User's Guide.

  \ref contents

*/
