/** \page h5mmain H5M File Format API
 *
 *\section Intro   Introduction
 *
 * MOAB's native file format is based on the HDF5 file format.  
 * The most common file extension used for such files is .h5m.
 * A .h5m file can be identified by the top-level \c tstt group
 * in the HDF5 file.
 *
 * The API implemented by this library is a wrapper on top of the
 * underlying HDF5 library.  It provides the following features:
 * - Enforces and hides MOAB's expected file layout 
 * - Provides a slightly higher-level API
 * - Provides some backwards compatibility for file layout changes
 *
 *
 *\section Overview   H5M File Layout
 *
 * The H5M file format relies on the use of a unique entity ID space for
 * all vertices, elements, and entity sets stored in the file.  This
 * ID space is defined by the application.  IDs must be unique over all
 * entity types (a vertex and an entity set may not have the same ID.)
 * The IDs must be positive (non-zero) integer values.  
 * There are no other requirements imposed by the format on the ID space.
 *
 * Elements, with the exception of polyhedra, are defined by a list of 
 * vertex IDs.  Polyhedra are defined by a list of face IDs.  Entity sets 
 * have a list of contained entity IDs, and lists of parent and child 
 * entity set IDs.  The set contents may include any valid entity ID,
 * including other sets.  The parent and child lists are expected to
 * contain only entity IDs corresponding to other entity sets.  A zero
 * entity ID may be used in some contexts (tag data with the mhdf_ENTITY_ID
 * property) to indicate a 'null' value,
 *
 * Element types are defined by the combination of a topology identifier (e.g. 
 * hexahedral topology) and the number of nodes in the element.  
 *
 *
 *\section Root   The tstt Group
 *
 * All file data is stored in the \c tstt group in the HDF5 root group.
 * The \c tstt group may have an optional scalar integer attribute 
 * named \c max_id .  This attribute, if present, should contain the
 * value of the largest entity ID used internally to the file.  It can
 * be used to verify that the code reading the file is using an integer
 * type of sufficient size to accommodate the entity IDs.
 *
 * The \c tstt group contains four sub-groups, a datatype object, and a 
 * dataset object.  The four sub-groups are: \c nodes, \c elements,
 * \c sets, and \c tags.  The dataset is named \c history .
 *
 * The \c elemtypes datatype is an enumeration of the elem topologies
 * used in the file.  The element topologies understood by MOAB are:
 * - \c Edge
 * - \c Tri
 * - \c Quad
 * - \c Polygon
 * - \c Tet
 * - \c Pyramid
 * - \c Prism
 * - \c Knife
 * - \c Hex
 * - \c Polyhedron
 * 
 *
 *\section History   The history DataSet
 *
 * The \c history DataSet is a list of variable-length strings with
 * application-defined meaning.  
 *
 *\section Nodes   The nodes Group
 *
 *
 * The \c nodes group contains a single DataSet and an optional
 * subgroup.  The \c tags subgroup is described in the 
 * \ref Dense "section on dense tag storage".  
 *
 * The \c coordinates
 * DataSet contains the coordinates of all vertices in the mesh.
 * The DataSet should contain floating point values and have a dimensions
 * \f$ n \times d \f$, where \c n is the number of vertices and \c d
 * is the number of coordinate values for each vertex.
 *
 * The \c coordinates DataSet must have an integer attribute named \c start_id .
 * The vertices are then defined to have IDs beginning with this value
 * and increasing sequentially in the order that they are defined in the
 * \c coordinates table.
 *
 *
 *\section Elements   The elements Group
 *
 * The \c elements group contains an application-defined number of 
 * subgroups.  Each subgroup defines one or more mesh elements that
 * have the same topology and length of connectivity (number of nodes
 * for any topology other than \c Polyhedron.)  The names of the subgroups
 * are application defined.  MOAB uses a combination of the element
 * topology name and connectivity length (e.g. "Hex8".).  
 *
 * Each subgroup must have an attribute named \c element_type that 
 * contains one of the enumerated element topology values defined 
 * in the \c elemtypes datatype described in a \ref Root "previous section".
 *
 * Each subgroup contains a single DataSet named \c connectivity and an 
 * optional subgroup named \c tags.  The \c tags subgroup is described in the 
 * \ref Dense "section on dense tag storage". 
 *
 * The \c connectivity DataSet is an \f$ n \times m \f$ array of integer
 * values.  The DataSet contains one row for each of the \c n contained
 * elements, where the connectivity of each element contains \c m IDs.  For
 * all element types supported by MOAB, with the exception of polyhedra,
 * the element connectivity list is expected to contain only IDs 
 * corresponding to nodes.  
 *
 * Each element \c connectivity DataSet must have an integer attribute 
 * named \c start_id .  The elements defined in the connectivity table
 * are defined to have IDs beginning with this value and increasing
 * sequentially in the order that they are defined in the table.
 *
 *
 *\section Sets   The sets Group
 *
 * The \c sets group contains the definitions of any entity sets stored
 * in the file.  It contains 1 to 4 DataSets and the optional \c tags 
 * subgroup.  The \c contents, \c parents, and \c children data sets
 * are one dimensional arrays containing the concatenation of the
 * corresponding lists for all of the sets represented in the file.
 *
 * The \c lists DataSet is a \f$ n \times 4 \f$ table, having one
 * row of four integer values for each set.  The first three values
 * for each set are the indices into the \c contents, \c children, 
 * and \c parents DataSets, respectively, at which the \em last value
 * for set is stored.  The contents, child, and parent lists for
 * sets are stored in the corresponding datasets in the same order as
 * the sets are listed in the \c lists DataSet, such that the index of
 * the first value in one of those tables is one greater than the 
 * corresponding end index in the \em previous row of the table.  The
 * number of content entries, parents, or children for a given set can
 * be calculated as the difference between the corresponding end index
 * entry for the current set and the same entry in the previous row 
 * of the table.  If the first set in the \c lists DataSet had no parent
 * sets, then the corresponding index in the third column of the table
 * would be \c -1.  If it had one parent, the index would be \c 0.  If it
 * had two parents, the index would be \c 1, as the first parent would be
 * stored at position 0 of the \c parents DataSet and the second at position
 * 1.
 *
 * The fourth column of the \c lists DataSet is a series of bit flags
 * defining some properties of the sets.  The four bit values currently
 * defined are:
 *  - 0x1 owner
 *  - 0x2 unique
 *  - 0x4 ordered
 *  - 0x8 range compressed
 *
 * The fourth (most significant) bit indicates that, in the \c contents 
 * data set, that the contents list for the corresponding set is stored
 * using a single range compression.  Rather than storing the IDs of the
 * contained entities individually, each ID \c i is followed by a count 
 * \c n indicating that the set contains the contiguous range of IDs
 * \f$ [i, i+n-1] \f$.
 *
 * The three least significant bits specify intended properties of the
 * set and are unrelated to how the set data is stored in the file.  These
 * properties, described briefly from least significant bit to most 
 * significant are: contained entities should track set membership;
 * the set should contain each entity only once (strict set); and
 * that the order of the entries in the set should be preserved.  
 *
 * Similar to the \c nodes/coordinates and \c elements/.../connectivity
 * DataSets, the \c lists DataSet must have an integer attribute 
 * named \c start_id .  IDs are assigned to to sets in the order that
 * they occur in the \c lists table, beginning with the attribute value.
 *
 * The \c sets group may contain a subgroup names \c tags.  The \c tags 
 * subgroup is described in the \ref Dense "section on dense tag storage". 
 *
 * 
 * \section Tags   The tags Group
 *
 * The \c tags group contains a sub-group for each tag defined
 * in the file.  These sub-groups contain the definition of the
 * tag and may contain some or all of the tag values associated with
 * entities in the file.  However, it should be noted that tag values
 * may also be stored in the "dense" format as described in the 
 * \ref Dense "section on dense tag storage".
 *
 * Each sub-group of the \c tags group contains the definition for
 * a single tag.  The name of each sub-group is the name of the 
 * corresponding tag.  Non-printable characters, characters
 * prohibited in group names in the HDF5 file format, and the
 * backslash ('\') character are encoded
 * in the name string by a backslash ('\') character followed by
 * the ASCII value of the character expressed as a pair of hexadecimal
 * digits.  Thus the backslash character would be represented as \c \5C .
 * Each tag group should also contain a comment which contains the
 * unencoded tag name.
 *
 * The tag sub-group may have any or all of the following four attributes:
 * \c default, \c global, \c is_handle, and \c variable_length.  
 * The \c default attribute, if present,
 * must contain a single tag value that is to be considered the 'default'
 * value of the tag.  The \c global attribute, if present, must contain a
 * single tag value that is the value of the tag as set on the mesh instance
 * (MOAB terminology) or root set (ITAPS terminology.)  The presence of the
 * \c is_handle attribute (the value, if any, is meaningless) indicates
 * that the tag values are to be considered entity IDs.  After reading the
 * file, the reader should map any such tag values to whatever mechanism
 * it uses to reference the corresponding entities read from the file.
 * The presence of the \c variable_length attribute indicates that each 
 * tag value is a variable-length array.  The reader should rely on the
 * presence of this attribute rather than the presence of the \c var_indices
 * DataSet discussed below because the file may contain the definition of
 * a variable length tag without containing any values for that tag.  In such
 * a case, the \c var_indices DataSet will not be present.
 *
 * Each tag sub-group will contain a committed type object named \c type .
 * This type must be the type instance used by the \c global and \c default
 * attributes discussed above and any tag value data sets.  For fixed-length
 * tag data, the tag types understood by MOAB are:
 *  - opaque data
 *  - a single floating point value
 *  - a single integer value
 *  - a bit field
 *  - an array of floating point values
 *  - an array of integer values
 * Any other data types will be treated as opaque data.
 * For Variable-length tag data, MOAB expects the \c type object to be
 * one of:
 *  - opaque data
 *  - a single floating point value
 *  - a single integer value
 *
 * For fixed-length tags, the tag sub-group may contain 'sparse' formatted
 * tag data, which is comprised of two data sets: \c id_list and \c values.
 * Both data sets must be 1-dimensional arrays of the same length.  The 
 * \c id_list data set contains a list of entity IDs and the \c values 
 * data set contains a list of corresponding tag values.  The data stored in
 * the \c values table must be of type \c type.  Fixed-length tag values
 * may also be stored in the "dense" format as described in the 
 * \ref Dense "section on dense tag storage".  A mixture of both sparse-
 * and dense-formatted tag values may be present for a single tag.
 *
 * For variable-length tags the tag values, if any, are always stored
 * in the tag sub-group of the \c tags group and are represented by
 * three one-dimensional data sets: \c id_list, \c var_indices, and \c values.  
 * Similar to the fixed-length sparse-formatted tag data, the \c id_list
 * contains the IDs of the entities for which tag values are defined.
 * The \c values dataset contains the concatenation of the tag values
 * for each of the entities referenced by ID in the \c id_list table, 
 * in the order that the entities are referenced in the \c id_list table.
 * The \c var_indices table contains an index into the \c values data set
 * for each entity in \c id_list.  The index indicates the position of
 * the \em last tag value for the entity in \c values.  The index of
 * the first value is one greater than the 
 * corresponding end index for the \em entry in \c var_indices.  The
 * number of tag values for a given entity can
 * be calculated as the difference between the corresponding end index
 * entry for the current entity and the previous value in the \c var_indices
 * dataset.  
 *
 *
 * \section Dense   The tags Sub-Groups
 *
 * Data for fixed-length tags may also be stored in the \c tags sub-group
 * of the \c nodes, \c sets, and subgroups of the \c elements group.  
 * Values for given tag are stored in a dataset within the \c tags sub-group
 * that has the following properties:
 *  - The name must be the same as that of the tag definition in the main 
 *      \c tags group
 *  - The type of the data set must be the committed type object stored
 *      as \c /tstt/tags/\c tagname/type .
 *  - The data set must have the same length as the data set in the
 *    parent group with the \c start_id attribute.  
 *
 * If dense-formatted data is specified for any entity in the group, then
 * it must be specified for every entity in the group.  The table is 
 * expected to contain one value for each entity in the corresponding 
 * primary definition table (\c /tstt/nodes/coordinates , 
 * \c /tstt/elements/\c name/connectivity , or \c /tstt/sets/list), in the
 * same order as the entities in that primary definition table.
 *
 *
 *\section mhdf_set mhdf Meshset data
 *
 * Meshset data is divided into three groups of data.  The set-list/meta-information table,
 * the set contents table and the set children table.  Each is written and read independently.
 *
 * The set list table contains one row for each set.  Each row contains four values:
 * {content list end index, child list end index, parent list end index, and flags}.  The flags 
 * value is a collection of bits with
 * values defined in \ref mhdf_set_flag .  The all the flags except \ref mhdf_SET_RANGE_BIT are
 * saved properties of the mesh data and are not relevant to the actual file in any way.  The
 * \ref mhdf_SET_RANGE_BIT flag is a toggle for how the meshset contents (not children) are saved.
 * It is an internal property of the file format and should not be passed on to the mesh database.
 * The content list end index and child list end index are the indices of the last entry for the
 * set in the contents and children tables respectively.  In the case where a set has either no
 * children or no contents, the last index of should be the same as the last index of the previous
 * set in the table, or -1 for the first set in the table.  Thus the first index is always one
 * greater than the last index of the previous set.  If the first index, calculated as one greater
 * that the last index of the previous set is greater than the last index of the current set, then
 * there are no values in the corresponding contents or children table for that set.
 *
 * The set contents table is a vector of integer global IDs that is the concatenation of the contents
 * data for all of the mesh sets.  The values are stored corresponding to the order of the sets
 * in the set list table.  Depending on the value of \ref mhdf_SET_RANGE_BIT in the flags field of
 * the set list table, the contents for a specific set may be stored in one of two formats.  If the
 * flag is set, the contents list is a list of pairs where each pair is a starting global Id and a 
 * count.  For each pair, the set contains the range of global Ids beginning at the start value. 
 * If the \ref mhdf_SET_RANGE_BIT flag is not set, the meshset contents are a simple list of global Ids.
 *
 * The meshset child table is a vector of integer global IDs.  It is a concatenation of the child
 * lists for all the mesh sets, in the order the sets occur in the meshset list table.  The values
 * are always simple lists.  The child table may never contain ranges of IDs.
 *
 *
 *\section mhdf_tag mhdf Tag data
 *
 * The data for each tag can be stored in two places/formats:  sparse and/or
 * dense.  The data may be stored in both, but there should not be redundant 
 * values for the same entity.  
 *
 * Dense tag data is stored as multiple tables of tag values, one for each
 * element group.  (Note:  special \ref mhdf_ElemHandle values are available
 * for accessing dense tag data on nodes or meshsets via the \ref mhdf_node_type_handle
 * and \ref mhdf_set_type_handle functions.)  Each dense tag table should contain
 * the same number of entries as the element connectivity table.  The tag values 
 * are associated with the corresponding element in the connectivity table.
 *
 * Sparse tag data is stored as a global table pair for each tag type.  The first
 * if the pair of tables is a list of Global IDs.  The second is the corresponding
 * tag value for each entity in the ID list.
 */

