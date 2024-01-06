-----------------------------
# Frequently Asked Questions
-----------------------------

## Where's the configure script?

  After installing the GNU autotools suite, execute the following
  command to generate the configure script (and other necessary
  generated files):
```  
    autoreconf -fi
```    
  If for some reason, the autoreconf command is not available,
  the following sequence of commands should have the same result:
```
    autoheader
    aclocal -I m4
    libtoolize -f
    autoconf 
    automake -a
```

1. Why aren't the configure script and other generated files in the Git repository?
  
    > Because they are generated files and often machine dependent.  Why save a version history for them?  Further, some of the above commands get re-run automatically when Makefile.am's or other files are changed.  This could lead to compatibility problems if some of the generated files in the Git repository are from a different version of the GNU autotools.

2. Aren't we requiring users to have GNU autotools installed in order to configure MOAB?

    > No.  Developers (or anyone else using source directly from the Git repository) must have the autotools installed.  When creating a tarball for distribution of MOAB, the commands below should be run. The resulting tarball will contain all necessary generated files, including the configure script.

3. What needs to be done to prepare MOAB for distribution as a tarball?
    - Check out a clean copy of MOAB.
    - Execute the following commands in the top-most directory:
        - `autoreconf -fi`
        - `./configure`
    - To create a distributable tarball from a working source directory, do `make dist`

-----------------------
Supported File Formats
-----------------------

Some of the file formats listed below may not be supported by a particular build of MOAB depending on the availability of external libraries.  An  up-to-date list of file formats supported in a particular build of MOAB can be obtained programatically using the MBReaderWriterSet API or as a simple list using the '`-l`' option of the mbconvert utility.

```
Format               Name     Read    Write   File name description
------------------  ------  -------- -------  ----------------
MOAB native (HDF5)  MOAB      yes      yes     h5m mhdf
Exodus II           EXODUS    yes      yes     exo exoII exo2 g gen
Climate NC          NC        yes      yes     nc
IDEAS Format        UNV       yes      no      unv
MCNP5 Format        MESHTAL   yes      no      meshtal
NASTRAN format      NAS       yes      no      nas bdf
Abaqus mesh format  ABAQUS    yes      no      abq
Kitware VTK         VTK       yes      yes     vtk
RPI SMS             SMS       yes      no      sms
Cubit               CUBIT     yes      no      cub
QSlim format        SMF       yes      yes     smf
SLAC                SLAC      no       yes     slac
GMV                 GMV       no       yes     gmv
Ansys               ANSYS     no       yes     ans
Gmsh                GMSH      yes      yes     msh gmsh
Stereo Lithography  STL       yes      yes     stl
TetGen mesh files   TETGEN    yes      no      node ele face edge
Template input      TEMPLATE  yes      yes     
```
Any of the values from the `Name` column may be passed as an option to the file I/O methods to request a particular file.  If no file  format is specified, the default is to choose the write format using the file extension and to try all file readers until one succeeds.

---------------
File IO Options
---------------

An options list as passed to MOAB file IO routines is a single C-style string containing the concatenation of a list of string options, where individual options are separated by a designated separator character. The default separator character is a semicolon (**`;`**).  To specify an alternate separator character, begin the options string with a semicolon followed by the desired separator.  Options are not case sensitive.

---------------
Common Options
---------------

`PRECISION=<N>`: Specify the precision to use when writing float and double values (such as node coordinates) to text-based file formats.

`CREATE`: Do not overwrite existing file.

`FACET_DISTANCE_TOLERANCE=<D>`: Max. distance deviation between facet and geometry, default:0.001.

`FACET_NORMAL_TOLERANCE=<N>`: Max. normal angle deviation (degrees) between facet and geometry, default:5.

`MAX_FACET_EDGE_LENGTH=<D>`: Max. facet edge length, default:0.

`CGM_ATTRIBS={yes|no}`: Actuation of all CGM attributes, default:no.

`DEBUG_IO=n`: Set threashold for debug messages from low-level (format-specific) reader/writer.  Default is `0` (none).

-------------------
Parallel IO Options
-------------------

MOAB must be built with parallel support enabled before these options can be used.

Parallel Read Options:
----------------------

`PARALLEL={NONE|BCAST|BCAST_DELETE|READ_DELETE|READ_PART}`: Set parallel read mode.

Options are:
  - **NONE**  - force serial read/write on each processor (default)
  - **BCAST** - read on one processor and broadcast a copy to all others
  - **BCAST_DELETE** - read file on one processor, broadcasting partitioned
                   data to other processors.
  - **READ_DELETE** - read file on all processors and partition by deleting
                  mesh from non-local partitions.
  - **READ_PART** - read only relevant subset of file on each processor,
                utilizing format-specific parallel I/O if available
                (this option is not supported for all file formats.)
  - **WRITE_PART** - write only relevant subset of file on each processor,
                 creating a single file utilizing format-specific parallel 
                 I/O if available (this option is not supported for all 
                 file formats.)
  - **FORMAT** - deprecated (use WRITE_PART)

`PARALLEL_RESOLVE_SHARED_ENTS`: Resolve which entities are shared between which processes, such that propogation of data across processes can be done.  This should
probably be the default behavior, as you almost certainly want this unless, perhaps, PARALLEL=BCAST.

  
`PARTITION` OR `PARTITION=<tag_name>`: Specify that mesh should be partitioned using partition IDs stored in a tag.  If the tag name is not specified, the default ("PARTITION") is used.


`PARTITION_VAL=<int_list>`: If a tag name is specified to the '`PARTITION`' option, then treat as partitions only those sets for which the tag *value* is a single integer
for which the value is one of the integers in the specified list.


`PARTITION_DISTRIBUTE`: Deprecated.  Implied by "`PARTITION`" option.

`PARTITION_BY_RANK`: Assume 1-1 mapping between MPI rank and part ID.  Assigning parts
to processors for which `rank == part id`.

`MPI_IO_RANK=<RANK>`: For IO modes in which a single processor handles all disk access, the MPI rank of the processor to use.  Default is **0**.

`PARALLEL_COMM=<id>`: Specify `moab::ParallelComm::id() == <id>` to use as communicator.

-----------------------
Format-specific options
-----------------------

Stereo Lithography (STL) files
------------------------------
  
`BINARY|ASCII`: Write binary or text STL file.  Default is text.

`BIG_ENDIAN|LITTLE_ENDIAN`: Force byte ordering of binary data.  Default is **BIG_ENDIAN** for writing and autodetection for reading (**BIG_ENDIAN** if autodetect fails).


MOAB Native (HDF5-based MHDF) format
------------------------------------

`ELEMENTS={EXPLICIT|NODES|SIDES}`: If reading only part of a file, specify which elements to read in addition to those contained in the specified set.  The options are:
>   - **EXPLICIT**    - read only explicitly designated elements
>   - **NODES**       - read any element for which all the nodes are being read.
>   - **SIDES**       - read explicilty specified elements and any elements that are sides of those elements.

Default is **SIDES** unless only nodes are explicitly specified, in which
case NODES will be used.

```
   CHILDREN = { NONE | SETS | CONTENTS }
   SETS     = { NONE | SETS | CONTENTS }
```
If reading only part of a file, specify whether or not child or contained sets (CHILDREN and SETS, respectively) of input sets are to be read. The options are:
>   + **NONE**     - Do not read sets because they are children of designated sets.
>   + **SETS**     - Read all child sets of designated input sets.
>   + **CONTENTS** - (Default).  Read all child sets and any entities contained in those sets.

`BUFFER_SIZE=<BYTES>`: Reduce buffer size for debugging purposes.

`KEEP`: For debugging purposes, retain partially written file if a failure occurs during the write process.

`BLOCKED_COORDINATE_IO={yes|no}`: During read of HDF5 file, read vertex coordinates in blocked format (read all X coordinates, followed by all Y coordinates, etc.)  Default is `'no'`.

`BCAST_SUMMARY={yes|no}`: During parallel read of HDF5 file, read file summary data on root process as serial IO and broadcast summary to other processes.  All processes then re-open file for parallel IO.  If 'no', file is opened only once by all processes for parallel IO and all processes read summary data.  Default is `'yes'`.

`BCAST_DUPLICATE_READS={yes|no}`: Work around MPI IO limitations when all processes read an indentical region of a file by reading the data only on the root process and using MPI_Bcast to send the data to other processes.  Default is `'yes'`.

`HYPERSLAB_OR`: During partial or parallel read, use documented (default) method for building HDF5 hyperslabs.  This option is deprecated.  Reader should automatically detect if `HYPERSLAB_APPEND` works and if not will default to `HYPERSLAB_OR`.  This option is retaind for debugging purposes only.

`HYPERSLAB_APPEND`: During partial or parallel read using a modified (i.e. hacked) HDF5 library, utilize hack for more efficient hyperslab selection construction. This option is deprecated.  Reader should automatically detect if `HYPERSLAB_APPEND` works and if not will default to `HYPERSLAB_OR`.

`HYPERSLAB_SELECT_LIMIT=n`: Set upper bound on the number of HDF5 hyperslabs that can be combined for a single read from an dataset.  The default is defined by `DEFAULT_HYPERSLAB_SELECT_LIMIT` in `ReadHDF5Dataset.cpp`.  If `HYPERSLAB_APPEND` is specified and this option is not, then the default is no limit.  This limit should be removed for future HDF5 releases (> 1.8.x) as said future releases will no longer require the `HYPERSLAB_APPEND` hack in order to avoid O(n^2) hyperslab selection behavior.

