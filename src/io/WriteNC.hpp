/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

#ifndef WRITENC_HPP_
#define WRITENC_HPP_

#ifndef IS_BUILDING_MB
#error "WriteNC.hpp isn't supposed to be included into an application"
#endif

#include <vector>
#include <map>
#include <set>
#include <string>

#include "moab/WriterIface.hpp"
#include "moab/ScdInterface.hpp"
#include "DebugOutput.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#endif

// Runtime dispatch — same rationale as ReadNC.hpp; see MBNcDispatch.hpp.
// The NCFUNC* macros below now route puts through mbnc_put_* wrappers
// (which select PNetCDF collective or standard NetCDF at runtime per the
// backend tag of the fileId returned by WriteNC's open/create path).
#include "MBNcDispatch.hpp"

//! Generic NC function dispatch
#define NCFUNC( func ) mbnc_##func

//! Collective I/O mode put (PNetCDF: collective; others: only mode)
#define NCFUNCAP( func ) mbnc_put##func

//! Independent I/O mode put. Routes to ncmpi_put_vara_* (no _all) on
//! PNetCDF; on non-PNetCDF backends maps to the same nc_put_vara_*
//! call as the collective variant. Caller is still responsible for
//! mbnc_begin_indep_data / mbnc_end_indep_data brackets when running
//! on PNetCDF (those are no-ops on other backends).
#define NCFUNCP( func ) mbnc_put##func##_indep

//! Nonblocking put (PNetCDF request aggregation; blocking on others)
#define NCFUNCREQP( func ) mbnc_iput##func

#define NCDF_SIZE size_t
#define NCDF_DIFF ptrdiff_t

namespace moab
{

class WriteUtilIface;
class NCWriteHelper;

/**
 * \brief Export NC files.
 */
class WriteNC : public WriterIface
{
    friend class NCWriteHelper;
    friend class ScdNCWriteHelper;
    friend class UcdNCWriteHelper;
    friend class NCWriteEuler;
    friend class NCWriteFV;
    friend class NCWriteHOMME;
    friend class NCWriteMPAS;
    friend class NCWriteGCRM;
    friend class NCWriteScrip;
    friend class NCWriteESMF;
    friend class NCWriteDomain;

  public:
    //! Factory method
    static WriterIface* factory( Interface* );

    //! Constructor
    WriteNC( Interface* impl = NULL );

    //! Destructor
    virtual ~WriteNC();

    //! Writes out a file
    ErrorCode write_file( const char* file_name,
                          const bool overwrite,
                          const FileOptions& opts,
                          const EntityHandle* output_list,
                          const int num_sets,
                          const std::vector< std::string >& qa_list,
                          const Tag* tag_list  = NULL,
                          int num_tags         = 0,
                          int export_dimension = 3 );

  private:
    //! ENTLOCNSEDGE for north/south edge
    //! ENTLOCWEEDGE for west/east edge
    enum EntityLocation
    {
        ENTLOCVERT = 0,
        ENTLOCNSEDGE,
        ENTLOCEWEDGE,
        ENTLOCFACE,
        ENTLOCSET,
        ENTLOCEDGE,
        ENTLOCREGION
    };

    class AttData
    {
      public:
        AttData() : attId( -1 ), attLen( 0 ), attVarId( -2 ), attDataType( NC_NAT ) {}
        int attId;
        NCDF_SIZE attLen;
        int attVarId;
        nc_type attDataType;
        std::string attValue;
    };

    class VarData
    {
      public:
        VarData() : varId( -1 ), numAtts( -1 ), entLoc( ENTLOCSET ), numLev( 0 ), sz( 0 ), has_tsteps( false ) {}
        int varId;
        int numAtts;
        nc_type varDataType;
        std::vector< int > varDims;  // The dimension indices making up this multi-dimensional variable
        std::map< std::string, AttData > varAtts;
        std::string varName;
        std::vector< Tag > varTags;            // Tags created for this variable, e.g. one tag per timestep
        std::vector< void* > memoryHogs;       // These will point to the real data; fill before writing the data
        std::vector< NCDF_SIZE > writeStarts;  // Starting index for writing data values along each dimension
        std::vector< NCDF_SIZE > writeCounts;  // Number of data values to be written along each dimension
        int entLoc;
        int numLev;
        int sz;
        bool has_tsteps;  // Indicate whether timestep numbers are appended to tag names
    };

    //! This info will be reconstructed from metadata stored on conventional fileSet tags
    //! Dimension names
    std::vector< std::string > dimNames;

    //! Dimension lengths
    std::vector< int > dimLens;

    //! Will collect used dimensions (coordinate variables)
    std::set< std::string > usedCoordinates;

    //! Dummy variables (for dimensions that have no corresponding coordinate variables)
    std::set< std::string > dummyVarNames;

    //! Global attribs
    std::map< std::string, AttData > globalAtts;

    //! Variable info
    std::map< std::string, VarData > varInfo;

    ErrorCode parse_options( const FileOptions& opts,
                             std::vector< std::string >& var_names,
                             std::vector< std::string >& desired_names,
                             std::vector< int >& tstep_nums,
                             std::vector< double >& tstep_vals );
    /*
     * Map out the header, from tags on file set; it is the inverse process from
     * ErrorCode NCHelper::create_conventional_tags
     */
    ErrorCode process_conventional_tags( EntityHandle fileSet );

    ErrorCode process_concatenated_attribute( const void* attPtr,
                                              int attSz,
                                              std::vector< int >& attLen,
                                              std::map< std::string, AttData >& attributes );

    //! Interface instance
    Interface* mbImpl;
    WriteUtilIface* mWriteIface;

    //! File var
    const char* fileName;

    //! File numbers assigned by (p)netcdf
    int fileId;

    //! Debug stuff
    DebugOutput dbgOut;

#ifdef MOAB_HAVE_MPI
    ParallelComm* myPcomm;
#endif

    //! Write options
    bool noMesh;
    bool noVars;
    bool append;

    //! Cached tags for writing. This will be important for ordering the data, in parallel
    Tag mGlobalIdTag;

    //! Are we writing in parallel? (probably in the future)
    bool isParallel;

    //! CAM Euler, etc,
    std::string grid_type;

    //! Helper class instance
    NCWriteHelper* myHelper;
};

}  // namespace moab

#endif
