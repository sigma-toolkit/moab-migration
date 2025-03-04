/*
 * NCHelperROMS.hpp
 *
 *  Created on: Mar. 4, 2025
 *      Author: iulian
 */

#ifndef NCHELPERROMS_HPP
#define NCHELPERROMS_HPP

#include "NCHelper.hpp"

namespace moab
{
#ifdef MOAB_HAVE_MPI
class ParallelComm;
#endif

class NCHelperROMS: public moab::ScdNCHelper {
public:
    NCHelperROMS( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
            : ScdNCHelper( readNC, fileId, opts, fileSet )
        {
        }
        static bool can_read_file( ReadNC* readNC, int fileId );

        ErrorCode create_mesh( Range& faces );

      private:
        virtual ErrorCode init_mesh_vals();
        virtual std::string get_mesh_type_name()
        {
            return "ROMS";
        }
        bool spherical;
        bool vertices_exist;
    };

}

#endif /* SRC_IO_NCHELPERROMS_HPP_ */
