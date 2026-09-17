/*
 * Intx2MeshOnSphere.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */

#ifndef INTX2MESHONSPHERE_HPP_
#define INTX2MESHONSPHERE_HPP_

#include "Intx2Mesh.hpp"

namespace moab
{

class Intx2MeshOnSphere : public moab::Intx2Mesh
{
  public:
    Intx2MeshOnSphere( Interface* mbimpl, IntxAreaUtils::AreaMethod amethod = IntxAreaUtils::DEFAULT_AREA_METHOD );

    virtual ~Intx2MeshOnSphere() override;

    void set_radius_source_mesh( double radius )
    {
        Rsrc = radius;
    }
    void set_radius_destination_mesh( double radius )
    {
        Rdest = radius;
    }

    double setup_tgt_cell( EntityHandle tgt, int& nsTgt ) override;

    // main method to intersect meshes on a sphere

    ErrorCode computeIntersectionBetweenTgtAndSrc( EntityHandle tgt,
                                                   EntityHandle src,
                                                   double* P,
                                                   int& nP,
                                                   double& area,
                                                   int markb[MAXEDGES],
                                                   int markr[MAXEDGES],
                                                   int& nsSrc,
                                                   int& nsTgt,
                                                   bool check_boxes_first = false ) override;

    ErrorCode findNodes( EntityHandle tgt, int nsTgt, EntityHandle src, int nsSrc, double* iP, int nP ) override;

    ErrorCode update_tracer_data( EntityHandle out_set, Tag& tagElem, Tag& tagArea );
#ifdef MOAB_HAVE_MPI
    virtual ErrorCode construct_covering_set( EntityHandle& initial_distributed_set,
                                              EntityHandle& covering_set,
                                              bool gnomonic       = true,
                                              int nb_ghost_layers = 0,
                                              bool include_edges  = false );

    virtual ErrorCode build_processor_euler_boxes( EntityHandle euler_set, Range& local_verts, bool gnomonic = true ) override;
#endif

    const IntxAreaUtils::AreaMethod areaMethod;

  private:
    int plane;           // current gnomonic plane
    double Rsrc, Rdest;  // radius of the sphere
};

} /* namespace moab */
#endif /* INTX2MESHONSPHERE_HPP_ */
