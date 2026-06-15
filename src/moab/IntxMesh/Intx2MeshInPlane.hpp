/*
 * Intx2MeshInPlane.hpp
 *
 *  Created on: Oct 24, 2012
 *      Author: iulian
 */

#ifndef INTX2MESHINPLANE_HPP_
#define INTX2MESHINPLANE_HPP_

#include "Intx2Mesh.hpp"
namespace moab
{

class Intx2MeshInPlane : public moab::Intx2Mesh
{
  public:
    Intx2MeshInPlane( Interface* mbimpl );

    virtual ~Intx2MeshInPlane() override;

    double setup_tgt_cell( EntityHandle tgt, int& nsTgt ) override;

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
};

}  // end namespace moab
#endif /* INTX2MESHINPLANE_HPP_ */
