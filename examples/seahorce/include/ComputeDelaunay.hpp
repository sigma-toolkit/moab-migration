#ifndef __compute_delaunay_hpp__
#define __compute_delaunay_hpp__

#include "RemapMPASROMS.hpp"


// moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& xyzd,
//                                             std::vector< double >& fd,
//                                             std::vector< double >& xyzi,
//                                             std::vector< double >& fi );


moab::ErrorCode ComputeDelaunayInterpolant( std::vector< double >& /*xyzd*/,
                                            std::vector< double >& /*fd*/,
                                            std::vector< double >& /*xyzi*/,
                                            std::vector< double >& /*fi*/ )
{
    // int nd = fd.size();
    // int ni = fi.size();

    return moab::MB_SUCCESS;
}


#endif  // __compute_delaunay_hpp__