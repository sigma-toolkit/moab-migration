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

#ifndef MB_UTIL_HPP
#define MB_UTIL_HPP

#include "moab/MOABConfig.h"
#include "moab/Forward.hpp"
#include "moab/CartVect.hpp"

#include <cmath>
#include <limits>
//#define moab_isfinite( f ) std::isfinite( f )
template<typename T>
inline bool moab_isfinite(T x) {
#if defined(_MSC_VER) && (_MSC_VER < 1800)
    // Old MSVC does not support std::isfinite
    return _finite(static_cast<double>(x)) != 0;
#else
    // Always use the standard C++11 version
    return std::isfinite(x);
#endif
}

namespace moab
{

/** \class Util
 *
 * \brief Utility functions for computational geometry and mathematical calculations
 */
class Util
{
  public:
    template < typename T >
    static bool is_finite( T value );

    static void normal( Interface* MB, EntityHandle handle, double& x, double& y, double& z );

    static void centroid( Interface* MB, EntityHandle handle, CartVect& coord );

    // static void edge_centers(Interface *MB, EntityHandle handle, std::vector<CartVect>
    // &coords_list);

    // static void face_centers(Interface *MB, EntityHandle handle, std::vector<CartVect>
    // &coords_list);

  private:
    Util() {}
};

template < typename T >
inline bool Util::is_finite( T value )
{
    return moab_isfinite( (double)value );
}

}  // namespace moab

#endif
