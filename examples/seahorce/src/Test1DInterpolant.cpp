/*
 * =====================================================================================
 *
 *       Filename:  test_pchip.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  07/17/2023 18:56:49
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *        Company:  Argonne National Lab
 *
 * =====================================================================================
 */

#include <iostream>
#include <vector>
#include <cmath>

#include "spline.h"
#include "HermiteCubicCurve.hpp"
#include "ComputeMBA.hpp"
#include "moab/Remapping/mlinterp.hpp"

typedef std::vector< double > Vector;

int main()
{
    // Example usage
    Vector xA    = { 1.0, 2.2, 3.7, 4.4, 5.6 };
    Vector yA    = { 2.0, 1.0, 3.0, 5.0, 4.0 };
    size_t nData = xA.size();

    Vector queryPoints = { 1.6, 5.0, 3.0, 0.2, 7.0, 2.2001 };
    size_t nInterp     = queryPoints.size();
    Vector queryMLIValue( nInterp );

    // Perform linear interpolation first
    mlinterp::interp( &nData, nInterp,                  // Number of points
                      yA.data(), queryMLIValue.data(),  // Output axis (y)
                      xA.data(), queryPoints.data()     // Input axis (x)
    );

    Vector queryTKC1Value( nInterp );
    tk::spline splSC1( xA, yA, tk::spline::cspline_hermite, true );
    for( size_t index = 0; index < nInterp; index++ )
        queryTKC1Value[index] = splSC1( queryPoints[index] );

    Vector queryTKC2Value( nInterp );
    tk::spline splSC2( xA, yA, tk::spline::cspline, false, tk::spline::first_deriv, 0.0, tk::spline::first_deriv, 0.0 );
    for( size_t index = 0; index < nInterp; index++ )
        queryTKC2Value[index] = splSC2( queryPoints[index] );

    Vector queryHCCValue( nInterp );
    HermiteCubicCurve< double > splHC_S;
    for( size_t index = 0; index < nData; index++ )
        splHC_S.add( xA[index], yA[index] );
    splHC_S.finish();
    for( size_t index = 0; index < nInterp; index++ )
        queryHCCValue[index] = splHC_S.at( queryPoints[index] );

    std::cout << "  Location     Value\n";
    std::cout << "-----------------------\n";
    for( size_t index = 0; index < nData; index++ )
        std::cout << "    " << xA[index] << "\t\t" << yA[index] << std::endl;
    std::cout << "-----------------------\n";
    for( size_t index = 0; index < nInterp; index++ )
        std::cout << "Location: " << queryPoints[index] << ", Linear: " << queryMLIValue[index]
                  << ", TKSplineC1: " << queryTKC1Value[index] << ", TKSplineC2: " << queryTKC2Value[index]
                  << ", HermiteCubicCurve: " << queryHCCValue[index] << std::endl;
}
