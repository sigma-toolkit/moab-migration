#include <iostream>
#include <vector>
#include <cmath>

#include "PCHIP.hpp"

int main()
{
    // Example usage
    Vector x = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    Vector y = { 2.0, 1.0, 3.0, 5.0, 4.0 };

    double query1 = 2.5, query2 = 7.0;

    // Perform PCHIP interpolation for interior nodes
    double interpolatedValue1 = pchipInterpolate( x, y, query1 );
    std::cout << "PCHIP Interpolated value at " << query1 << ": " << interpolatedValue1 << std::endl;

    // Perform PCHIP interpolation with extrapolation based on nearest-neighbor for extrapolations
    double interpolatedValue2 = pchipInterpolate( x, y, query2 );
    std::cout << "PCHIP Interpolated value at " << query2 << ": " << interpolatedValue2 << std::endl;

    // Perform monotone linear interpolation with extrapolation based on nearest-neighbor
    double linearInterpolatedValue1 = linearInterpolate( x, y, query1 );
    double linearInterpolatedValue2 = linearInterpolate( x, y, query2 );

    std::cout << "Linear interpolated value at " << query1 << ": " << linearInterpolatedValue1 << std::endl;
    std::cout << "Linear interpolated value at " << query2 << ": " << linearInterpolatedValue2 << std::endl;

    return 0;
}
