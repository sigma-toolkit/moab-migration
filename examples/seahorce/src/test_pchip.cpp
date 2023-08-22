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

typedef std::vector<double> Vector;

// Helper function for finding the index of the nearest neighbor
int findNearestNeighbor(const Vector& x, double query)
{
    int n = x.size();
    int idx = 0;
    double minDist = std::abs(x[0] - query);

    for (int i = 1; i < n; ++i)
    {
        double dist = std::abs(x[i] - query);
        if (dist < minDist)
        {
            minDist = dist;
            idx = i;
        }
    }

    return idx;
}

// Helper function for computing the slope using finite differences
double computeSlope(double x1, double y1, double x2, double y2)
{
    return (y2 - y1) / (x2 - x1);
}

// Function to perform Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) interpolation with extrapolation based on nearest-neighbor
double pchipInterpolate(const Vector& x, const Vector& y, double query)
{
    int n = x.size();

    // Handle cases of extrapolation
    if (query <= x[0])
        return y[0];
    else if (query >= x[n - 1])
        return y[n - 1];

    // Find the interval containing the query point
    int idx = findNearestNeighbor(x, query);

    // Compute slopes for the neighboring points
    double h0, h1, m0, m1;
    if (idx > 0 && idx < n - 1)
    {
        h0 = x[idx] - x[idx - 1];
        h1 = x[idx + 1] - x[idx];
        m0 = computeSlope(x[idx - 1], y[idx - 1], x[idx], y[idx]);
        m1 = computeSlope(x[idx], y[idx], x[idx + 1], y[idx + 1]);
    }
    else if (idx == 0)
    {
        h0 = x[1] - x[0];
        h1 = x[2] - x[1];
        m0 = computeSlope(x[0], y[0], x[1], y[1]);
        m1 = computeSlope(x[1], y[1], x[2], y[2]);
    }
    else
    {
        h0 = x[n - 2] - x[n - 3];
        h1 = x[n - 1] - x[n - 2];
        m0 = computeSlope(x[n - 3], y[n - 3], x[n - 2], y[n - 2]);
        m1 = computeSlope(x[n - 2], y[n - 2], x[n - 1], y[n - 1]);
    }

    // Compute the interpolation polynomial coefficients
    double t = (query - x[idx]) / h1;
    double t2 = t * t;
    double t3 = t * t2;
    double c0 = 2 * t3 - 3 * t2 + 1;
    double c1 = t3 - 2 * t2 + t;
    double c2 = -2 * t3 + 3 * t2;
    double c3 = t3 - t2;

    // Evaluate the interpolated value
    double interpolatedValue = c0 * y[idx] + c1 * m1 * h1 + c2 * y[idx + 1] + c3 * m0 * h1;

    return interpolatedValue;
}

int main()
{
    // Example usage
    Vector x = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector y = {2.0, 1.0, 3.0, 5.0, 4.0};

    double query = 1.5;

    // Perform PCHIP interpolation with extrapolation based on nearest-neighbor
    double interpolatedValue = pchipInterpolate(x, y, query);

    std::cout << "Interpolated value at " << query << ": " << interpolatedValue << std::endl;

    return 0;
}

