### Scattered data interpolation with multilevel B-Splines

This library provides the adaptive MBA algorithm from [1] implemented in C++.
This is a fast algorithm for scattered N-dimensional data interpolation and
approximation. Python bindings are also provided.

Example of 2D interpolation in C++:
~~~{.cpp}
#include <mba.hpp>

typedef boost::array<double, 2> point2d;

int main() {
    // Coordinates of data points.
    std::vector<point2d> coo = {
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.4, 0.4},
        {0.6, 0.6}
    };

    // Data values.
    std::vector<double> val = {
        0.2, 0.0, 0.0, -0.2, -1.0, 1.0
    };

    // Bounding box containing the data points.
    point2d lo = {-0.1, -0.1};
    point2d hi = { 1.1,  1.1};

    // Initial grid size.
    boost::array<size_t, 2> grid = {3, 3};

    // Algorithm setup.
    mba::MBA<2> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    point2d p = {0.3, 0.7};
    double w = interp(p);
}
~~~

Same example in python:
~~~{.py}
from mba import *

interp = mba2(lo=[-0.1,-0.1], hi=[1.1,1.1], grid=[3,3],
              coo=[ [0.0, 0.0], [0.0, 1.0], [1.0, 0.0],
                    [1.0, 1.0], [0.4, 0.4], [0.6, 0.6] ],
              val=[0.2, 0.0, 0.0, -0.2, -1.0, 1.0]
              )

w = interp([[0.3, 0.7]])
~~~

### References

1. S. Lee, G. Wolberg, and S. Y. Shin. Scattered data interpolation with
   multilevel B-Splines. IEEE Transactions on Visualization and
   Computer Graphics, 3:228–244, 1997.

