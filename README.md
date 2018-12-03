### Scattered data interpolation with multilevel B-Splines

[![Build Status](https://travis-ci.org/ddemidov/mba.svg?branch=master)](https://travis-ci.org/ddemidov/mba)

This library provides the adaptive MBA algorithm from [1] implemented in C++11.
This is a fast algorithm for scattered N-dimensional data interpolation and
approximation. Python bindings are also provided.

Example of 2D interpolation in C++:
```cpp
#include <mba.hpp>

int main() {
    // Coordinates of data points.
    std::vector<mba::point<2>> coo = {
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
    mba::point<2> lo = {-0.1, -0.1};
    mba::point<2> hi = { 1.1,  1.1};

    // Initial grid size.
    mba::index<2> grid = {3, 3};

    // Algorithm setup.
    mba::MBA<2> interp(lo, hi, grid, coo, val);

    // Get interpolated value at arbitrary location.
    double w = interp(mba::point<2>{0.3, 0.7});
}
```

Same example in python:
```.py
from mba import *

interp = mba2(lo=[-0.1,-0.1], hi=[1.1,1.1], grid=[3,3],
              coo=[ [0.0, 0.0], [0.0, 1.0], [1.0, 0.0],
                    [1.0, 1.0], [0.4, 0.4], [0.6, 0.6] ],
              val=[0.2, 0.0, 0.0, -0.2, -1.0, 1.0]
              )

w = interp([[0.3, 0.7]])
```

Also see
[python/example.ipynb](http://nbviewer.ipython.org/github/ddemidov/mba/blob/master/python/example.ipynb),
[python/layered.ipynb](http://nbviewer.ipython.org/github/ddemidov/mba/blob/master/python/layered.ipynb).

### References

1. S. Lee, G. Wolberg, and S. Y. Shin. Scattered data interpolation with
   multilevel B-Splines. IEEE Transactions on Visualization and
   Computer Graphics, 3:228–244, 1997, [doi:10.1109/2945.620490](https://doi.org/10.1109/2945.620490).

