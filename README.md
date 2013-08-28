### Scattered data interpolation with multilevel B-Splines

This is an implementation of the MBA algorithm from [1]. This is a fast
algorithm for scattered N-dimensional data interpolation and approximation.

Example of 2D interpolation:
~~~{.cpp}
#include <mba.hpp>

int main() {
    // Coordinates of data points.
    std::vector< std::array<double,2> > p = {
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.4, 0.4},
        {0.6, 0.6}
    };

    // Data values.
    std::vector<double> v = {
        0.2, 0.0, 0.0, -0.2, -1.0, 1.0
    };

    // Algorithm setup.
    mba::cloud<2> c(
            mba::make_array<double>(-0.01, -0.01),
            mba::make_array<double>( 1.01,  1.01),
            p, v, mba::make_array<size_t>(5, 5)
            );

    // Get interpolated value at arbitrary location.
    double w = c(0.3, 0.7);
}
~~~

### References

1. S. Lee, G. Wolberg, and S. Y. Shin. Scattered data interpolation with
   multilevel B-Splines. IEEE Transactions on Visualization and
   Computer Graphics, 3:228–244, 1997.

