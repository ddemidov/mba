#include <iostream>
#include <fstream>

#include <mba.hpp>

int main() {
    std::vector< std::array<double,2> > p = {
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.4, 0.4},
        {0.6, 0.6}
    };

    std::vector<double> v = {
        0.2, 0.0, 0.0, -0.2, -1.0, 1.0
    };

    mba::cloud<2> c(
            mba::make_array<double>(-0.01, -0.01),
            mba::make_array<double>( 1.01,  1.01),
            p, v, mba::make_array<size_t>(5, 5)
            );

    std::ofstream dat("c.dat");
    dat << std::scientific;
    for(double y = 0; y <= 1; y += 0.01) {
        for(double x = 0; x <= 1; x += 0.01)
            dat << c(x, y) << " ";
        dat << std::endl;
    }
}
