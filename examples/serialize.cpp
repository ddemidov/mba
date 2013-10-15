#include <iostream>
#include <fstream>
#include <mba.hpp>

void store(const char *fname) {
    std::vector< std::array<double,2> > p = {
        mba::make_array<double>(0.0, 0.0),
        mba::make_array<double>(0.0, 1.0),
        mba::make_array<double>(1.0, 0.0),
        mba::make_array<double>(1.0, 1.0),
        mba::make_array<double>(0.4, 0.4),
        mba::make_array<double>(0.6, 0.6)
    };

    std::vector<double> v = {
        0.2, 0.0, 0.0, -0.2, -1.0, 1.0
    };

    std::array<double, 2> xmin = {-0.01, -0.01};
    std::array<double, 2> xmax = { 1.01,  1.01};

    mba::cloud<2> c(
            xmin, xmax, p, v, mba::default_grid(xmin, xmax)
            );

    std::ofstream f(fname, std::ios::binary);
    c.save(f);
}

mba::cloud<2> load(const char *fname) {
    std::ifstream f(fname, std::ios::binary);
    mba::cloud<2> c(f);

    return c;
}

int main() {
    // Initialize and store the cloud.
    store("cloud.dat");

    // Load the cloud.
    mba::cloud<2> c = load("cloud.dat");

    // Use the cloud.
    std::ofstream dat("c.dat");
    dat << std::scientific;
    for(double y = 0; y <= 1; y += 0.01) {
        for(double x = 0; x <= 1; x += 0.01)
            dat << c(x, y) << " ";
        dat << std::endl;
    }
}
