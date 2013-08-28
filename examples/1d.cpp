#include <iostream>
#include <fstream>
#include <cmath>

#include <mba.hpp>

int main() {
    std::vector< std::array<double,1> > p;
    std::vector< double > v;

    const size_t n = 8;
    for(size_t i = 0; i < n; ++i) {
        double x = i * M_PI / (n - 1);
        p.push_back( mba::make_array<double>(x) );
        v.push_back( sin(x) );
    }

    mba::cloud<1> c(
            mba::make_array<double>(-0.01), mba::make_array<double>(M_PI + 0.01),
            p, v, mba::make_array<size_t>(2)
            );

    std::ofstream dat("c.dat");
    dat << std::scientific;
    for(double x = 0; x < M_PI; x += 0.01)
        dat << x << " " << c(x) << std::endl;
}
