#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include <mba.hpp>

int main(int argc, char *argv[]) {
    const size_t n = argc < 2 ? 1024 * 1024 : std::stoi(argv[1]);

    std::default_random_engine rng(0);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);

    std::vector< std::array<double,2> > p(n);
    std::generate(p.begin(), p.end(), [&]() { return mba::make_array<double>(rnd(rng), rnd(rng)); });

    std::vector<double> v(n);
    std::transform(p.begin(), p.end(), v.begin(), [](const std::array<double,2> &c) {
            double x = c[0] - 0.5;
            double y = c[1] - 0.5;
            return x * x + y * y;
            });

    std::array<double, 2> xmin = {-0.01, -0.01};
    std::array<double, 2> xmax = { 1.01,  1.01};

    mba::cloud<2> c(xmin, xmax, p, v, mba::default_grid(xmin, xmax));

    std::ofstream dat("c.dat");
    dat << std::scientific;
    for(double y = 0; y <= 1; y += 0.01) {
        for(double x = 0; x <= 1; x += 0.01)
            dat << c(x, y) << " ";
        dat << std::endl;
    }
}
