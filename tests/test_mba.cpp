#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <mba/mba.hpp>

TEST_CASE( "Grid iterator" ) {
    mba::index<2> dim = {2, 3};
    mba::detail::grid_iterator<2> g(dim);

    for(size_t j = 0; j < 2; ++j) {
        for(size_t i = 0; i < 3; ++i, ++g) {
            REQUIRE(static_cast<bool>(g));
            REQUIRE(g[0] == j);
            REQUIRE(g[1] == i);
        }
    }

    REQUIRE(!static_cast<bool>(g));
}

TEST_CASE( "Control lattice" ) {
    mba::point<2> lo = {-0.1, -0.1};
    mba::point<2> hi = { 1.1,  1.1};

    mba::index<2> grid = {9, 9};

    std::array<mba::point<2>, 2> coo = {0.0, 0.0, 1.0, 1.0};
    std::array<double, 2> val = {1.0, 0.0};

    SECTION("Dense") {
        mba::detail::control_lattice_dense<2> phi(
                lo, hi, grid,
                std::begin(coo), std::end(coo), std::begin(val)
                );

        REQUIRE(val[0] == phi(coo[0]));
        REQUIRE(val[1] == phi(coo[1]));
    }

    SECTION("Sparse") {
        mba::detail::control_lattice_sparse<2> phi(
                lo, hi, grid,
                std::begin(coo), std::end(coo), std::begin(val)
                );

        REQUIRE(val[0] == phi(coo[0]));
        REQUIRE(val[1] == phi(coo[1]));
    }
}

TEST_CASE( "MBA" ) {
    mba::point<2> lo = {-0.1, -0.1};
    mba::point<2> hi = { 1.1,  1.1};

    mba::index<2> grid = {2, 2};

    SECTION( "few points" ) {
        const int n = 2;
        std::array<mba::point<2>, n> coo = {0, 0, 1, 1};
        std::array<double, n> val = {1.0, 0.0};

        mba::MBA<2> phi(lo, hi, grid, coo, val, 8, 1e-8);

        for(size_t i = 0; i < n; ++i) {
            REQUIRE(std::abs(val[i] - phi(coo[i])) < 1e-8);
        }
    }

    SECTION( "enough points" ) {
        const int n = 4;
        std::array<mba::point<2>, n> coo = {0, 0, 0, 1, 1, 0, 1, 1};
        std::array<double, n> val = {1.0, 0.5, 0.75, 0.0};

        mba::MBA<2> phi(lo, hi, grid, coo, val, 8, 1e-8);

        for(size_t i = 0; i < n; ++i) {
            REQUIRE(std::abs(val[i] - phi(coo[i])) < 1e-8);
        }
    }

    SECTION("small test function") {
        std::vector<mba::point<2>> coo;
        std::vector<double> val;
        for (unsigned int y = 0; y < 10; y++) {
            for (unsigned int x = 0; x < 10; x++) {
                if (x % 2 == 0 && y % 2 == 0) {
                    const double x_coord = x / 10.;
                    const double y_coord = y / 10.;
                    const double current_value = sin(x * 8.) * sin(y * 8.);
                    coo.push_back({x_coord, y_coord});
                    val.push_back(current_value);
                }
            }
        }

        const mba::MBA<2> phi(lo, hi, grid, coo, val, 8, 1e-8);

        for (unsigned int y = 0; y < 10; y++) {
            for (unsigned int x = 0; x < 10; x++) {
                const double x_coord = x / 10.;
                const double y_coord = y / 10.;
                const double true_value = sin(x * 8.) * sin(y * 8.);
                const double inter_value = phi({x_coord, y_coord});
                REQUIRE(std::abs(true_value - inter_value) < 1.5);
            }
        }
    }
}
