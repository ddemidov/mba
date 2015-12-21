#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <mba/mba.hpp>

TEST_CASE( "Grid iterator" ) {
    boost::array<size_t, 2> dim = {2, 3};
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
    boost::array<double, 2> lo = {-0.1, -0.1};
    boost::array<double, 2> hi = { 1.1,  1.1};

    boost::array<size_t, 2> grid = {9, 9};

    boost::array<boost::array<double, 2>, 2> coo = {0.0, 0.0, 1.0, 1.0};
    boost::array<double, 2> val = {1.0, 0.0};

    SECTION("Dense") {
        mba::detail::control_lattice_dense<2> phi(
                lo, hi, grid,
                boost::begin(coo), boost::end(coo), boost::begin(val)
                );

        REQUIRE(val[0] == phi(coo[0]));
        REQUIRE(val[1] == phi(coo[1]));
    }

    SECTION("Sparse") {
        mba::detail::control_lattice_sparse<2> phi(
                lo, hi, grid,
                boost::begin(coo), boost::end(coo), boost::begin(val)
                );

        REQUIRE(val[0] == phi(coo[0]));
        REQUIRE(val[1] == phi(coo[1]));
    }
}

TEST_CASE( "MBA" ) {
    boost::array<double, 2> lo = {-0.1, -0.1};
    boost::array<double, 2> hi = { 1.1,  1.1};

    boost::array<size_t, 2> grid = {2, 2};

    boost::array<boost::array<double, 2>, 2> coo = {0.0, 0.0, 1.0, 1.0};
    boost::array<double, 2> val = {1.0, 0.0};

    mba::MBA<2> phi(lo, hi, grid, coo, val, 8, 1e-8);

    REQUIRE(std::abs(val[0] - phi(coo[0])) < 1e-8);
    REQUIRE(std::abs(val[1] - phi(coo[1])) < 1e-8);
}
