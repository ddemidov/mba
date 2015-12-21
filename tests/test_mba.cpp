#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <mba/mba.hpp>

TEST_CASE( "Extent generator" ) {
    REQUIRE(mba::extents[6][7].dims[0]     == 6);
    REQUIRE(mba::extents[6][7].dims[1]     == 7);
    REQUIRE(mba::extents[6][7].dims.size() == 2);
}

TEST_CASE( "Grid iterator" ) {
    mba::detail::grid_iterator<2> g(mba::extents[2][3]);

    for(size_t j = 0; j < 2; ++j) {
        for(size_t i = 0; i < 3; ++i, ++g) {
            REQUIRE(static_cast<bool>(g));
            REQUIRE(g[0] == j);
            REQUIRE(g[1] == i);
        }
    }

    REQUIRE(!static_cast<bool>(g));
}

TEST_CASE( "Array operations" ) {
    using namespace mba::detail;

    boost::array<double, 4> a = {2, 4, 6, 8};
    boost::array<double, 4> b = {1, 2, 3, 4};

    SECTION( "addition" ) {
        boost::array<double, 4> c = a   + b;
        boost::array<double, 4> d = a   + 2.0;
        boost::array<double, 4> e = 2.0 + a;

        for(int i = 0; i < 4; ++i) {
            REQUIRE(c[i] == a[i] + b[i]);
            REQUIRE(d[i] == a[i] + 2.0);
            REQUIRE(e[i] == 2.0  + a[i]);
        }
    }

    SECTION( "subtraction" ) {
        boost::array<double, 4> c = a   - b;
        boost::array<double, 4> d = a   - 2.0;
        boost::array<double, 4> e = 2.0 - a;

        for(int i = 0; i < 4; ++i) {
            REQUIRE(c[i] == a[i] - b[i]);
            REQUIRE(d[i] == a[i] - 2.0);
            REQUIRE(e[i] == 2.0  - a[i]);
        }
    }

    SECTION( "multiplication" ) {
        boost::array<double, 4> c = a   * b;
        boost::array<double, 4> d = a   * 2.0;
        boost::array<double, 4> e = 2.0 * a;

        for(int i = 0; i < 4; ++i) {
            REQUIRE(c[i] == a[i] * b[i]);
            REQUIRE(d[i] == a[i] * 2.0);
            REQUIRE(e[i] == 2.0  * a[i]);
        }
    }

    SECTION( "division" ) {
        boost::array<double, 4> c = a   / b;
        boost::array<double, 4> d = a   / 2.0;
        boost::array<double, 4> e = 2.0 / a;

        for(int i = 0; i < 4; ++i) {
            REQUIRE(c[i] == a[i] / b[i]);
            REQUIRE(d[i] == a[i] / 2.0);
            REQUIRE(e[i] == 2.0  / a[i]);
        }
    }
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
