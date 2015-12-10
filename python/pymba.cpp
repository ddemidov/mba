#include <memory>

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/stl_iterator.hpp>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy_boost_python.hpp"

#include <mba.hpp>

void precondition(bool cond, std::string error_message) {
    if (!cond) throw std::runtime_error(error_message);
}

template <unsigned NDIM>
struct python_mba {
    python_mba(
        const numpy_boost<double, 1> &lo,
        const numpy_boost<double, 1> &hi,
        const numpy_boost<double, 2> &coo,
        const numpy_boost<double, 1> &val,
        const numpy_boost<long,   1> &grd,
        int levels, double tol
        )
    {
        precondition(lo.shape()[0]  == NDIM, "wrong lo shape");
        precondition(hi.shape()[0]  == NDIM, "wrong lo shape");
        precondition(coo.shape()[1] == NDIM, "wrong coo shape");
        precondition(grd.shape()[0] == NDIM, "wrong grid shape");

        precondition(coo.shape()[0] == val.shape()[0], "inconsistent coo and val shapes");

        int n = coo.shape()[0];

        std::vector<std::array<double,NDIM>> p;
        std::vector<double> v;

        p.reserve(n);
        v.reserve(n);

        for (int i = 0; i < n; ++i) {
            std::array<double, NDIM> c;
            for (int j = 0; j < NDIM; ++j)
                c[j] = coo[i][j];
            p.push_back(c);
            v.push_back(val[i]);
        }

        std::array<double, NDIM> xmin;
        std::array<double, NDIM> xmax;
        std::array<size_t, NDIM> grid;

        for (int i = 0; i < NDIM; ++i) {
            xmin[i] = lo[i];
            xmax[i] = hi[i];
            grid[i] = grd[i];
        }

        cloud = std::make_shared<mba::cloud<NDIM>>(xmin, xmax, p, v, grid, levels, tol);
    }

    double get(const numpy_boost<double, 1> &p) const {
        std::array<double, NDIM> x;
        for (int i = 0; i < NDIM; ++i) x[i] = p[i];
        return (*cloud)(x);
    }

    std::shared_ptr<mba::cloud<NDIM>> cloud;
};

#if PY_MAJOR_VERSION >= 3
void*
#else
void
#endif
call_import_array() {
    import_array();
    return NUMPY_IMPORT_ARRAY_RETVAL;
}

//---------------------------------------------------------------------------
BOOST_PYTHON_MODULE(pymba)
{
    using namespace boost::python;
    docstring_options docopts(true, true, false);

    call_import_array();

    numpy_boost_python_register_type<double, 1>();
    numpy_boost_python_register_type<double, 2>();
    numpy_boost_python_register_type<long,   1>();

    class_< python_mba<1> >(
            "cloud1d", "Multilevel B-Spline in 1D",
            init<
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 2>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<long,   1>&,
                int,
                double
                >( args("lo", "hi", "coo", "val", "grid", "levels", "tol") )
            )
        .def("__call__", &python_mba<1>::get, "get approximation")
        ;

    class_< python_mba<2> >(
            "cloud2d", "Multilevel B-Spline in 2D",
            init<
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 2>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<long,   1>&,
                int,
                double
                >( args("lo", "hi", "coo", "val", "grid", "levels", "tol") )
            )
        .def("__call__", &python_mba<2>::get, "get approximation")
        ;

    class_< python_mba<3> >(
            "cloud3d", "Multilevel B-Spline in 3D",
            init<
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<double, 2>&,
                const numpy_boost<double, 1>&,
                const numpy_boost<long,   1>&,
                int,
                double
                >( args("lo", "hi", "coo", "val", "grid", "levels", "tol") )
            )
        .def("__call__", &python_mba<3>::get, "get approximation")
        ;
}
