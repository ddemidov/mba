#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <functional>

#include <mba/mba.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void precondition(bool cond, std::string error_message) {
    if (!cond) throw std::runtime_error(error_message);
}

template <unsigned NDim>
struct python_mba {
    python_mba(
            py::array_t<double> lo,
            py::array_t<double> hi,
            py::array_t<int>    grid,
            py::array_t<double> coo,
            py::array_t<double> val,
            std::function<double(std::vector<double>)> f,
            int max_levels = 8,
            double tol = 1e-8,
            double min_fill = 0.5
            )
    {
        auto initial = [f](const std::array<double, NDim> &a) {
                return f(std::vector<double>(std::begin(a), std::end(a)));
            };

        init(lo, hi, grid, coo, val, initial, max_levels, tol, min_fill);
    }

    python_mba(
            py::array_t<double> lo,
            py::array_t<double> hi,
            py::array_t<int>    grid,
            py::array_t<double> coo,
            py::array_t<double> val,
            int max_levels = 8,
            double tol = 1e-8,
            double min_fill = 0.5
            )
    {
        std::function<double(const std::array<double,NDim>&)> none;
        init(lo, hi, grid, coo, val, none, max_levels, tol, min_fill);
    }

    py::array_t<double> apply(py::array_t<double> _coo) const {
        typedef std::array<double, NDim> point;

        py::buffer_info coo = _coo.request();

        const size_t ndim = coo.ndim;

        precondition(ndim >= 2 && coo.shape[ndim-1] == NDim,
                "coo should be a n x " + std::to_string(NDim) + " matrix"
                );

        std::vector<size_t> strides(ndim-1);
        std::vector<size_t> shape(ndim-1);

        shape.back() = coo.shape[ndim-2];
        strides.back() = sizeof(double);
        for(int i = ndim-2; i > 0; --i) {
            shape[i-1] = coo.shape[i-1];
            strides[i-1] = strides[i] * shape[i];
        }

        const size_t n = coo.size / NDim;

        const point *coo_begin = static_cast<const point*>(coo.ptr);
        const point *coo_end   = coo_begin + n;

        std::vector<double> val(n);

        std::transform(coo_begin, coo_end, val.begin(), std::ref(*m));

        return py::array(py::buffer_info(val.data(), sizeof(double),
                    py::format_descriptor<double>::value,
                    ndim-1, shape, strides));
    }

    std::shared_ptr< mba::MBA<NDim> > m;

    void init(
            py::array_t<double> _lo,
            py::array_t<double> _hi,
            py::array_t<int>    _grid,
            py::array_t<double> _coo,
            py::array_t<double> _val,
            std::function<double(const std::array<double,NDim>&)> initial,
            int max_levels = 8,
            double tol = 1e-8,
            double min_fill = 0.5
            )
    {
        typedef std::array<size_t, NDim> index;
        typedef std::array<double, NDim> point;

        py::buffer_info lo   = _lo.request();
        py::buffer_info hi   = _hi.request();
        py::buffer_info grid = _grid.request();
        py::buffer_info coo  = _coo.request();
        py::buffer_info val  = _val.request();

        precondition(lo.ndim == 1 && lo.shape[0] == NDim,
                "lo should be a vector of size " + std::to_string(NDim)
                );
        precondition(hi.ndim == 1 && hi.shape[0] == NDim,
                "hi should be a vector of size " + std::to_string(NDim)
                );
        precondition(grid.ndim == 1 && grid.shape[0] == NDim,
                "grid should be a vector of size " + std::to_string(NDim)
                );
        precondition(coo.ndim == 2 && coo.shape[1] == NDim,
                "coo should be a n x " + std::to_string(NDim) + " matrix"
                );
        precondition(val.ndim == 1 && val.shape[0] == coo.shape[0],
                "coo and val dimensions disagree"
                );

        const size_t n = coo.shape[0];

        const point *coo_begin = static_cast<const point*>(coo.ptr);
        const point *coo_end   = coo_begin + n;

        const double *val_begin = static_cast<const double*>(val.ptr);

        point cmin, cmax;
        index grid_size;

        std::copy_n(static_cast<const double*>(lo.ptr), NDim, std::begin(cmin));
        std::copy_n(static_cast<const double*>(hi.ptr), NDim, std::begin(cmax));
        std::copy_n(static_cast<const int*>(grid.ptr),  NDim, std::begin(grid_size));

        if (!initial) {
            initial = mba::linear_approximation<NDim>(coo_begin, coo_end, val_begin);
        }

        m = std::make_shared< mba::MBA<NDim> >(
                cmin, cmax, grid_size,
                coo_begin, coo_end, val_begin,
                max_levels, tol, min_fill, initial
                );
    }


};


template <unsigned NDim>
void register_mba(py::module &m) {
    std::string name = "mba" + std::to_string(NDim);
    std::string desc = "Multilevel B-Spline in " + std::to_string(NDim) + "D";

    py::class_< python_mba<NDim> >(m, name.c_str(), desc.c_str())
        .def(py::init<
                    py::array_t<double>,
                    py::array_t<double>,
                    py::array_t<int>,
                    py::array_t<double>,
                    py::array_t<double>,
                    int, double, double
                >(), "Constructor",
                py::arg("lo"),
                py::arg("hi"),
                py::arg("grid"),
                py::arg("coo"),
                py::arg("val"),
                py::arg("max_levels") = 8,
                py::arg("tol") = 1e-8,
                py::arg("min_fill") = 0.5
            )
        .def(py::init<
                    py::array_t<double>,
                    py::array_t<double>,
                    py::array_t<int>,
                    py::array_t<double>,
                    py::array_t<double>,
                    std::function<double(std::vector<double>)>,
                    int, double, double
                >(), "Constructor",
                py::arg("lo"),
                py::arg("hi"),
                py::arg("grid"),
                py::arg("coo"),
                py::arg("val"),
                py::arg("initial"),
                py::arg("max_levels") = 8,
                py::arg("tol") = 1e-8,
                py::arg("min_fill") = 0.5
            )
        .def("__call__", &python_mba<NDim>::apply)
        .def("__repr__", [](const python_mba<NDim> &m){
                std::ostringstream s;
                s << *m.m;
                return s.str();
                })
        ;

}

PYBIND11_MODULE(mba, m) {
    register_mba<1>(m);
    register_mba<2>(m);
    register_mba<3>(m);
}
