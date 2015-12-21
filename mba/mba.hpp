#ifndef MBA_MBA_HPP
#define MBA_MBA_HPP

/*
The MIT License

Copyright (c) 2015 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   mba/mba.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Multilevel B-spline interpolation.
 */

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/type_traits.hpp>
#include <boost/scoped_ptr.hpp>

namespace mba {
namespace detail {

/// Extent generator
template <unsigned NDim>
struct extent_gen {
    boost::array<size_t, NDim> dims;

    extent_gen<NDim+1> operator[](size_t n) const {
        extent_gen<NDim+1> e;
        boost::copy(dims, boost::begin(e.dims));
        e.dims[NDim] = n;
        return e;
    }
};

} // namespace detail

const detail::extent_gen<0> extents;

namespace detail {

template <size_t N, size_t M>
struct power : boost::integral_constant<size_t, N * power<N, M-1>::value> {};

template <size_t N>
struct power<N, 0> : boost::integral_constant<size_t, 1> {};

/// N-dimensional grid iterator (nested loop with variable depth).
template <unsigned NDim>
class grid_iterator {
    public:
        typedef boost::array<size_t, NDim> index;

        explicit grid_iterator(const extent_gen<NDim> &ext)
            : N(ext.dims)
        {
            boost::fill(i, 0);
            done = (i == N);
        }

        explicit grid_iterator(const boost::array<size_t, NDim> &dims)
            : N(dims)
        {
            boost::fill(i, 0);
            done = (i == N);
        }

        explicit grid_iterator(size_t dim) {
            boost::fill(N, dim);
            boost::fill(i, 0);
            done = (0 == dim);
        }

        size_t operator[](size_t d) const {
            return i[d];
        }

        const index& operator*() const {
            return i;
        }

        grid_iterator& operator++() {
            done = true;
            for(size_t d = NDim; d--; ) {
                if (++i[d] < N[d]) {
                    done = false;
                    break;
                }
                i[d] = 0;
            }
        }

        operator bool() const { return !done; }

    private:
        index N, i;
        bool  done;
};

template <typename T, size_t N>
boost::array<T, N> operator+(boost::array<T, N> a, const boost::array<T, N> &b) {
    boost::transform(a, b, boost::begin(a), std::plus<T>());
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator+(T a, boost::array<T, N> b) {
    boost::transform(b, boost::begin(b), std::bind1st(std::plus<T>(), a));
    return b;
}

template <typename T, size_t N>
boost::array<T, N> operator+(boost::array<T, N> a, T b) {
    boost::transform(a, boost::begin(a), std::bind2nd(std::plus<T>(), b));
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator-(boost::array<T, N> a, const boost::array<T, N> &b) {
    boost::transform(a, b, boost::begin(a), std::minus<T>());
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator-(T a, boost::array<T, N> b) {
    boost::transform(b, boost::begin(b), std::bind1st(std::minus<T>(), a));
    return b;
}

template <typename T, size_t N>
boost::array<T, N> operator-(boost::array<T, N> a, T b) {
    boost::transform(a, boost::begin(a), std::bind2nd(std::minus<T>(), b));
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator*(boost::array<T, N> a, const boost::array<T, N> &b) {
    boost::transform(a, b, boost::begin(a), std::multiplies<T>());
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator*(T a, boost::array<T, N> b) {
    boost::transform(b, boost::begin(b), std::bind1st(std::multiplies<T>(), a));
    return b;
}

template <typename T, size_t N>
boost::array<T, N> operator*(boost::array<T, N> a, T b) {
    boost::transform(a, boost::begin(a), std::bind2nd(std::multiplies<T>(), b));
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator/(boost::array<T, N> a, const boost::array<T, N> &b) {
    boost::transform(a, b, boost::begin(a), std::divides<T>());
    return a;
}

template <typename T, size_t N>
boost::array<T, N> operator/(T a, boost::array<T, N> b) {
    boost::transform(b, boost::begin(b), std::bind1st(std::divides<T>(), a));
    return b;
}

template <typename T, size_t N>
boost::array<T, N> operator/(boost::array<T, N> a, T b) {
    boost::transform(a, boost::begin(a), std::bind2nd(std::divides<T>(), b));
    return a;
}

template <unsigned NDim>
class control_lattice {
    public:
        typedef boost::array<size_t, NDim> index;
        typedef boost::array<double, NDim> point;

        template <class CooIter, class ValIter>
        control_lattice(
                const point &coo_min, const point &coo_max, index grid_size,
                CooIter coo_begin, CooIter coo_end, ValIter val_begin
                ) : cmin(coo_min), cmax(coo_max), grid(grid_size)
        {
            for(unsigned i = 0; i < NDim; ++i) {
                hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
                cmin[i] -= 1 / hinv[i];
                grid[i] += 2;
            }

            boost::multi_array<double, NDim> delta(grid);
            boost::multi_array<double, NDim> omega(grid);

            std::fill(delta.data(), delta.data() + delta.num_elements(), 0.0);
            std::fill(omega.data(), omega.data() + omega.num_elements(), 0.0);

            CooIter p = coo_begin;
            ValIter v = val_begin;

            for(; p != coo_end; ++p, ++v) {
                if (!boxed(coo_min, *p, coo_max)) continue;

                index i;
                point s;

                for(unsigned d = 0; d < NDim; ++d) {
                    double u = ((*p)[d] - cmin[d]) * hinv[d];
                    i[d] = floor(u) - 1;
                    s[d] = u - floor(u);
                }

                boost::array< double, power<4, NDim>::value > w;
                double sum_w2 = 0.0;

                {
                    size_t idx = 0;
                    for(grid_iterator<NDim> d(4); d; ++d, ++idx) {
                        double prod = 1.0;
                        for(unsigned k = 0; k < NDim; ++k) prod *= B(d[k], s[k]);

                        w[idx] = prod;
                        sum_w2 += prod * prod;
                    }
                }

                {
                    size_t idx = 0;
                    for(grid_iterator<NDim> d(4); d; ++d, ++idx) {
                        double w1  = w[idx];
                        double w2  = w1 * w1;
                        double phi = (*v) * w1 / sum_w2;

                        index j = i + (*d);

                        delta(j) += w2 * phi;
                        omega(j) += w2;
                    }
                }
            }

            phi.resize(grid);

            std::transform(
                    delta.data(), delta.data() + delta.num_elements(),
                    omega.data(), phi.data(), safe_divide
                    );
        }

        double operator()(const point &p) const {
            index i;
            point s;

            for(unsigned d = 0; d < NDim; ++d) {
                double u = (p[d] - cmin[d]) * hinv[d];
                i[d] = floor(u) - 1;
                s[d] = u - floor(u);
            }

            double f = 0;

            for(grid_iterator<NDim> d(4); d; ++d) {
                double w = 1.0;
                for(unsigned k = 0; k < NDim; ++k) w *= B(d[k], s[k]);

                f += w * phi(i + (*d));
            }

            return f;
        }

        template <class CooIter, class ValIter>
        double residual(CooIter coo_begin, CooIter coo_end, ValIter val_begin) const {
            double res = 0.0;

            CooIter p = coo_begin;
            ValIter v = val_begin;

            for(; p != coo_end; ++p, ++v) {
                (*v) -= (*this)(*p);
                res += (*v) * (*v);
            }

            return res;
        }

        void append_refined(const control_lattice &r) {
            static const boost::array<double, 5> s = {
                0.125, 0.500, 0.750, 0.500, 0.125
            };

            for(grid_iterator<NDim> i(r.grid); i; ++i) {
                double f = r.phi(*i);

                for(grid_iterator<NDim> d(5); d; ++d) {
                    index j;
                    bool skip = false;
                    for(unsigned k = 0; k < NDim; ++k) {
                        j[k] = 2 * i[k] + d[k] - 3;
                        if (j[k] >= grid[k]) {
                            skip = true;
                            break;
                        }
                    }

                    if (skip) continue;

                    double c = 1.0;
                    for(unsigned k = 0; k < NDim; ++k) c *= s[d[k]];

                    phi(j) += f * c;
                }
            }
        }

    private:
        point cmin, cmax, hinv;
        index grid;

        boost::multi_array<double, NDim> phi;

        // Value of k-th B-Spline at t.
        static double B(size_t k, double t) {
            assert(0 <= t && t < 1);
            assert(k < 4);

            switch (k) {
                case 0:
                    return (t * (t * (-t + 3) - 3) + 1) / 6;
                case 1:
                    return (t * t * (3 * t - 6) + 4) / 6;
                case 2:
                    return (t * (t * (-3 * t + 3) + 3) + 1) / 6;
                case 3:
                    return t * t * t / 6;
                default:
                    return 0;
            }
        }

        static bool boxed(const point &lo, const point &p, const point &hi) {
            for(unsigned i = 0; i < NDim; ++i) {
                if (p[i] < lo[i] || p[i] > hi[i]) return false;
            }
            return true;
        }

        static double safe_divide(double a, double b) {
            return b == 0.0 ? 0.0 : a / b;
        }
};

} // namespace detail

template <unsigned NDim>
class MBA {
    public:
        typedef boost::array<size_t, NDim> index;
        typedef boost::array<double, NDim> point;

        template <class CooIter, class ValIter>
        MBA(
                const point &coo_min, const point &coo_max, index grid,
                CooIter coo_begin, CooIter coo_end, ValIter val_begin,
                unsigned max_levels = 8, double tol = 1e-8
           )
        {
            init(coo_min, coo_max, grid, coo_begin, coo_end, val_begin, max_levels, tol);
        }

        template <class CooRange, class ValRange>
        MBA(
                const point &coo_min, const point &coo_max, index grid,
                CooRange coo, ValRange val,
                unsigned max_levels = 8, double tol = 1e-8
           )
        {
            init(coo_min, coo_max, grid, boost::begin(coo), boost::end(coo), boost::begin(val), max_levels, tol);
        }

        double operator()(const point &p) const {
            return (*psi)(p);
        }
    private:
        typedef detail::control_lattice<NDim> control_lattice;
        boost::scoped_ptr<control_lattice> psi;

        template <class CooIter, class ValIter>
        void init(
                const point &cmin, const point &cmax, index grid,
                CooIter coo_begin, CooIter coo_end, ValIter val_begin,
                unsigned max_levels, double tol
                )
        {
            using namespace mba::detail;

            const ptrdiff_t n = std::distance(coo_begin, coo_end);
            std::vector<double> val(val_begin, val_begin + n);

            psi.reset(new control_lattice(cmin, cmax, grid, coo_begin, coo_end, val.begin()));

            double eps = tol * boost::inner_product(val, val, 0.0);
            double res = psi->residual(coo_begin, coo_end, val.begin());

            for(size_t k = 1; (k < max_levels) && (res > eps); ++k) {
                grid = 2ul * grid - 1ul;

                boost::scoped_ptr<control_lattice> f(
                        new control_lattice(cmin, cmax, grid, coo_begin, coo_end, val.begin())
                        );

                res = f->residual(coo_begin, coo_end, val.begin());
                f->append_refined(*psi);
                psi.swap(f);
            }
        }
};

} // namespace mba

#endif
