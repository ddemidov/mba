#ifndef MBA_HPP
#define MBA_HPP

/*
The MIT License

Copyright (c) 2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   mba.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Scattered data interpolation with multilevel B-Splines.
 */

#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <cassert>

#ifdef MBA_VERBOSE
#  include <iostream>
#  include <iomanip>
#endif

namespace mba {

/// Helper function for std::array creation.
template <class T, class... X>
inline std::array<T, sizeof...(X)> make_array(X... x) {
    std::array<T, sizeof...(X)> p = {static_cast<T>(x)...};
    return p;
}

/// Default initial grid size.
template <size_t N>
inline std::array<size_t, N> default_grid(
        const std::array<double, N> &xmin,
        const std::array<double, N> &xmax
        )
{
    size_t imin = 0;
    std::array<double, N> delta = {xmax[0] - xmin[0]};

    for(size_t k = 1; k < N; ++k) {
        delta[k] = xmax[k] - xmin[k];
        if (delta[k] < delta[imin])
            imin = k;
    }

    std::array<size_t, N> grid;

    for(size_t k = 0; k < N; ++k)
        grid[k] = static_cast<size_t>(2 * delta[k] / delta[imin] + 0.5);

    return grid;
}

/// Scattered data interpolation with multilevel B-Splines.
/**
 * This is an implementation of the MBA algorithm from [1]. This is a fast
 * algorithm for scattered N-dimensional data interpolation and approximation.
 *
 * [1] S. Lee, G. Wolberg, and S. Y. Shin. Scattered data interpolation with
 *     multilevel B-Splines. IEEE Transactions on Visualization and
 *     Computer Graphics, 3:228–244, 1997.
 */
template <size_t NDIM>
class cloud {
    public:
        typedef std::array<double, NDIM> point;
        typedef std::array<size_t, NDIM> index;

    private:
        class clattice;
        std::unique_ptr<clattice> psi;  // Control lattice.

    public:
        /**
         * \param cmin   corner of bounding box with smallest coordinates.
         * \param cmax   corner of bounding box with largest coordinates.
         * \param coo    coordinates of data points.
         * \param val    values of data points.
         * \param grid   initial control lattice size (excluding boundary points).
         * \param levels number of levels in hierarchy.
         * \param tol    stop if residual is less than this.
         */
        cloud(
                const point &cmin, const point &cmax,
                const std::vector<point> &coo, std::vector<double> val,
                std::array<size_t, NDIM> grid, size_t levels = 8, double tol = 1e-8
             )
        {
#ifndef NDEBUG
            assert(coo.size() == val.size());

            for(size_t k = 0; k < NDIM; ++k)
                assert(grid[k] > 1);
#endif

            double res0 = std::accumulate(val.begin(), val.end(), 0.0,
                    [](double sum, double v) { return sum + v * v; });

            psi.reset( new clattice(cmin, cmax, grid, coo, val) );
            double res = psi->update_data(coo, val);
#ifdef MBA_VERBOSE
            std::cout << "level  0: res = " << std::scientific << res << std::endl;
#endif

            for (size_t k = 1; (res > res0 * tol) && (k < levels); ++k) {
                for(size_t d = 0; d < NDIM; ++d) grid[d] = 2 * grid[d] - 1;

                clattice *f = new clattice(cmin, cmax, grid, coo, val);
                res = f->update_data(coo, val);

#ifdef MBA_VERBOSE
                std::cout << "level " << std::setw(2) << k << std::scientific << ": res = " << res << std::endl;
#endif

                f->append_refined(*psi);
                psi.reset(f);
            }

        }

        /// Get interpolated value at given position.
        double operator()(const point &p) const {
            return (*psi)(p);
        }

        /// Get interpolated value at given position.
        template <class... X>
        typename std::enable_if< sizeof...(X) == NDIM, double >::type
        operator()(X... x) const {
            return (*psi)(make_array<double>(x...));
        }
    private:
        /// Control lattice
        class clattice {
            private:
                point cmin, hinv;
                index n, stride;
                std::vector<double> phi;

            public:
                // Control lattice initialization.
                clattice(
                        const point &c0, const point &cmax, std::array<size_t, NDIM> grid,
                        const std::vector<point> &coo, const std::vector<double> &val
                        ) : cmin(c0), n(grid)
                {
                    for(size_t d = 0; d < NDIM; ++d) {
                        hinv[d] = (grid[d] - 1) / (cmax[d] - cmin[d]);
                        cmin[d] -= 1 / hinv[d];
                        n[d]    += 2;
                    }

                    stride[NDIM - 1] = 1;
                    for(size_t d = NDIM - 1; d--; )
                        stride[d] = stride[d + 1] * n[d + 1];

                    std::vector<double> delta(n[0] * stride[0], 0.0);
                    std::vector<double> omega(n[0] * stride[0], 0.0);

                    auto p = coo.begin();
                    auto v = val.begin();
                    for(; p != coo.end(); ++p, ++v) {
                        if (!contained(c0, cmax, *p)) continue;

                        index i;
                        point s;

                        bool valid = true;

                        for(size_t d = 0; d < NDIM; ++d) {
                            double u = ((*p)[d] - cmin[d]) * hinv[d];
                            i[d] = floor(u) - 1;
                            s[d] = u - floor(u);
                        }

                        std::array<double, power<4, NDIM>::value> w;
                        double sw2 = 0;

                        for(scounter<4, NDIM> d; d.valid(); ++d) {
                            double buf = 1;
                            for(size_t k = 0; k < NDIM; ++k)
                                buf *= B(d[k], s[k]);

                            w[d] = buf;
                            sw2 += buf * buf;
                        }

                        for(scounter<4, NDIM> d; d.valid(); ++d) {
                            double phi = (*v) * w[d] / sw2;

                            size_t idx = 0;
                            for(size_t k = 0; k < NDIM; ++k) {
                                assert(i[k] + d[k] < n[k]);

                                idx += (i[k] + d[k]) * stride[k];
                            }

                            double w2 = w[d] * w[d];

                            assert(idx < delta.size());

                            delta[idx] += w2 * phi;
                            omega[idx] += w2;
                        }
                    }

                    phi.resize(omega.size());

                    for(auto w = omega.begin(), d = delta.begin(), f = phi.begin();
                            w != omega.end();
                            ++w, ++d, ++f
                       )
                    {
                        if (fabs(*w) < 1e-32)
                            *f = 0;
                        else
                            *f = (*d) / (*w);
                    }
                }

                // Get interpolated value at given position.
                double operator()(std::array<double, NDIM> p) const {
                    index i;
                    point s;

                    for(size_t d = 0; d < NDIM; ++d) {
                        double u = (p[d] - cmin[d]) * hinv[d];
                        i[d] = floor(u) - 1;
                        s[d] = u - floor(u);
                    }

                    double f = 0;

                    for(scounter<4, NDIM> d; d.valid(); ++d) {
                        double w = 1;
                        for(size_t k = 0; k < NDIM; ++k)
                            w *= B(d[k], s[k]);

                        f += w * get(i, d);
                    }

                    return f;
                }

                // Subtract interpolated values from data points.
                double update_data(const std::vector<point> &coo, std::vector<double> &val) const {
                    auto c = coo.begin();
                    auto v = val.begin();

                    double res = 0;

                    for(; c != coo.end(); ++c, ++v) {
                        *v -= (*this)(*c);

                        res += (*v) * (*v);
                    }

                    return res;
                }

                // Refine r and append it to the current control lattice.
                void append_refined(const clattice &r) {
                    static const std::array<double, 5> s = {
                        0.125, 0.500, 0.750, 0.500, 0.125
                    };

                    for(dcounter<NDIM> i(r.n); i.valid(); ++i) {
                        double f = r.phi[i];
                        for(scounter<5, NDIM> d; d.valid(); ++d) {
                            index j;
                            bool skip = false;
                            size_t idx = 0;
                            for(size_t k = 0; k < NDIM; ++k) {
                                j[k] = 2 * i[k] + d[k] - 3;
                                if (j[k] >= n[k]) { skip = true; break; }

                                idx += j[k] * stride[k];
                            }

                            if (skip) continue;

                            double c = 1;
                            for(size_t k = 0; k < NDIM; ++k) c *= s[d[k]];

                            phi[idx] += f * c;
                        }
                    }
                }
            private:
                clattice() {}

                // Compile time value of N^M.
                template <size_t N, size_t M>
                struct power : std::integral_constant<size_t, N * power<N, M-1>::value> {};

                template <size_t N>
                struct power<N, 0> : std::integral_constant<size_t, 1> {};

                // Nested loop counter of compile-time size (M loops of size N).
                template <size_t N, size_t M>
                class scounter {
                    public:
                        scounter() : idx(0) {
                            std::fill(i.begin(), i.end(), static_cast<size_t>(0));
                        }

                        size_t operator[](size_t d) const {
                            return i[d];
                        }

                        scounter& operator++() {
                            for(size_t d = M; d--; ) {
                                if (++i[d] < N) break;
                                i[d] = 0;
                            }

                            ++idx;

                            return *this;
                        }

                        operator size_t() const {
                            return idx;
                        }

                        bool valid() const {
                            return idx < power<N, M>::value;
                        }
                    private:
                        size_t idx;
                        std::array<size_t, M> i;
                };

                // Nested loop counter of run-time size (M loops of given sizes).
                template <size_t M>
                class dcounter {
                    public:
                        dcounter(const std::array<size_t, M> &N)
                            : idx(0), N(N),
                              size(std::accumulate(N.begin(), N.end(),
                                        static_cast<size_t>(1), std::multiplies<size_t>()))
                        {
                            std::fill(i.begin(), i.end(), static_cast<size_t>(0));
                        }

                        size_t operator[](size_t d) const {
                            return i[d];
                        }

                        dcounter& operator++() {
                            for(size_t d = M; d--; ) {
                                if (++i[d] < N[d]) break;
                                i[d] = 0;
                            }

                            ++idx;

                            return *this;
                        }

                        operator size_t() const {
                            return idx;
                        }

                        bool valid() const {
                            return idx < size;
                        }
                    private:
                        size_t idx, size;
                        std::array<size_t, M> N, i;
                };

                // Value of k-th B-Spline at t.
                static inline double B(size_t k, double t) {
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

                // x is within [xmin, xmax].
                static bool contained(
                        const point &xmin, const point &xmax, const point &x)
                {
                    for(size_t d = 0; d < NDIM; ++d) {
                        static const double eps = 1e-16;

                        if (x[d] - eps <  xmin[d]) return false;
                        if (x[d] + eps >= xmax[d]) return false;
                    }

                    return true;
                }

                // Get value of phi at index (i + d).
                template <class Shift>
                inline double get(const index &i, const Shift &d) const {
                    size_t idx = 0;

                    for(size_t k = 0; k < NDIM; ++k) {
                        size_t j = i[k] + d[k];

                        if (j >= n[k]) return 0;
                        idx += j * stride[k];
                    }

                    return phi[idx];
                }
        };
};

} //namespace mba

#endif
