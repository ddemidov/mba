#include "golang_mba.h"
#include <array>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <cstring>
#include "mba.hpp"

mba::MBA<2> *GetMba_2(void *ptr) {
    return reinterpret_cast<mba::MBA<2> *>(ptr);
}

void* Mba_2_New(double *lo,
            double *hi,
            int64_t *grid,
            int coo_len,
            double *x,
            double *y,
            double *sval,
            int max_levels,
            double tol,
            double min_fill
) {

    std::vector<mba::point<2> > coo;
    coo.reserve(coo_len);

    for (int i = 0; i < coo_len; i++) {
        coo.push_back({x[i], y[i]});
    }

    mba::point<2> cmin, cmax;
    std::copy_n(lo, 2, std::begin(cmin));
    std::copy_n(hi, 2, std::begin(cmax));

    mba::index<2> grid_int;
    grid_int[0] = size_t(grid[0]);
    grid_int[1] = size_t(grid[1]);

    auto initial = mba::linear_approximation<2>(coo.begin(), coo.end(), sval);

    auto m = new mba::MBA<2>(
                cmin, cmax, grid_int,
                coo.begin(), coo.end(), sval,
                max_levels, tol, min_fill, initial
                );
    return m;
}

void Mba_2_Destroy(void* mba) {
    delete GetMba_2(mba);
    return;
}

char *Mba_2_String(void *mba) {
    // encode to stream
    std::ostringstream s;
    s << *GetMba_2(mba);

    // get string and length
    auto str = s.str();
    auto n = str.length();

    // copy into newly allocated c string
    char *cstr = (char *)malloc(n + 1);
    std::memcpy(cstr, str.data(), n);
    cstr[n] = '\0';

    // caller should call free!
    return cstr;
}

void Mba_2_Apply(void *mba, int coo_len, double *x, double *y, double *res) {
    auto m = GetMba_2(mba);
    for (int i = 0; i < coo_len; i++) {
        res[i] = (*m)(mba::point<2>{x[i], y[i]});
    }
    return;
}
