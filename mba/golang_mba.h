#pragma once
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>

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
);
void Mba_2_Destroy(void* mba);

// Returns the string repr of the MBA.
// It is the caller's responsibility to free the string
char *Mba_2_String(void *mba);

void Mba_2_Apply(void *mba, int coo_len, double *x, double *y, double *res);

#ifdef __cplusplus
}  // extern "C"
#endif