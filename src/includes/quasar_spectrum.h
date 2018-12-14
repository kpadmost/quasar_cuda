#ifndef QUASAR_SPECTRUM_H
#define QUASAR_SPECTRUM_H
#include <builtin_types.h>

#define ASTRO_OBJ_SIZE 4096
extern "C"
void generateWavelengths(
    double* h_output,
    const double4* params,
    const size_t spectrum_number
);

extern "C"
void singleInterpolation(
    double *h_matrix_x,
    double *h_matrix_y,
    size_t *h_sizes_x,
    const size_t size, // number of quasars(coln)
    double *h_matrix_s,
    double *h_matrix_t,
    const size_t size_s,
    double *h_output
);

#endif
