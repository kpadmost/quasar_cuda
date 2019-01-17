#ifndef TOOLS_H
#define TOOLS_H

#include "../kernels/cuda_rcpp_common.h"

extern "C"
void calculateReglin(
    const double *h_x,
    const double *h_y,
    const size_t *h_col_sizes,
    double8 *h_reglin_results,
    const size_t width,
    const size_t height
);

extern "C"
void calculateChisq(
    const double *h_x,
    const double *h_y,
    const double *h_e,
    const size_t *h_col_sizes,
    double *h_results,
    const size_t width,
    const size_t height
);

#endif