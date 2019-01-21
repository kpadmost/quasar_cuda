#ifndef GAUSSIAN_H
#define GAUSSIAN_H


#define MAX_FITGAUSSIAN_ITER 500

extern "C"
void fitGaussian(
    const double *h_x,
    const double *h_y,
    const size_t *h_col_sizes,
    double4 *results,
    const size_t width,
    const size_t height,
    const int max_iterations
);

extern "C"
void calculateFwhm(
    const double4 *h_fit_params,
    double *h_fwhm,
    const uint width
);

extern "C"
void calculateGaussianChisq(
    const double *h_x,
    const double *h_y,
    const double *h_e,
    const double4 *h_gparams,
    const uint *h_col_sizes,
    double *h_chisq,
    const uint width,
    const uint height
);

extern "C"
void calculateGaussian(
    const double *h_x,
    double *h_y,
    const double4 *h_gaussian,
    const uint *h_sizes,
    const uint width,
    const uint height
);

#endif