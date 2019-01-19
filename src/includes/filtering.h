#ifndef FILTERING_H
#define FILTERING_H
#include <builtin_types.h>
extern "C"
void countNotInf(
  const double *h_input,
  size_t *h_output,
  const size_t width,
  const size_t height
);

extern "C"
void copyIfNotInf(
    const double *h_input,
    double *h_output,
    const size_t width,
    const size_t height,
    const size_t newHeight
);

extern "C"
void filterNonpositive(
    double *h_spectrums,
    double *h_wavelength,
    double *h_errors,
    const size_t *real_sizes,
    const size_t width
);

extern "C"
void filterWithWavelengthsWindows(
    double *h_wavelengths,
    double *h_spectrums,
    double *h_errors,
    const size_t *h_sizes,
    const double2 *h_windows,
    const size_t windows_number,
    const size_t width,
    const size_t height
);

extern "C"
void filterWithParagon(
    const double *h_paragon,
    double *h_a,
    double *h_b,
    size_t *h_sizes,
    const size_t width,
    const size_t height
);

#endif