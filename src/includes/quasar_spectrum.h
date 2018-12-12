#ifndef QUASAR_SPECTRUM_H
#define QUASAR_SPECTRUM_H
#include <builtin_types.h>

#define ASTRO_OBJ_SIZE 4096
extern "C"
void generateWavelengths(
    double* h_output,
    double4* params,
    const size_t spectrum_number
);

#endif
